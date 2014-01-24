/*
 * nodecutunrooted.h
 *
 *  Created on: 28-apr-2013
 *      Author: M. El-Kebir
 */

#ifndef NODECUTUNROOTED_H
#define NODECUTUNROOTED_H

#include <ilcplex/ilocplex.h>
#include <ilcplex/ilocplexi.h>
#include <lemon/bfs.h>
#include <lemon/adaptors.h>
#include <lemon/preflow.h>
#include <lemon/tolerance.h>
#include <lemon/time_measure.h>
#include <lemon/edmonds_karp.h>
#include <lemon/smart_graph.h>
#include <vector>
#include <set>
#include <queue>
#include <stack>
#include "bk_alg.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class NodeCutUnrootedCallback : public IloCplex::LazyConstraintCallbackI
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;

  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  typedef lemon::SmartDigraph Digraph;
  typedef typename Digraph::Arc DiArc;
  typedef typename Digraph::Node DiNode;
  typedef typename Digraph::NodeIt DiNodeIt;
  typedef typename Digraph::ArcIt DiArcIt;
  typedef typename Digraph::InArcIt DiInArcIt;
  typedef typename Digraph::OutArcIt DiOutArcIt;

  typedef typename Digraph::template NodeMap<Node> DiNodeNodeMap;
  typedef typename Graph::template NodeMap<DiNode> NodeDiNodeMap;
  typedef typename Graph::template NodeMap<DiArc> NodeDiArcMap;
  typedef typename Digraph::ArcMap<double> CapacityMap;

  NodeCutUnrootedCallback(IloEnv env,
                          IloBoolVarArray x,
                          IloBoolVarArray y,
                          const Graph& g,
                          const WeightNodeMap& weight,
                          const IntNodeMap& nodeMap,
                          const IntArcMap& arcMap,
                          int n,
                          int m,
                          int maxNumberOfCuts,
                          const IntNodeMap& comp)
    : IloCplex::LazyConstraintCallbackI(env)
    , _x(x)
    , _y(y)
    , _g(g)
    , _weight(weight)
    , _nodeMap(nodeMap)
    , _arcMap(arcMap)
    , _n(n)
    , _m(m)
    , _maxNumberOfCuts(maxNumberOfCuts)
    , _comp(comp)
    , _tol()
    , _h()
    , _cap(_h)
    , _tmpCap(_h)
    , _g2h1(_g)
    , _g2h2(_g)
    , _h2g(_h)
    , _g2hRootArc(_g)
    , _root(lemon::INVALID)
  {
    init();
  }

  NodeCutUnrootedCallback(const NodeCutUnrootedCallback& other)
    : IloCplex::LazyConstraintCallbackI(other)
    , _x(other._x)
    , _y(other._y)
    , _g(other._g)
    , _weight(other._weight)
    , _nodeMap(other._nodeMap)
    , _arcMap(other._arcMap)
    , _n(other._n)
    , _m(other._m)
    , _maxNumberOfCuts(other._maxNumberOfCuts)
    , _comp(other._comp)
    , _tol(other._tol)
    , _h()
    , _cap(_h)
    , _tmpCap(_h)
    , _g2h1(_g)
    , _g2h2(_g)
    , _h2g(_h)
    , _g2hRootArc(_g)
    , _root(lemon::INVALID)
  {
    typename Digraph::template NodeMap<DiNode> nodeMap(other._h);
    typename Digraph::template ArcMap<DiArc> arcMap(other._h);

    lemon::digraphCopy(other._h, _h)
        .nodeRef(nodeMap)
        .arcRef(arcMap)
        .nodeMap(other._h2g, _h2g)
        .arcMap(other._cap, _cap)
        .arcMap(other._tmpCap, _tmpCap)
        .run();

    for (NodeIt v(_g); v != lemon::INVALID; ++v)
    {
      DiNode v1 = other._g2h1[v];
      DiNode v2 = other._g2h2[v];
      DiArc root_arc_v = other._g2hRootArc[v];

      _g2h1[v] = nodeMap[v1];
      _g2h2[v] = nodeMap[v2];
      _g2hRootArc[v] = arcMap[root_arc_v];
    }

    _root = nodeMap[other._root];
  }

  virtual ~NodeCutUnrootedCallback()
  {
  }

protected:
  virtual void main();
  virtual IloCplex::CallbackI* duplicateCallback() const
  {
    return (new (getEnv()) NodeCutUnrootedCallback(*this));
  }

private:
  IloBoolVarArray _x;
  IloBoolVarArray _y;
  const Graph& _g;
  const WeightNodeMap& _weight;
  const IntNodeMap& _nodeMap;
  const IntArcMap& _arcMap;
  const int _n;
  const int _m;
  const int _maxNumberOfCuts;
  const IntNodeMap& _comp;
  const lemon::Tolerance<double> _tol;

  Digraph _h;
  CapacityMap _cap;
  CapacityMap _tmpCap;
  NodeDiNodeMap _g2h1;
  NodeDiNodeMap _g2h2;
  DiNodeNodeMap _h2g;
  NodeDiArcMap _g2hRootArc;
  DiNode _root;

  typedef lemon::ReverseDigraph<const Digraph> ReverseDigraph;
  typedef lemon::Preflow<ReverseDigraph, CapacityMap> RevPreFlowAlg;
  typedef lemon::Preflow<Digraph, CapacityMap> PreFlowAlg;

  typedef std::set<Node> NodeSet;
  typedef typename NodeSet::const_iterator NodeSetIt;
  typedef std::set<DiArc> DiArcSet;
  typedef typename DiArcSet::const_iterator DiArcSetIt;

  typedef std::queue<Node> NodeQueue;
  typedef std::queue<DiNode> DiNodeQueue;
  typedef std::set<DiNode> DiNodeSet;
  //typedef typename lemon::ResidualDigraph<const Graph, const DoubleArcMap, const DoubleArcMap, typename PreFlowAlg::Tolerance> ResidualGraphType;
  //typedef typename lemon::ReverseDigraph<ResidualGraphType> ReverseResidualGraphType;

  void init()
  {
    // we initialize _h:
    // - for every node i, there will be two nodes i1 and i2
    //   connected by an arc from i1 to i2
    // - for every edge (i,j) there are two arcs
    //   in h: (i2,j1) and (j2,i1) with capacties 1
    _h.clear();
    _root = _h.addNode();
    for (NodeIt i(_g); i != lemon::INVALID; ++i)
    {
      DiNode i1 = _h.addNode();
      DiNode i2 = _h.addNode();
      _g2h1[i] = i1;
      _g2h2[i] = i2;
      _h2g[i1] = i;
      _h2g[i2] = i;

      DiArc i1i2 = _h.addArc(i1, i2);
      _cap[i1i2] = 0;

      DiArc ri1 = _h.addArc(_root, i1);
      _g2hRootArc[i] = ri1;
      _cap[ri1] = 1;
    }

    for (EdgeIt e(_g); e != lemon::INVALID; ++e)
    {
      Node i = _g.u(e);
      Node j = _g.v(e);
      DiNode i1 = _g2h1[i];
      DiNode i2 = _g2h2[i];
      DiNode j1 = _g2h1[j];
      DiNode j2 = _g2h2[j];

      DiArc i2j1 = _h.addArc(i2, j1);
      DiArc j2i1 = _h.addArc(j2, i1);
      _cap[i2j1] = _cap[j2i1] = 1;
    }
  }

  void computeCapacities(CapacityMap& capacity)
  {
    // cap((i,j)) = x_i
    for (NodeIt v(_g); v != lemon::INVALID; ++v)
    {
      double val = 1e-9 + getValue(_x[_nodeMap[v]]);
      DiNode v1 = _g2h1[v];
      capacity[DiOutArcIt(_h, v1)] = val;

      // cap((r,i)) = y_i
      val = 1e-9 + getValue(_y[_nodeMap[v]]);
      capacity[_g2hRootArc[v]] = val;
    }
  }

  void addViolatedConstraint(Node target, const NodeSet& dS, const NodeSet& S)
  {
    IloExpr expr(getEnv());

    for (NodeSetIt nodeIt = dS.begin(); nodeIt != dS.end(); nodeIt++)
    {
      expr += _x[_nodeMap[*nodeIt]];
    }

    for (NodeSetIt nodeIt = S.begin(); nodeIt != S.end(); nodeIt++)
    {
      expr += _y[_nodeMap[*nodeIt]];
    }

    add(_x[_nodeMap[target]] <= expr);
  }

  template<class DGR, class T>
  void determineCutSet(const DGR& h,
                       const T& pf,
                       NodeSet& dS,
                       NodeSet& S,
                       bool reverse)
  {
    for (NodeIt v(_g); v != lemon::INVALID; ++v)
    {
      DiNode v1 = _g2h1[v];
      DiNode v2 = _g2h2[v];
      if (pf.minCut(v1) != pf.minCut(v2))
      {
        dS.insert(v);
      }
    }

    if (reverse)
    {
      for (typename DGR::InArcIt a(h, _root); a != lemon::INVALID; ++a)
      {
        if (pf.minCut(h.source(a)) != pf.minCut(h.target(a)))
        {
          S.insert(_h2g[h.source(a)]);
        }
      }
    }
    else
    {
      for (typename DGR::OutArcIt a(h, _root); a != lemon::INVALID; ++a)
      {
        if (pf.minCut(h.source(a)) != pf.minCut(h.target(a)))
        {
          S.insert(_h2g[h.target(a)]);
        }
      }
    }
  }
};

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void NodeCutUnrootedCallback<GR, NWGHT, NLBL, EWGHT>::main()
{
  //lemon::Timer t;
  computeCapacities(_cap);
  lemon::mapCopy(_h, _cap, _tmpCap);

  ReverseDigraph revH = lemon::reverseDigraph(_h);
  PreFlowAlg pf(_h, _tmpCap, _root, _root);
  RevPreFlowAlg revPf(revH, _tmpCap, _root, _root);

  int nCuts = 0;
  for (NodeIt i(_g); i != lemon::INVALID && (nCuts <= _maxNumberOfCuts || _maxNumberOfCuts == -1); ++i)
  {
    double x_i_value = getValue(_x[_nodeMap[i]]);
    if (!_tol.nonZero(x_i_value))
      continue;

    pf.target(_g2h2[i]);

    bool done = false;
    bool update = false;
    while (!done)
    {
      pf.runMinCut();

      // let's see if there's a violated constraint
      // idea: we could use addLocal here...
      double minCutValue = pf.flowValue();
      if (_tol.less(minCutValue, x_i_value))
      {
        // determine N (forward)
        NodeSet fwdS, fwdDS;
        determineCutSet(_h, pf, fwdDS, fwdS, false);

        revPf.target(_root);
        revPf.source(_g2h2[i]);
        revPf.runMinCut();
        NodeSet bwdS, bwdDS;
        determineCutSet(revH, revPf, bwdDS, bwdS, false);

        // add violated constraints
        addViolatedConstraint(i, fwdDS, fwdS);
        if (fwdDS.size() != bwdDS.size() || fwdS.size() != bwdS.size() ||
            fwdDS != bwdDS || bwdS != fwdS)
        {
          addViolatedConstraint(i, bwdDS, bwdS);
          nCuts++;
        }

        // generate nested-cuts
        for (NodeSetIt nodeIt = fwdDS.begin(); nodeIt != fwdDS.end(); nodeIt++)
        {
          // update the capactity to generate nested-cuts
          _tmpCap[DiOutArcIt(_h, _g2h1[*nodeIt])] = 1;
        }

        nCuts++;
        update = true;
      }
      else
      {
        if (update)
          lemon::mapCopy(_h, _cap, _tmpCap);
        done = true;
      }
    }
  }

  //std::cerr << "Generated cuts: " << nCuts << std::endl;
  //std::cerr << "Time: " << t.realTime() << "s" << std::endl;
}

} // namespace mwcs
} // namespace nina

#endif // NODECUTUNROOTED_H
