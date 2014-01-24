/*
 * nodecut.h
 *
 *  Created on: 19-apr-2013
 *      Author: M. El-Kebir
 */

#ifndef NODECUT_H
#define NODECUT_H

#include <ilcplex/ilocplex.h>
#include <ilcplex/ilocplexi.h>
#include <lemon/bfs.h>
#include <lemon/adaptors.h>
#include <lemon/preflow.h>
#include <lemon/tolerance.h>
#include <lemon/time_measure.h>
#include <lemon/smart_graph.h>
#include <vector>
#include <set>
#include <queue>
#include <stack>

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class NodeCutCallback : public IloCplex::LazyConstraintCallbackI
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

  typedef typename Graph::template NodeMap<DiNode> NodeNodeMap;
  typedef typename Digraph::ArcMap<double> CapacityMap;

  NodeCutCallback(IloEnv env,
                  IloBoolVarArray x,
                  const Graph& g,
                  const WeightNodeMap& weight,
                  Node root,
                  const IntNodeMap& nodeMap,
                  const IntArcMap& arcMap,
                  int n,
                  int m,
                  int maxNumberOfCuts,
                  const IntNodeMap& comp)
    : IloCplex::LazyConstraintCallbackI(env)
    , _x(x)
    , _g(g)
    , _weight(weight)
    , _root(root)
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
  {
    init();
  }

  NodeCutCallback(const NodeCutCallback& other)
    : IloCplex::LazyConstraintCallbackI(other)
    , _x(other._x)
    , _g(other._g)
    , _weight(other._weight)
    , _root(other._root)
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
  {
    typename Digraph::template NodeMap<DiNode> nodeMap(other._h);
    lemon::digraphCopy(other._h, _h)
        .nodeRef(nodeMap)
        .arcMap(other._cap, _cap)
        .arcMap(other._tmpCap, _tmpCap)
        .run();

    for (NodeIt v(_g); v != lemon::INVALID; ++v)
    {
      DiNode v1 = other._g2h1[v];
      DiNode v2 = other._g2h2[v];

      _g2h1[v] = nodeMap[v1];
      _g2h2[v] = nodeMap[v2];
    }
  }

  virtual ~NodeCutCallback()
  {
  }

protected:
  virtual void main();
  virtual IloCplex::CallbackI* duplicateCallback() const
  {
    return (new (getEnv()) NodeCutCallback(*this));
  }

private:
  IloBoolVarArray _x;
  const Graph& _g;
  const WeightNodeMap& _weight;
  const Node _root;
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
  NodeNodeMap _g2h1;
  NodeNodeMap _g2h2;

  typedef lemon::ReverseDigraph<const Digraph> ReverseDigraph;
  typedef lemon::Preflow<ReverseDigraph, CapacityMap> RevPreFlowAlg;
  typedef lemon::Preflow<Digraph, CapacityMap> PreFlowAlg;
  typedef std::set<Node> NodeSet;
  typedef typename NodeSet::const_iterator NodeSetIt;
  typedef std::set<DiArc> DiArcSet;
  typedef typename DiArcSet::const_iterator DiArcSetIt;
  //typedef typename lemon::ResidualDigraph<const Graph, const DoubleArcMap, const DoubleArcMap, typename PreFlowAlg::Tolerance> ResidualGraphType;
  //typedef typename lemon::ReverseDigraph<ResidualGraphType> ReverseResidualGraphType;

  void init()
  {
    // we initialize _h:
    // - for every non-root node i, there will be two nodes i1 and i2
    //   connected by an arc from i1 to i2
    // - the root r only has one counterpart in h, i.e. r'
    // - for every edge (i,j) where i != r and j != r there are two arcs
    //   in h: (i2,j1) and (j2,i1) with capacties 1
    // - for every edge (r, i) there is an arc (r,i1) in h with capacity 1
    _h.clear();
    for (NodeIt i(_g); i != lemon::INVALID; ++i)
    {
      if (i == _root)
      {
        DiNode cpy_i = _h.addNode();
        _g2h1[i] = _g2h2[i] = cpy_i;
      }
      else
      {
        DiNode i1 = _h.addNode();
        DiNode i2 = _h.addNode();
        _g2h1[i] = i1;
        _g2h2[i] = i2;

        DiArc a = _h.addArc(i1, i2);
        _cap[a] = 0;
      }
    }

    for (EdgeIt e(_g); e != lemon::INVALID; ++e)
    {
      Node i = _g.u(e);
      Node j = _g.v(e);
      if (i == _root)
      {
        DiArc a = _h.addArc(_g2h2[i], _g2h1[j]);
        _cap[a] = 1;
      }
      else if (j == _root)
      {
        DiArc a = _h.addArc(_g2h2[j], _g2h1[i]);
        _cap[a] = 1;
      }
      else
      {
        DiNode i1 = _g2h1[i];
        DiNode i2 = _g2h2[i];
        DiNode j1 = _g2h1[j];
        DiNode j2 = _g2h2[j];

        DiArc i2j1 = _h.addArc(i2, j1);
        DiArc j2i1 = _h.addArc(j2, i1);
        _cap[i2j1] = _cap[j2i1] = 1;
      }
    }
  }

  void computeCapacities(CapacityMap& capacity) const
  {
    // cap((i,j)) = x_i
    for (NodeIt v(_g); v != lemon::INVALID; ++v)
    {
      if (v != _root)
      {
        DiNode v1 = _g2h1[v];
        double val = getValue(_x[_nodeMap[v]]);
        capacity[DiOutArcIt(_h, v1)] = 1e-9 + val;
      }
    }
  }

  void addViolatedConstraint(Node target, const NodeSet& N)
  {
    IloExpr expr(getEnv());

    for (NodeSetIt nodeIt = N.begin(); nodeIt != N.end(); nodeIt++)
    {
      expr += _x[_nodeMap[*nodeIt]];
    }

    add(_x[_nodeMap[target]] <= expr);
  }

  template<class T>
  NodeSet determineCutSet(const T& pf)
  {
    NodeSet result;
    for (NodeIt v(_g); v != lemon::INVALID; ++v)
    {
      if (v != _root)
      {
        DiNode v1 = _g2h1[v];
        DiNode v2 = _g2h2[v];
        if (pf.minCut(v1) != pf.minCut(v2))
        {
          result.insert(v);
        }
      }
    }
    return result;
  }
};

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void NodeCutCallback<GR, NWGHT, NLBL, EWGHT>::main()
{
  lemon::Timer t;

  computeCapacities(_cap);
  lemon::mapCopy(_h, _cap, _tmpCap);

  ReverseDigraph revH = lemon::reverseDigraph(_h);
  PreFlowAlg pf(_h, _tmpCap, _g2h1[_root], _g2h1[_root]);
  RevPreFlowAlg revPf(revH, _tmpCap, _g2h1[_root], _g2h1[_root]);

  const int rootComp = _comp[_root];
  int nCuts = 0;
  for (NodeIt i(_g); i != lemon::INVALID && (nCuts <= _maxNumberOfCuts || _maxNumberOfCuts == -1); ++i)
  {
    double x_i_value = getValue(_x[_nodeMap[i]]);
    if (i == _root || !_tol.nonZero(x_i_value) || _comp[i] != rootComp)
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
        NodeSet fwdN = determineCutSet(pf);

        revPf.target(_g2h1[_root]);
        revPf.source(_g2h2[i]);
        revPf.runMinCut();
        NodeSet bwdN = determineCutSet(revPf);

        // add violated constraints
        addViolatedConstraint(i, fwdN);
        if (bwdN.size() != fwdN.size() || bwdN != fwdN)
        {
          addViolatedConstraint(i, bwdN);
          nCuts++;
        }

        // generate nested-cuts
        for (NodeSetIt nodeIt = fwdN.begin(); nodeIt != fwdN.end(); nodeIt++)
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

  std::cerr << "Generated cuts: " << nCuts << std::endl;
  std::cerr << "Time: " << t.realTime() << "s" << std::endl;
}

} // namespace mwcs
} // namespace nina

#endif // NODECUT_H
