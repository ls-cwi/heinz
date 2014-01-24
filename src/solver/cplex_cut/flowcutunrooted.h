/*
 * flowcutunrooted.h
 *
 *  Created on: 26-apr-2013
 *      Author: M. El-Kebir
 */

#ifndef FLOWCUTUNROOTED_H
#define FLOWCUTUNROOTED_H

#include <ilcplex/ilocplex.h>
#include <ilcplex/ilocplexi.h>
#include <lemon/adaptors.h>
#include <lemon/preflow.h>
#include <lemon/tolerance.h>
#include <lemon/time_measure.h>
#include <lemon/smart_graph.h>
#include <vector>
#include <set>

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class FlowCutUnrootedCallback : public IloCplex::LazyConstraintCallbackI
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
  typedef std::set<DiNode> DiNodeSet;
  typedef std::set<DiArc> DiArcSet;
  typedef typename DiArcSet::const_iterator DiArcSetIt;

  FlowCutUnrootedCallback(IloEnv env,
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
    , _g2h(_g)
    , _h2g(_h)
    , _g2hRootArc(_g)
    , _root(lemon::INVALID)
  {
    init();
  }

  FlowCutUnrootedCallback(const FlowCutUnrootedCallback& other)
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
    , _g2h(_g)
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
      DiNode v_h = other._g2h[v];
      DiArc root_arc_v = other._g2hRootArc[v];

      _g2h[v] = nodeMap[v_h];
      _g2hRootArc[v] = arcMap[root_arc_v];
    }

    _root = nodeMap[other._root];
  }

  virtual ~FlowCutUnrootedCallback()
  {
  }

protected:
  virtual void main();
  virtual IloCplex::CallbackI* duplicateCallback() const
  {
    return (new (getEnv()) FlowCutUnrootedCallback(*this));
  }

protected:
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
  NodeDiNodeMap _g2h;
  DiNodeNodeMap _h2g;
  NodeDiArcMap _g2hRootArc;
  DiNode _root;

  typedef lemon::ReverseDigraph<const Digraph> ReverseDigraph;
  typedef lemon::Preflow<ReverseDigraph, CapacityMap> RevPreFlowAlg;
  typedef lemon::Preflow<Digraph, CapacityMap> PreFlowAlg;

  void init()
  {
    _h.clear();

    _root = _h.addNode();
    _h2g[_root] = lemon::INVALID;

    for (NodeIt v(_g); v != lemon::INVALID; ++v)
    {
      DiNode v_h = _h.addNode();
      _g2h[v] = v_h;
      _h2g[v_h] = v;

      DiArc root_arc = _h.addArc(_root, v_h);
      _cap[root_arc] = 0;
      _g2hRootArc[v] = root_arc;
    }

    for (EdgeIt e(_g); e != lemon::INVALID; ++e)
    {
      Node i = _g.u(e);
      Node j = _g.v(e);

      DiNode i_h = _g2h[i];
      DiNode j_h = _g2h[j];

      _cap[_h.addArc(i_h, j_h)] = 0;
      _cap[_h.addArc(j_h, i_h)] = 0;
    }
  }

  void computeCapacities(CapacityMap& capacity) const
  {
    // cap((i,j)) = x_i
    for (NodeIt v(_g); v != lemon::INVALID; ++v)
    {
      double val = 1e-9 + getValue(_x[_nodeMap[v]]);

      DiNode v_h = _g2h[v];
      for (DiOutArcIt a(_h, v_h); a != lemon::INVALID; ++a)
      {
        capacity[a] = val;
      }

      capacity[_g2hRootArc[v]] = 1e-9 + getValue(_y[_nodeMap[v]]);
    }
  }

  template<class DGR, class T>
  DiArcSet determineCutSet(const DGR& h, const T& pf, bool reverse)
  {
    DiArcSet arcDeltaS;
    for (NodeIt v(_g); v != lemon::INVALID; ++v)
    {
      DiNode v_h = _g2h[v];

      if (pf.minCut(v_h))
      {
        for (typename DGR::OutArcIt a(h, v_h); a != lemon::INVALID; ++a)
        {
          DiNode w_h = h.target(a);
          if (!pf.minCut(w_h))
          {
            arcDeltaS.insert(a);
          }
        }
      }
    }

    if (reverse)
    {
      for (typename DGR::InArcIt a(h, _root); a != lemon::INVALID; ++a)
      {
        if (pf.minCut(h.source(a)) != pf.minCut(h.target(a)))
        {
          arcDeltaS.insert(a);
        }
      }
    }
    else
    {
      for (typename DGR::OutArcIt a(h, _root); a != lemon::INVALID; ++a)
      {
        if (pf.minCut(h.source(a)) != pf.minCut(h.target(a)))
        {
          arcDeltaS.insert(a);
        }
      }
    }

    return arcDeltaS;
  }

  template<class T>
  void addViolatedConstraint(const T& h, Node target, const DiArcSet& deltaS)
  {
    IloExpr expr(getEnv());
    DiNodeSet xSet;

    //bool first = true;
    for (DiArcSetIt arcIt = deltaS.begin(); arcIt != deltaS.end(); arcIt++)
    {
      DiNode u_h = h.source(*arcIt);
      DiNode v_h = h.target(*arcIt);

      if (u_h == _root)
      {
        expr += _y[_nodeMap[_h2g[v_h]]];
        //if (!first)
        //{
        //  std::cerr << " + ";
        //}
        //else
        //{
        //  first = false;
        //}
        //std::cerr << _y[_nodeMap[_h2g[v_h]]].getName();
      }
      else if (xSet.find(u_h) == xSet.end())
      {
        expr += _x[_nodeMap[_h2g[u_h]]];
        xSet.insert(u_h);
        //if (!first)
        //{
        //  std::cerr << " + ";
        //}
        //else
        //{
        //  first = false;
        //}
        //std::cerr << _x[_nodeMap[_h2g[u_h]]].getName();
      }
    }

    add(_x[_nodeMap[target]] <= expr);
    //std::cerr << " >= " << _x[_nodeMap[target]].getName() << std::endl;
  }
};

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void FlowCutUnrootedCallback<GR, NWGHT, NLBL, EWGHT>::main()
{
  lemon::Timer t;

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

    pf.target(_g2h[i]);

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
        // determine deltaS (forward)
        DiArcSet fwdArcDeltaS = determineCutSet(_h, pf, false);

        revPf.target(_root);
        revPf.source(_g2h[i]);
        revPf.runMinCut();
        DiArcSet bwdArcDeltaS = determineCutSet(revH, revPf, true);

        // add violated constraints
        addViolatedConstraint(_h, i, fwdArcDeltaS);
        if (bwdArcDeltaS.size() != fwdArcDeltaS.size() || bwdArcDeltaS != fwdArcDeltaS)
        {
          addViolatedConstraint(revH, i, bwdArcDeltaS);
          nCuts++;
        }

        // generate nested-cuts
        for (DiArcSetIt arcIt = fwdArcDeltaS.begin(); arcIt != fwdArcDeltaS.end(); arcIt++)
        {
          // update the capactity to generate nested-cuts (skip ze root node)
          if (_h.source(*arcIt) != _root)
            _tmpCap[*arcIt] = 1;
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

#endif // FLOWCUTUNROOTED_H
