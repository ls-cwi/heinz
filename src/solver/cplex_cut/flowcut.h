/*
 * flowcut.h
 *
 *  Created on: 19-apr-2013
 *      Author: M. El-Kebir
 */

#ifndef FLOWCUT_H
#define FLOWCUT_H

#include <ilcplex/ilocplex.h>
#include <ilcplex/ilocplexi.h>
#include <lemon/bfs.h>
#include <lemon/adaptors.h>
#include <lemon/preflow.h>
#include <lemon/tolerance.h>
#include <lemon/time_measure.h>
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
class FlowCutCallback : public IloCplex::LazyConstraintCallbackI
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;

  TEMPLATE_GRAPH_TYPEDEFS(Graph);

  FlowCutCallback(IloEnv env,
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
    , _capacity(_g)
    , _tmpCapacity(_g)
  {
  }

  FlowCutCallback(const FlowCutCallback& other)
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
    , _capacity(_g)
    , _tmpCapacity(_g)
  {
    // todo: is dit wel nodig?
    lemon::mapCopy(_g, other._capacity, _capacity);
    lemon::mapCopy(_g, other._tmpCapacity, _tmpCapacity);
  }

  virtual ~FlowCutCallback()
  {
  }

protected:
  virtual void main();
  virtual IloCplex::CallbackI* duplicateCallback() const
  {
    return (new (getEnv()) FlowCutCallback(*this));
  }

protected:
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
  DoubleArcMap _capacity;
  DoubleArcMap _tmpCapacity;

  typedef lemon::ReverseDigraph<const Graph> ReverseGraph;
  typedef lemon::Preflow<ReverseGraph, DoubleArcMap> RevPreFlowAlg;
  //typedef std::queue<Node> NodeQueue;
  typedef std::stack<Node> NodeQueue;
  typedef lemon::Preflow<Graph, DoubleArcMap> PreFlowAlg;
  typedef std::set<Node> NodeSet;
  typedef typename NodeSet::const_iterator NodeSetIt;
  typedef std::set<Arc> ArcSet;
  typedef typename ArcSet::const_iterator ArcSetIt;
  typedef typename lemon::ResidualDigraph<const Graph, const DoubleArcMap, const DoubleArcMap, typename PreFlowAlg::Tolerance> ResidualGraphType;
  typedef typename lemon::ReverseDigraph<ResidualGraphType> ReverseResidualGraphType;

  virtual void computeCapacities(DoubleArcMap& capacity) const
  {
    // cap((i,j)) = x_i
    for (ArcIt a(_g); a != lemon::INVALID; ++a)
      capacity[a] = 1e-9 + IloNum(getValue(_x[_nodeMap[_g.source(a)]]));
      //capacity[a] = IloNum(getValue(_x[_nodeMap[_g.source(a)]]));
  }

  virtual void addViolatedConstraint(Node target, const ArcSet& deltaS)
  {
    IloExpr expr(getEnv());

    for (ArcSetIt arcIt = deltaS.begin(); arcIt != deltaS.end(); arcIt++)
    {
      expr += _x[_nodeMap[_g.source(*arcIt)]];
    }

    add(_x[_nodeMap[target]] <= expr);
  }

  template<class T>
  ArcSet determineCutSet(const T& pf)
  {
    ArcSet arcDeltaS;
    for (NodeIt v(_g); v != lemon::INVALID; ++v)
    {
      if (pf.minCut(v))
      {
        for (OutArcIt a(_g, v); a != lemon::INVALID; ++a)
        {
          Node w = _g.target(a);
          if (!pf.minCut(w))
          {
            arcDeltaS.insert(a);
          }
        }
      }
    }

    return arcDeltaS;
  }
};

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void FlowCutCallback<GR, NWGHT, NLBL, EWGHT>::main()
{
  lemon::Timer t;

  computeCapacities(_capacity);
  lemon::mapCopy(_g, _capacity, _tmpCapacity);

  ReverseGraph revG = lemon::reverseDigraph(_g);
  PreFlowAlg pf(_g, _tmpCapacity, _root, _root);
  RevPreFlowAlg revPf(revG, _tmpCapacity, _root, _root);

  const int rootComp = _comp[_root];
  int nCuts = 0;
  for (NodeIt i(_g); i != lemon::INVALID && (nCuts <= _maxNumberOfCuts || _maxNumberOfCuts == -1); ++i)
  {
    double x_i_value = getValue(_x[_nodeMap[i]]);
    if (i == _root || !_tol.nonZero(x_i_value) || _comp[i] != rootComp)
      continue;

    pf.target(i);

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
        ArcSet fwdArcDeltaS = determineCutSet(pf);

        revPf.target(_root);
        revPf.source(i);
        revPf.runMinCut();
        ArcSet bwdArcDeltaS = determineCutSet(revPf);

        // add violated constraints
        addViolatedConstraint(i, fwdArcDeltaS);
        if (bwdArcDeltaS.size() != fwdArcDeltaS.size() || bwdArcDeltaS != fwdArcDeltaS)
        {
          addViolatedConstraint(i, bwdArcDeltaS);
          nCuts++;
        }

        // generate nested-cuts
        for (ArcSetIt arcIt = fwdArcDeltaS.begin(); arcIt != fwdArcDeltaS.end(); arcIt++)
        {
          // update the capactity to generate nested-cuts
          _tmpCapacity[*arcIt] = 1;
        }

        nCuts++;
        update = true;
      }
      else
      {
        if (update)
          lemon::mapCopy(_g, _capacity, _tmpCapacity);
        done = true;
      }
    }
  }

  std::cerr << "Generated cuts: " << nCuts << std::endl;
  std::cerr << "Time: " << t.realTime() << "s" << std::endl;
}

} // namespace mwcs
} // namespace nina

#endif // FLOWCUT_H
