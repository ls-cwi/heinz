#ifndef FLOWMINCUT_H
#define FLOWMINCUT_H

#include <ilcplex/ilocplex.h>
#include <ilcplex/ilocplexi.h>
#include <lemon/bfs.h>
#include <lemon/adaptors.h>
#include <lemon/preflow.h>
#include <lemon/tolerance.h>
#include <lemon/time_measure.h>
#include <vector>
#include <set>
#include "flowcut.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class FlowCutMinCallback : public FlowCutCallback<GR, NWGHT, NLBL, EWGHT>
{
public:
  typedef FlowCutCallback<GR, NWGHT, NLBL, EWGHT> Parent;
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;

  TEMPLATE_GRAPH_TYPEDEFS(Graph);

  FlowCutMinCallback(IloEnv env,
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
    : Parent(env, x, g, weight, root, nodeMap, arcMap, n, m, maxNumberOfCuts, comp)
  {
  }

  virtual ~FlowCutMinCallback()
  {
  }

  using Parent::getValue;
  using Parent::getEnv;
  using Parent::add;
  using Parent::determineCutSet;
  using Parent::_x;
  using Parent::_g;
  using Parent::_weight;
  using Parent::_root;
  using Parent::_nodeMap;
  using Parent::_arcMap;
  using Parent::_n;
  using Parent::_m;
  using Parent::_maxNumberOfCuts;
  using Parent::_comp;
  using Parent::_tol;
  using Parent::_capacity;
  using Parent::_tmpCapacity;

protected:
  virtual void main();
  virtual IloCplex::CallbackI* duplicateCallback() const
  {
    return (new (getEnv()) FlowCutMinCallback(*this));
  }

  typedef typename Parent::PreFlowAlg PreFlowAlg;
  typedef typename Parent::NodeSet NodeSet;
  typedef typename Parent::NodeSetIt NodeSetIt;
  typedef typename Parent::ArcSet ArcSet;
  typedef typename Parent::ArcSetIt ArcSetIt;
  typedef typename Parent::ResidualGraphType ResidualGraphType;
  typedef typename Parent::ReverseResidualGraphType ReverseResidualGraphType;

  virtual void computeCapacities(DoubleArcMap& capacity) const
  {
    // cap((i,j)) = min{x_i, x_j}
    for (EdgeIt e(_g); e != lemon::INVALID; ++e)
    {
      Arc a1 = _g.direct(e, true);
      Arc a2 = _g.direct(e, false);
      capacity[a1] = capacity[a2] = 1e-9 + IloNum(std::min(getValue(_x[_nodeMap[_g.u(e)]]),
                                                           getValue(_x[_nodeMap[_g.v(e)]])));
    }
  }

  virtual void addViolatedConstraint(Node target, const ArcSet& deltaS)
  {
    IloExpr expr(getEnv());

    for (ArcSetIt arcIt = deltaS.begin(); arcIt != deltaS.end(); arcIt++)
    {
      if (getValue(_x[_nodeMap[_g.source(*arcIt)]]) < getValue(_x[_nodeMap[_g.target(*arcIt)]]))
        expr += _x[_nodeMap[_g.source(*arcIt)]];
      else
        expr += _x[_nodeMap[_g.target(*arcIt)]];
    }

    add(_x[_nodeMap[target]] <= expr);
  }
};

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void FlowCutMinCallback<GR, NWGHT, NLBL, EWGHT>::main()
{
  lemon::Timer t;

  computeCapacities(_capacity);
  lemon::mapCopy(_g, _capacity, _tmpCapacity);

  PreFlowAlg pf(_g, _tmpCapacity, _root, _root);

  const int rootComp = _comp[_root];
  int nCuts = 0;
  for (NodeIt i(_g); i != lemon::INVALID && (nCuts <= _maxNumberOfCuts || _maxNumberOfCuts == -1); ++i)
  {
    double x_i_value = getValue(_x[_nodeMap[i]]);
    if (i == _root || !_tol.nonZero(x_i_value) || _comp[i] != rootComp)
      continue;

    bool done = false;
    bool update = false;
    while (!done)
    {
      pf.source(_root);
      pf.target(i);
      pf.runMinCut();

      // let's see if there's a violated constraint
      // idea: we could use addLocal here...
      double minCutValue = pf.flowValue();
      if (_tol.less(minCutValue, x_i_value))
      {
        // determine deltaS (forward)
        ArcSet fwdArcDeltaS = determineCutSet(pf);

        pf.source(i);
        pf.target(_root);
        pf.runMinCut();

        ArcSet bwdArcDeltaS = determineCutSet(pf);

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

#endif // FLOWMINCUT_H
