/*
 * nodecutlazy.h
 *
 *  Created on: 24-feb-2014
 *      Author: M. El-Kebir
 */

#ifndef NODECUTLAZY_H
#define NODECUTLAZY_H

#include <ilcplex/ilocplex.h>
#include <ilcplex/ilocplexi.h>
#include <ilconcert/ilothread.h>
#include <lemon/tolerance.h>
#include <lemon/smart_graph.h>
#include <vector>
#include <set>
#include <queue>
#include <list>
#include "nodecut.h"

namespace nina {
namespace mwcs {
  
template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class NodeCutLazy : public IloCplex::LazyConstraintCallbackI,
                    public NodeCut<GR, NWGHT, NLBL, EWGHT>
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;
  typedef NodeCut<GR, NWGHT, NLBL, EWGHT> Parent;
  
  using Parent::_x;
  using Parent::_y;
  using Parent::_g;
  using Parent::_weight;
  using Parent::_root;
  using Parent::_nodeMap;
  using Parent::_n;
  using Parent::_m;
  using Parent::_maxNumberOfCuts;
  using Parent::_tol;
  using Parent::_pNodeBoolMap;
  using Parent::_pSubG;
  using Parent::_pComp;
  using Parent::_pMutex;
  using Parent::_epsilon;
  
  using Parent::lock;
  using Parent::unlock;
  using Parent::determineConnectedComponents;
  using Parent::separateConnectedComponent;
  
protected:
  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  
public:
  NodeCutLazy(IloEnv env,
              IloBoolVarArray x,
              IloBoolVarArray y,
              const Graph& g,
              const WeightNodeMap& weight,
              const IntNodeMap& nodeMap,
              int n,
              int m,
              int maxNumberOfCuts,
              IloFastMutex* pMutex)
    : IloCplex::LazyConstraintCallbackI(env)
    , Parent(x, y, g, weight, nodeMap, n, m, maxNumberOfCuts, pMutex)
  {
  }
  
  NodeCutLazy(const NodeCutLazy& other)
    : IloCplex::LazyConstraintCallbackI(other)
    , Parent(other)
  {
  }
  
  virtual ~NodeCutLazy()
  {
  }
};

} // namespace mwcs
} // namespace nina


#endif // NODECUTLAZY_H
