/*
 * edgecutlazy.h
 *
 *  Created on: 11-sep-2014
 *      Author: M. El-Kebir
 */

#ifndef EDGECUTLAZY_H
#define EDGECUTLAZY_H

#include <ilcplex/ilocplex.h>
#include <ilcplex/ilocplexi.h>
#include <ilconcert/ilothread.h>
#include <lemon/tolerance.h>
#include <lemon/smart_graph.h>
#include <vector>
#include <set>
#include <queue>
#include <list>
#include "edgecut.h"

namespace nina {
namespace mwcs {
  
template<typename DGR,
         typename NWGHT = typename DGR::template NodeMap<double> >
class EdgeCutLazy : public IloCplex::LazyConstraintCallbackI,
                    public EdgeCut<DGR, NWGHT>
{
public:
  typedef DGR Graph;
  typedef NWGHT WeightNodeMap;
  typedef EdgeCut<DGR, NWGHT> Parent;
  
  using Parent::_x;
  using Parent::_z;
  using Parent::_d;
  using Parent::_weight;
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
  using Parent::determineStronglyConnectedComponents;
  using Parent::separateStronglyConnectedComponent;
  
protected:
  TEMPLATE_DIGRAPH_TYPEDEFS(Graph);
  
public:
  EdgeCutLazy(IloEnv env,
              IloBoolVarArray x,
              IloBoolVarArray z,
              const Graph& d,
              const WeightNodeMap& weight,
              const IntNodeMap& nodeMap,
              const IntArcMap& arcMap,
              int n,
              int m,
              int maxNumberOfCuts,
              IloFastMutex* pMutex)
    : IloCplex::LazyConstraintCallbackI(env)
    , Parent(x, z, d, weight, nodeMap, arcMap, n, m, maxNumberOfCuts, pMutex)
  {
  }
  
  EdgeCutLazy(const EdgeCutLazy& other)
    : IloCplex::LazyConstraintCallbackI(other)
    , Parent(other)
  {
  }
  
  virtual ~EdgeCutLazy()
  {
  }
};

} // namespace mwcs
} // namespace nina


#endif // EDGECUTLAZY_H
