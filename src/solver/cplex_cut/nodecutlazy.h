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
#include <lemon/adaptors.h>
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
class NodeCutLazy : public NodeCut<GR, NWGHT, NLBL, EWGHT>,
                    public IloCplex::LazyConstraintCallbackI
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
  using Parent::_pMutex;
  using Parent::_epsilon;
  
  using Parent::lock;
  using Parent::unlock;
  
protected:
  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  
  typedef lemon::FilterNodes<const Graph, const BoolNodeMap> SubGraph;
  typedef typename SubGraph::NodeIt SubNodeIt;
  
protected:
  const SubGraph* _pSubG;
  IntNodeMap* _pComp;
  
public:
  NodeCutLazy(IloEnv env,
              IloBoolVarArray x,
              IloBoolVarArray y,
              const Graph& g,
              const WeightNodeMap& weight,
              Node root,
              const IntNodeMap& nodeMap,
              int n,
              int m,
              int maxNumberOfCuts,
              IloFastMutex* pMutex)
    : Parent(x, y, g, weight, root, nodeMap, n, m, maxNumberOfCuts, pMutex)
    , IloCplex::LazyConstraintCallbackI(env)
    , _pSubG(NULL)
    , _pComp(NULL)
  {
    lock();
    _pSubG = new SubGraph(_g, *_pNodeBoolMap);
    _pComp = new IntNodeMap(_g);
    unlock();
  }
  
  NodeCutLazy(const NodeCutLazy& other)
    : Parent(other)
    , IloCplex::LazyConstraintCallbackI(other)
    , _pSubG(NULL)
    , _pComp(NULL)
  {
    lock();
    _pSubG = new SubGraph(_g, *_pNodeBoolMap);
    _pComp = new IntNodeMap(_g);
    unlock();
  }
  
  virtual ~NodeCutLazy()
  {
    lock();
    delete _pSubG;
    delete _pComp;
    unlock();
  }
  
protected:
  virtual void main()
  {
    separate();
  }
  
  virtual void separate() = 0;
};

} // namespace mwcs
} // namespace nina


#endif // NODECUTLAZY_H
