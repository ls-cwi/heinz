/*
 * edgecutuser.h
 *
 *  Created on: 11-sep-2014
 *      Author: M. El-Kebir
 */

#ifndef EDGECUTUSER_H
#define EDGECUTUSER_H

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
#include "backoff.h"
#include "bk_alg.h"

namespace nina {
namespace mwcs {

template<typename DGR,
         typename NWGHT = typename DGR::template NodeMap<double> >
class EdgeCutUser : public IloCplex::UserCutCallbackI,
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
  using Parent::_pMutex;
  using Parent::_epsilon;
  using Parent::lock;
  using Parent::unlock;
  
protected:
  TEMPLATE_DIGRAPH_TYPEDEFS(Graph);
  
  typedef typename Parent::NodeSet NodeSet;
  typedef typename Parent::NodeSetIt NodeSetIt;
  typedef typename Parent::NodeSetVector NodeSetVector;
  typedef typename Parent::NodeSetVectorIt NodeSetVectorIt;
  typedef typename Parent::SubGraph SubGraph;
  typedef typename Parent::SubNodeIt SubNodeIt;
  typedef typename Parent::SubArcIt SubArcIt;
  typedef typename Parent::NodeQueue NodeQueue;
  typedef typename Parent::NodeList NodeList;
  typedef typename Parent::NodeListIt NodeListIt;
  typedef DoubleArcMap CapacityMap;
  
  typedef nina::BkFlowAlg<Graph> BkAlg;
  
protected:
  CapacityMap* _pCap;
  BkAlg* _pBK;
  BoolNodeMap* _pMarked;
  NodeSet _rootSet;
  
  int _cutCount;
  int _nodeNumber;
  
protected:
  static const double _cutEpsilon = 0.00001 * _epsilon;
  const lemon::Tolerance<double> _cutTol;
  BackOff _backOff;
  bool _makeAttempt;
  
public:
  EdgeCutUser(IloEnv env,
              IloBoolVarArray x,
              IloBoolVarArray z,
              const Graph& d,
              const WeightNodeMap& weight,
              const IntNodeMap& nodeMap,
              const IntArcMap& arcMap,
              int n,
              int m,
              int maxNumberOfCuts,
              IloFastMutex* pMutex,
              const BackOff& backOff)
    : IloCplex::UserCutCallbackI(env)
    , Parent(x, z, d, weight, nodeMap, arcMap, n, m, maxNumberOfCuts, pMutex)
    , _pCap(NULL)
    , _pBK(NULL)
    , _pMarked(NULL)
    , _cutCount(0)
    , _nodeNumber(0)
    , _cutTol(_cutEpsilon)
    , _backOff(backOff)
    , _makeAttempt(true)
  {
    lock();
    _pCap = new CapacityMap(_d);
    _pMarked = new BoolNodeMap(_d, false);
    unlock();
  }
  
  EdgeCutUser(const EdgeCutUser& other)
    : IloCplex::UserCutCallbackI(other)
    , Parent(other)
    , _pBK(NULL)
    , _pMarked(NULL)
    , _rootSet(other._rootSet)
    , _cutCount(0)
    , _nodeNumber(0)
    , _cutTol(other._cutTol)
    , _backOff(other._backOff)
    , _makeAttempt(other._makeAttempt)
  {
    // TODO: to what values should I set cutCount and nodeNumber??
    lock();
    _pCap = new CapacityMap(_d);
    _pMarked = new BoolNodeMap(_d, false);
    unlock();
  }
  
  virtual ~EdgeCutUser()
  {
    delete _pBK;
    
    lock();
    delete _pCap;
    delete _pMarked;
    unlock();
  }
  
protected:
  virtual void main()
  {
//    if (!isAfterCutLoop())
//    {
//      return;
//    }
    
    if (_nodeNumber != getNnodes())
    {
      _nodeNumber = getNnodes();
      _cutCount = 0;
      _makeAttempt = _backOff.makeAttempt();
    }
    
    if (_makeAttempt && (_cutCount < _maxNumberOfCuts || _cutCount == -1 || (_nodeNumber == 0 && _cutCount < 50)))
    {
      separate();
      ++_cutCount;
    }
  }
  
  virtual void separate() = 0;
  
  void determineFwdCutSet(const Graph& h,
                          const BkAlg& bk,
                          const Node diRoot,
                          BoolNodeMap& marked,
                          NodeList& diS)
  {
    // we do a BFS on the *residual network* starting from _diRoot
    // and only following arcs that have nonzero residual capacity
    
    lemon::mapFill(h, marked, false);
    
    NodeQueue queue;
    queue.push(diRoot);
    marked[diRoot] = true;
    
    while (!queue.empty())
    {
      Node v = queue.front();
      queue.pop();
      diS.push_back(v);
      
      for (OutArcIt a(h, v); a != lemon::INVALID; ++a)
      {
        Node w = h.target(a);
        
        if (!marked[w] && _cutTol.nonZero(bk.resCap(a)))
        {
          queue.push(w);
          marked[w] = true;
        }
      }
      
      for (InArcIt a(h, v); a != lemon::INVALID; ++a)
      {
        Node w = h.source(a);
        
        if (!marked[w] && _cutTol.nonZero(bk.revResCap(a)))
        {
          queue.push(w);
          marked[w] = true;
        }
      }
    }
    
    assert(marked[diRoot] != marked[bk.getTarget()]);
  }
  
  void determineBwdCutSet(const Graph& h,
                          const BkAlg& bk,
                          const Node diRoot,
                          const Node target,
                          BoolNodeMap& marked,
                          NodeList& diS)
  {
    // we do a BFS on the reversed *residual network* starting from target
    // and only following arcs that have nonzero residual capacity
    
    lemon::mapFill(h, marked, false);
    
    NodeQueue queue;
    queue.push(target);
    marked[target] = true;
    
    while (!queue.empty())
    {
      Node v = queue.front();
      queue.pop();
      diS.push_back(v);
      
      for (InArcIt a(h, v); a != lemon::INVALID; ++a)
      {
        Node u = h.source(a);
        
        if (!marked[u] && _cutTol.nonZero(bk.resCap(a)))
        {
          queue.push(u);
          marked[u] = true;
        }
      }
      
      for (OutArcIt a(h, v); a != lemon::INVALID; ++a)
      {
        Node u = h.target(a);
        
        if (!marked[u] && _cutTol.nonZero(bk.revResCap(a)))
        {
          queue.push(u);
          marked[u] = true;
        }
      }
    }
    
    assert(marked[diRoot] != marked[target]);
  }
};
  
} // namespace mwcs
} // namespace nina


#endif // EDGECUTUSER_H
