/*
 * nodecutuser.h
 *
 *  Created on: 24-feb-2014
 *      Author: M. El-Kebir
 */

#ifndef NODECUTUSER_H
#define NODECUTUSER_H

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

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class NodeCutUser : public IloCplex::UserCutCallbackI,
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
  using Parent::_nodeMap;
  using Parent::_n;
  using Parent::_maxNumberOfCuts;
  using Parent::_tol;
  using Parent::_pNodeBoolMap;
  using Parent::_pMutex;
  using Parent::_epsilon;
  using Parent::lock;
  using Parent::unlock;
  
protected:
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
  
  typedef typename Parent::NodeSet NodeSet;
  typedef typename Parent::NodeSetIt NodeSetIt;
  typedef typename Parent::NodeSetVector NodeSetVector;
  typedef typename Parent::NodeSetVectorIt NodeSetVectorIt;
  typedef typename Parent::SubGraph SubGraph;
  typedef typename Parent::SubNodeIt SubNodeIt;
  typedef typename Parent::SubEdgeIt SubEdgeIt;
  
  typedef std::queue<Node> NodeQueue;
  typedef std::queue<DiNode> DiNodeQueue;
  typedef std::set<DiNode> DiNodeSet;
  typedef std::list<DiNode> DiNodeList;
  typedef DiNodeList::const_iterator DiNodeListIt;
  
  typedef nina::BkFlowAlg<Digraph> BkAlg;
  typedef typename Digraph::NodeMap<bool> DiBoolNodeMap;
  
protected:
  Digraph _h;
  CapacityMap _cap;
  NodeDiNodeMap* _pG2h1;
  NodeDiNodeMap* _pG2h2;
  NodeDiArcMap* _pG2hRootArc;
  DiNodeNodeMap _h2g;
  DiNode _diRoot;
  BkAlg* _pBK;
  DiBoolNodeMap _marked;
  
  int _cutCount;
  int _nodeNumber;
  
protected:
  static const double _cutEpsilon = 0.00001 * _epsilon;
  const lemon::Tolerance<double> _cutTol;
  BackOff _backOff;
  bool _makeAttempt;
  
public:
  NodeCutUser(IloEnv env,
              IloBoolVarArray x,
              IloBoolVarArray y,
              const Graph& g,
              const WeightNodeMap& weight,
              const IntNodeMap& nodeMap,
              int n,
              int maxNumberOfCuts,
              IloFastMutex* pMutex,
              const BackOff& backOff)
    : IloCplex::UserCutCallbackI(env)
    , Parent(x, y, g, weight, nodeMap, n, maxNumberOfCuts, pMutex)
    , _h()
    , _cap(_h)
    , _pG2h1(NULL)
    , _pG2h2(NULL)
    , _pG2hRootArc(NULL)
    , _h2g(_h)
    , _diRoot(lemon::INVALID)
    , _pBK(NULL)
    , _marked(_h, false)
    , _cutCount(0)
    , _nodeNumber(0)
    , _cutTol(_cutEpsilon)
    , _backOff(backOff)
    , _makeAttempt(true)
  {
    lock();
    _pG2h1 = new NodeDiNodeMap(_g);
    _pG2h2 = new NodeDiNodeMap(_g);
    unlock();
  }
  
  NodeCutUser(const NodeCutUser& other)
    : IloCplex::UserCutCallbackI(other)
    , Parent(other)
    , _h()
    , _cap(_h)
    , _pG2h1(NULL)
    , _pG2h2(NULL)
    , _pG2hRootArc(NULL)
    , _h2g(_h)
    , _diRoot(lemon::INVALID)
    , _pBK(NULL)
    , _marked(_h, false)
    , _cutCount(0)
    , _nodeNumber(0)
    , _cutTol(other._cutTol)
    , _backOff(other._backOff)
    , _makeAttempt(other._makeAttempt)
  {
    // TODO: to what values should I set cutCount and nodeNumber??
  }
  
  virtual ~NodeCutUser()
  {
    delete _pBK;
    
    lock();
    delete _pG2h1;
    delete _pG2h2;
    delete _pG2hRootArc;
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
    
    if (_makeAttempt && (_cutCount < _maxNumberOfCuts || _cutCount == -1 || (_nodeNumber == 0 && _cutCount < 250)))
    {
      separate();
      ++_cutCount;
    }
  }
  
  virtual void separate() = 0;
  
  void determineFwdCutSet(const Digraph& h,
                          const BkAlg& bk,
                          const DiNode diRoot,
                          DiBoolNodeMap& marked,
                          DiNodeList& diS)
  {
    // we do a DFS on the *residual network* starting from _diRoot
    // and only following arcs that have nonzero residual capacity
    
    lemon::mapFill(h, marked, false);
    
    DiNodeQueue queue;
    queue.push(diRoot);
    marked[diRoot] = true;
    
    while (!queue.empty())
    {
      DiNode v = queue.front();
      queue.pop();
      diS.push_back(v);
      
      for (DiOutArcIt a(h, v); a != lemon::INVALID; ++a)
      {
        DiNode w = h.target(a);
        
        if (!marked[w] && _cutTol.nonZero(bk.resCap(a)))
        {
          queue.push(w);
          marked[w] = true;
        }
      }
      
      for (DiInArcIt a(h, v); a != lemon::INVALID; ++a)
      {
        DiNode w = h.source(a);
        
        if (!marked[w] && _cutTol.nonZero(bk.revResCap(a)))
        {
          queue.push(w);
          marked[w] = true;
        }
      }
    }
    
    assert(_marked[_diRoot] != _marked[bk.getTarget()]);
  }
  
  void determineBwdCutSet(const Digraph& h,
                          const BkAlg& bk,
                          const DiNode diRoot,
                          const DiNode target,
                          DiBoolNodeMap& marked,
                          DiNodeList& diS)
  {
    // we do a DFS on the reversed *residual network* starting from target
    // and only following arcs that have nonzero residual capacity
    
    lemon::mapFill(h, marked, false);
    
    DiNodeQueue queue;
    queue.push(target);
    marked[target] = true;
    
    while (!queue.empty())
    {
      DiNode v = queue.front();
      queue.pop();
      diS.push_back(v);
      
      for (DiInArcIt a(h, v); a != lemon::INVALID; ++a)
      {
        DiNode u = h.source(a);
        
        if (!marked[u] && _cutTol.nonZero(bk.resCap(a)))
        {
          queue.push(u);
          marked[u] = true;
        }
      }
      
      for (DiOutArcIt a(h, v); a != lemon::INVALID; ++a)
      {
        DiNode u = h.target(a);
        
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


#endif // NODECUTUSER_H
