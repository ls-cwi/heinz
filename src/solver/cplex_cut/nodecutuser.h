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
  
  typedef typename Parent::NodeVector NodeVector;
  typedef typename Parent::NodeVectorIt NodeVectorIt;
  typedef typename Parent::NodeMatrix NodeMatrix;
  typedef typename Parent::NodeSet NodeSet;
  typedef typename Parent::NodeSetIt NodeSetIt;
  
  typedef std::queue<Node> NodeQueue;
  typedef std::queue<DiNode> DiNodeQueue;
  typedef std::set<DiNode> DiNodeSet;
  typedef std::list<DiNode> DiNodeList;
  typedef DiNodeList::const_iterator DiNodeListIt;
  
  typedef lemon::FilterNodes<const Graph, const BoolNodeMap> SubGraph;
  typedef typename SubGraph::NodeIt SubNodeIt;
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
  
public:
  NodeCutUser(IloEnv env,
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
    : IloCplex::UserCutCallbackI(env)
    , Parent(x, y, g, weight, root, nodeMap, n, m, maxNumberOfCuts, pMutex)
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
    if (!isAfterCutLoop())
    {
      return;
    }
    
    if (_nodeNumber != getNnodes())
    {
      _nodeNumber = getNnodes();
      _cutCount = 0;
    }
    
    if (_nodeNumber == 0 || _cutCount < _maxNumberOfCuts)
    {
      separate();
      ++_cutCount;
    }
  }
  
  virtual void separate() = 0;
  
  void determineFwdCutSet(const Digraph& h,
                          const BkAlg& bk,
                          DiNodeList& diS)
  {
    // we do a DFS on the *residual network* starting from _diRoot
    // and only following arcs that have nonzero residual capacity
    
    lemon::mapFill(_h, _marked, false);
    
    DiNodeQueue queue;
    queue.push(_diRoot);
    _marked[_diRoot] = true;
    
    while (!queue.empty())
    {
      DiNode v = queue.front();
      queue.pop();
      diS.push_back(v);
      
      for (DiOutArcIt a(h, v); a != lemon::INVALID; ++a)
      {
        DiNode w = _h.target(a);
        
        if (!_marked[w] && _cutTol.nonZero(bk.resCap(a)))
        {
          queue.push(w);
          _marked[w] = true;
        }
      }
      
      for (DiInArcIt a(h, v); a != lemon::INVALID; ++a)
      {
        DiNode w = _h.source(a);
        
        if (!_marked[w] && _cutTol.nonZero(bk.revResCap(a)))
        {
          queue.push(w);
          _marked[w] = true;
        }
      }
    }
    
    assert(_marked[_diRoot] != _marked[bk.getTarget()]);
  }
  
  void determineBwdCutSet(const Digraph& h,
                          const BkAlg& bk,
                          const DiNode target,
                          DiNodeList& diS)
  {
    // we do a DFS on the reversed *residual network* starting from target
    // and only following arcs that have nonzero residual capacity
    
    lemon::mapFill(_h, _marked, false);
    
    DiNodeQueue queue;
    queue.push(target);
    _marked[target] = true;
    
    while (!queue.empty())
    {
      DiNode v = queue.front();
      queue.pop();
      diS.push_back(v);
      
      for (DiInArcIt a(h, v); a != lemon::INVALID; ++a)
      {
        DiNode u = _h.source(a);
        
        if (!_marked[u] && _cutTol.nonZero(bk.resCap(a)))
        {
          queue.push(u);
          _marked[u] = true;
        }
      }
      
      for (DiOutArcIt a(h, v); a != lemon::INVALID; ++a)
      {
        DiNode u = _h.target(a);
        
        if (!_marked[u] && _cutTol.nonZero(bk.revResCap(a)))
        {
          queue.push(u);
          _marked[u] = true;
        }
      }
    }
    
    assert(_marked[_diRoot] != _marked[target]);
  }
};
  
} // namespace mwcs
} // namespace nina


#endif // NODECUTUSER_H
