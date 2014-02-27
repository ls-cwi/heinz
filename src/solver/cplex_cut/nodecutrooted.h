/*
 * nodecutrooted.h
 *
 *  Created on: 18-feb-2014
 *      Author: M. El-Kebir
 */

#ifndef NODECUTROOTED_H
#define NODECUTROOTED_H

#include "nodecutlazy.h"
#include "nodecutuser.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class NodeCutRootedLazyConstraint : public NodeCutLazy<GR, NWGHT, NLBL, EWGHT>
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;
  typedef NodeCutLazy<GR, NWGHT, NLBL, EWGHT> Parent;

protected:
  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  typedef typename Parent::NodeVector NodeVector;
  typedef typename Parent::NodeVectorIt NodeVectorIt;
  typedef typename Parent::NodeMatrix NodeMatrix;
  typedef typename Parent::NodeSet NodeSet;
  typedef typename Parent::NodeSetIt NodeSetIt;
  typedef typename Parent::SubGraph SubGraph;
  typedef typename Parent::SubNodeIt SubNodeIt;
  
  using Parent::_x;
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
  using Parent::_y;
  using Parent::_pSubG;
  using Parent::_pComp;
  
  using Parent::lock;
  using Parent::unlock;
  using Parent::addViolatedConstraint;
  using Parent::getEnv;
  using Parent::getValues;
  
  friend class NodeCut<GR, NWGHT, NLBL, EWGHT>;

public:
  NodeCutRootedLazyConstraint(IloEnv env,
                              IloBoolVarArray x,
                              const Graph& g,
                              const WeightNodeMap& weight,
                              Node root,
                              const IntNodeMap& nodeMap,
                              int n,
                              int m,
                              int maxNumberOfCuts,
                              IloFastMutex* pMutex)
    : Parent(env, x, IloBoolVarArray(), g, weight, root, nodeMap, n, m, maxNumberOfCuts, pMutex)
  {
  }

  NodeCutRootedLazyConstraint(const NodeCutRootedLazyConstraint& other)
    : Parent(other)
  {
  }

  virtual ~NodeCutRootedLazyConstraint()
  {
  }

protected:
  virtual void main()
  {
    separate();
  }

  virtual IloCplex::CallbackI* duplicateCallback() const
  {
    return (new (getEnv()) NodeCutRootedLazyConstraint(*this));
  }
  
  void separate()
  {
    IloNumArray x_values(getEnv(), _n);
    getValues(x_values, _x);
    
    typedef std::vector<NodeSet> NodeSetVector;
    typedef typename NodeSetVector::const_iterator NodeSetVectorIt;
    
    for (NodeIt i(_g); i != lemon::INVALID; ++i)
    {
      int idx_i = _nodeMap[i];
      _pNodeBoolMap->set(i, _tol.nonZero(x_values[idx_i]));
    }
    
    const int nComp = lemon::connectedComponents(*_pSubG, *_pComp);
    if (nComp > 1)
    {
      NodeSetVector compMatrix(nComp, NodeSet());
      for (SubNodeIt i(*_pSubG); i != lemon::INVALID; ++i)
      {
        int compIdx = (*_pComp)[i];
        compMatrix[compIdx].insert(i);
      }
      
      int nCuts = 0;
      for (int compIdx = 0; compIdx < nComp; compIdx++)
      {
        const NodeSet& S = compMatrix[compIdx];
        
        if (S.find(_root) != S.end())
        {
          // nothing to separate as root is in S
          continue;
        }
        
        // determine dS
        NodeSet dS;
        for (NodeSetIt it = S.begin(); it != S.end(); ++it)
        {
          const Node i = *it;
          for (OutArcIt a(_g, i); a != lemon::INVALID; ++a)
          {
            const Node j = _g.target(a);
            if (S.find(j) == S.end())
            {
              dS.insert(j);
            }
          }
        }
        
        for (NodeSetIt it = S.begin(); it != S.end(); ++it)
        {
          const Node i = *it;
          addViolatedConstraint(*this, i, dS);
          ++nCuts;
        }
      }
      
//      std::cerr << "Generated " << nCuts << " lazy cuts" << std::endl;
//      std::cerr << "[";
//      for (int idx = 0; idx < nComp; idx++)
//      {
//        std::cerr << " " << compMatrix[idx].size();
//      }
//      std::cerr << " ]" << std::endl;
    }
    
    x_values.end();
  }
};

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class NodeCutRootedUserCut : public NodeCutUser<GR, NWGHT, NLBL, EWGHT>
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;
  typedef NodeCutUser<GR, NWGHT, NLBL, EWGHT> Parent;

protected:
  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  typedef typename Parent::Digraph Digraph;
  typedef typename Parent::DiArc DiArc;
  typedef typename Parent::DiNode DiNode;
  typedef typename Parent::DiNodeIt DiNodeIt;
  typedef typename Parent::DiArcIt DiArcIt;
  typedef typename Parent::DiInArcIt DiInArcIt;
  typedef typename Parent::DiOutArcIt DiOutArcIt;
  typedef typename Parent::DiNodeNodeMap DiNodeNodeMap;
  typedef typename Parent::NodeDiNodeMap NodeDiNodeMap;
  typedef typename Parent::NodeDiArcMap NodeDiArcMap;
  typedef typename Parent::CapacityMap CapacityMap;
  typedef typename Parent::NodeVector NodeVector;
  typedef typename Parent::NodeVectorIt NodeVectorIt;
  typedef typename Parent::NodeMatrix NodeMatrix;
  typedef typename Parent::NodeSet NodeSet;
  typedef typename Parent::NodeSetIt NodeSetIt;
  typedef typename Parent::NodeQueue NodeQueue;
  typedef typename Parent::DiNodeQueue DiNodeQueue;
  typedef typename Parent::DiNodeSet DiNodeSet;
  typedef typename Parent::DiNodeList DiNodeList;
  typedef typename Parent::DiNodeListIt DiNodeListIt;
  typedef typename Parent::SubGraph SubGraph;
  typedef typename Parent::SubNodeIt SubNodeIt;
  typedef typename Parent::BkAlg BkAlg;
  typedef typename Parent::DiBoolNodeMap DiBoolNodeMap;
  
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
  using Parent::_cutEpsilon;
  using Parent::_h;
  using Parent::_cap;
  using Parent::_pG2h1;
  using Parent::_pG2h2;
  using Parent::_pG2hRootArc;
  using Parent::_h2g;
  using Parent::_diRoot;
  using Parent::_pBK;
  using Parent::_marked;
  using Parent::_cutCount;
  
  using Parent::lock;
  using Parent::unlock;
  using Parent::determineFwdCutSet;
  using Parent::determineBwdCutSet;
  using Parent::addViolatedConstraint;
  using Parent::getEnv;
  using Parent::getValues;
  using Parent::printNodeSet;
  
  friend class NodeCut<GR, NWGHT, NLBL, EWGHT>;

public:
  NodeCutRootedUserCut(IloEnv env,
                       IloBoolVarArray x,
                       const Graph& g,
                       const WeightNodeMap& weight,
                       Node root,
                       const IntNodeMap& nodeMap,
                       int n,
                       int m,
                       int maxNumberOfCuts,
                       IloFastMutex* pMutex)
    : Parent(env, x, IloBoolVarArray(), g, weight, root, nodeMap, n, m, maxNumberOfCuts, pMutex)
  {
    init();
    _pBK = new BkAlg(_h, _cap);
  }

  NodeCutRootedUserCut(const NodeCutRootedUserCut& other)
    : Parent(other)
  {
    typename Digraph::template NodeMap<DiNode> nodeMap(other._h);
    typename Digraph::template ArcMap<DiArc> arcMap(other._h);
    
    lemon::digraphCopy(other._h, _h)
    .nodeRef(nodeMap)
    .arcRef(arcMap)
    .nodeMap(other._h2g, _h2g)
    .arcMap(other._cap, _cap)
    .run();
    
    lock();
    _pG2h1 = new NodeDiNodeMap(_g);
    _pG2h2 = new NodeDiNodeMap(_g);
    unlock();
    
    for (NodeIt v(_g); v != lemon::INVALID; ++v)
    {
      DiNode v1 = (*other._pG2h1)[v];
      DiNode v2 = (*other._pG2h2)[v];
      
      _pG2h1->set(v, nodeMap[v1]);
      if (v2 != lemon::INVALID)
        _pG2h2->set(v, nodeMap[v2]);
    }
    
    _diRoot = nodeMap[other._diRoot];
    _pBK = new BkAlg(_h, _cap);
  }

  virtual ~NodeCutRootedUserCut()
  {
  }

protected:
  virtual IloCplex::CallbackI* duplicateCallback() const
  {
    return (new (getEnv()) NodeCutRootedUserCut(*this));
  }

  void separate()
  {
    IloNumArray x_values(getEnv(), _n);
    getValues(x_values, _x);
    
    computeCapacities(_cap, x_values);
    
    int nCuts = 0;
    int nBackCuts = 0;
    int nNestedCuts = 0;

    _pBK->setSource(_diRoot);
    _pNodeBoolMap->set(_root, false);
    for (NodeIt i(_g); i != lemon::INVALID; ++i)
    {
      // skip if node was already considered or its x-value is 0
      if (!(*_pNodeBoolMap)[i]) continue;
      
      const double x_i_value = x_values[_nodeMap[i]];
      
      _pBK->setTarget((*_pG2h2)[i]);
      _pBK->setCap(_cap);

      bool nestedCut = false;
      bool first = true;
      while (true)
      {
        if (first)
        {
          _pBK->run();
          first = false;
        }
        else
        {
          _pBK->run(true);
        }
        
        // let's see if there's a violated constraint
        double minCutValue = _pBK->maxFlow();
        if (_tol.less(minCutValue, x_i_value))
        {
          // determine N (forward)
          NodeSet fwdDS;
          determineFwdCutSet(_h, *_pBK, fwdDS);
          
          double m = 0;
          for (NodeSetIt it = fwdDS.begin(); it != fwdDS.end(); ++it)
          {
            m += _cap[DiOutArcIt(_h, (*_pG2h1)[*it])];
          }
          if (_tol.different(m, minCutValue))
          {
            printNodeSet(fwdDS, _x, x_values);
            assert(false);
          }
          
          // numerical instability may cause minCutValue < x_i_value
          // even though there is nothing to cut ( || fwdDS.find(i) != fwdDS.end())
          if (fwdDS.empty()) break;
          
          // determine N (forward)
          NodeSet bwdDS;
          determineBwdCutSet(_h, *_pBK, bwdDS);
          
          bool backCuts = fwdDS.size() != bwdDS.size() || fwdDS != bwdDS;
          
          // add violated constraints
          _pNodeBoolMap->set(i, false);
          addViolatedConstraint(*this, i, fwdDS);
          ++nCuts;
          if (nestedCut) ++nNestedCuts;
          
          if (backCuts)
          {
            addViolatedConstraint(*this, i, bwdDS);
            ++nCuts;
            ++nBackCuts;
          }
          
          // generate nested-cuts
          for (NodeSetIt nodeIt = fwdDS.begin(); nodeIt != fwdDS.end(); nodeIt++)
          {
            nestedCut = true;
            // update the capactity to generate nested-cuts
            _pBK->incCap(DiOutArcIt(_h, (*_pG2h1)[*nodeIt]), 1);
          }
        }
        else
        {
          break;
        }
      }

//      std::cerr <<  "#" << _cutCount << ": generated " << nCuts
//                << " user cuts of which " << nBackCuts << " are back-cuts and "
//                << nNestedCuts << " are nested cuts" << std::endl;
    }
    
    x_values.end();
  }
  
  void computeCapacities(CapacityMap& capacity,
                         IloNumArray x_values)
  {
    // cap((i,j)) = x_i
    for (NodeIt v(_g); v != lemon::INVALID; ++v)
    {
      double val = x_values[_nodeMap[v]];
      if (!_tol.nonZero(val))
      {
        _pNodeBoolMap->set(v, false);
        val = 10* _cutEpsilon;
      }
      else
      {
        _pNodeBoolMap->set(v, true);
      }
      
      DiNode v1 = (*_pG2h1)[v];
      DiOutArcIt a(_h, v1);
      
      if (a != lemon::INVALID)
      {
        capacity[a] = val;
      }
      else
      {
        assert(v == _root);
      }
    }
  }
  
  void determineBwdCutSet(const Digraph& h,
                          const BkAlg& bk,
                          NodeSet& dS)
  {
    DiNode target = bk.getTarget();
    DiNodeList diS;
    determineBwdCutSet(h, bk, target, diS);
    
    for (DiNodeListIt nodeIt = diS.begin(); nodeIt != diS.end(); nodeIt++)
    {
      DiNode v = *nodeIt;
      assert(_marked[v]);
      if (v == target) continue;
      
      for (DiInArcIt a(h, v); a != lemon::INVALID; ++a)
      {
        DiNode u = h.source(a);
        if (!_marked[u])
        {
          //std::cout << _h.id(u) << " -> "
          //          << _h.id(v) << " "
          //          << bk.flow(a) << "/" << bk.cap(a) << std::endl;
          dS.insert(_h2g[v]);
        }
      }
    }
  }
  
  void determineFwdCutSet(const Digraph& h,
                          const BkAlg& bk,
                          NodeSet& dS)
  {
    DiNode target = bk.getTarget();
    DiNodeList diS;
    determineFwdCutSet(h, bk, diS);
    
    for (DiNodeListIt nodeIt = diS.begin(); nodeIt != diS.end(); nodeIt++)
    {
      DiNode v = *nodeIt;
      assert(_marked[v]);
      
      for (DiOutArcIt a(h, v); a != lemon::INVALID; ++a)
      {
        DiNode w = h.target(a);
        if (!_marked[w] && w != target)
        {
          assert(v != _diRoot);
          //std::cout << _h.id(v) << " -> "
          //          << _h.id(w) << " "
          //          << bk.flow(a) << "/" << bk.cap(a) << " : " << bk.resCap(a) << std::endl;
          dS.insert(_h2g[w]);
        }
      }
    }
  }
  
  void init()
  {
    // we initialize _h:
    // - for every node i, there will be two nodes i1 and i2
    //   connected by an arc from i1 to i2
    // - for every edge (i,j) there are two arcs
    //   in h: (i2,j1) and (j2,i1) with capacties 1
    _h.clear();
    _diRoot = _h.addNode();
    
    for (NodeIt i(_g); i != lemon::INVALID; ++i)
    {
      if (i != _root)
      {
        DiNode i1 = _h.addNode();
        DiNode i2 = _h.addNode();
        _pG2h1->set(i, i1);
        _pG2h2->set(i, i2);
        _h2g[i1] = i;
        _h2g[i2] = i;
        
        DiArc i1i2 = _h.addArc(i1, i2);
        _cap[i1i2] = 0;
      }
      else
      {
        _pG2h1->set(i, _diRoot);
        _pG2h2->set(i, _diRoot);
        _h2g[_diRoot] = _root;
      }
    }
    
    for (EdgeIt e(_g); e != lemon::INVALID; ++e)
    {
      Node i = _g.u(e);
      Node j = _g.v(e);
      DiNode i1 = (*_pG2h1)[i];
      DiNode i2 = (*_pG2h2)[i];
      DiNode j1 = (*_pG2h1)[j];
      DiNode j2 = (*_pG2h2)[j];
      
      if (i != _root && j != _root)
      {
        DiArc i2j1 = _h.addArc(i2, j1);
        DiArc j2i1 = _h.addArc(j2, i1);
        _cap[i2j1] = _cap[j2i1] = 1;
      }
      else if (i == _root)
      {
        DiArc ij1 = _h.addArc(_diRoot, j1);
        _cap[ij1] = 1;
      }
      else if (j == _root)
      {
        DiArc ji1 = _h.addArc(_diRoot, i1);
        _cap[ji1] = 1;
      }
    }
  }
};

} // namespace mwcs
} // namespace nina

#endif // NODECUTROOTED_H
