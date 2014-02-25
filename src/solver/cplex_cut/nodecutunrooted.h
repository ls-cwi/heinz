/*
 * nodecutunrooted.h
 *
 *  Created on: 18-feb-2014
 *      Author: M. El-Kebir
 */

#ifndef NODECUTUNROOTED_H
#define NODECUTUNROOTED_H

#include "nodecutlazy.h"
#include "nodecutuser.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class NodeCutUnrootedLazyConstraint : public NodeCutLazy<GR, NWGHT, NLBL, EWGHT>
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
  using Parent::_pSubG;
  using Parent::_pComp;
  
  using Parent::lock;
  using Parent::unlock;
  using Parent::addViolatedConstraint;
  using Parent::getEnv;
  using Parent::getValues;
  
  friend class NodeCut<GR, NWGHT, NLBL, EWGHT>;

public:
  NodeCutUnrootedLazyConstraint(IloEnv env,
                                IloBoolVarArray x,
                                IloBoolVarArray y,
                                const Graph& g,
                                const WeightNodeMap& weight,
                                const IntNodeMap& nodeMap,
                                int n,
                                int m,
                                int maxNumberOfCuts,
                                IloFastMutex* pMutex)
    : Parent(env, x, y, g, weight, lemon::INVALID, nodeMap, n, m, maxNumberOfCuts, pMutex)
  {
  }

  NodeCutUnrootedLazyConstraint(const NodeCutUnrootedLazyConstraint& other)
    : Parent(other)
  {
  }

//  virtual ~NodeCutUnrootedLazyConstraint()
//  {
//  }

protected:
  virtual void main()
  {
    separate();
  }

  virtual IloCplex::CallbackI* duplicateCallback() const
  {
    std::cout << "Duplicate!" << std::endl;
    return (new (getEnv()) NodeCutUnrootedLazyConstraint(*this));
  }

  void separate()
  {
    IloNumArray x_values(getEnv(), _n);
    getValues(x_values, _x);
    
    IloNumArray y_values(getEnv(), _n);
    getValues(y_values, _y);
    
    // determine non-zero y-vars
    NodeSet Y;
    for (NodeIt i(_g); i != lemon::INVALID; ++i)
    {
      int idx_i = _nodeMap[i];
      if (_tol.nonZero(y_values[idx_i]))
      {
        Y.insert(i);
      }
      _pNodeBoolMap->set(i, _tol.nonZero(x_values[idx_i]));
    }
    assert(Y.size() == 1);

    typedef std::vector<NodeSet> NodeSetVector;
    typedef typename NodeSetVector::const_iterator NodeSetVectorIt;
    
    const int nComp = lemon::connectedComponents(*_pSubG, *_pComp);
    int nCuts = 0;
    if (nComp != 1)
    {
      NodeSetVector compMatrix(nComp, NodeSet());
      for (SubNodeIt i(*_pSubG); i != lemon::INVALID; ++i)
      {
        int compIdx = (*_pComp)[i];
        compMatrix[compIdx].insert(i);
      }

      for (int compIdx = 0; compIdx < nComp; compIdx++)
      {
        const NodeSet& S = compMatrix[compIdx];

        NodeSet S_cap_Y;
        std::set_intersection(S.begin(), S.end(),
                              Y.begin(), Y.end(),
                              std::inserter(S_cap_Y, S_cap_Y.begin()));

        if (!S_cap_Y.empty())
        {
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
          addViolatedConstraint(*this, i, dS, S);
          ++nCuts;
        }
      }

      std::cerr << "[";
      for (int idx = 0; idx < nComp; idx++)
      {
        std::cerr << " " << compMatrix[idx].size();
      }
      std::cerr << " ]" << std::endl;
    }
    
    x_values.end();
    y_values.end();
    
    std::cerr << "Generated " << nCuts << " lazy cuts" << std::endl;
  }
};

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class NodeCutUnrootedUserCut : public NodeCutUser<GR, NWGHT, NLBL, EWGHT>
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
  using Parent::_h;
  using Parent::_cap;
  using Parent::_pG2h1;
  using Parent::_pG2h2;
  using Parent::_pG2hRootArc;
  using Parent::_h2g;
  using Parent::_diRoot;
  using Parent::_pBK;
  using Parent::_marked;

  using Parent::lock;
  using Parent::unlock;
  using Parent::determineFwdCutSet;
  using Parent::determineBwdCutSet;
  using Parent::addViolatedConstraint;
  using Parent::getEnv;
  using Parent::getValues;
  
  friend class NodeCut<GR, NWGHT, NLBL, EWGHT>;

public:
  NodeCutUnrootedUserCut(IloEnv env,
                         IloBoolVarArray x,
                         IloBoolVarArray y,
                         const Graph& g,
                         const WeightNodeMap& weight,
                         const IntNodeMap& nodeMap,
                         int n,
                         int m,
                         int maxNumberOfCuts,
                         IloFastMutex* pMutex)
    : Parent(env, x, y, g, weight, lemon::INVALID, nodeMap, n, m, maxNumberOfCuts, pMutex)
  {
    lock();
    _pG2hRootArc = new NodeDiArcMap(_g);
    unlock();
    
    init();
    _pBK = new BkAlg(_h, _cap);
  }

  NodeCutUnrootedUserCut(const NodeCutUnrootedUserCut& other)
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
    _pG2hRootArc = new NodeDiArcMap(_g);
    unlock();
    
    for (NodeIt v(_g); v != lemon::INVALID; ++v)
    {
      DiNode v1 = (*other._pG2h1)[v];
      DiNode v2 = (*other._pG2h2)[v];
      DiArc root_arc_v = (*other._pG2hRootArc)[v];
      
      _pG2h1->set(v, nodeMap[v1]);
      _pG2h2->set(v, nodeMap[v2]);
      _pG2hRootArc->set(v, arcMap[root_arc_v]);
    }
    
    _diRoot = nodeMap[other._diRoot];
    _pBK = new BkAlg(_h, _cap);
  }

  virtual ~NodeCutUnrootedUserCut()
  {
  }

protected:
  virtual IloCplex::CallbackI* duplicateCallback() const
  {
    return (new (getEnv()) NodeCutUnrootedUserCut(*this));
  }

  void separate()
  {
    IloNumArray x_values(getEnv(), _n);
    getValues(x_values, _x);
    
    IloNumArray y_values(getEnv(), _n);
    getValues(y_values, _y);
    
    computeCapacities(_cap, x_values, y_values);
    
    int nCuts = 0;
    int nBackCuts = 0;
    int nNestedCuts = 0;
    
    //printNonZeroVars(cbk, _y, y_values);
    
    for (NodeIt i(_g); i != lemon::INVALID; ++i)
    {
      // skip if node was already considered or its x-value is 0
      if (!(*_pNodeBoolMap)[i]) continue;

      const double x_i_value = x_values[_nodeMap[i]];
      
      _pBK->setSource(_diRoot);
      _pBK->setTarget((*_pG2h2)[i]);
      _pBK->setCap(_cap);
      
      bool nestedCut = false;
      bool first = true;
      //        while (true)
      {
        if (first)
        {
          _pBK->run();
          first = false;
        }
        else
        {
          _pBK->run(false);
        }
          
        // let's see if there's a violated constraint
        double minCutValue = _pBK->maxFlow();
        if (_tol.less(minCutValue, x_i_value))
        {
          // determine N (forward)
          NodeSet fwdS, fwdDS;
          determineFwdCutSet(_h, *_pBK, fwdDS, fwdS);
          
          // numerical instability may cause minCutValue < x_i_value
          // even though there is nothing to cut
          if (fwdS.empty() && fwdDS.empty()) break;
          
          // determine N (backward)
          NodeSet bwdS, bwdDS;
          determineBwdCutSet(_h, *_pBK, bwdDS, bwdS);
          
          bool backCuts = fwdDS.size() != bwdDS.size() ||
          fwdS.size() != bwdS.size() ||
          fwdDS != bwdDS || bwdS != fwdS;
          
          // add violated constraints for all nodes j in fwdS with x_j >= x_i
          for (NodeSetIt it2 = fwdS.begin(); it2 != fwdS.end(); ++it2)
          {
            const Node j = *it2;
            const double x_j_value = x_values[_nodeMap[j]];
            
            if (_tol.less(minCutValue, x_j_value))
            {
              _pNodeBoolMap->set(j, false);
              addViolatedConstraint(*this, j, fwdDS, fwdS);
              ++nCuts;
              if (nestedCut) ++nNestedCuts;
            }
          }
          
          if (backCuts)
          {
            // add violated constraints for all nodes j in bwdS with x_j >= x_i
            for (NodeSetIt it2 = bwdS.begin(); it2 != bwdS.end(); ++it2)
            {
              const Node j = *it2;
              const double x_j_value = x_values[_nodeMap[j]];
              
              if (_tol.less(minCutValue, x_j_value))
              {
                _pNodeBoolMap->set(j, false);
                addViolatedConstraint(*this, j, bwdDS, bwdS);
                ++nCuts;
                ++nBackCuts;
              }
            }
          }
          
          // generate nested-cuts
          //            for (NodeSetIt nodeIt = fwdDS.begin(); nodeIt != fwdDS.end(); nodeIt++)
          //            {
          //              nestedCut = true;
          //              // update the capactity to generate nested-cuts
          //              _pBK->incCap(DiOutArcIt(_h, (*_pG2h1)[*nodeIt]), 1);
          //            }
        }
        //          else
        //          {
        //            break;
        //          }
      }
    }
    
    x_values.end();
    y_values.end();
    
    std::cerr << "Generated " << nCuts
              << " user cuts of which " << nBackCuts << " are back-cuts and "
              << nNestedCuts << " are nested cuts" << std::endl;
    
    //    std::cerr << "[";
    //    for (int idx = 0; idx < nComp; idx++)
    //    {
    //      std::cerr << " " << compMatrix[idx].size();
    //    }
    //    std::cerr << " ]" << std::endl;
    
    //std::cerr << "Time: " << t.realTime() << "s" << std::endl;
  }
  
  
  void determineBwdCutSet(const Digraph& h,
                          const BkAlg& bk,
                          NodeSet& dS,
                          NodeSet& S)
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
        if (u == _diRoot)
        {
          assert(!_marked[u]);
          //std::cout << _h.id(u) << " -> "
          //          << _h.id(v) << " "
          //          << bk.flow(a) << "/" << bk.cap(a) << std::endl;
          S.insert(_h2g[v]);
        }
        else if (!_marked[u])
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
                          NodeSet& dS,
                          NodeSet& S)
  {
    DiNode target = bk.getTarget();
    DiNodeList diS;
    determineFwdCutSet(h, bk, diS);
    
    for (DiNodeListIt nodeIt = diS.begin(); nodeIt != diS.end(); nodeIt++)
    {
      DiNode v = *nodeIt;
      assert(_marked[v]);
      
      if (v == _diRoot)
      {
        for (DiOutArcIt a(h, v); a != lemon::INVALID; ++a)
        {
          DiNode w = h.target(a);
          if (!_marked[w])
          {
            //std::cout << _h.id(v) << " -> "
            //          << _h.id(w) << " "
            //          << bk.flow(a) << "/" << bk.cap(a) << std::endl;
            S.insert(_h2g[w]);
          }
        }
      }
      else
      {
        for (DiOutArcIt a(h, v); a != lemon::INVALID; ++a)
        {
          DiNode w = h.target(a);
          if (!_marked[w] && w != target)
          {
            //std::cout << _h.id(v) << " -> "
            //          << _h.id(w) << " "
            //          << bk.flow(a) << "/" << bk.cap(a) << std::endl;
            dS.insert(_h2g[w]);
          }
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
      DiNode i1 = _h.addNode();
      DiNode i2 = _h.addNode();
      _pG2h1->set(i, i1);
      _pG2h2->set(i, i2);
      _h2g[i1] = i;
      _h2g[i2] = i;
      
      DiArc i1i2 = _h.addArc(i1, i2);
      _cap[i1i2] = 0;
      
      DiArc ri1 = _h.addArc(_diRoot, i1);
      _pG2hRootArc->set(i, ri1);
      _cap[ri1] = 1;
    }
    
    for (EdgeIt e(_g); e != lemon::INVALID; ++e)
    {
      Node i = _g.u(e);
      Node j = _g.v(e);
      DiNode i1 = (*_pG2h1)[i];
      DiNode i2 = (*_pG2h2)[i];
      DiNode j1 = (*_pG2h1)[j];
      DiNode j2 = (*_pG2h2)[j];
      
      DiArc i2j1 = _h.addArc(i2, j1);
      DiArc j2i1 = _h.addArc(j2, i1);
      _cap[i2j1] = _cap[j2i1] = 1;
    }
  }
  
  NodeSet computeCapacities(CapacityMap& capacity,
                            IloNumArray x_values,
                            IloNumArray y_values)
  {
    NodeSet Y;
    
    for (NodeIt v(_g); v != lemon::INVALID; ++v)
    {
      // cap((i,j)) = x_i
      double val = x_values[_nodeMap[v]];
      if (!_tol.nonZero(val))
      {
        _pNodeBoolMap->set(v, false);
        val = 10 * _epsilon;
      }
      else
      {
        _pNodeBoolMap->set(v, true);
      }
      DiNode v1 = (*_pG2h1)[v];
      capacity[DiOutArcIt(_h, v1)] = val;
      
      // cap((r,i)) = y_i
      val = y_values[_nodeMap[v]];
      if (!_tol.nonZero(val))
      {
        val = 10 * _epsilon;
      }
      else
      {
        Y.insert(v);
      }
      capacity[(*_pG2hRootArc)[v]] = val;
    }
    
    return Y;
  }
};

} // namespace mwcs
} // namespace nina

#endif // NODECUTUNROOTED_H
