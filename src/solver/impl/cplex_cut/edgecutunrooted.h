/*
 * edgecutunrooted.h
 *
 *  Created on: 11-sep-2014
 *      Author: M. El-Kebir
 */

#ifndef EDGECUTUNROOTED_H
#define EDGECUTUNROOTED_H

#include "edgecutlazy.h"
#include "edgecutuser.h"
#include <algorithm>

namespace nina {
namespace mwcs {

template<typename DGR,
         typename NWGHT = typename DGR::template NodeMap<double> >
class EdgeCutUnrootedLazyConstraint : public EdgeCutLazy<DGR, NWGHT>
{
public:
  typedef DGR Graph;
  typedef NWGHT WeightNodeMap;
  typedef EdgeCutLazy<DGR, NWGHT> Parent;

protected:
  TEMPLATE_DIGRAPH_TYPEDEFS(Graph);
  
  typedef typename Parent::NodeSet NodeSet;
  typedef typename Parent::NodeSetIt NodeSetIt;
  typedef typename Parent::NodeSetVector NodeSetVector;
  typedef typename Parent::NodeSetVectorIt NodeSetVectorIt;
  typedef typename Parent::SubGraph SubGraph;
  typedef typename Parent::SubNodeIt SubNodeIt;

  using Parent::_x;
  using Parent::_z;
  using Parent::_d;
  using Parent::_weight;
  using Parent::_nodeMap;
  using Parent::_arcMap;
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
  using Parent::getEnv;
  using Parent::getValues;
  using Parent::constructRHS;
  using Parent::add;
  using Parent::isValid;
  using Parent::determineStronglyConnectedComponents;
  using Parent::separateStronglyConnectedComponent;
  using Parent::printNonZeroVars;
  
  friend class EdgeCut<DGR, NWGHT>;

public:
  EdgeCutUnrootedLazyConstraint(IloEnv env,
                                IloBoolVarArray x,
                                IloBoolVarArray z,
                                const Graph& d,
                                const Node root,
                                const WeightNodeMap& weight,
                                const IntNodeMap& nodeMap,
                                const IntArcMap& arcMap,
                                int n,
                                int m,
                                int maxNumberOfCuts,
                                IloFastMutex* pMutex)
    : Parent(env, x, z, d, weight, nodeMap, arcMap, n, m, maxNumberOfCuts, pMutex)
    , _root(root)
  {
  }

  EdgeCutUnrootedLazyConstraint(const EdgeCutUnrootedLazyConstraint& other)
    : Parent(other)
  {
  }

  virtual ~EdgeCutUnrootedLazyConstraint()
  {
  }

protected:
  Node _root;
  
  virtual void main()
  {
    separate();
  }

  virtual IloCplex::CallbackI* duplicateCallback() const
  {
    return (new (getEnv()) EdgeCutUnrootedLazyConstraint(*this));
  }

  void separate()
  {
    IloNumArray x_values(getEnv(), _n);
    getValues(x_values, _x);

    IloNumArray z_values(getEnv(), _m);
    getValues(z_values, _z);
    
//    printNonZeroVars(*this, _z, z_values);
//    printNonZeroVars(*this, _x, x_values);
    
    // determine non-zero y-vars
    NodeSet rootNodes;
    for (OutArcIt a(_d, _root); a != lemon::INVALID; ++a)
    {
      int idx_a = _arcMap[a];
      if (_tol.nonZero(z_values[idx_a]))
      {
        rootNodes.insert(_d.target(a));
      }
    }
    
    assert(rootNodes.size() == 1);

    // determine connected components
    NodeSetVector nonZeroComponents = determineStronglyConnectedComponents(x_values);

    int nCuts = 0;
    for (NodeSetVectorIt it = nonZeroComponents.begin(); it != nonZeroComponents.end(); ++it)
    {
      const NodeSet& nonZeroComponent = *it;
      if (nonZeroComponent.find(_root) != nonZeroComponent.end())
      {
        assert(nonZeroComponent.size() == 1);
        continue;
      }
      
      NodeSet intersection;
      std::set_intersection(rootNodes.begin(), rootNodes.end(),
                            nonZeroComponent.begin(), nonZeroComponent.end(),
                            std::inserter(intersection, intersection.begin()));
      
      if (intersection.size() == 0)
      {
        separateStronglyConnectedComponent(*it, rootNodes, x_values, z_values, *this, nCuts);
      }
    }

    x_values.end();
    z_values.end();
    
    std::cerr << "#comps: " << nonZeroComponents.size() << ", generated " << nCuts << " lazy cuts" << std::endl;
  }
};

template<typename DGR,
         typename NWGHT = typename DGR::template NodeMap<double> >
class EdgeCutUnrootedUserCut : public EdgeCutUser<DGR, NWGHT>
{
public:
  typedef DGR Graph;
  typedef NWGHT WeightNodeMap;
  typedef EdgeCutUser<DGR, NWGHT> Parent;

protected:
  TEMPLATE_DIGRAPH_TYPEDEFS(Graph);
  typedef typename Parent::CapacityMap CapacityMap;
  typedef typename Parent::NodeSet NodeSet;
  typedef typename Parent::NodeSetIt NodeSetIt;
  typedef typename Parent::ArcSet ArcSet;
  typedef typename Parent::ArcSetIt ArcSetIt;
  typedef typename Parent::NodeSetVector NodeSetVector;
  typedef typename Parent::NodeSetVectorIt NodeSetVectorIt;
  typedef typename Parent::SubGraph SubGraph;
  typedef typename Parent::SubNodeIt SubNodeIt;
  typedef typename Parent::SubArcIt SubArcIt;
  typedef typename Parent::NodeQueue NodeQueue;
  typedef typename Parent::BkAlg BkAlg;
  typedef typename Parent::NodeList NodeList;
  typedef typename Parent::NodeListIt NodeListIt;

  using Parent::_x;
  using Parent::_z;
  using Parent::_d;
  using Parent::_weight;
  using Parent::_nodeMap;
  using Parent::_arcMap;
  using Parent::_n;
  using Parent::_m;
  using Parent::_maxNumberOfCuts;
  using Parent::_tol;
  using Parent::_pNodeBoolMap;
  using Parent::_pMutex;
  using Parent::_epsilon;
  using Parent::_cutEpsilon;
  using Parent::_rootSet;
  using Parent::_pBK;
  using Parent::_pCap;
  using Parent::_pMarked;
  using Parent::_cutCount;
  using Parent::_nodeNumber;
  using Parent::_pSubG;
  using Parent::_pComp;

  using Parent::lock;
  using Parent::unlock;
  using Parent::determineFwdCutSet;
  using Parent::determineBwdCutSet;
  using Parent::add;
  using Parent::getEnv;
  using Parent::getValues;
  using Parent::getNnodes;
  using Parent::constructRHS;
  using Parent::isValid;
  using Parent::determineStronglyConnectedComponents;
  using Parent::separateStronglyConnectedComponent;
  
  friend class EdgeCut<DGR, NWGHT>;

public:
  EdgeCutUnrootedUserCut(IloEnv env,
                         IloBoolVarArray x,
                         IloBoolVarArray z,
                         const Graph& d,
                         const Node root,
                         const WeightNodeMap& weight,
                         const IntNodeMap& nodeMap,
                         const IntArcMap& arcMap,
                         int n,
                         int m,
                         int maxNumberOfCuts,
                         IloFastMutex* pMutex,
                         BackOff backOff)
    : Parent(env, x, z, d, weight, nodeMap, arcMap, n, m, maxNumberOfCuts, pMutex, backOff)
    , _root(root)
  {
    _pBK = new BkAlg(_d, *_pCap);
    _rootSet.insert(root);
  }

  EdgeCutUnrootedUserCut(const EdgeCutUnrootedUserCut& other)
    : Parent(other)
  {
    _pBK = new BkAlg(_d, *_pCap);
  }

  virtual ~EdgeCutUnrootedUserCut()
  {
  }

protected:
  Node _root;
  
  virtual IloCplex::CallbackI* duplicateCallback() const
  {
    return (new (getEnv()) EdgeCutUnrootedUserCut(*this));
  }
  
  void separateMinCut(const NodeSet& nonZeroComponent,
                      const IloNumArray& x_values,
                      const IloNumArray& z_values,
                      int& nCuts, int& nBackCuts, int& nNestedCuts)
  {
    IloExpr rhs(getEnv());

    _pBK->setSource(_root);
    _pNodeBoolMap->set(_root, false);
    for (NodeSetIt it = nonZeroComponent.begin(); it != nonZeroComponent.end(); ++it)
    {
      Node i = *it;
      // skip if node was already considered or its x-value is 0
      if (!(*_pNodeBoolMap)[i]) continue;
      
      const double x_i_value = x_values[_nodeMap[i]];
      
      _pBK->setTarget(i);
      _pBK->setCap(*_pCap);
      
      bool nestedCut = false;
      bool first = true;
//      while (true)
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
          NodeSet fwdS;
          ArcSet fwdDS;
          determineFwdCutSet(_d, *_pBK, _root, *_pMarked, fwdDS, fwdS);
          
          // numerical instability may cause minCutValue < x_i_value
          // even though there is nothing to cut
          if (fwdDS.empty()) break;
          
          // determine N (backward)
          ArcSet bwdDS;
//          determineBwdCutSet(_h, *_pBK, diRoot, _h2g, _marked, bwdDS, bwdS);
          
          bool backCuts = fwdDS.size() != bwdDS.size() || fwdDS != bwdDS;
          
          // add violated constraints for all nodes j in fwdS with x_j >= x_i
          constructRHS(rhs, fwdDS);
          for (NodeSetIt it2 = fwdS.begin(); it2 != fwdS.end(); ++it2)
          {
            const Node j = *it2;
            if (j == _root) continue;
            const double x_j_value = x_values[_nodeMap[j]];

            if (_tol.less(minCutValue, x_j_value))
            {
              //              std::cerr << std::endl << "nonZeroComp" << std::endl;
              //              this->printNodeSet(nonZeroComponent, _y, y_values);
              //
              //              std::cerr << std::endl << "S" << std::endl;
              //              this->printNodeSet(fwdS, _y, y_values);
              //
              //              std::cerr << std::endl << "dS" << std::endl;
              //              this->printNodeSet(fwdDS, _x, x_values);

              assert(isValid(j, fwdDS));
              //              std::cerr << x_j_value - minCutValue << std::endl;
              
              _pNodeBoolMap->set(j, false);
              add(_x[_nodeMap[j]] <= rhs, IloCplex::UseCutPurge).end();
              
              ++nCuts;
              if (nestedCut) ++nNestedCuts;
            }
          }
          
          if (backCuts)
          {
            // add violated constraints for all nodes j in bwdS with x_j >= x_i
//            constructRHS(rhs, bwdDS, bwdS);
//            for (NodeSetIt it2 = bwdS.begin(); it2 != bwdS.end(); ++it2)
//            {
//              const Node j = *it2;
//              const double x_j_value = x_values[_nodeMap[j]];
//              
//              if (_tol.less(minCutValue, x_j_value))
//              {
////                assert(isValid(j, bwdDS, bwdS));
//                
//                _pNodeBoolMap->set(j, false);
//                add(_x[_nodeMap[j]] <= rhs, IloCplex::UseCutPurge).end();
//                
//                ++nCuts;
//                ++nBackCuts;
//              }
//            }
          }
          
          // generate nested-cuts
//          for (NodeSetIt nodeIt = fwdS.begin(); nodeIt != fwdS.end(); ++nodeIt)
//          {
//            nestedCut = true;
//            // update the capactity to generate nested-cuts
////            _pBK->incCap(DiOutArcIt(_h, (*_pG2h1)[*nodeIt]), 1);
//            _pBK->incCap((*_pG2hRootArc)[*nodeIt], 1);
//          }
        }
//        else
//        {
//          break;
//        }
      }
    }
    
    rhs.end();
  }
  
  void separate()
  {
    IloNumArray x_values(getEnv(), _n);
    getValues(x_values, _x);
    
    IloNumArray z_values(getEnv(), _m);
    getValues(z_values, _z);
    
    // determine connected components
    NodeSet rootNodes = computeCapacities(*_pCap, x_values, z_values);
    NodeSetVector nonZeroComponents = determineStronglyConnectedComponents(x_values);
    
    int nCuts = 0;
    int nBackCuts = 0;
    int nNestedCuts = 0;

    for (NodeSetVectorIt it = nonZeroComponents.begin(); it != nonZeroComponents.end(); ++it)
    {
      const NodeSet& nonZeroComponent = *it;
      
      NodeSet intersection;
      std::set_intersection(rootNodes.begin(), rootNodes.end(),
                            nonZeroComponent.begin(), nonZeroComponent.end(),
                            std::inserter(intersection, intersection.begin()));
      
      if (_nodeNumber == 0 || intersection.size() > 0)
      {
        separateMinCut(nonZeroComponent, x_values, z_values, nCuts, nBackCuts, nNestedCuts);
      }
      else
      {
        separateStronglyConnectedComponent(nonZeroComponent, rootNodes, x_values, z_values, *this, nCuts);
      }
    }
    
    x_values.end();
    z_values.end();
    
//    std::cerr << "[";
//    for (NodeSetVectorIt it = nonZeroComponents.begin(); it != nonZeroComponents.end(); ++it)
//    {
//      NodeSet intersection;
//      std::set_intersection(rootNodes.begin(), rootNodes.end(),
//                            it->begin(), it->end(),
//                            std::inserter(intersection, intersection.begin()));
//      
//      std::cerr << " " << it->size();
//      if (intersection.size() > 0) std::cerr << "*";
//    }
//    std::cerr << " ]" << std::endl;
    
//    if (nCuts != 0 )
//    {
//      std::cerr <<  "#" << _cutCount << ", #comp = " << nonZeroComponents.size()
//                << ": generated " << nCuts
//                << " user cuts of which " << nBackCuts << " are back-cuts and "
//                << nNestedCuts << " are nested cuts" << std::endl;
//    }
    //std::cerr << "Time: " << t.realTime() << "s" << std::endl;
  }
  
  
//  void determineBwdCutSet(const Digraph& h,
//                          const BkAlg& bk,
//                          const DiNode diRoot,
//                          const DiNodeNodeMap& h2g,
//                          DiBoolNodeMap& marked,
//                          NodeSet& dS,
//                          NodeSet& S)
//  {
//    DiNode target = bk.getTarget();
//    DiNodeList diS;
//    determineBwdCutSet(h, bk, diRoot, target, marked, diS);
//    
//    for (DiNodeListIt nodeIt = diS.begin(); nodeIt != diS.end(); nodeIt++)
//    {
//      DiNode v = *nodeIt;
//      assert(marked[v]);
//      if (v == target) continue;
//      
//      for (DiInArcIt a(h, v); a != lemon::INVALID; ++a)
//      {
//        DiNode u = h.source(a);
//        if (u == diRoot)
//        {
//          assert(!marked[u]);
//          //std::cout << _h.id(u) << " -> "
//          //          << _h.id(v) << " "
//          //          << bk.flow(a) << "/" << bk.cap(a) << std::endl;
//          S.insert(h2g[v]);
//        }
//        else if (!marked[u])
//        {
//          //std::cout << _h.id(u) << " -> "
//          //          << _h.id(v) << " "
//          //          << bk.flow(a) << "/" << bk.cap(a) << std::endl;
//          dS.insert(h2g[v]);
//        }
//      }
//    }
//  }
//  
  void determineFwdCutSet(const Graph& h,
                          const BkAlg& bk,
                          const Node root,
                          BoolNodeMap& marked,
                          ArcSet& dS,
                          NodeSet& S)
  {
    Node target = bk.getTarget();
    NodeList diS;
    determineFwdCutSet(h, bk, root, marked, diS);
    
    for (NodeListIt nodeIt = diS.begin(); nodeIt != diS.end(); nodeIt++)
    {
      Node v = *nodeIt;
      assert(marked[v]);
      S.insert(v);
      
      for (OutArcIt a(h, v); a != lemon::INVALID; ++a)
      {
        Node w = h.target(a);
        if (!marked[w])
        {
          //std::cout << _h.id(v) << " -> "
          //          << _h.id(w) << " "
          //          << bk.flow(a) << "/" << bk.cap(a) << std::endl;
          dS.insert(a);
        }
      }
    }
  }
  
  NodeSet computeCapacities(CapacityMap& capacity,
                            IloNumArray x_values,
                            IloNumArray z_values)
  {
    NodeSet rootNodes;
    
    for (NodeIt v(_d); v != lemon::INVALID; ++v)
    {
      double val = x_values[_nodeMap[v]];
      _pNodeBoolMap->set(v, _tol.nonZero(val));
    }
    
    for (ArcIt a(_d); a != lemon::INVALID; ++a)
    {
      // cap((r,i)) = y_i
      double val = z_values[_arcMap[a]];
      if (!_tol.nonZero(val))
      {
        capacity[a] = val;
      }
      else
      {
        capacity[a] = 10 * _cutEpsilon;
      }
      
      if (_d.source(a) == _root && _tol.nonZero(val))
      {
        rootNodes.insert(_d.target(a));
      }
      capacity[a] = val;
    }
    
    return rootNodes;
  }
};

} // namespace mwcs
} // namespace nina

#endif // EDGECUTUNROOTED_H
