/*
 * nodecutunrooted.h
 *
 *  Created on: 18-feb-2014
 *      Author: M. El-Kebir, G. W. Klau
 */

#ifndef NODECUTUNROOTED_H
#define NODECUTUNROOTED_H

#include "nodecutlazy.h"
#include "nodecutuser.h"
#include <algorithm>

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
  
  typedef typename Parent::NodeSet NodeSet;
  typedef typename Parent::NodeSetIt NodeSetIt;
  typedef typename Parent::NodeSetVector NodeSetVector;
  typedef typename Parent::NodeSetVectorIt NodeSetVectorIt;
  typedef typename Parent::SubGraph SubGraph;
  typedef typename Parent::SubNodeIt SubNodeIt;

  using Parent::_x;
  using Parent::_y;
  using Parent::_g;
  using Parent::_weight;
  using Parent::_nodeMap;
  using Parent::_n;
  using Parent::_k;
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
  using Parent::determineConnectedComponents;
  using Parent::separateConnectedComponent;
  
  friend class NodeCut<GR, NWGHT, NLBL, EWGHT>;

public:
  NodeCutUnrootedLazyConstraint(IloEnv env,
                                IloBoolVarArray x,
                                IloBoolVarArray y,
                                const Graph& g,
                                const WeightNodeMap& weight,
                                const IntNodeMap& nodeMap,
                                int n,
                                int k,
                                int maxNumberOfCuts,
                                IloFastMutex* pMutex)
    : Parent(env, x, y, g, weight, nodeMap, n, k, maxNumberOfCuts, pMutex)
  {
  }

  NodeCutUnrootedLazyConstraint(const NodeCutUnrootedLazyConstraint& other)
    : Parent(other)
  {
  }

  virtual ~NodeCutUnrootedLazyConstraint()
  {
  }

protected:
  virtual void main()
  {
    separate();
  }

  virtual IloCplex::CallbackI* duplicateCallback() const
  {
    return (new (getEnv()) NodeCutUnrootedLazyConstraint(*this));
  }

  void separate()
  {
    IloNumArray x_values(getEnv(), _n);
    getValues(x_values, _x);
    
    IloNumArray y_values(getEnv(), _n);
    getValues(y_values, _y);
    
    // determine non-zero y-vars
    NodeSet rootNodes;
    for (NodeIt i(_g); i != lemon::INVALID; ++i)
    {
      int idx_i = _nodeMap[i];
      if (_tol.nonZero(y_values[idx_i]))
      {
        rootNodes.insert(i);
      }
    }
    
    cout << "[in lazy callback] k = " << _k << endl;
    
    assert(rootNodes.size() <= _k);
    
    // determine connected components
    NodeSetVector nonZeroComponents = determineConnectedComponents(x_values);

    int nCuts = 0;
    for (NodeSetVectorIt it = nonZeroComponents.begin(); it != nonZeroComponents.end(); ++it)
    {
      const NodeSet& nonZeroComponent = *it;
      
      NodeSet intersection;
      std::set_intersection(rootNodes.begin(), rootNodes.end(),
                            nonZeroComponent.begin(), nonZeroComponent.end(),
                            std::inserter(intersection, intersection.begin()));
      
      if (intersection.size() == 0)
      {
        separateConnectedComponent(*it, rootNodes, x_values, y_values, *this, nCuts);
      }
    }
    
    x_values.end();
    y_values.end();
    
//    std::cerr << "#comps: " << nonZeroComponents.size() << ", generated " << nCuts << " lazy cuts" << std::endl;
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
  typedef typename Parent::NodeSet NodeSet;
  typedef typename Parent::NodeSetIt NodeSetIt;
  typedef typename Parent::NodeSetVector NodeSetVector;
  typedef typename Parent::NodeSetVectorIt NodeSetVectorIt;
  typedef typename Parent::SubGraph SubGraph;
  typedef typename Parent::SubNodeIt SubNodeIt;
  typedef typename Parent::SubEdgeIt SubEdgeIt;
  typedef typename Parent::NodeQueue NodeQueue;
  typedef typename Parent::DiNodeQueue DiNodeQueue;
  typedef typename Parent::DiNodeSet DiNodeSet;
  typedef typename Parent::DiNodeSetIt DiNodeSetIt;
  typedef typename Parent::DiNodeList DiNodeList;
  typedef typename Parent::DiNodeListIt DiNodeListIt;
  typedef typename Parent::BkAlg BkAlg;
  typedef typename Parent::DiBoolNodeMap DiBoolNodeMap;

  using Parent::_x;
  using Parent::_y;
  using Parent::_g;
  using Parent::_weight;
  using Parent::_nodeMap;
  using Parent::_n;
  using Parent::_k;
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
  using Parent::_diRootSet;
  using Parent::_pBK;
  using Parent::_marked;
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
  using Parent::determineConnectedComponents;
  using Parent::separateConnectedComponent;
  
  friend class NodeCut<GR, NWGHT, NLBL, EWGHT>;

public:
  NodeCutUnrootedUserCut(IloEnv env,
                         IloBoolVarArray x,
                         IloBoolVarArray y,
                         const Graph& g,
                         const WeightNodeMap& weight,
                         const IntNodeMap& nodeMap,
                         int n,
                         int k,
                         int maxNumberOfCuts,
                         IloFastMutex* pMutex,
                         BackOff backOff)
    : Parent(env, x, y, g, weight, nodeMap, n, k, maxNumberOfCuts, pMutex, backOff)
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
    
    for (DiNodeSetIt diRootIt = other._diRootSet.begin();
         diRootIt != other._diRootSet.end(); ++diRootIt)
    {
      _diRootSet.insert(nodeMap[*diRootIt]);
    }
    
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
  
  void separateMinCut(const NodeSet& nonZeroComponent,
                      const IloNumArray& x_values,
                      const IloNumArray& y_values,
                      int& nCuts, int& nBackCuts, int& nNestedCuts)
  {
    IloExpr rhs(getEnv());
    
    assert(_diRootSet.size() == 1);
    DiNode diRoot = *_diRootSet.begin();

    _pBK->setSource(diRoot);
    for (NodeSetIt it = nonZeroComponent.begin(); it != nonZeroComponent.end(); ++it)
    {
      Node i = *it;
      // skip if node was already considered or its x-value is 0
      if (!(*_pNodeBoolMap)[i]) continue;
      
      const double x_i_value = x_values[_nodeMap[i]];
      
      _pBK->setTarget((*_pG2h2)[i]);
      _pBK->setCap(_cap);
      
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
          NodeSet fwdS, fwdDS;
          determineFwdCutSet(_h, *_pBK, diRoot, _h2g, _marked, fwdDS, fwdS);
          
          // numerical instability may cause minCutValue < x_i_value
          // even though there is nothing to cut
          if (fwdS.empty() && fwdDS.empty()) break;
          
          // determine N (backward)
          NodeSet bwdS, bwdDS;
          determineBwdCutSet(_h, *_pBK, diRoot, _h2g, _marked, bwdDS, bwdS);
          
          bool backCuts = fwdDS.size() != bwdDS.size() ||
          fwdS.size() != bwdS.size() ||
          fwdDS != bwdDS || bwdS != fwdS;
          
          // add violated constraints for all nodes j in fwdS with x_j >= x_i
          constructRHS(rhs, fwdDS, fwdS);
          for (NodeSetIt it2 = fwdS.begin(); it2 != fwdS.end(); ++it2)
          {
            const Node j = *it2;
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
              
              assert(isValid(j, fwdDS, fwdS));
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
            constructRHS(rhs, bwdDS, bwdS);
            for (NodeSetIt it2 = bwdS.begin(); it2 != bwdS.end(); ++it2)
            {
              const Node j = *it2;
              const double x_j_value = x_values[_nodeMap[j]];
              
              if (_tol.less(minCutValue, x_j_value))
              {
                assert(isValid(j, bwdDS, bwdS));
                
                _pNodeBoolMap->set(j, false);
                add(_x[_nodeMap[j]] <= rhs, IloCplex::UseCutPurge).end();
                
                ++nCuts;
                ++nBackCuts;
              }
            }
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
  
//  void separateMinCutLocal(const NodeSet& nonZeroComponent,
//                           const Node root,
//                           const IloNumArray& x_values,
//                           const IloNumArray& y_values,
//                           int& nCuts, int& nBackCuts, int& nNestedCuts)
//  {
//    IloExpr rhs(getEnv());
//    
//    Digraph h;
//    DiNodeNodeMap h2g(h);
//    DiNode diRoot;
//    CapacityMap cap(h);
//    
//    init(nonZeroComponent, root, x_values, y_values, h, h2g, diRoot, cap);
//    
//    BkAlg bk(h, cap);
//
//    bk.setSource(diRoot);
//    _pNodeBoolMap->set(root, false);
//    for (NodeSetIt it = nonZeroComponent.begin(); it != nonZeroComponent.end(); ++it)
//    {
//      Node i = *it;
//      // skip if node was already considered or its x-value is 0
//      if (!(*_pNodeBoolMap)[i]) continue;
//      
//      const double x_i_value = x_values[_nodeMap[i]];
//
//      bk.setTarget((*_pG2h2)[i]);
//      bk.setCap(cap);
//
//      bool nestedCut = false;
//      bool first = true;
//      while (true)
//      {
//        if (first)
//        {
//          bk.run();
//          first = false;
//        }
//        else
//        {
//          bk.run(false);
//        }
//        
//        // let's see if there's a violated constraint
//        double minCutValue = bk.maxFlow();
//        if (_tol.less(minCutValue, x_i_value))
//        {
//          // determine N (forward)
//          NodeSet fwdS, fwdDS;
//          determineFwdCutSet(h, bk, diRoot, h2g, _marked, fwdDS, fwdS);
//          
//          //          std::cout << x_i_value << "\t"
//          //                    << minCutValue << "\t"
//          //                    << fwdS.size() << "\t"
//          //                    << fwdDS.size() << std::endl;
//          
//          // numerical instability may cause minCutValue < x_i_value
//          // even though there is nothing to cut
//          if (fwdS.empty() && fwdDS.empty()) break;
//          
//          // determine N (backward)
//          NodeSet bwdS, bwdDS;
//          determineBwdCutSet(h, bk, diRoot, h2g, _marked, bwdDS, bwdS);
//          
//          bool backCuts = fwdDS.size() != bwdDS.size() ||
//          fwdS.size() != bwdS.size() ||
//          fwdDS != bwdDS || bwdS != fwdS;
//          
//          // add violated constraints for all nodes j in fwdS with x_j >= x_i
//          constructRHS(rhs, fwdDS, fwdS);
//          for (NodeSetIt it2 = fwdS.begin(); it2 != fwdS.end(); ++it2)
//          {
//            const Node j = *it2;
//            const double x_j_value = x_values[_nodeMap[j]];
//            
//            if (_tol.less(minCutValue, x_j_value))
//            {
////              std::cerr << std::endl << "nonZeroComp" << std::endl;
////              this->printNodeSet(nonZeroComponent, _y, y_values);
////
////              std::cerr << std::endl << "S" << std::endl;
////              this->printNodeSet(fwdS, _y, y_values);
////              
////              std::cerr << std::endl << "dS" << std::endl;
////              this->printNodeSet(fwdDS, _x, x_values);
//              
//              assert(isValid(j, fwdDS, fwdS));
////              std::cerr << x_j_value - minCutValue << std::endl;
//              
//              _pNodeBoolMap->set(j, false);
//              add(_x[_nodeMap[j]] <= rhs, IloCplex::UseCutPurge).end();
//              
//              ++nCuts;
//              if (nestedCut) ++nNestedCuts;
//            }
//          }
//          
//          if (backCuts)
//          {
//            // add violated constraints for all nodes j in bwdS with x_j >= x_i
//            constructRHS(rhs, bwdDS, bwdS);
//            for (NodeSetIt it2 = bwdS.begin(); it2 != bwdS.end(); ++it2)
//            {
//              const Node j = *it2;
//              const double x_j_value = x_values[_nodeMap[j]];
//              
//              if (_tol.less(minCutValue, x_j_value))
//              {
//                assert(isValid(j, bwdDS, bwdS));
//                
//                _pNodeBoolMap->set(j, false);
//                add(_x[_nodeMap[j]] <= rhs, IloCplex::UseCutPurge).end();
//                
//                ++nCuts;
//                ++nBackCuts;
//              }
//            }
//          }
//          
//          // generate nested-cuts
//          for (NodeSetIt nodeIt = fwdDS.begin(); nodeIt != fwdDS.end(); ++nodeIt)
//          {
//            nestedCut = true;
//            // update the capactity to generate nested-cuts
//            bk.incCap(DiOutArcIt(h, (*_pG2h1)[*nodeIt]), 1);
//          }
//        }
//        else
//        {
//          break;
//        }
//      }
//    }
//    
//    rhs.end();
//  }
  
  void separate()
  {
    cout << "*** [construction area] in user callback" << endl;
    IloNumArray x_values(getEnv(), _n);
    getValues(x_values, _x);
    
    IloNumArray y_values(getEnv(), _n);
    getValues(y_values, _y);
    
    // determine connected components
    NodeSet rootNodes = computeCapacities(_cap, x_values, y_values);
    NodeSetVector nonZeroComponents = determineConnectedComponents(x_values);
    
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
        separateMinCut(nonZeroComponent, x_values, y_values, nCuts, nBackCuts, nNestedCuts);
      }
      else
      {
        separateConnectedComponent(nonZeroComponent, rootNodes, x_values, y_values, *this, nCuts);
      }
    }
    
    x_values.end();
    y_values.end();
    
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
  
  
  void determineBwdCutSet(const Digraph& h,
                          const BkAlg& bk,
                          const DiNode diRoot,
                          const DiNodeNodeMap& h2g,
                          DiBoolNodeMap& marked,
                          NodeSet& dS,
                          NodeSet& S)
  {
    DiNode target = bk.getTarget();
    DiNodeList diS;
    determineBwdCutSet(h, bk, diRoot, target, marked, diS);
    
    for (DiNodeListIt nodeIt = diS.begin(); nodeIt != diS.end(); nodeIt++)
    {
      DiNode v = *nodeIt;
      assert(marked[v]);
      if (v == target) continue;
      
      for (DiInArcIt a(h, v); a != lemon::INVALID; ++a)
      {
        DiNode u = h.source(a);
        if (u == diRoot)
        {
          assert(!marked[u]);
          //std::cout << _h.id(u) << " -> "
          //          << _h.id(v) << " "
          //          << bk.flow(a) << "/" << bk.cap(a) << std::endl;
          S.insert(h2g[v]);
        }
        else if (!marked[u])
        {
          //std::cout << _h.id(u) << " -> "
          //          << _h.id(v) << " "
          //          << bk.flow(a) << "/" << bk.cap(a) << std::endl;
          dS.insert(h2g[v]);
        }
      }
    }
  }
  
  void determineFwdCutSet(const Digraph& h,
                          const BkAlg& bk,
                          const DiNode diRoot,
                          const DiNodeNodeMap& h2g,
                          DiBoolNodeMap& marked,
                          NodeSet& dS,
                          NodeSet& S)
  {
    DiNode target = bk.getTarget();
    DiNodeList diS;
    determineFwdCutSet(h, bk, diRoot, marked, diS);
    
    for (DiNodeListIt nodeIt = diS.begin(); nodeIt != diS.end(); nodeIt++)
    {
      DiNode v = *nodeIt;
      assert(marked[v]);
      
      if (v == diRoot)
      {
        for (DiOutArcIt a(h, v); a != lemon::INVALID; ++a)
        {
          DiNode w = h.target(a);
          if (!marked[w])
          {
            //std::cout << _h.id(v) << " -> "
            //          << _h.id(w) << " "
            //          << bk.flow(a) << "/" << bk.cap(a) << std::endl;
            S.insert(h2g[w]);
          }
        }
      }
      else
      {
        for (DiOutArcIt a(h, v); a != lemon::INVALID; ++a)
        {
          DiNode w = h.target(a);
          if (!marked[w] && w != target)
          {
            //std::cout << _h.id(v) << " -> "
            //          << _h.id(w) << " "
            //          << bk.flow(a) << "/" << bk.cap(a) << std::endl;
            dS.insert(h2g[w]);
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
    
    _diRootSet.clear();
    DiNode diRoot = _h.addNode();
    _diRootSet.insert(diRoot);
    
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
      
      DiArc ri1 = _h.addArc(diRoot, i1);
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
  
  void init(const NodeSet& nonZeroComponent,
            const Node root,
            const IloNumArray& x_values,
            const IloNumArray& y_values,
            Digraph& h,
            DiNodeNodeMap& h2g,
            DiNode& diRoot,
            CapacityMap& cap)
  {
    NodeSet shell;
    
    lemon::mapFill(_g, *_pG2h1, lemon::INVALID);
    lemon::mapFill(_g, *_pG2h2, lemon::INVALID);
    
    h.clear();
    diRoot = h.addNode();
    for (NodeSetIt it = nonZeroComponent.begin(); it != nonZeroComponent.end(); ++it)
    {
      const Node i = *it;
      
      double x_i_value = x_values[_nodeMap[i]];
      double y_i_value = y_values[_nodeMap[i]];
      
      DiNode i1 = h.addNode();
      DiNode i2 = h.addNode();
      _pG2h1->set(i, i1);
      _pG2h2->set(i, i2);
      h2g[i1] = i;
      h2g[i2] = i;
      
      DiArc i1i2 = h.addArc(i1, i2);
      cap[i1i2] = _tol.nonZero(x_i_value) ? x_i_value : 10 * _cutEpsilon;
      
      DiArc ri1 = h.addArc(diRoot, i1);
      cap[ri1] = _tol.nonZero(y_i_value) ? y_i_value : 10 * _cutEpsilon;
      
      for (IncEdgeIt e(_g, i); e != lemon::INVALID; ++e)
      {
        Node j = _g.oppositeNode(i, e);
        if ((*_pG2h1)[j] == lemon::INVALID && nonZeroComponent.find(j) == nonZeroComponent.end())
        {
          shell.insert(j);
          
          double x_j_value = x_values[_nodeMap[j]];
          double y_j_value = y_values[_nodeMap[j]];
          
          assert(!_tol.nonZero(x_j_value));
          
          DiNode j1 = h.addNode();
          DiNode j2 = h.addNode();
          _pG2h1->set(j, j1);
          _pG2h2->set(j, j2);
          h2g[j1] = j;
          h2g[j2] = j;
          
          DiArc j1j2 = h.addArc(j1, j2);
          cap[j1j2] = _tol.nonZero(x_j_value) ? x_j_value : 10 * _cutEpsilon;
          
          DiArc rj1 = h.addArc(diRoot, j1);
          // shell things may never end up in S unless they have a nonzero y
          cap[rj1] = _tol.nonZero(y_j_value) ? y_j_value : 1;
        }
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
      
      if (i1 != lemon::INVALID && j1 != lemon::INVALID)
      {
        assert(i2 != lemon::INVALID);
        assert(j2 != lemon::INVALID);
        
        DiArc i2j1 = h.addArc(i2, j1);
        DiArc j2i1 = h.addArc(j2, i1);
        cap[i2j1] = cap[j2i1] = 1;
      }
    }
  }
  
  NodeSet computeCapacities(CapacityMap& capacity,
                            IloNumArray x_values,
                            IloNumArray y_values)
  {
    NodeSet rootNodes;
    
    for (NodeIt v(_g); v != lemon::INVALID; ++v)
    {
      // cap((i,j)) = x_i
      double val = x_values[_nodeMap[v]];
      if (!_tol.nonZero(val))
      {
        _pNodeBoolMap->set(v, false);
        val = 10 * _cutEpsilon;
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
        val = 10 * _cutEpsilon;
      }
      else
      {
        rootNodes.insert(v);
      }
      capacity[(*_pG2hRootArc)[v]] = val;
    }
    
    return rootNodes;
  }
};

} // namespace mwcs
} // namespace nina

#endif // NODECUTUNROOTED_H
