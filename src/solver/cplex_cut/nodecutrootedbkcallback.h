/*
 * nodecutrootedbkcallback.h
 *
 *  Created on: 18-feb-2014
 *      Author: M. El-Kebir
 */

#ifndef NODECUTROOTEDBKCALLBACK_H
#define NODECUTROOTEDBKCALLBACK_H

#include "nodecutrootedbk.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class NodeCutRootedBkLazyCallback : public IloCplex::LazyConstraintCallbackI,
                                    public NodeCutRootedBk<GR, NWGHT, NLBL, EWGHT>
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;
  typedef NodeCutRootedBk<GR, NWGHT, NLBL, EWGHT> Parent;

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
  typedef typename Parent::DoubleVector DoubleVector;
  typedef typename Parent::NodeSet NodeSet;
  typedef typename Parent::NodeSetIt NodeSetIt;
  typedef typename Parent::DiArcSet DiArcSet;
  typedef typename Parent::DiArcSetIt DiArcSetIt;
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
  using Parent::_g;
  using Parent::_weight;
  using Parent::_root;
  using Parent::_nodeMap;
  using Parent::_n;
  using Parent::_m;
  using Parent::_maxNumberOfCuts;
  using Parent::_comp;
  using Parent::_tol;
  using Parent::_h;
  using Parent::_cap;
  using Parent::_pG2h1;
  using Parent::_pG2h2;
  using Parent::_h2g;
  using Parent::_diRoot;
  using Parent::_pBK;
  using Parent::_marked;
  using Parent::_pMutex;
  using Parent::_pNodeFilterMap;
  using Parent::_pSubG;
  using Parent::_pComp;
  using Parent::lock;
  using Parent::unlock;
  using Parent::separate;
  using Parent::isUser;

public:
  NodeCutRootedBkLazyCallback(IloEnv env,
                              IloBoolVarArray x,
                              const Graph& g,
                              const WeightNodeMap& weight,
                              Node root,
                              const IntNodeMap& nodeMap,
                              int n,
                              int m,
                              int maxNumberOfCuts,
                              const IntNodeMap& comp,
                              IloFastMutex* pMutex)
    : IloCplex::LazyConstraintCallbackI(env)
    , Parent(x, g, weight, root, nodeMap, n, m, maxNumberOfCuts, comp, pMutex)
  {
  }

  NodeCutRootedBkLazyCallback(const NodeCutRootedBkLazyCallback& other)
    : IloCplex::LazyConstraintCallbackI(other)
    , Parent(other)
  {
  }

  virtual ~NodeCutRootedBkLazyCallback()
  {
  }

protected:
  virtual void main()
  {
    //std::cerr << "L: " << getNnodes() << std::endl;
    separate(*this);
  }

  virtual IloCplex::CallbackI* duplicateCallback() const
  {
    return (new (getEnv()) NodeCutRootedBkLazyCallback(*this));
  }

  virtual void addConstraint(IloConstraint con)
  {
    add(con);
  }

  virtual bool isUser() const
  {
    return false;
  }
};

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class NodeCutRootedBkUserCallback : public IloCplex::UserCutCallbackI,
                                    public NodeCutRootedBk<GR, NWGHT, NLBL, EWGHT>
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;
  typedef NodeCutRootedBk<GR, NWGHT, NLBL, EWGHT> Parent;

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
  typedef typename Parent::DoubleVector DoubleVector;
  typedef typename Parent::NodeSet NodeSet;
  typedef typename Parent::NodeSetIt NodeSetIt;
  typedef typename Parent::DiArcSet DiArcSet;
  typedef typename Parent::DiArcSetIt DiArcSetIt;
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
  using Parent::_g;
  using Parent::_weight;
  using Parent::_root;
  using Parent::_nodeMap;
  using Parent::_n;
  using Parent::_m;
  using Parent::_maxNumberOfCuts;
  using Parent::_comp;
  using Parent::_tol;
  using Parent::_h;
  using Parent::_cap;
  using Parent::_pG2h1;
  using Parent::_pG2h2;
  using Parent::_h2g;
  using Parent::_diRoot;
  using Parent::_pBK;
  using Parent::_marked;
  using Parent::_pMutex;
  using Parent::_pNodeFilterMap;
  using Parent::_pSubG;
  using Parent::_pComp;
  using Parent::lock;
  using Parent::unlock;
  using Parent::separate;
  using Parent::isUser;

public:
  NodeCutRootedBkUserCallback(IloEnv env,
                              IloBoolVarArray x,
                              const Graph& g,
                              const WeightNodeMap& weight,
                              Node root,
                              const IntNodeMap& nodeMap,
                              int n,
                              int m,
                              int maxNumberOfCuts,
                              const IntNodeMap& comp,
                              IloFastMutex* pMutex)
    : IloCplex::UserCutCallbackI(env)
    , Parent(x, g, weight, root, nodeMap, n, m, maxNumberOfCuts, comp, pMutex)
  {
  }

  NodeCutRootedBkUserCallback(const NodeCutRootedBkUserCallback& other)
    : IloCplex::UserCutCallbackI(other)
    , Parent(other)
  {
  }

  virtual ~NodeCutRootedBkUserCallback()
  {
  }

protected:
  virtual void main()
  {
    //std::cerr << "U: " << getNnodes() << std::endl;
    separate(*this);
  }

  virtual IloCplex::CallbackI* duplicateCallback() const
  {
    return (new (getEnv()) NodeCutRootedBkUserCallback(*this));
  }

  virtual void addConstraint(IloConstraint con)
  {
    addLocal(con);
  }

  virtual bool isUser() const
  {
    return true;
  }
};

} // namespace mwcs
} // namespace nina

#endif // NODECUTROOTEDBKCALLBACK_H
