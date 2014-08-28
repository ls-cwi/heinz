/*
 * heuristicrooted.h
 *
 *  Created on: 27-feb-2014
 *      Author: M. El-Kebir
 */

#ifndef HEURISTICROOTED_H
#define HEURISTICROOTED_H

#include <ilcplex/ilocplex.h>
#include <ilcplex/ilocplex.h>
#include <ilconcert/ilothread.h>
#include <lemon/adaptors.h>
#include <lemon/bfs.h>
#include <lemon/kruskal.h>
#include <set>
#include "solver/impl/treesolverrootedimpl.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class HeuristicRooted : public IloCplex::HeuristicCallbackI
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;
  
protected:
  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  typedef lemon::FilterEdges<const Graph, const BoolEdgeMap> SubGraphType;
  typedef MwcsGraph<const SubGraphType, const WeightNodeMap, LabelNodeMap, DoubleEdgeMap> MwcsSubGraphType;
  typedef TreeSolverRootedImpl<const SubGraphType, const WeightNodeMap, LabelNodeMap, DoubleEdgeMap> TreeSolverRootedImplType;
  typedef typename MwcsSubGraphType::BoolNodeMap SubBoolNodeMap;
  typedef typename std::set<Node> NodeSet;
  typedef typename NodeSet::const_iterator NodeSetIt;
  
public:
  HeuristicRooted(IloEnv env,
                  IloBoolVarArray x,
//                  IloBoolVarArray z,
                  const Graph& g,
                  const WeightNodeMap& weight,
                  NodeSet rootNodes,
                  const IntNodeMap& nodeMap,
//                  const IntEdgeMap& edgeMap,
                  int n,
                  int m,
                  IloFastMutex* pMutex)
    : IloCplex::HeuristicCallbackI(env)
    , _x(x)
//    , _z(z)
    , _g(g)
    , _weight(weight)
    , _rootNodes(rootNodes)
    , _nodeMap(nodeMap)
//    , _edgeMap(edgeMap)
    , _n(n)
    , _m(m)
    , _pEdgeCost(NULL)
    , _pEdgeFilterMap(NULL)
    , _pSubG(NULL)
    , _pSubSolutionMap(NULL)
    , _pMwcsSubGraph(NULL)
    , _pMwcsSubTreeSolver(NULL)
    , _pMutex(pMutex)
  {
    lock();
    _pEdgeCost = new DoubleEdgeMap(_g);
    _pEdgeFilterMap = new BoolEdgeMap(_g, false);
    _pSubG = new SubGraphType(_g, *_pEdgeFilterMap);
    _pSubSolutionMap = new SubBoolNodeMap(*_pSubG, false);
    _pMwcsSubGraph = new MwcsSubGraphType();
    _pMwcsSubGraph->init(_pSubG, NULL, &_weight, NULL);
    _pMwcsSubTreeSolver = new TreeSolverRootedImplType();
    unlock();
  }
  
  HeuristicRooted(const HeuristicRooted& other)
    : IloCplex::HeuristicCallbackI(other._env)
    , _x(other._x)
//    , _z(other._z)
    , _g(other._g)
    , _weight(other._weight)
    , _rootNodes(other._rootNodes)
    , _nodeMap(other._nodeMap)
//    , _edgeMap(other._edgeMap)
    , _n(other._n)
    , _m(other._m)
    , _pEdgeCost(NULL)
    , _pEdgeFilterMap(NULL)
    , _pSubG(NULL)
    , _pSubSolutionMap(NULL)
    , _pMwcsSubGraph(NULL)
    , _pMwcsSubTreeSolver(NULL)
    , _pMutex(other._pMutex)
  {
    lock();
    _pEdgeCost = new DoubleEdgeMap(_g);
    _pEdgeFilterMap = new BoolEdgeMap(_g, false);
    _pSubG = new SubGraphType(_g, *_pEdgeFilterMap);
    _pSubSolutionMap = new SubBoolNodeMap(*_pSubG, false);
    _pMwcsSubGraph = new MwcsSubGraphType();
    _pMwcsSubGraph->init(_pSubG, NULL, &_weight, NULL);
    _pMwcsSubTreeSolver = new TreeSolverRootedImplType();
    unlock();
  }
  
  ~HeuristicRooted()
  {
    lock();
    delete _pMwcsSubTreeSolver;
    delete _pMwcsSubGraph;
    delete _pEdgeCost;
    delete _pEdgeFilterMap;
    delete _pSubSolutionMap;
    delete _pSubG;
    unlock();
  }
  
protected:
  virtual void main()
  {
    computeEdgeWeights();
    computeMinimumCostSpanningTree();
    
    IloBoolVarArray solutionVar(_env, 0);
    IloNumArray solution(_env, 0);
    double solutionWeight = hasIncumbent() ? getIncumbentObjValue() : -1;
    
    if (computeMaxWeightConnectedSubtree(solutionVar, solution, solutionWeight))
    {
      setCplexSolution(solutionVar, solution, solutionWeight);
    }
    
    solutionVar.end();
    solution.end();
  }

  virtual IloCplex::CallbackI* duplicateCallback() const
  {
    return (new (_env) HeuristicRooted(*this));
  }
  
  void lock()
  {
    if (_pMutex)
      _pMutex->lock();
  }
  
  void unlock()
  {
    if (_pMutex)
      _pMutex->unlock();
  }
  
  void setCplexSolution(IloBoolVarArray solutionVar, IloNumArray solution, double solutionWeight)
  {
    setSolution(solutionVar, solution, solutionWeight);
  }
  
  void computeEdgeWeights()
  {
    IloNumArray x_values(_env, _n);
    getValues(x_values, _x);
    
    for (EdgeIt e(_g); e != lemon::INVALID; ++e)
    {
      Node u = _g.u(e);
      Node v = _g.v(e);
      
      double x_u = x_values[_nodeMap[u]];
      double x_v = x_values[_nodeMap[v]];
      
      _pEdgeCost->set(e, 2 - (x_u + x_v));
    }
    
    x_values.end();
  }
  
  void computeMinimumCostSpanningTree()
  {
    lock();
    lemon::kruskal(_g, *_pEdgeCost, *_pEdgeFilterMap);
    unlock();
  }
  
  bool computeMaxWeightConnectedSubtree(IloBoolVarArray& solutionVar,
                                        IloNumArray& solution,
                                        double& solutionWeight)
  {
    _pMwcsSubTreeSolver->init(*_pMwcsSubGraph, _rootNodes);
    
    NodeSet solutionSet;
    double score;
    
    _pMwcsSubTreeSolver->solve(score, *_pSubSolutionMap, solutionSet);
    
    if (score > solutionWeight)
    {
      solutionWeight = score;
      solutionVar.add(_x);
      solution.add(_x.getSize(), 0);

//      solutionVar.add(_z);
//      solution.add(_z.getSize(), 0);
      
      for (NodeSetIt it = solutionSet.begin(); it != solutionSet.end(); ++it)
      {
        solution[_nodeMap[*it]] = 1;
      }
      
//      for (EdgeIt e(_g); e != lemon::INVALID; ++e)
//      {
//        if (module.find(_g.u(e)) != module.end() && module.find(_g.v(e)) != module.end())
//        {
//          solution[_x.getSize() + _edgeMap[e]] = 1;
//        }
//      }
      
      return true;
    }
    
    return false;
  }
  
protected:
  IloBoolVarArray _x;
//  IloBoolVarArray _z;
  const Graph& _g;
  const WeightNodeMap& _weight;
  NodeSet _rootNodes;
  const IntNodeMap& _nodeMap;
//  const IntEdgeMap& _edgeMap;
  const int _n;
  const int _m;
  DoubleEdgeMap* _pEdgeCost;
  BoolEdgeMap* _pEdgeFilterMap;
  const SubGraphType* _pSubG;
  SubBoolNodeMap* _pSubSolutionMap;
  MwcsSubGraphType* _pMwcsSubGraph;
  TreeSolverRootedImplType* _pMwcsSubTreeSolver;
  IloFastMutex* _pMutex;
};
  
} // namespace mwcs
} // namespace nina

#endif // HEURISTICROOTED_H
