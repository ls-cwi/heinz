/*
 * heuristicunrooted.h
 *
 *  Created on: 27-feb-2014
 *      Author: M. El-Kebir
 */

#ifndef HEURISTICUNROOTED_H
#define HEURISTICUNROOTED_H

#include <ilcplex/ilocplex.h>
#include <ilcplex/ilocplex.h>
#include <ilconcert/ilothread.h>
#include <lemon/adaptors.h>
#include <lemon/bfs.h>
#include <lemon/tolerance.h>
#include <set>
#include "heuristicrooted.h"

namespace nina {
namespace mwcs {
  
template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class HeuristicUnrooted : public HeuristicRooted<GR, NWGHT, NLBL, EWGHT>
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;
  typedef HeuristicRooted<GR, NWGHT, NLBL, EWGHT> Parent;
  
  using Parent::_env;
  using Parent::_x;
  using Parent::_z;
  using Parent::_g;
  using Parent::_weight;
  using Parent::_root;
  using Parent::_nodeMap;
  using Parent::_edgeMap;
  using Parent::_n;
  using Parent::_m;
  using Parent::_pMutex;
  using Parent::hasIncumbent;
  using Parent::getIncumbentObjValue;
  using Parent::computeMaxWeightConnectedSubtree;
  using Parent::computeEdgeWeights;
  using Parent::computeMinimumCostSpanningTree;
  using Parent::getValue;
  using Parent::getValues;
  using Parent::setCplexSolution;
  
protected:
  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  typedef typename Parent::SubGraphType SubGraphType;
  typedef typename Parent::MwcsSubGraphType MwcsSubGraphType;
  typedef typename Parent::MwcsSubTreeSolver MwcsSubTreeSolver;
  typedef typename Parent::NodeSet NodeSet;
  typedef typename Parent::NodeSetIt NodeSetIt;
  
public:
  HeuristicUnrooted(IloEnv env,
                    IloBoolVarArray x,
                    IloBoolVarArray y,
                    IloBoolVarArray z,
                    const Graph& g,
                    const WeightNodeMap& weight,
                    const IntNodeMap& nodeMap,
                    const IntEdgeMap& edgeMap,
                    int n,
                    int m,
                    IloFastMutex* pMutex)
    : Parent(env, x, z, g, weight, lemon::INVALID, nodeMap, edgeMap, n, m, pMutex)
    , _y(y)
    , _tol(_epsilon)
  {
  }
    
  HeuristicUnrooted(const HeuristicUnrooted& other)
    : Parent(other)
    , _y(other._y)
    , _tol(other._tol)
  {
  }
    
  ~HeuristicUnrooted()
  {
  }
    
protected:
  virtual void main()
  {
    computeEdgeWeights();
    computeMinimumCostSpanningTree();
    determineRootNode();
    
    IloBoolVarArray solutionVar(_env, 0);
    IloNumArray solution(_env, 0);
    double solutionWeight = hasIncumbent() ? getIncumbentObjValue() : -1;
    
    if (computeMaxWeightConnectedSubtree(solutionVar, solution, solutionWeight))
    {
      solutionVar.add(_y);
      solution.add(_y.getSize(), 0);
      solution[_x.getSize() + _z.getSize() + _nodeMap[_root]] = 1;
      
      std::cout << "Root: " << _y[_nodeMap[_root]].getName()
                << " = " << getValue(_y[_nodeMap[_root]])
                << " and found solution with weight: " << solutionWeight << std::endl;
      setCplexSolution(solutionVar, solution, solutionWeight);
    }
    
    solutionVar.end();
    solution.end();
  }
    
  virtual IloCplex::CallbackI* duplicateCallback() const
  {
    return (new (_env) HeuristicUnrooted(*this));
  }
  
  void determineRootNode()
  {
    IloNumArray y_values(_env, _n);
    getValues(y_values, _y);
    
    for (NodeIt i(_g); i != lemon::INVALID; ++i)
    {
      if (_tol.nonZero(y_values[_nodeMap[i]]))
      {
        _root = i;
        break;
      }
    }
    
    y_values.end();
  }
  
protected:
  IloBoolVarArray _y;
  const lemon::Tolerance<double> _tol;
  
  // 1e-5 is the epsilon that CPLEX uses (for deciding integrality),
  // i.e. if |x| < 1e-5 it's considered to be 0 by CPLEX.
  // if 1 - |x| < 1e-5 it's considered to be 1 by CPLEX.
  // we use the same epsilon to separate violated cuts.
  static const double _epsilon = 1e-5;
};
    
} // namespace mwcs
} // namespace nina

#endif // HEURISTICUNROOTED_H
