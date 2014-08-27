/*
 * treeheuristicsolverunrootedimpl.h
 *
 *  Created on: 11-jan-2013
 *      Author: M. El-Kebir
 */

#ifndef TREEHEURISTICSOLVERUNROOTEDIMPL_H
#define TREEHEURISTICSOLVERUNROOTEDIMPL_H

#include "solverunrootedimpl.h"
#include "treeheuristicsolverimpl.h"
#include "treesolverunrootedimpl.h"
#include <assert.h>
#include <lemon/kruskal.h>

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class TreeHeuristicSolverUnrootedImpl : public SolverUnrootedImpl<GR, NWGHT, NLBL, EWGHT>,
                                        public TreeHeuristicSolverImpl<GR, NWGHT, NLBL, EWGHT>
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;

  typedef SolverUnrootedImpl<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> Parent1;
  typedef TreeHeuristicSolverImpl<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> Parent2;
  
  typedef typename Parent1::NodeSet NodeSet;
  typedef typename Parent1::NodeSetIt NodeSetIt;
  typedef typename Parent1::MwcsGraphType MwcsGraphType;
  
  typedef typename Parent2::Options Options;
  typedef typename Parent2::EdgeHeuristic EdgeHeuristic;

  TEMPLATE_GRAPH_TYPEDEFS(Graph);

  using Parent1::_pMwcsGraph;
  
  using Parent2::_options;
  using Parent2::_pFilterMap;
  using Parent2::_pSubG;
  using Parent2::_pMwcsSubGraph;
  using Parent2::_pEdgeCost;
  using Parent2::solveMonteCarlo;

protected:
  typedef typename Parent2::SubGraphType SubGraphType;
  typedef typename Parent2::MwcsSubGraphType MwcsSubGraphType;
  typedef typename Parent2::MwcsAnalyzeType MwcsAnalyzeType;
  typedef typename Parent2::SubBoolNodeMap SubBoolNodeMap;

  typedef TreeSolverUnrootedImpl<const SubGraphType, const WeightNodeMap, LabelNodeMap, DoubleEdgeMap> MwcsUnrootedSubTreeSolverType;

public:
  TreeHeuristicSolverUnrootedImpl(const Options& options)
    : Parent1()
    , Parent2(options)
    , _mwcsUnrootedSubTreeSolver()
  {
  }

  void init(const MwcsGraphType& mwcsGraph)
  {
    Parent1::init(mwcsGraph);
    Parent2::init(mwcsGraph);
  }
  
  bool solve(double& score, BoolNodeMap& solutionMap, NodeSet& solutionSet)
  {
    return solveMonteCarlo(*_pMwcsGraph, _mwcsUnrootedSubTreeSolver, score, solutionMap, solutionSet);
  }
  
protected:
  void initSubTreeSolver()
  {
    _mwcsUnrootedSubTreeSolver.init(*_pMwcsSubGraph);
  }
  
  bool solveSubTree(double& score, SubBoolNodeMap& solutionMap, NodeSet& solutionSet)
  {
    return _mwcsUnrootedSubTreeSolver.solve(score, solutionMap, solutionSet);
  }
  
private:
  MwcsUnrootedSubTreeSolverType _mwcsUnrootedSubTreeSolver;
};

} // namespace mwcs
} // namespace nina

#endif // TREEHEURISTICSOLVERUNROOTED_H
