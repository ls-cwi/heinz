/*
 * treeheurisiticsolverrooted.h
 *
 *  Created on: 11-jan-2013
 *      Author: M. El-Kebir
 */

#ifndef TREEHEURISTICSOLVERROOTED_H
#define TREEHEURISTICSOLVERROOTED_H

#include "solverrooted.h"
#include "treeheuristicsolver.h"
#include <set>
#include <vector>
#include <map>
#include <limits>
#include <assert.h>
#include <lemon/bfs.h>
#include <lemon/connectivity.h>
#include <lemon/adaptors.h>
#include <lemon/kruskal.h>
#include <lemon/random.h>
#include "treesolver.h"
#include "analysis.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class TreeHeuristicSolverRooted : public SolverRooted<GR, NWGHT, NLBL, EWGHT>,
                                  public TreeHeuristicSolver<GR, NWGHT, NLBL, EWGHT>
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;

  typedef Solver<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> GrandParent;
  typedef SolverRooted<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> Parent1;
  typedef TreeHeuristicSolver<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> Parent2;

  typedef typename Parent1::NodeSet NodeSet;
  typedef typename Parent1::NodeSetIt NodeSetIt;
  typedef typename Parent1::MwcsGraphType MwcsGraphType;
  typedef typename Parent2::EdgeHeuristic EdgeHeuristic;
  typedef typename Parent2::SubGraphType SubGraphType;
  typedef typename Parent2::MwcsSubGraphType MwcsSubGraphType;
  typedef typename Parent2::MwcsSubTreeSolverType MwcsSubTreeSolverType;
  typedef typename Parent2::MwcsAnalyzeType MwcsAnalyzeType;
  
  TEMPLATE_GRAPH_TYPEDEFS(Graph);

  using GrandParent::_mwcsGraph;
  using Parent1::_rootNodes;
  using Parent1::_score;
  using Parent1::_solutionMap;
  using Parent1::_solutionSet;
  using Parent1::printSolution;
  using Parent1::getSolutionWeight;
  using Parent1::getSolutionNodeMap;
  using Parent1::getSolutionModule;
  using Parent1::isNodeInSolution;
  using Parent1::init;
  using Parent2::_filterMap;
  using Parent2::_subG;
  using Parent2::_pMwcsSubGraph;
  using Parent2::_pMwcsSubSolver;
  using Parent2::_pEdgeCost;
  using Parent2::setEdgeCosts;

public:
  TreeHeuristicSolverRooted(const MwcsGraphType& mwcsGraph,
                            const NodeSet& rootNodes);
  TreeHeuristicSolverRooted(const MwcsGraphType& mwcsGraph,
                            const NodeSet& rootNodes,
                            EdgeHeuristic heuristic,
                            MwcsAnalyzeType* pAnalyze);

  virtual void init();
  virtual bool solve();
  virtual void setLowerBound(double LB) {}
};

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline TreeHeuristicSolverRooted<GR, NWGHT, NLBL, EWGHT>::TreeHeuristicSolverRooted(const MwcsGraphType& mwcsGraph,
                                                                                    const NodeSet& rootNodes)
  : GrandParent(mwcsGraph)
  , Parent1(mwcsGraph, rootNodes)
  , Parent2(mwcsGraph)
{
}
  
template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline TreeHeuristicSolverRooted<GR, NWGHT, NLBL, EWGHT>::TreeHeuristicSolverRooted(const MwcsGraphType& mwcsGraph,
                                                                                    const NodeSet& rootNodes,
                                                                                    EdgeHeuristic heuristic,
                                                                                    MwcsAnalyzeType* pAnalyze)
  : GrandParent(mwcsGraph)
  , Parent1(mwcsGraph, rootNodes)
  , Parent2(mwcsGraph, heuristic, pAnalyze)
{
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void TreeHeuristicSolverRooted<GR, NWGHT, NLBL, EWGHT>::init()
{
  // construct mwcssubgraph instance
  delete _pMwcsSubGraph;
  _pMwcsSubGraph = new MwcsSubGraphType();
  _pMwcsSubGraph->init(&_subG, NULL, &_mwcsGraph.getScores(), NULL);

  // construct mwcssubsolver instance
  delete _pMwcsSubSolver;
  _pMwcsSubSolver = new MwcsSubTreeSolverType(*_pMwcsSubGraph, _rootNodes);
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline bool TreeHeuristicSolverRooted<GR, NWGHT, NLBL, EWGHT>::solve()
{
  const Graph& g = _mwcsGraph.getGraph();

  // compute minimum spanning tree
  lemon::kruskal(g, *_pEdgeCost, _filterMap);

  _pMwcsSubSolver->init();
  bool res = _pMwcsSubSolver->solve();
  _score = _pMwcsSubSolver->getSolutionWeight();
  lemon::mapCopy(g, _pMwcsSubSolver->getSolutionNodeMap(), _solutionMap);
  _solutionSet = _pMwcsSubSolver->getSolutionModule();

  return res;
}

} // namespace mwcs
} // namespace nina

#endif // TREEHEURISTICSOLVERROOTED_H
