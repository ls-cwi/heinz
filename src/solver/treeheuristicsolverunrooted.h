/*
 * treeheurisiticsolverunrooted.h
 *
 *  Created on: 11-jan-2013
 *      Author: M. El-Kebir
 */

#ifndef TREEHEURISTICSOLVERUNROOTED_H
#define TREEHEURISTICSOLVERUNROOTED_H

#include "solverunrooted.h"
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
#include "analysis.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class TreeHeuristicSolverUnrooted : public SolverUnrooted<GR, NWGHT, NLBL, EWGHT>,
                                    public TreeHeuristicSolver<GR, NWGHT, NLBL, EWGHT>
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;

  typedef Solver<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> GrandParent;
  typedef SolverUnrooted<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> Parent1;
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
  TreeHeuristicSolverUnrooted(const MwcsGraphType& mwcsGraph);
  TreeHeuristicSolverUnrooted(const MwcsGraphType& mwcsGraph,
                              EdgeHeuristic heuristic,
                              MwcsAnalyzeType* pAnalyze);

  virtual void init();
  virtual bool solve();
};

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline TreeHeuristicSolverUnrooted<GR, NWGHT, NLBL, EWGHT>::TreeHeuristicSolverUnrooted(const MwcsGraphType& mwcsGraph)
  : GrandParent(mwcsGraph)
  , Parent1(mwcsGraph)
  , Parent2(mwcsGraph)
{
}
  
template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline TreeHeuristicSolverUnrooted<GR, NWGHT, NLBL, EWGHT>::TreeHeuristicSolverUnrooted(const MwcsGraphType& mwcsGraph,
                                                                                        EdgeHeuristic heuristic,
                                                                                        MwcsAnalyzeType* pAnalyze)
  : GrandParent(mwcsGraph)
  , Parent1(mwcsGraph)
  , Parent2(mwcsGraph, heuristic, pAnalyze)
{
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void TreeHeuristicSolverUnrooted<GR, NWGHT, NLBL, EWGHT>::init()
{
  // construct mwcssubgraph instance
  delete _pMwcsSubGraph;
  _pMwcsSubGraph = new MwcsSubGraphType();
  _pMwcsSubGraph->init(&_subG, NULL, &_mwcsGraph.getScores(), NULL);

  // construct mwcssubsolver instance
  delete _pMwcsSubSolver;
  _pMwcsSubSolver = new MwcsSubTreeSolverType(*_pMwcsSubGraph, NodeSet());
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline bool TreeHeuristicSolverUnrooted<GR, NWGHT, NLBL, EWGHT>::solve()
{
  const Graph& g = _mwcsGraph.getGraph();

  // compute minimum spanning tree
  lemon::kruskal(g, *_pEdgeCost, _filterMap);

  bool res = false;

  for (NodeIt v(g); v != lemon::INVALID; ++v)
  {
    if (_mwcsGraph.getScore(v) > 0)
    {
      _pMwcsSubSolver->setRoot(v);
      _pMwcsSubSolver->init();
      res = _pMwcsSubSolver->solve();
      if (_score < _pMwcsSubSolver->getSolutionWeight())
      {
        _score = _pMwcsSubSolver->getSolutionWeight();
        lemon::mapCopy(g, _pMwcsSubSolver->getSolutionNodeMap(), _solutionMap);
        _solutionSet = _pMwcsSubSolver->getSolutionModule();
      }
    }
  }

  return res;
}

} // namespace mwcs
} // namespace nina

#endif // TREEHEURISTICSOLVERUNROOTED_H
