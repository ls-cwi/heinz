/*
 * treesolverrooted.h
 *
 *  Created on: 11-jan-2013
 *      Author: M. El-Kebir
 */

#ifndef TREESOLVERROOTED_H
#define TREESOLVERROOTED_H

#include "solverrooted.h"
#include "treesolver.h"
#include <set>
#include <map>
#include <limits>
#include <assert.h>
#include <lemon/bfs.h>
#include <lemon/connectivity.h>

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class TreeSolverRooted : public SolverRooted<GR, NWGHT, NLBL, EWGHT>,
                         public TreeSolver<GR, NWGHT, NLBL, EWGHT>
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;
  
  typedef Solver<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> GrandParent;
  typedef SolverRooted<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> Parent1;
  typedef TreeSolver<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> Parent2;

  typedef typename Parent1::MwcsGraphType MwcsGraphType;
  typedef typename Parent1::NodeSet NodeSet;
  typedef typename Parent1::NodeSetIt NodeSetIt;
  typedef typename Parent1::NodeVector NodeVector;
  typedef typename Parent1::NodeVectorIt NodeVectorIt;
  
  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  
  using Parent1::_mwcsGraph;
  using Parent1::_score;
  using Parent1::_solutionMap;
  using Parent1::_solutionSet;
  using Parent1::_rootNodes;
  using Parent1::printSolution;
  using Parent1::getSolutionWeight;
  using Parent1::getSolutionNodeMap;
  using Parent1::getSolutionModule;
  using Parent1::isNodeInSolution;
  
  using Parent2::init;
  using Parent2::_dpMap;
  using Parent2::_nodesPerLevel;
  using Parent2::_pred;
  using Parent2::_level;
  using Parent2::_bfs;
  using Parent2::_root;
  
protected:
  typedef typename Parent2::PredMap PredMap;
  typedef typename Parent2::DistMap DistMap;
  typedef typename Parent2::NodeSetVector NodeSetVector;
  typedef typename Parent2::NodeSetVectorIt NodeSetVectorIt;
  typedef typename Parent2::NodeSetVectorRevIt NodeSetVectorRevIt;
  typedef typename Parent2::DpEntry DpEntry;
  typedef typename Parent2::DpEntryMap DpEntryMap;
  typedef typename Parent2::BfsType BfsType;
  
public:
  TreeSolverRooted(const MwcsGraphType& mwcsGraph, const NodeSet& rootNodes);
  virtual bool solve();
  
  void setRoot(Node root)
  {
    assert(_rootNodes.find(root) != _rootNodes.end());
    _root = root;
  }
};

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline TreeSolverRooted<GR, NWGHT, NLBL, EWGHT>::TreeSolverRooted(const MwcsGraphType& mwcsGraph,
                                                              const NodeSet& rootNodes)
  : GrandParent(mwcsGraph)
  , Parent1(mwcsGraph, rootNodes)
  , Parent2(mwcsGraph)
{
  _root = *_rootNodes.begin();
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline bool TreeSolverRooted<GR, NWGHT, NLBL, EWGHT>::solve()
{
  const WeightNodeMap& weight = _mwcsGraph.getScores();
  
  // work bottom-up
  for (NodeSetVectorRevIt nodeSetIt = _nodesPerLevel.rbegin();
       nodeSetIt != _nodesPerLevel.rend(); nodeSetIt++)
  {
    for (NodeSetIt nodeIt = nodeSetIt->begin(); nodeIt != nodeSetIt->end(); nodeIt++)
    {
      Node node = *nodeIt;
      _dpMap[node]._solution.insert(node);
      _dpMap[node]._weight = weight[node];
      
      const NodeSet& children = _dpMap[node]._children;
      for (NodeSetIt childIt = children.begin(); childIt != children.end(); childIt++)
      {
        double childWeight = _dpMap[*childIt]._weight;
        if (childWeight > 0 || _rootNodes.find(*childIt) != _rootNodes.end())
        {
          _dpMap[node]._weight += childWeight;
          _dpMap[node]._solution.insert(
                                        _dpMap[*childIt]._solution.begin(), _dpMap[*childIt]._solution.end());
        }
      }
    }
  }
  
  // construct the solution
  _score = _dpMap[_root]._weight;
  _solutionSet = _dpMap[_root]._solution;
  for (NodeSetIt nodeIt = _solutionSet.begin(); nodeIt != _solutionSet.end(); nodeIt++)
  {
    _solutionMap[*nodeIt] = true;
  }
  
  return true;
}
  
} // namespace mwcs
} // namespace nina

#endif // TREESOLVERROOTED_H
