/*
 * treesolverrootedimpl.h
 *
 *  Created on: 11-jan-2013
 *      Author: M. El-Kebir
 */

#ifndef TREESOLVERROOTEDIMPL_H
#define TREESOLVERROOTEDIMPL_H

#include "solverrootedimpl.h"
#include "treesolverimpl.h"
#include <set>
#include <map>
#include <limits>
#include <assert.h>

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class TreeSolverRootedImpl : public SolverRootedImpl<GR, NWGHT, NLBL, EWGHT>,
                             public TreeSolverImpl<GR, NWGHT, NLBL, EWGHT>
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;
  
  typedef SolverRootedImpl<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> Parent1;
  typedef TreeSolverImpl<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> Parent2;

  typedef typename Parent1::MwcsGraphType MwcsGraphType;
  typedef typename Parent1::NodeSet NodeSet;
  typedef typename Parent1::NodeSetIt NodeSetIt;
  
  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  
  using Parent1::_pMwcsGraph;
  using Parent1::_rootNodes;
  
  using Parent2::_pDpMap;
  using Parent2::_nodesPerLevel;
  using Parent2::_pPred;
  using Parent2::_pLevel;
  using Parent2::_pBfs;
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
  TreeSolverRootedImpl()
    : Parent1()
    , Parent2()
  {
  }
  
  virtual bool solve(double& score, BoolNodeMap& solutionMap, NodeSet& solutionSet);
  
  virtual void init(const MwcsGraphType& mwcsGraph,
                    const NodeSet& rootNodes)
  {
    Parent1::init(mwcsGraph, rootNodes);
    Parent2::init(mwcsGraph, *_rootNodes.begin());
  }
};


template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline bool TreeSolverRootedImpl<GR, NWGHT, NLBL, EWGHT>::solve(double& score,
                                                                BoolNodeMap& solutionMap,
                                                                NodeSet& solutionSet)
{
  const WeightNodeMap& weight = _pMwcsGraph->getScores();
  
  // work bottom-up
  for (NodeSetVectorRevIt nodeSetIt = _nodesPerLevel.rbegin();
       nodeSetIt != _nodesPerLevel.rend(); nodeSetIt++)
  {
    for (NodeSetIt nodeIt = nodeSetIt->begin(); nodeIt != nodeSetIt->end(); nodeIt++)
    {
      Node node = *nodeIt;
      (*_pDpMap)[node]._solution.insert(node);
      (*_pDpMap)[node]._weight = weight[node];
      
      const NodeSet& children = (*_pDpMap)[node]._children;
      for (NodeSetIt childIt = children.begin(); childIt != children.end(); childIt++)
      {
        double childWeight = (*_pDpMap)[*childIt]._weight;
        if (childWeight > 0 || _rootNodes.find(*childIt) != _rootNodes.end())
        {
          (*_pDpMap)[node]._weight += childWeight;
          (*_pDpMap)[node]._solution.insert(
          (*_pDpMap)[*childIt]._solution.begin(), (*_pDpMap)[*childIt]._solution.end());
        }
      }
    }
  }
  
  // construct the solution
  score = (*_pDpMap)[_root]._weight;
  solutionSet = (*_pDpMap)[_root]._solution;
  for (NodeSetIt nodeIt = solutionSet.begin(); nodeIt != solutionSet.end(); nodeIt++)
  {
    solutionMap[*nodeIt] = true;
  }
  
  return true;
}
  
} // namespace mwcs
} // namespace nina

#endif // TREESOLVERROOTEDIMPL_H
