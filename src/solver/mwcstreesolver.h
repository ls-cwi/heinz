/*
 * mwcstreesolver.h
 *
 *  Created on: 11-jan-2013
 *      Author: M. El-Kebir
 */

#ifndef MWCSTREESOLVER_H
#define MWCSTREESOLVER_H

#include "mwcssolver.h"
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
class MwcsTreeSolver : public MwcsSolver<GR, NWGHT, NLBL, EWGHT>
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;

  typedef MwcsSolver<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> Parent;
  typedef typename Parent::MwcsGraphType MwcsGraphType;

  TEMPLATE_GRAPH_TYPEDEFS(Graph);

  using Parent::_mwcsGraph;
  using Parent::_root;
  using Parent::_score;
  using Parent::_solutionMap;
  using Parent::_solutionSet;
  using Parent::printSolution;
  using Parent::getSolutionWeight;
  using Parent::getSolutionNodeMap;
  using Parent::getSolutionModule;
  using Parent::isNodeInSolution;
  using Parent::init;

private:
  typedef std::set<Node> NodeSet;
  typedef typename NodeSet::const_iterator NodeSetIt;

  typedef typename Graph::template NodeMap<NodeSet> NodeSetMap;
  typedef typename lemon::Bfs<Graph>::PredMap PredMap;
  typedef typename lemon::Bfs<Graph>::DistMap DistMap;

  typedef std::vector<NodeSet> NodeSetVector;
  typedef typename NodeSetVector::const_iterator NodeSetVectorIt;
  typedef typename NodeSetVector::const_reverse_iterator NodeSetVectorRevIt;

  struct DpEntry {
    double _weight;
    NodeSet _solution;
    NodeSet _children;

    DpEntry()
      : _weight(-std::numeric_limits<double>::max())
      , _solution()
      , _children()
    {
    }

    void clear()
    {
      _weight = -std::numeric_limits<double>::max();
      _solution.clear();
      _children.clear();
    }
  };
  typedef typename Graph::template NodeMap<DpEntry> DpEntryMap;
  typedef lemon::Bfs<Graph> BfsType;

public:
  MwcsTreeSolver(const MwcsGraphType& mwcsGraph);
  virtual void init(Node root);
  virtual bool solve();
  virtual void setLowerBound(double LB) {}

private:
  DpEntryMap _dpMap;
  NodeSetVector _nodesPerLevel;
  PredMap _pred;
  DistMap _level;
  BfsType _bfs;
};

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline MwcsTreeSolver<GR, NWGHT, NLBL, EWGHT>::MwcsTreeSolver(const MwcsGraphType& mwcsGraph)
  : Parent(mwcsGraph)
  , _dpMap(mwcsGraph.getGraph())
  , _nodesPerLevel()
  , _pred(mwcsGraph.getGraph())
  , _level(mwcsGraph.getGraph())
  , _bfs(mwcsGraph.getGraph())
{
  _bfs.predMap(_pred);
  _bfs.distMap(_level);
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsTreeSolver<GR, NWGHT, NLBL, EWGHT>::init(Node root)
{
  const Graph& g = _mwcsGraph.getGraph();

  _nodesPerLevel.clear();
  _score = 0;
  _solutionSet.clear();
  for (NodeIt v(g); v != lemon::INVALID; ++v)
  {
    _solutionMap[v] = false;
    _dpMap[v].clear();
  }

  _root = root;

  //assert(lemon::acyclic(g));
  assert(root != lemon::INVALID);

  // start by doing a bfs to determine levels
  _bfs.init();
  _bfs.addSource(root);

  while (!_bfs.emptyQueue())
  {
    Node node = _bfs.processNextNode();
    int bfsLevel = _level[node];

    if (bfsLevel == static_cast<int>(_nodesPerLevel.size()))
      _nodesPerLevel.push_back(NodeSet());

    _nodesPerLevel[bfsLevel].insert(node);

    Node parent = _pred[node] != lemon::INVALID ?
          g.oppositeNode(node, _pred[node]) : lemon::INVALID;

    if (parent != lemon::INVALID)
    {
      //std::cout << "Child: " << _mwcsGraph.getLabel(node)
      //          << " parent: " << _mwcsGraph.getLabel(parent) << std::endl;
      _dpMap[parent]._children.insert(node);
    }
  }
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline bool MwcsTreeSolver<GR, NWGHT, NLBL, EWGHT>::solve()
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
        if (childWeight > 0)
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

#endif // MWCSTREESOLVER_H
