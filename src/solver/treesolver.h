/*
 * treesolver.h
 *
 *  Created on: 11-jan-2013
 *      Author: M. El-Kebir
 */

#ifndef TREESOLVER_H
#define TREESOLVER_H

#include "solver.h"
#include <set>
#include <vector>
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
class TreeSolver : public virtual Solver<GR, NWGHT, NLBL, EWGHT>
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;

  typedef Solver<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> Parent;
  typedef typename Parent::MwcsGraphType MwcsGraphType;
  
  typedef typename Parent::NodeSet NodeSet;
  typedef typename Parent::NodeSetIt NodeSetIt;
  typedef typename Parent::NodeSetNonConstIt NodeSetNonConstIt;
  typedef typename Parent::NodeVector NodeVector;
  typedef typename Parent::NodeVectorIt NodeVectorIt;
  typedef typename Parent::NodeVectorNonConstIt NodeVectorNonConstIt;

  TEMPLATE_GRAPH_TYPEDEFS(Graph);

  using Parent::_mwcsGraph;
  using Parent::_score;
  using Parent::_solutionMap;
  using Parent::_solutionSet;
  using Parent::printSolution;
  using Parent::getSolutionWeight;
  using Parent::getSolutionNodeMap;
  using Parent::getSolutionModule;
  using Parent::isNodeInSolution;
  using Parent::init;

protected:
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
  TreeSolver(const MwcsGraphType& mwcsGraph);
  virtual void init();

  virtual void setRoot(Node root)
  {
    assert(root != lemon::INVALID);
    _root = root;
  }

protected:
  DpEntryMap _dpMap;
  NodeSetVector _nodesPerLevel;
  PredMap _pred;
  DistMap _level;
  BfsType _bfs;
  Node _root;
};

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline TreeSolver<GR, NWGHT, NLBL, EWGHT>::TreeSolver(const MwcsGraphType& mwcsGraph)
  : Parent(mwcsGraph)
  , _dpMap(mwcsGraph.getGraph())
  , _nodesPerLevel()
  , _pred(mwcsGraph.getGraph())
  , _level(mwcsGraph.getGraph())
  , _bfs(mwcsGraph.getGraph())
  , _root(lemon::INVALID)
{
  _bfs.predMap(_pred);
  _bfs.distMap(_level);
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void TreeSolver<GR, NWGHT, NLBL, EWGHT>::init()
{
  const Graph& g = _mwcsGraph.getGraph();
  assert(lemon::acyclic(g));
  assert(_root != lemon::INVALID);

  // initialize data structures
  _nodesPerLevel.clear();
  _score = 0;
  _solutionSet.clear();
  for (NodeIt v(g); v != lemon::INVALID; ++v)
  {
    _solutionMap[v] = false;
    _dpMap[v].clear();
  }
  
  // start by doing a bfs to determine levels
  _bfs.init();
  _bfs.addSource(_root);

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

} // namespace mwcs
} // namespace nina

#endif // TREESOLVER_H
