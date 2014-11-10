/*
 * treesolverimpl.h
 *
 *  Created on: 11-jan-2013
 *      Author: M. El-Kebir
 */

#ifndef TREESOLVERIMPL_H
#define TREESOLVERIMPL_H

#include "../solver.h"
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
class TreeSolverImpl
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;

  typedef MwcsGraph<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> MwcsGraphType;
  
  TEMPLATE_GRAPH_TYPEDEFS(Graph);

protected:
  typedef typename lemon::Bfs<Graph>::PredMap PredMap;
  typedef typename lemon::Bfs<Graph>::DistMap DistMap;

  typedef std::set<Node> NodeSet;
  typedef typename NodeSet::const_iterator NodeSetIt;
  
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

protected:
  TreeSolverImpl()
    : _pDpMap(NULL)
    , _nodesPerLevel()
    , _pPred(NULL)
    , _pLevel(NULL)
    , _pBfs(NULL)
    , _root(lemon::INVALID)
  {
  }
  
  virtual ~TreeSolverImpl()
  {
    delete _pDpMap;
    delete _pPred;
    delete _pLevel;
    delete _pBfs;
  }

public:
  void init(const MwcsGraphType& mwcsGraph, Node root);
 
  virtual bool solve(double& score,
                     double& scoreUB,
                     BoolNodeMap& solutionMap,
                     NodeSet& solutionSet) = 0;

protected:
  DpEntryMap* _pDpMap;
  NodeSetVector _nodesPerLevel;
  PredMap* _pPred;
  DistMap* _pLevel;
  BfsType* _pBfs;
  Node _root;
};

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void TreeSolverImpl<GR, NWGHT, NLBL, EWGHT>::init(const MwcsGraphType& mwcsGraph,
                                                         Node root)
{
  const Graph& g = mwcsGraph.getGraph();

  // initialize data structures
  delete _pDpMap;
  _nodesPerLevel.clear();
  delete _pPred;
  delete _pLevel;
  delete _pBfs;
  
  _pDpMap = new DpEntryMap(g);
  _pPred = new PredMap(g, lemon::INVALID);
  _pLevel = new DistMap(g);
  _pBfs = new BfsType(g);
  _pBfs->predMap(*_pPred);
  _pBfs->distMap(*_pLevel);
  _root = root;
  
  assert(lemon::acyclic(g));
  assert(_root != lemon::INVALID);
  
  // start by doing a bfs to determine levels
  _pBfs->init();
  _pBfs->addSource(_root);

  while (!_pBfs->emptyQueue())
  {
    Node node = _pBfs->processNextNode();
    int bfsLevel = (*_pLevel)[node];

    if (bfsLevel == static_cast<int>(_nodesPerLevel.size()))
      _nodesPerLevel.push_back(NodeSet());

    _nodesPerLevel[bfsLevel].insert(node);

    Node parent = (*_pPred)[node] != lemon::INVALID ?
          g.oppositeNode(node, (*_pPred)[node]) : lemon::INVALID;

    if (parent != lemon::INVALID)
    {
      //std::cout << "Child: " << _mwcsGraph.getLabel(node)
      //          << " parent: " << _mwcsGraph.getLabel(parent) << std::endl;
      (*_pDpMap)[parent]._children.insert(node);
    }
  }
}

} // namespace mwcs
} // namespace nina

#endif // TREESOLVERIMPL_H
