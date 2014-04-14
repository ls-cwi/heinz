/*
 * blocktree.h
 *
 *  Created on: 11-apr-2014
 *      Author: M. El-Kebir
 */

#ifndef BLOCKTREE_H
#define BLOCKTREE_H

#include <lemon/core.h>
#include <lemon/smart_graph.h>
#include <lemon/connectivity.h>
#include <set>
#include <vector>
#include <algorithm>
#include <ostream>

namespace nina {

template<typename GR>
class BlockTree
{
public:
  typedef GR Graph;
  typedef lemon::SmartGraph Tree;
  
  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  
  typedef std::vector<Edge> EdgeVector;
  typedef typename EdgeVector::const_iterator EdgeVectorIt;
  
  typedef typename Graph::template EdgeMap<Tree::Node> ToTreeNodeMap;
  typedef Tree::template EdgeMap<Node> ArticulationMap;
  typedef Tree::template NodeMap<EdgeVector> RealEdgesMap;
  
  BlockTree(const Graph& g)
    : _G(g)
    , _toTreeNode(g, lemon::INVALID)
    , _T()
    , _articulationPoint(_T)
    , _realEdges(_T)
    , _numTreeNodes(0)
  {
  }
  
  bool run();
  
  const Tree& getBlockTree() const
  {
    return _T;
  }
  
  const Node getArticulationPoint(Tree::Edge e) const
  {
    return _articulationPoint[e];
  }
  
  const EdgeVector& getRealEdges(Tree::Node n) const
  {
    return _realEdges[n];
  }
  
  Tree::Node toBlockTreeNode(Edge e) const
  {
    return _toTreeNode[e];
  }
  
  int getNumBlockTreeNodes() const
  {
    return _numTreeNodes;
  }
  
  void printNodes(std::ostream& out) const
  {
    for (Tree::NodeIt v(_T); v != lemon::INVALID; ++v)
    {
      out << _T.id(v) << " " << _realEdges[v].size() << std::endl;
    }
  }
  
  void printEdges(std::ostream& out) const
  {
    for (Tree::EdgeIt e(_T); e != lemon::INVALID; ++e)
    {
      out << _T.id(_T.u(e)) << " " << _T.id(_T.v(e)) << std::endl;
    }
  }
  
private:
  const Graph& _G;
  ToTreeNodeMap _toTreeNode;
  Tree _T;
  ArticulationMap _articulationPoint;
  RealEdgesMap _realEdges;
  int _numTreeNodes;
};

template<typename GR>
inline bool BlockTree<GR>::run()
{
  IntEdgeMap biCompMap(_G, -1);
  
  _numTreeNodes = lemon::biNodeConnectedComponents(_G, biCompMap);
  _T.clear();
  _T.reserveNode(_numTreeNodes);
  
  // construct mapping from int to Tree::Node
  std::vector<Tree::Node> treeNodes;
  treeNodes.reserve(_numTreeNodes);
  for (int i = 0; i < _numTreeNodes; ++i)
  {
    treeNodes.push_back(_T.addNode());
  }
  
  typedef std::set<int> IntSet;
  typedef std::set<Node> NodeSet;
  typedef Tree::template NodeMap<NodeSet> NodeSetTreeMap;
  NodeSetTreeMap articulationPoints(_T);
  
  // construct _toTreeNode mapping: from edge to Tree::Node
  // and identify articulation nodes
  for (NodeIt v(_G); v != lemon::INVALID; ++v)
  {
    IntSet biCompIndices;
    for (IncEdgeIt e(_G, v); e != lemon::INVALID; ++e)
    {
      int biCompMap_e = biCompMap[e];
      biCompIndices.insert(biCompMap_e);
      _toTreeNode[e] = treeNodes[biCompMap_e];
      _realEdges[treeNodes[biCompMap_e]].push_back(e);
    }
    
    if (biCompIndices.size() > 1)
    {
      for (typename IntSet::const_iterator it = biCompIndices.begin(); it != biCompIndices.end(); ++it)
      {
        articulationPoints[treeNodes[*it]].insert(v);
      }
    }
  }
  
  // finish constructing _T
  for (Tree::NodeIt u(_T); u != lemon::INVALID; ++u)
  {
    const NodeSet& set_u = articulationPoints[u];
    for (Tree::NodeIt v = u; v != lemon::INVALID; ++v)
    {
      if (u == v) continue;
      const NodeSet& set_v = articulationPoints[v];
  
      // let's compute the intersection
      NodeSet set_uv;
      std::set_intersection(set_u.begin(), set_u.end(),
                            set_v.begin(), set_v.end(),
                            std::insert_iterator<NodeSet>(set_uv, set_uv.begin()));
      
      assert(set_uv.size() == 0 || set_uv.size() == 1);
      if (set_uv.size() == 1)
      {
        _articulationPoint[_T.addEdge(u, v)] = *set_uv.begin();
      }
    }
  }
  
  return true;
}
  
}

#endif /* BLOCKTREE_H */