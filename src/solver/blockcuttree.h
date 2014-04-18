/*
 * blockcuttree.h
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
class BlockCutTree
{
public:
  typedef GR Graph;
  typedef lemon::SmartBpGraph Tree;
  
  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  
  typedef std::vector<Edge> EdgeVector;
  typedef typename EdgeVector::const_iterator EdgeVectorIt;
  
  typedef typename Graph::template EdgeMap<Tree::BlueNode> ToBlockNodeMap;
  typedef Tree::template RedNodeMap<Node> ArticulationMap;
  typedef Tree::template BlueNodeMap<EdgeVector> RealEdgesMap;
  typedef Tree::template NodeMap<int> IntTreeNodeMap;
  
  BlockCutTree(const Graph& g)
    : _G(g)
    , _toBlockNode(g, lemon::INVALID)
    , _T()
    , _articulationPoint(_T)
    , _realEdges(_T)
    , _deg(_T)
    , _numTreeNodes(0)
  {
  }
  
  bool run();
  
  const Tree& getBlockCutTree() const
  {
    return _T;
  }
  
  const Node getArticulationPoint(Tree::RedNode n) const
  {
    return _articulationPoint[n];
  }
  
  const EdgeVector& getRealEdges(Tree::BlueNode n) const
  {
    return _realEdges[n];
  }
  
  Tree::Node toBlockNode(Edge e) const
  {
    return _toBlockNode[e];
  }
  
  int getNumBlockTreeNodes() const
  {
    return _numTreeNodes;
  }
  
  int getDegree(Tree::Node n) const
  {
    return _deg[n];
  }
  
  void printNodes(std::ostream& out) const
  {
    for (Tree::NodeIt v(_T); v != lemon::INVALID; ++v)
    {
      if (_T.blue(v))
      {
        out << "Block node: " << _T.id(v) << ", deg: " << _deg[v]
            << ", #edges: " << _realEdges[_T.asBlueNode(v)].size() << std::endl;
      }
      else
      {
        out << "Cut node: " << _T.id(v) << ", deg: " << _deg[v] << std::endl;
      }
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
  ToBlockNodeMap _toBlockNode;
  Tree _T;
  ArticulationMap _articulationPoint;
  RealEdgesMap _realEdges;
  IntTreeNodeMap _deg;
  int _numTreeNodes;
};

template<typename GR>
inline bool BlockCutTree<GR>::run()
{
  IntEdgeMap biCompMap(_G, -1);
  
  _numTreeNodes = lemon::biNodeConnectedComponents(_G, biCompMap);
  _T.clear();
  _T.reserveNode(2 * _numTreeNodes);
  
  // construct mapping from int to Tree::Node
  std::vector<Tree::BlueNode> blockNodes;
  blockNodes.reserve(_numTreeNodes);
  for (int i = 0; i < _numTreeNodes; ++i)
  {
    blockNodes.push_back(_T.addBlueNode());
    _deg[blockNodes.back()] = 0;
  }
  
  typedef std::set<int> IntSet;
  typedef std::set<Node> NodeSet;
  typedef Tree::template NodeMap<NodeSet> NodeSetTreeMap;
  typedef typename Graph::template NodeMap<Tree::RedNode> RevArticulationMap;
  
  NodeSetTreeMap articulationPoints(_T);
  RevArticulationMap revArticulationMap(_G, lemon::INVALID);
  
  // construct _toTreeNode mapping: from edge to Tree::BlueNode
  // and identify articulation nodes
  for (NodeIt v(_G); v != lemon::INVALID; ++v)
  {
    IntSet biCompIndices;
    for (IncEdgeIt e(_G, v); e != lemon::INVALID; ++e)
    {
      int biCompMap_e = biCompMap[e];
      biCompIndices.insert(biCompMap_e);
      _toBlockNode[e] = blockNodes[biCompMap_e];
      _realEdges[blockNodes[biCompMap_e]].push_back(e);
    }
    
    if (biCompIndices.size() > 1)
    {
      typename Tree::RedNode cut_v = revArticulationMap[v];
      if (cut_v == lemon::INVALID)
      {
        cut_v = revArticulationMap[v] = _T.addRedNode();
        _deg[cut_v] = 0;
        _articulationPoint[cut_v] = v;
      }
      
      for (typename IntSet::const_iterator it = biCompIndices.begin(); it != biCompIndices.end(); ++it)
      {
        _T.addEdge(cut_v, blockNodes[*it]);
        ++_deg[cut_v];
        ++_deg[blockNodes[*it]];
      }
    }
  }
  
  return true;
}
  
}

#endif /* BLOCKCUTTREE_H */