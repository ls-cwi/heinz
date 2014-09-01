/*
 * blockcuttree.h
 *
 *  Created on: 11-apr-2014
 *      Author: M. El-Kebir
 */

#ifndef BLOCKTREE_H
#define BLOCKTREE_H

#include <lemon/core.h>
#include <lemon/list_graph.h>
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
  typedef lemon::ListBpGraph Tree;
  
  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  
  typedef std::vector<Edge> EdgeVector;
  typedef typename EdgeVector::const_iterator EdgeVectorIt;
  
  typedef std::set<Node> NodeSet;
  typedef typename NodeSet::const_iterator NodeSetIt;
  
  typedef Tree::BlueNode BlockNode;
  typedef Tree::BlueNodeIt BlockNodeIt;
  typedef Tree::RedNode CutNode;
  typedef Tree::RedNodeIt CutNodeIt;
  
  typedef typename Graph::template EdgeMap<BlockNode> ToBlockNodeMap;
  typedef Tree::template RedNodeMap<Node> ArticulationMap;
  typedef Tree::template BlueNodeMap<EdgeVector> RealEdgesMap;
  typedef Tree::template NodeMap<int> IntTreeNodeMap;
  typedef Tree::template BlueNodeMap<bool> BoolTreeBlockNodeMap;
  
  typedef std::set<BlockNode> BlockNodeSet;
  typedef typename BlockNodeSet::const_iterator BlockNodeSetIt;
  typedef std::vector<BlockNodeSet> BlockNodeSetVector;
  typedef typename BlockNodeSetVector::const_iterator BlockNodeSetVectorIt;
  
  BlockCutTree(const Graph& g)
    : _G(g)
    , _toBlockNode(g, lemon::INVALID)
    , _T()
    , _articulationPoint(_T)
    , _realEdges(_T)
    , _deg(_T)
    , _numBlockNodes(0)
    , _blockNodesByDegree()
  {
  }
  
  bool run();
  
  bool removeBlockNode(Tree::BlueNode n)
  {
    if (getDegree(n) == 0)
    {
      _blockNodesByDegree[0].erase(n);
      _T.erase(n);
      --_numBlockNodes;
      
      return true;
    }
    else if (getDegree(n) == 1)
    {
      CutNode cut = _T.redNode(Tree::IncEdgeIt(_T, n));
      
      _blockNodesByDegree[1].erase(n);
      _T.erase(n);
      --_numBlockNodes;
      
      // remove cut node if its degree becomes 1
      if (--_deg[cut] == 1)
      {
        BlockNode b2 = _T.blueNode(Tree::IncEdgeIt(_T, cut));
        _blockNodesByDegree[_deg[b2]--].erase(b2);
        _blockNodesByDegree[_deg[b2]].insert(b2);
        _T.erase(cut);
      }
      
      return true;
    }
    else
    {
      return false;
    }
  }
  
  const Tree& getBlockCutTree() const
  {
    return _T;
  }
  
  const Node getArticulationPoint(CutNode n) const
  {
    return _articulationPoint[n];
  }
  
  const EdgeVector& getRealEdges(BlockNode n) const
  {
    return _realEdges[n];
  }
  
  void getRealNodes(BlockNode n, NodeSet& nodeSet) const
  {
    const EdgeVector& realEdges = getRealEdges(n);
    for (EdgeVectorIt edgeIt = realEdges.begin(); edgeIt != realEdges.end(); ++edgeIt)
    {
      nodeSet.insert(_G.u(*edgeIt));
      nodeSet.insert(_G.v(*edgeIt));
    }
  }
  
  Tree::Node toBlockNode(Edge e) const
  {
    return _toBlockNode[e];
  }
  
  int getNumBlockTreeNodes() const
  {
    return _numBlockNodes;
  }
  
  int getDegree(Tree::Node n) const
  {
    return _deg[n];
  }
  
  const BlockNodeSet& getBlockNodeSetByDegree(int deg) const
  {
    assert(0 <= deg && static_cast<size_t>(deg) < _blockNodesByDegree.size());
    return _blockNodesByDegree[deg];
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
  int _numBlockNodes;
  BlockNodeSetVector _blockNodesByDegree;
};

template<typename GR>
inline bool BlockCutTree<GR>::run()
{
  IntEdgeMap biCompMap(_G, -1);
  
  _numBlockNodes = lemon::biNodeConnectedComponents(_G, biCompMap);
  _T.clear();
  _T.reserveNode(2 * _numBlockNodes);
  
  // construct mapping from int to Tree::Node
  std::vector<BlockNode> blockNodes;
  blockNodes.reserve(_numBlockNodes);
  for (int i = 0; i < _numBlockNodes; ++i)
  {
    blockNodes.push_back(_T.addBlueNode());
    _deg[blockNodes.back()] = 0;
  }
  
  typedef std::set<int> IntSet;
  typedef std::set<Node> NodeSet;
  typedef Tree::template NodeMap<NodeSet> NodeSetTreeMap;
  typedef typename Graph::template NodeMap<CutNode> RevArticulationMap;
  
  NodeSetTreeMap articulationPoints(_T);
  RevArticulationMap revArticulationMap(_G, lemon::INVALID);
  
  // construct _toTreeNode mapping: from edge to BlockNode
  // and identify articulation nodes
  for (NodeIt v(_G); v != lemon::INVALID; ++v)
  {
    IntSet biCompIndices;
    for (IncEdgeIt e(_G, v); e != lemon::INVALID; ++e)
    {
      int biCompMap_e = biCompMap[e];
      biCompIndices.insert(biCompMap_e);
      _toBlockNode[e] = blockNodes[biCompMap_e];
    }
    
    if (biCompIndices.size() > 1)
    {
      CutNode cut_v = revArticulationMap[v];
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
  
  for (EdgeIt e(_G); e != lemon::INVALID; ++e)
  {
    int biCompMap_e = biCompMap[e];
    _realEdges[blockNodes[biCompMap_e]].push_back(e);
  }
  
  // determine maxBlockDegree
  int maxBlockDegree = -1;
  for (BlockNodeIt b(_T); b != lemon::INVALID; ++b)
  {
    maxBlockDegree = std::max(_deg[b], maxBlockDegree);
  }
  
  // construct _blockNodeByDegree
  _blockNodesByDegree = BlockNodeSetVector(maxBlockDegree + 1);
  for (BlockNodeIt b(_T); b != lemon::INVALID; ++b)
  {
    _blockNodesByDegree[_deg[b]].insert(b);
  }
  
  return true;
}
  
}

#endif /* BLOCKCUTTREE_H */