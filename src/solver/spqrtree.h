/*
 * spqrtree.h
 *
 *  Created on: 30-mar-2014
 *      Author: M. El-Kebir
 */

#ifndef SPQRTREE_H
#define SPQRTREE_H

#include <vector>
#include <algorithm>
#include <lemon/core.h>
#include <lemon/smart_graph.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/NodeArray.h>
#include <ogdf/basic/EdgeArray.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/decomposition/StaticSPQRTree.h>

namespace nina {
  
template<typename GR>
class SpqrTree
{
public:
  typedef GR Graph;
  typedef lemon::SmartGraph Tree;
  
  typedef enum { SPQR_S, SPQR_P, SPQR_R } SpqrNodeType;
  
  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  
  typedef typename Graph::template NodeMap<ogdf::node> ToOgdfNodeMap;
  typedef typename Graph::template EdgeMap<ogdf::edge> ToOgdfEdgeMap;
  
  typedef ogdf::NodeArray<Node> ToLemonNodeMap;
  typedef ogdf::EdgeArray<Edge> ToLemonEdgeMap;

  typedef std::pair<Node, Node> NodePair;
  typedef std::vector<Edge> EdgeVector;
  
  typedef Tree::template NodeMap<SpqrNodeType> SpqrNodeTypeMap;
  typedef typename Graph::template EdgeMap<Tree::Node> ToSpqrNodeMap;
  typedef Tree::template EdgeMap<NodePair> SpqrCutPairMap;
  typedef Tree::template NodeMap<EdgeVector> SpqrRealEdgesMap;
  
  SpqrTree(const Graph& g)
    : _lemonG(g)
    , _ogdfG()
    , _toOgdfNode(_lemonG)
    , _toLemonNode(_ogdfG)
    , _toOgdfEdge(_lemonG)
    , _toLemonEdge(_ogdfG)
    , _toSpqrNode(_lemonG)
    , _T()
    , _spqrNodeType(_T)
    , _spqrCutPair(_T)
    , _spqrRealEdges(_T)
    , _numSpqrNodes(0)
  {
    // construct graph
    construct();
  }
  
  bool run();
  
  const Tree& getSpqrTree() const
  {
    return _T;
  }
  
  const NodePair& getCutPair(Tree::Edge e) const
  {
    return _spqrCutPair[e];
  }
  
  const EdgeVector& getRealEdges(Tree::Node n) const
  {
    return _spqrRealEdges[n];
  }
  
  Tree::Node toSpqrNode(Edge e) const
  {
    return _toSpqrNode[e];
  }
  
  int getNumSpqrNodes() const
  {
    return _numSpqrNodes;
  }
  
  SpqrNodeType getSpqrNodeType(Tree::Node n) const
  {
    return _spqrNodeType[n];
  }
   
  SpqrNodeType getSpqrNodeType(Edge e) const
  {
    return _spqrNodeType[_toSpqrNode[e]];
  }
  
  static char toChar(SpqrNodeType type)
  {
    switch (type)
    {
      case SPQR_S:
        return 'S';
      case SPQR_P:
        return 'P';
      case SPQR_R:
        return 'R';
    }
  }
  
private:
  void construct();
  
private:
  // input LEMON graph
  const Graph& _lemonG;
  
  // input OGDF graph
  ogdf::Graph _ogdfG;
  
  // mappings
  ToOgdfNodeMap _toOgdfNode;
  ToLemonNodeMap _toLemonNode;
  ToOgdfEdgeMap _toOgdfEdge;
  ToLemonEdgeMap _toLemonEdge;
  ToSpqrNodeMap _toSpqrNode;
  
  // SPQR tree
  Tree _T;
  SpqrNodeTypeMap _spqrNodeType;
  SpqrCutPairMap _spqrCutPair;
  SpqrRealEdgesMap _spqrRealEdges;
  int _numSpqrNodes;
};

template<typename GR>
inline void SpqrTree<GR>::construct()
{
  for (NodeIt v(_lemonG); v != lemon::INVALID; ++v)
  {
    ogdf::node new_vs = _ogdfG.newNode();
    _toOgdfNode[v] = new_vs;
    _toLemonNode[new_vs] = v;
  }
  
  for (EdgeIt e(_lemonG); e != lemon::INVALID; ++e)
  {
    Node u = _lemonG.u(e);
    Node v = _lemonG.v(e);
    ogdf::node new_u = _toOgdfNode[u];
    ogdf::node new_v = _toOgdfNode[v];
    
    ogdf::edge new_e = _ogdfG.newEdge(new_u, new_v);
    _toOgdfEdge[e] = new_e;
    _toLemonEdge[new_e] = e;
  }
}

template<typename GR>
inline bool SpqrTree<GR>::run()
{
  if (!ogdf::isBiconnected(_ogdfG) ||
      _ogdfG.numberOfEdges() <= 2 ||
      !ogdf::isLoopFree(_ogdfG))
  {
    std::cerr << "Graph is not a valid input for SPQR-tree decomposition!" << std::endl;
    return false;
  }
  
  ogdf::StaticSPQRTree spqr(_ogdfG);
  const ogdf::Graph & ogdfT = spqr.tree();

  ogdf::node ogdf_n;
  ogdf::edge ogdf_e;
 
  _numSpqrNodes = ogdfT.numberOfNodes();
  _T.clear();
  _T.reserveNode(_numSpqrNodes);
  _T.reserveEdge(ogdfT.numberOfEdges());
  
  ogdf::NodeArray<Tree::Node> toSpqrLemonNode(ogdfT, lemon::INVALID);
  
  forall_nodes(ogdf_n, ogdfT)
  {
    const ogdf::Skeleton& Sn = spqr.skeleton(ogdf_n);
    const ogdf::Graph& Gn = Sn.getGraph();

    const Tree::Node lemon_n = toSpqrLemonNode[ogdf_n] = _T.addNode();
    
    SpqrNodeType nodeType;
    switch (spqr.typeOf(ogdf_n))
    {
      case ogdf::SPQRTree::SNode:
        nodeType = SPQR_S;
        break;
      case ogdf::SPQRTree::PNode:
        nodeType = SPQR_P;
        break;
      case ogdf::SPQRTree::RNode:
        nodeType = SPQR_R;
        break;
    }
    
    _spqrNodeType[lemon_n] = nodeType;

    // determine real edges
    forall_edges(ogdf_e, Gn)
    {
      if (!Sn.isVirtual(ogdf_e))
      {
        Edge lemon_e = _toLemonEdge[Sn.realEdge(ogdf_e)];
        _toSpqrNode[lemon_e] = lemon_n;
        _spqrRealEdges[lemon_n].push_back(lemon_e);
      }
    }
  }
  
  return true;
}

} // namespace nina

#endif /* SPQRTREE_H */
