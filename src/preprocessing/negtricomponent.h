/*
 * negtricomponent.h
 *
 *  Created on: 30-mar-2014
 *      Author: M. El-Kebir
 */

#ifndef NEGTRICOMPONENT_H
#define NEGTRICOMPONENT_H

#include <lemon/adaptors.h>
#include <lemon/core.h>
#include <lemon/connectivity.h>
#include <lemon/dijkstra.h>
#include <string>
#include <vector>
#include <list>
#include <set>
#include <limits>
#include "unrootedrule.h"
#include "solver/blockcuttree.h"
#include "solver/spqrtree.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename WGHT = typename GR::template NodeMap<double> >
class NegTriComponent : public UnrootedRule<GR, WGHT>
{
public:
  typedef GR Graph;
  typedef WGHT WeightNodeMap;
  typedef UnrootedRule<GR, WGHT> Parent;
  typedef typename Parent::NodeMap NodeMap;
  typedef typename Parent::NodeSet NodeSet;
  typedef typename Parent::NodeSetIt NodeSetIt;
  typedef typename Parent::NodeSetMap NodeSetMap;
  typedef typename Parent::DegreeNodeMap DegreeNodeMap;
  typedef typename Parent::DegreeNodeSetVector DegreeNodeSetVector;
  typedef typename Parent::LabelNodeMap LabelNodeMap;
  typedef typename Parent::ArcLookUpType ArcLookUpType;

  TEMPLATE_GRAPH_TYPEDEFS(Graph);

  using Parent::remove;
  using Parent::merge;

  NegTriComponent();
  virtual ~NegTriComponent() {}
  virtual int apply(Graph& g,
                    const ArcLookUpType& arcLookUp,
                    LabelNodeMap& label,
                    WeightNodeMap& score,
                    NodeMap& mapToPre,
                    NodeSetMap& preOrigNodes,
                    NodeSetMap& neighbors,
                    int& nNodes,
                    int& nArcs,
                    int& nEdges,
                    DegreeNodeMap& degree,
                    DegreeNodeSetVector& degreeVector,
                    double& LB);

  virtual std::string name() const { return "NegTriComponent"; }
  
private:
  typedef BlockCutTree<Graph> BlockCutTreeType;
  typedef lemon::SubGraph<Graph> SubGraphType;
  typedef typename nina::SpqrTree<SubGraphType> SpqrType;
  
  typedef typename BlockCutTreeType::Tree BlockT;
  typedef typename std::pair<Node, Node> NodePair;
  typedef typename std::vector<Edge> EdgeVector;
  typedef typename std::list<Node> NodeList;
  typedef typename NodeList::const_iterator NodeListIt;
  typedef typename EdgeVector::const_iterator EdgeVectorIt;
  typedef typename SpqrType::Tree SpqrT;
  typedef typename SpqrT::template EdgeMap<bool> BoolSpqrTreeEdgeMap;
  typedef typename SpqrT::template NodeMap<bool> BoolSpqrTreeNodeMap;
  typedef typename SpqrT::template NodeMap<int> IntSpqrTreeNodeMap;
  typedef typename SpqrT::template NodeMap<NodeSet> NodeSetTreeNodeMap;
  typedef typename SpqrT::template ArcMap<int> IntSpqrTreeArcMap;
  typedef typename lemon::Orienter<const SpqrT> RootedSpqrT;
  typedef typename SpqrT::Arc SpqrTreeArc;
  typedef typename std::vector<typename SpqrT::Node> SpqrTreeNodeVector;
  
  void root(const Graph& g,
            const SpqrType& spqr,
            const BoolSpqrTreeNodeMap& neg,
            const typename SpqrT::Node u,
            const typename SpqrT::Edge uv,
            const typename SpqrT::Node v,
            RootedSpqrT& rootedT,
            BoolSpqrTreeEdgeMap& dir,
            IntSpqrTreeNodeMap& depth,
            IntSpqrTreeNodeMap& size,
            IntSpqrTreeNodeMap& realEdgesSize,
            BoolSpqrTreeNodeMap& negSubTree,
            NodeSetTreeNodeMap& orgNodesInSubTree);
  
  void identifyRoots(const RootedSpqrT& rootedT,
                     const typename SpqrT::Node v,
                     const BoolSpqrTreeNodeMap& negSubTree,
                     SpqrTreeNodeVector& roots);
  
  void shortCircuit(const Graph& g,
                    const ArcLookUpType& arcLookUp,
                    const WeightNodeMap& score,
                    const SpqrType& spqr,
                    const RootedSpqrT& rootedT,
                    const typename SpqrT::Node v,
                    NodeList& path);
  
  double induce(const SpqrType& spqr,
                const RootedSpqrT& rootedT,
                const WeightNodeMap& score,
                const typename SpqrT::Node v,
                const Graph& g,
                BoolNodeMap& presence);
};
  
template<typename GR, typename WGHT>
inline NegTriComponent<GR, WGHT>::NegTriComponent()
  : Parent()
{
}
  
  
template<typename GR, typename WGHT>
inline double NegTriComponent<GR, WGHT>::induce(const SpqrType& spqr,
                                                const RootedSpqrT& rootedT,
                                                const WeightNodeMap& score,
                                                const typename SpqrT::Node v,
                                                const Graph& g,
                                                BoolNodeMap& presence)
{
  double minScore = std::numeric_limits<double>::max();
  
  const EdgeVector& edges = spqr.getRealEdges(v);
  for (EdgeVectorIt edgeIt = edges.begin(); edgeIt != edges.end(); ++edgeIt)
  {
    Node u = g.u(*edgeIt);
    Node w = g.v(*edgeIt);
    presence[u] = presence[w] = true;
    
    if (score[u] < minScore)
    {
      minScore = score[u];
    }
    if (score[w] < minScore)
    {
      minScore = score[w];
    }
  }
  
  for (typename RootedSpqrT::OutArcIt a(rootedT, v); a != lemon::INVALID; ++a)
  {
    const typename SpqrT::Node w = rootedT.target(a);
    minScore = std::min(minScore, induce(spqr, rootedT, score, w, g, presence));
  }
  
  return minScore;
}

template<typename GR, typename WGHT>
inline void NegTriComponent<GR, WGHT>::shortCircuit(const Graph& g,
                                                    const ArcLookUpType& arcLookUp,
                                                    const WeightNodeMap& score,
                                                    const SpqrType& spqr,
                                                    const RootedSpqrT& rootedT,
                                                    const typename SpqrT::Node v,
                                                    NodeList& path)
{
  typedef typename lemon::FilterNodes<const Graph> InducedSubGraph;
  
  // identify cut pair
  typename RootedSpqrT::InArcIt inArc(rootedT, v);
  
  if (inArc == lemon::INVALID)
  {
    // degenerate case: complete biconnected component is negative
    return;
  }
  
  const NodePair& cutPair = spqr.getCutPair(inArc);
  path.clear();
  
  // check: is there an edge between the cut pair?
  if (arcLookUp(cutPair.first, cutPair.second) == lemon::INVALID)
  {
    // find the shortest path
    BoolNodeMap filter(g, false);
    double minScore = induce(spqr, rootedT, score, v, g, filter);
    
    // let's construct the subgraph
    IntArcMap cost(g, 0);
    InducedSubGraph subG(g, filter);
    for (typename InducedSubGraph::ArcIt a(subG); a != lemon::INVALID; ++a)
    {
      cost[a] = -minScore - score[subG.target(a)];
      assert(cost[a] >= 0);
    }
    
    // and now let's do a Dijkstra
    lemon::Dijkstra<InducedSubGraph, IntArcMap> dijkstra(subG, cost);
    dijkstra.run(cutPair.first);
    
    typename InducedSubGraph::Node u = cutPair.second;
    while (dijkstra.predArc(u) != lemon::INVALID)
    {
      u = dijkstra.predNode(u);
      if (u != cutPair.first)
      {
        path.push_front(u);
      }
    }
  }
}
  
template<typename GR, typename WGHT>
inline void NegTriComponent<GR, WGHT>::identifyRoots(const RootedSpqrT& rootedT,
                                                     const typename SpqrT::Node v,
                                                     const BoolSpqrTreeNodeMap& negSubTree,
                                                     SpqrTreeNodeVector& roots)
{
  if (negSubTree[v])
  {
    roots.push_back(v);
  }
  else
  {
    for (typename RootedSpqrT::OutArcIt a(rootedT, v); a != lemon::INVALID; ++a)
    {
      identifyRoots(rootedT, rootedT.target(a), negSubTree, roots);
    }
  }
}

template<typename GR, typename WGHT>
inline void NegTriComponent<GR, WGHT>::root(const Graph& g,
                                            const SpqrType& spqr,
                                            const BoolSpqrTreeNodeMap& neg,
                                            const typename SpqrT::Node u,
                                            const typename SpqrT::Edge uv,
                                            const typename SpqrT::Node v,
                                            RootedSpqrT& rootedT,
                                            BoolSpqrTreeEdgeMap& dir,
                                            IntSpqrTreeNodeMap& depth,
                                            IntSpqrTreeNodeMap& size,
                                            IntSpqrTreeNodeMap& realEdgesSize,
                                            BoolSpqrTreeNodeMap& negSubTree,
                                            NodeSetTreeNodeMap& orgNodesInSubTree)
{
  const SpqrT& T = spqr.getSpqrTree();
  
  if (u == lemon::INVALID)
  {
    depth[v] = 0;
  }
  else
  {
    // direct arc from u to v
    dir[uv] = T.direction(T.direct(uv, u));
    depth[v] = depth[u] + 1;
  }
  
  negSubTree[v] = neg[v];
  size[v] = 1;
  
  const EdgeVector& realEdges = spqr.getRealEdges(v);
  realEdgesSize[v] = static_cast<int>(realEdges.size());
  for (EdgeVectorIt edgeIt = realEdges.begin(); edgeIt != realEdges.end(); ++edgeIt)
  {
    orgNodesInSubTree[v].insert(g.u(*edgeIt));
    orgNodesInSubTree[v].insert(g.v(*edgeIt));
  }
  
  for (typename SpqrT::IncEdgeIt e(T, v); e != lemon::INVALID; ++e)
  {
    typename SpqrT::Node w = T.oppositeNode(v, e);
    if (w != u)
    {
      root(g, spqr, neg, v, e, w, rootedT, dir, depth, size, realEdgesSize, negSubTree, orgNodesInSubTree);
      negSubTree[v] = negSubTree[v] && negSubTree[w];
      size[v] += size[w];
      realEdgesSize[v] += realEdgesSize[w];
      orgNodesInSubTree[v].insert(orgNodesInSubTree[w].begin(), orgNodesInSubTree[w].end());
    }
  }
}
  
template<typename GR, typename WGHT>
inline int NegTriComponent<GR, WGHT>::apply(Graph& g,
                                            const ArcLookUpType& arcLookUp,
                                            LabelNodeMap& label,
                                            WeightNodeMap& score,
                                            NodeMap& mapToPre,
                                            NodeSetMap& preOrigNodes,
                                            NodeSetMap& neighbors,
                                            int& nNodes,
                                            int& nArcs,
                                            int& nEdges,
                                            DegreeNodeMap& degree,
                                            DegreeNodeSetVector& degreeVector,
                                            double& LB)
{
  BlockCutTreeType blockCutTree(g);
  blockCutTree.run();
  
  BoolNodeMap filterNodes(g, false);
  BoolEdgeMap filterEdges(g, false);
  SubGraphType subG(g, filterNodes, filterEdges);
  
  const typename BlockCutTreeType::Tree& blockT = blockCutTree.getBlockCutTree();
  
  int res = 0;
  
  for (typename BlockCutTreeType::Tree::BlueNodeIt block(blockT);
       block != lemon::INVALID; ++block)
  {
    lemon::mapFill(g, filterNodes, false);
    lemon::mapFill(g, filterEdges, false);
    
    // construct subgraph for biconnected component
    const EdgeVector& realBlockEdges = blockCutTree.getRealEdges(block);
    for (EdgeVectorIt edgeIt = realBlockEdges.begin(); edgeIt != realBlockEdges.end(); ++edgeIt)
    {
      filterEdges[*edgeIt] = true;
      filterNodes[g.u(*edgeIt)] = true;
      filterNodes[g.v(*edgeIt)] = true;
    }
    
    // now construct the spqr tree
    SpqrType spqr(subG);
    spqr.run();
    
    const SpqrT& spqrT = spqr.getSpqrTree();
    BoolSpqrTreeNodeMap negative(spqrT, true);
    
    typename SpqrT::Node rootNode = lemon::INVALID;
    for (typename SpqrT::NodeIt v(spqrT); v != lemon::INVALID; ++v)
    {
      const EdgeVector& realSpqrEdges = spqr.getRealEdges(v);
      for (EdgeVectorIt edgeIt = realSpqrEdges.begin(); edgeIt != realSpqrEdges.end(); ++edgeIt)
      {
        negative[v] = negative[v] && (score[g.u(*edgeIt)] <= 0 && score[g.v(*edgeIt)] <= 0);
      }
      
      if (rootNode == lemon::INVALID || spqr.getDegree(v) > spqr.getDegree(rootNode))
      {
        rootNode = v;
      }
    }
    
    // root at max degree node and determine depth
    IntSpqrTreeNodeMap depth(spqrT, -1), size(spqrT, 1), realEdgesSize(spqrT, 0);
    BoolSpqrTreeEdgeMap dir(spqrT, true);
    RootedSpqrT rootedSpqrT(spqrT, dir);
    BoolSpqrTreeNodeMap negativeSubTree(spqrT, false);
    NodeSetTreeNodeMap orgNodesInSubTree(spqrT);
    root(g, spqr, negative, lemon::INVALID, lemon::INVALID,
         rootNode, rootedSpqrT, dir, depth, size,
         realEdgesSize, negativeSubTree, orgNodesInSubTree);
    
    // determine roots of negative subtrees
    SpqrTreeNodeVector negativeRoots;
    identifyRoots(rootedSpqrT, rootNode, negativeSubTree, negativeRoots);
    
    for (typename SpqrTreeNodeVector::const_iterator it = negativeRoots.begin(); it != negativeRoots.end(); ++it)
    {
      if (realEdgesSize[*it] > 2)
      {
        NodeList path;
        shortCircuit(g, arcLookUp, score, spqr, rootedSpqrT, *it, path);
        
        std::cerr << "Node " << spqrT.id(*it)
        << ": size " << size[*it]
        << ": edges size " << realEdgesSize[*it]
        << ", degree " << spqr.getDegree(*it)
        << std::endl;
        
        const NodePair& cutPair = spqr.getCutPair(typename RootedSpqrT::InArcIt(rootedSpqrT, *it));
        
        NodeSet nodesToRemove = orgNodesInSubTree[*it];
        nodesToRemove.erase(cutPair.first);
        nodesToRemove.erase(cutPair.second);
        for (NodeListIt it2 = path.begin(); it2 != path.end(); ++it2)
        {
          nodesToRemove.erase(*it2);
        }
        
        for (NodeSetIt it2 = nodesToRemove.begin(); it2 != nodesToRemove.end(); ++it2)
        {
          remove(g, mapToPre, preOrigNodes, neighbors,
                 nNodes, nArcs, nEdges,
                 degree, degreeVector, *it2);
          ++res;
        }

        bool first = true;
        Node v = cutPair.first;
        for (NodeListIt it2 = path.begin(); it2 != path.end(); ++it2)
        {
          v = merge(g, arcLookUp, label, score,
                    mapToPre, preOrigNodes, neighbors,
                    nNodes, nArcs, nEdges,
                    degree, degreeVector, v, *it2, LB);
          ++res;
          
          if (first) first = false;
          else std::cerr << " -- ";
          std::cerr << g.id(*it2) << " (" << score[*it2] << ")";
        }
        std::cerr << std::endl;
      }
    }
  }
  
  return res;
}

} // namespace mwcs
} // namespace nina

#endif // NEGTRICOMPONENT_H
