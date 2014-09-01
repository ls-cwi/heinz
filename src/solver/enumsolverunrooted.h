/*
 * enumsolverunrooted.h
 *
 *  Created on: 30-jan-2013
 *      Author: M. El-Kebir
 */

#ifndef ENUMSOLVERUNROOTED_H
#define ENUMSOLVERUNROOTED_H

#include <algorithm>
#include <set>
#include <vector>
#include <list>
#include <assert.h>
#include <ostream>

#include "mwcs.h"
#include "mwcsgraph.h"
#include "mwcspreprocessedgraph.h"

#include "solver/solverunrooted.h"
#include "solver/impl/solverrootedimpl.h"

#include "blockcuttree.h"
#include "spqrtree.h"

#include <lemon/adaptors.h>

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class EnumSolverUnrooted : public SolverUnrooted<GR, NWGHT, NLBL, EWGHT>
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;
  
  typedef SolverUnrooted<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> Parent;
  typedef typename Parent::SolverUnrootedImplType SolverUnrootedImplType;
  typedef typename Parent::MwcsGraphType MwcsGraphType;
  typedef typename Parent::NodeSet NodeSet;
  typedef typename Parent::NodeSetIt NodeSetIt;
  typedef typename Parent::NodeVector NodeVector;
  typedef typename Parent::NodeVectorIt NodeVectorIt;

  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  
  typedef SolverRootedImpl<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> SolverRootedImplType;
  
  using Parent::_score;
  using Parent::_pSolutionMap;
  using Parent::_solutionSet;
  using Parent::_pImpl;

public:
  EnumSolverUnrooted(SolverUnrootedImplType* pUnrootedImpl,
                     SolverRootedImplType* pRootedImpl,
                     bool preprocess)
    : Parent(pUnrootedImpl)
    , _pRootedImpl(pRootedImpl)
    , _preprocess(preprocess)
  {
  }
  
  ~EnumSolverUnrooted()
  {
    delete _pRootedImpl;
  }
  
  bool solve(const MwcsGraphType& mwcsGraph);

protected:
  typedef MwcsPreprocessedGraph<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> MwcsPreGraphType;
  typedef lemon::FilterNodes<const Graph, BoolNodeMap> SubGraph;
  typedef typename SubGraph::Node SubNode;
  typedef typename SubGraph::NodeIt SubNodeIt;
  typedef typename SubGraph::ArcIt SubArcIt;
  typedef typename Graph::template NodeMap<Node> NodeMap;
  
  typedef std::vector<Edge> EdgeVector;
  typedef typename EdgeVector::const_iterator EdgeVectorIt;
  
  typedef BlockCutTree<Graph> BlockCutTreeType;
  typedef typename BlockCutTreeType::Tree BcTree;
  typedef typename BlockCutTreeType::BoolTreeBlockNodeMap BcTreeBoolBlockNodeMap;
  typedef typename BlockCutTreeType::BlockNode BcTreeBlockNode;
  typedef typename BlockCutTreeType::BlockNodeIt BcTreeBlockNodeIt;
  typedef typename BlockCutTreeType::CutNode BcTreeCutNode;
  typedef typename BlockCutTreeType::CutNodeIt BcTreeCutNodeIt;
  typedef typename BlockCutTreeType::BlockNodeSet BcTreeBlockNodeSet;
  typedef typename BlockCutTreeType::BlockNodeSetIt BcTreeBlockNodeSetIt;
  
  typedef typename nina::SpqrTree<SubGraph> SpqrType;
  typedef typename SpqrType::Tree SpqrTree;
  typedef typename SpqrTree::template NodeMap<bool> SpqrTreeBoolNodeMap;
  typedef typename SpqrTree::template NodeMap<int> SpqrTreeIntNodeMap;
  typedef typename SpqrTree::template NodeMap<NodeSet> SpqrTreeNodeSetNodeMap;
  typedef typename SpqrTree::template EdgeMap<bool> SpqrTreeBoolEdgeMap;
  typedef typename SpqrTree::Node SpqrTreeNode;
  typedef typename SpqrTree::NodeIt SpqrTreeNodeIt;
  typedef typename SpqrTree::Edge SpqrTreeEdge;
  typedef typename SpqrTree::EdgeIt SpqrTreeEdgeIt;
  typedef typename SpqrTree::IncEdgeIt SpqrTreeIncEdgeIt;
  typedef typename lemon::Orienter<const SpqrTree> RootedSpqrTree;
  typedef typename RootedSpqrTree::InArcIt RootedSpqrTreeInArcIt;
  typedef typename RootedSpqrTree::OutArcIt RootedSpqrTreeOutArcIt;
  typedef std::vector<SpqrTreeNode> SpqrTreeNodeVector;
  typedef typename SpqrTreeNodeVector::const_iterator SpqrTreeNodeVectorIt;
  typedef std::pair<Node, Node> NodePair;
  typedef std::list<Node> NodeList;
  typedef typename NodeList::const_iterator NodeListIt;
  
private:
  SolverRootedImplType* _pRootedImpl;
  bool _preprocess;
  
  bool solveComponent(MwcsPreGraphType& mwcsGraph,
                      NodeSet& solutionSet,
                      double& solutionScore);
  
  bool processBlock(MwcsPreGraphType& mwcsGraph,
                    BlockCutTreeType& bcTree,
                    BcTreeBlockNode b,
                    BoolNodeMap& sameBlock);
  
  bool solveBlock(MwcsPreGraphType& mwcsGraph,
                  BlockCutTreeType& bcTree,
                  const BcTreeBoolBlockNodeMap& bcTreeNeg,
                  BcTreeBlockNode b,
                  int blockIndex,
                  int nBlocks);
  
  bool solveBlockUnrooted(MwcsPreGraphType& mwcsGraph,
                          NodeSet& solutionSet,
                          double& solutionScore);

  bool solveBlockRooted(MwcsPreGraphType& mwcsGraph,
                        Node rootNode,
                        NodeSet& solutionSet,
                        double& solutionScore);
  
  void map(const MwcsGraphType& mwcsGraph,
           const NodeMap& m,
           const NodeSet& source,
           NodeSet& target) const
  {
    for (NodeSetIt nodeIt = source.begin(); nodeIt != source.end(); nodeIt++)
    {
      const NodeSet& orgNodes = mwcsGraph.getOrgNodes(*nodeIt);
      for (NodeSetIt orgNodeIt = orgNodes.begin(); orgNodeIt != orgNodes.end(); orgNodeIt++)
      {
        target.insert(m[*orgNodeIt]);
      }
    }
  }
  
  void initLocalGraph(const Graph& g,
                      const WeightNodeMap& scoreG,
                      const LabelNodeMap& labelG,
                      BoolNodeMap& filterG,
                      Graph& subG,
                      DoubleNodeMap& weightSubG,
                      LabelNodeMap& labelSubG,
                      NodeMap& mapToG,
                      MwcsPreGraphType& mwcsSubGraph)
  {
    SubGraph subTmpSameCompG(g, filterG);
    lemon::graphCopy(subTmpSameCompG, subG)
      .nodeMap(scoreG, weightSubG)
      .nodeMap(labelG, labelSubG)
      .nodeCrossRef(mapToG)
      .run();
    
    mwcsSubGraph.init(&subG, &labelSubG, &weightSubG, NULL);
  }
  
  void initLocalGraph(const Graph& g,
                      const WeightNodeMap& scoreG,
                      const LabelNodeMap& labelG,
                      BoolNodeMap& filterG,
                      Graph& subG,
                      DoubleNodeMap& weightSubG,
                      LabelNodeMap& labelSubG,
                      NodeMap& mapToG,
                      NodeMap& mapToSubG,
                      MwcsPreGraphType& mwcsSubGraph)
  {
    SubGraph subTmpSameCompG(g, filterG);
    lemon::graphCopy(subTmpSameCompG, subG)
      .nodeMap(scoreG, weightSubG)
      .nodeMap(labelG, labelSubG)
      .nodeCrossRef(mapToG)
      .nodeRef(mapToSubG)
      .run();
    
    mwcsSubGraph.init(&subG, &labelSubG, &weightSubG, NULL);
  }
  
  void printNodeSet(const MwcsGraphType& mwcsGraph,
                    const NodeSet& nodeSet) const
  {
    for (NodeSetIt nodeIt = nodeSet.begin(); nodeIt != nodeSet.end(); ++nodeIt)
    {
      Node v = *nodeIt;
      std::cout << mwcsGraph.getLabel(v) << " : " << mwcsGraph.getScore(v) << std::endl;
    }
  }
  
  void rootSpqrTree(const SubGraph& g,
                    const WeightNodeMap& score,
                    const SpqrType& spqr,
                    const SpqrTreeBoolNodeMap& neg,
                    const SpqrTreeNode u,
                    const SpqrTreeEdge uv,
                    const SpqrTreeNode v,
                    RootedSpqrTree& rootedT,
                    SpqrTreeBoolEdgeMap& dir,
                    SpqrTreeIntNodeMap& depth,
                    SpqrTreeIntNodeMap& size,
                    SpqrTreeIntNodeMap& realEdgesSize,
                    SpqrTreeBoolNodeMap& negSubTree,
                    SpqrTreeNodeSetNodeMap& orgNodesInSubTree,
                    SpqrTreeNodeSetNodeMap& orgPosNodesInSubTree);
  
  void identifyNegativeRoots(const SpqrType& spqr,
                             const RootedSpqrTree& rootedT,
                             const SpqrTreeNode v,
                             const SpqrTreeNodeSetNodeMap& orgPosNodesInSubTree,
                             SpqrTreeNodeVector& roots);
  
  double shortCircuit(const Graph& orgG,
                      const SubGraph& g,
                      const WeightNodeMap& score,
                      const SpqrType& spqr,
                      const RootedSpqrTree& rootedT,
                      const SpqrTreeNode v,
                      const SpqrTreeNodeSetNodeMap& orgNodesInSubTree,
                      const bool reverse,
                      BoolNodeMap& filter,
                      DoubleArcMap& arcCost,
                      NodeList& path);

};

template<typename GR, typename WGHT, typename NLBL, typename EWGHT>
inline bool EnumSolverUnrooted<GR, WGHT, NLBL, EWGHT>::solve(const MwcsGraphType& mwcsGraph)
{
  const Graph& g = mwcsGraph.getGraph();
  BoolNodeMap allowedNodesSameComp(g);
  
  // 1. iterate over the components
  int nComponents = mwcsGraph.getComponentCount();
  const IntNodeMap& comp = mwcsGraph.getComponentMap();
  
  for (int compIdx = 0; compIdx < nComponents; ++compIdx)
  {
    lemon::mapFill(g, allowedNodesSameComp, false);
    for (NodeIt node(g); node != lemon::INVALID; ++node)
    {
      allowedNodesSameComp[node] = comp[node] == compIdx;
    }
    
    Graph subG;
    DoubleNodeMap weightSubG(subG);
    LabelNodeMap labelSubG(subG);
    NodeMap mapToG(subG);
    MwcsPreGraphType mwcsSubGraph;
    
    // 2b. create and preprocess subgraph
    initLocalGraph(g,
                   mwcsGraph.getScores(),
                   mwcsGraph.getLabels(),
                   allowedNodesSameComp,
                   subG,
                   weightSubG,
                   labelSubG,
                   mapToG,
                   mwcsSubGraph);
    
    if (g_verbosity >= VERBOSE_ESSENTIAL)
    {
      std::cout << std::endl;
      std::cout << "// Considering component " << compIdx + 1 << "/" << nComponents
                << ": contains " << mwcsSubGraph.getNodeCount() << " nodes and "
                << mwcsSubGraph.getEdgeCount() << " edges" << std::endl;
    }
    
    // 3. solve
    double solutionScore;
    NodeSet solutionSet;
    if (!solveComponent(mwcsSubGraph, solutionSet, solutionScore))
    {
      return false;
    }
    
    if (solutionScore > _score)
    {
      _score = solutionScore;
      _solutionSet.clear();
      map(mwcsSubGraph, mapToG, solutionSet, _solutionSet);
    }
  }
  
  _pSolutionMap = new BoolNodeMap(g, false);
  if (_solutionSet.size() > 0)
  {
    for (NodeSetIt nodeIt = _solutionSet.begin();
         nodeIt != _solutionSet.end(); ++nodeIt)
    {
      _pSolutionMap->set(*nodeIt, true);
    }
  }

  return _solutionSet.size() > 0;
}
  
template<typename GR, typename WGHT, typename NLBL, typename EWGHT>
inline bool EnumSolverUnrooted<GR, WGHT, NLBL, EWGHT>::solveComponent(MwcsPreGraphType& mwcsGraph,
                                                                      NodeSet& solutionSet,
                                                                      double& solutionScore)
{
  // mwcsGraph corresponds to a single component
  const Graph& g = mwcsGraph.getGraph();
  const WeightNodeMap& score = mwcsGraph.getScores();
  
  if (_preprocess)
  {
    // preprocess the graph
    mwcsGraph.preprocess(NodeSet());
  }
  
  // generate block-cut vertex tree
  BlockCutTreeType bcTree(g);
  bcTree.run();
  
  // determine block negativity
  const BcTree& T = bcTree.getBlockCutTree();
  BcTreeBoolBlockNodeMap treeNeg(T, true);
  for (BcTreeBlockNodeIt v(T); v != lemon::INVALID; ++v)
  {
    const EdgeVector& realEdges = bcTree.getRealEdges(v);
    for (EdgeVectorIt edgeIt = realEdges.begin();
         edgeIt != realEdges.end(); ++edgeIt)
    {
      treeNeg[v] = treeNeg[v] && (score[g.u(*edgeIt)] <= 0 && score[g.v(*edgeIt)] <= 0);
    }
  }
  
  bcTree.printNodes(std::cout);
  bcTree.printEdges(std::cout);
  
  int nBlocks = bcTree.getNumBlockTreeNodes();
  int blockIndex = 0;
  
  for (int blockDegree = 1; blockDegree >= 0; --blockDegree)
  {
    const BcTreeBlockNodeSet& leaves = bcTree.getBlockNodeSetByDegree(blockDegree);
    while (!leaves.empty())
    {
      BcTreeBlockNode b = *leaves.begin();
      if (!solveBlock(mwcsGraph,
                      bcTree, treeNeg, b,
                      blockIndex, nBlocks))
      {
        return false;
      }
      
      // update block-cut vertex tree
      bcTree.removeBlockNode(b);
      ++blockIndex;
    }
  }
  
  // solve
  mwcsGraph.print(std::cout);

  BoolNodeMap solutionMap(g, false);
  _pImpl->init(mwcsGraph);
  return _pImpl->solve(solutionScore, solutionMap, solutionSet);

}
  
template<typename GR, typename WGHT, typename NLBL, typename EWGHT>
inline double EnumSolverUnrooted<GR, WGHT, NLBL, EWGHT>::shortCircuit(const Graph& orgG,
                                                                      const SubGraph& g,
                                                                      const WeightNodeMap& score,
                                                                      const SpqrType& spqr,
                                                                      const RootedSpqrTree& rootedT,
                                                                      const SpqrTreeNode v,
                                                                      const SpqrTreeNodeSetNodeMap& orgNodesInSubTree,
                                                                      const bool reverse,
                                                                      BoolNodeMap& filter,
                                                                      DoubleArcMap& arcCost,
                                                                      NodeList& path)
{
  double pathLength = 0;
  
  // identify cut pair
  const RootedSpqrTreeInArcIt cutArc(rootedT, v);
  
  if (cutArc == lemon::INVALID)
  {
    // degenerate case (reverse == false): complete biconnected component is negative, i.e. v is the root
    return pathLength;
  }
  
  const NodePair& cutPair = spqr.getCutPair(cutArc);
  path.clear();
  
  // find the shortest path
  const NodeSet& orgNodes = orgNodesInSubTree[v];
  lemon::mapFill(g, filter, reverse);
  for (NodeSetIt nodeIt = orgNodes.begin(); nodeIt != orgNodes.end(); ++nodeIt)
  {
    filter[*nodeIt] = !reverse;
  }
  if (reverse)
  {
    filter[cutPair.first] = filter[cutPair.second] = true;
  }
  assert(filter[cutPair.first] && filter[cutPair.second]);
  
  // let's construct the subgraph
  lemon::mapFill(g, arcCost, 0);
  SubGraph subG(orgG, filter);
  
  for (SubArcIt a(subG); a != lemon::INVALID; ++a)
  {
    arcCost[a] = score[subG.target(a)] > 0 ? 0 : -score[subG.target(a)];
    assert(arcCost[a] >= 0);
  }
  
  // and now let's do a Dijkstra
  lemon::Dijkstra<SubGraph, DoubleArcMap> dijkstra(subG, arcCost);
  dijkstra.run(cutPair.first);
  
  SubNode u = cutPair.second;
  while (dijkstra.predArc(u) != lemon::INVALID)
  {
    u = dijkstra.predNode(u);
    if (u != cutPair.first)
    {
      path.push_front(u);
      if (score[u] < 0)
        pathLength += score[u];
    }
  }
  
  return pathLength;
}


template<typename GR, typename WGHT, typename NLBL, typename EWGHT>
inline bool EnumSolverUnrooted<GR, WGHT, NLBL, EWGHT>::processBlock(MwcsPreGraphType& mwcsGraph,
                                                                    BlockCutTreeType& bcTree,
                                                                    BcTreeBlockNode b,
                                                                    BoolNodeMap& sameBlock)
{
  const WeightNodeMap& score = mwcsGraph.getScores();
  const Graph& g = mwcsGraph.getGraph();
  
  SubGraph subG(g, sameBlock);
  
  if (lemon::countEdges(subG) <= 2)
  {
    return false;
  }
  
  bool result = false;
  
  SpqrType spqr(subG);
  spqr.run();
  
  const SpqrTree& T = spqr.getSpqrTree();

  // determine root node (max deg node) and negative nodes
  SpqrTreeBoolNodeMap negative(T, true);
  SpqrTreeNode rootNode = lemon::INVALID;
  for (SpqrTreeNodeIt v(T); v != lemon::INVALID; ++v)
  {
    const EdgeVector& realSpqrEdges = spqr.getRealEdges(v);
    for (EdgeVectorIt edgeIt = realSpqrEdges.begin(); edgeIt != realSpqrEdges.end(); ++edgeIt)
    {
      negative[v] = negative[v] && (score[subG.u(*edgeIt)] <= 0 && score[subG.v(*edgeIt)] <= 0);
    }
    
    if (rootNode == lemon::INVALID || spqr.getDegree(v) > spqr.getDegree(rootNode))
    {
      rootNode = v;
    }
  }
  
  // root and determine depth
  SpqrTreeIntNodeMap depth(T, -1), size(T, 1), realEdgesSize(T, 0);
  SpqrTreeBoolNodeMap negativeSubTree(T, false);
  SpqrTreeBoolEdgeMap dir(T, true);
  RootedSpqrTree rootedT(T, dir);
  SpqrTreeNodeSetNodeMap orgNodesInSubTree(T), orgPosNodesInSubTree(T);
  
  rootSpqrTree(subG, score, spqr, negative,
               lemon::INVALID, lemon::INVALID,
               rootNode, rootedT, dir, depth, size,
               realEdgesSize, negativeSubTree,
               orgNodesInSubTree, orgPosNodesInSubTree);
  
  // determine roots of negative subtrees
  SpqrTreeNodeVector negativeRoots;
  identifyNegativeRoots(spqr, rootedT, rootNode, orgPosNodesInSubTree, negativeRoots);
  
  // first let's replace the negative tri components by single nodes and if possible remove them
  int nNodesRemoved = 0;
  BoolNodeMap filter(g, false);
  DoubleArcMap arcCost(g);
  for (SpqrTreeNodeVectorIt it = negativeRoots.begin(); it != negativeRoots.end(); ++it)
  {
    const NodePair& cutPair = spqr.getCutPair(RootedSpqrTreeInArcIt(rootedT, *it));
    
    // first reduce to a single node (if necessary)
    Node singleNode = lemon::INVALID;
    if (realEdgesSize[*it] > 2)
    {
      NodeList path;
      shortCircuit(g, subG, score, spqr, rootedT, *it, orgNodesInSubTree, false, filter, arcCost, path);
      
  //    std::cout << mwcsGraph.getLabel(cutPair.first) << " (" << mwcsGraph.getScore(cutPair.first) << ")";
  //    for (NodeListIt nodeIt = path.begin(); nodeIt != path.end(); ++nodeIt)
  //    {
  //      std::cout << " -- " << mwcsGraph.getLabel(*nodeIt) << " (" << mwcsGraph.getScore(*nodeIt) << ")";
  //    }
  //    std::cout << " -- " << mwcsGraph.getLabel(cutPair.second) << " (" << mwcsGraph.getScore(cutPair.second) << ")" << std::endl;
      
      NodeSet nodesToMerge(path.begin(), path.end());
      singleNode = mwcsGraph.merge(nodesToMerge);
      
      NodeSet nodesToRemove = orgNodesInSubTree[*it];
      for (NodeListIt nodeIt = path.begin(); nodeIt != path.end(); ++nodeIt)
      {
        nodesToRemove.erase(*nodeIt);
      }
      nodesToRemove.erase(cutPair.first);
      nodesToRemove.erase(cutPair.second);
      
      // update block cut tree
      NodeSet nodesToRemoveFromB = nodesToRemove;
      nodesToRemoveFromB.insert(nodesToMerge.begin(), nodesToMerge.end());
      nodesToRemoveFromB.erase(singleNode);
//      bcTree.removeFromBlockNode(b, nodesToRemoveFromB);
      
      mwcsGraph.remove(nodesToRemove);
      
      if (g_verbosity >= VERBOSE_ESSENTIAL)
      {
        std::cout << "// Replacing negative triconnected component of "
                  << orgNodesInSubTree[*it].size() << " nodes"
                  << " and " << realEdgesSize[*it] << " edges with a single node"
                  << std::endl;
      }
      
      result = true;
    }
    else
    {
      assert(realEdgesSize[*it] == 2);
      const Edge& e = spqr.getRealEdges(*it).front();
      if (g.u(e) == cutPair.first || g.u(e) == cutPair.second)
      {
        singleNode = g.v(e);
      }
      else
      {
        singleNode = g.u(e);
      }
    }
    
    // now let's check whether we can remove the isolated node
    NodeList path;
    double pathLength = shortCircuit(g, subG, score, spqr, rootedT, *it, orgNodesInSubTree, true, filter, arcCost, path);
    
    if (pathLength >= score[singleNode])
    {
      ++nNodesRemoved;

      // update block cut tree
      NodeSet nodesToRemoveFromB;
      nodesToRemoveFromB.insert(singleNode);
//      bcTree.removeFromBlockNode(b, nodesToRemoveFromB);
      mwcsGraph.remove(nodesToRemoveFromB);
      
      result = true;
    }
  }
  
  if (g_verbosity >= VERBOSE_ESSENTIAL)
  {
    std::cout << "// Removed " << nNodesRemoved
              << " negative triconnected component of 1 node and 2 edges" << std::endl;
  }
  
//  mwcsGraph.preprocess(NodeSet());
  bcTree.recomputeRealNodesAndEdges();
  
  // Next, look for subtrees with more than 4 nodes
//  for (RootedSpqrTreeOutArcIt a(rootedT, rootNode); a != lemon::INVALID; ++a)
//  {
//    const SpqrTreeNode child = T.v(a);
//
//    if (!negativeSubTree[child] && orgNodesInSubTree[child].size() > 4)
//    {
//      const NodePair& cutPair = spqr.getCutPair(a);
//      std::cout << orgNodesInSubTree[child].size() << std::endl;
//    }
//  }
  
  return result;
  
//  {
//    std::ofstream eOut("spqr.edges.txt");
//    // let's print the spqr tree
//    for (SpqrTreeEdgeIt e(T); e != lemon::INVALID; ++e)
//    {
//      eOut << T.id(T.u(e)) << " (pp) " << T.id(T.v(e))
//           << "\t" << subG.id(spqr.getCutPair(e).first)
//           << "\t" << subG.id(spqr.getCutPair(e).second)
//           << std::endl;
//    }
//    eOut.close();
//    
//    std::ofstream nOut("spqr.nodes.txt");
//    // let's print the spqr tree
//    for (SpqrTreeNodeIt v(T); v != lemon::INVALID; ++v)
//    {
//      nOut << T.id(v) << "\t" << spqr.toChar(spqr.getSpqrNodeType(v))
//           << "\t" << negative[v] << "\t"
//           << realEdgesSize[v] << "\t"
//           << orgNodesInSubTree[v].size() << "\t"
//           << orgPosNodesInSubTree[v].size() << std::endl;
//    }
//    nOut.close();
//  }
}
  
template<typename GR, typename WGHT, typename NLBL, typename EWGHT>
inline void EnumSolverUnrooted<GR, WGHT, NLBL, EWGHT>::identifyNegativeRoots(const SpqrType& spqr,
                                                                             const RootedSpqrTree& rootedT,
                                                                             const SpqrTreeNode v,
                                                                             const SpqrTreeNodeSetNodeMap&  orgPosNodesInSubTree,
                                                                             SpqrTreeNodeVector& roots)
{
  // identify cut pair
  const typename RootedSpqrTree::InArcIt cutArc(rootedT, v);
  
  NodeSet posNodes = orgPosNodesInSubTree[v];
  if (cutArc != lemon::INVALID)
  {
    const NodePair& cutPair = spqr.getCutPair(cutArc);
    posNodes.erase(cutPair.first);
    posNodes.erase(cutPair.second);
  }
  
  if (posNodes.size() == 0)
  {
    roots.push_back(v);
  }
  else
  {
    for (typename RootedSpqrTree::OutArcIt a(rootedT, v); a != lemon::INVALID; ++a)
    {
      identifyNegativeRoots(spqr, rootedT, rootedT.target(a), orgPosNodesInSubTree, roots);
    }
  }
}
  
template<typename GR, typename WGHT, typename NLBL, typename EWGHT>
inline void EnumSolverUnrooted<GR, WGHT, NLBL, EWGHT>::rootSpqrTree(const SubGraph& g,
                                                                    const WeightNodeMap& score,
                                                                    const SpqrType& spqr,
                                                                    const SpqrTreeBoolNodeMap& neg,
                                                                    const SpqrTreeNode u,
                                                                    const SpqrTreeEdge uv,
                                                                    const SpqrTreeNode v,
                                                                    RootedSpqrTree& rootedT,
                                                                    SpqrTreeBoolEdgeMap& dir,
                                                                    SpqrTreeIntNodeMap& depth,
                                                                    SpqrTreeIntNodeMap& size,
                                                                    SpqrTreeIntNodeMap& realEdgesSize,
                                                                    SpqrTreeBoolNodeMap& negSubTree,
                                                                    SpqrTreeNodeSetNodeMap& orgNodesInSubTree,
                                                                    SpqrTreeNodeSetNodeMap& orgPosNodesInSubTree)
{
  const SpqrTree& T = spqr.getSpqrTree();
  
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
    Node org_u = g.u(*edgeIt);
    Node org_v = g.v(*edgeIt);
    
    orgNodesInSubTree[v].insert(org_u);
    orgNodesInSubTree[v].insert(org_v);
    
    if (score[org_u] > 0)
    {
      orgPosNodesInSubTree[v].insert(org_u);
    }
    if (score[org_v] > 0)
    {
      orgPosNodesInSubTree[v].insert(org_v);
    }
  }
  
  for (SpqrTreeIncEdgeIt e(T, v); e != lemon::INVALID; ++e)
  {
    SpqrTreeNode w = T.oppositeNode(v, e);
    if (w != u)
    {
      rootSpqrTree(g, score, spqr, neg, v, e, w,
                   rootedT, dir, depth, size,
                   realEdgesSize, negSubTree,
                   orgNodesInSubTree, orgPosNodesInSubTree);
      negSubTree[v] = negSubTree[v] && negSubTree[w];
      size[v] += size[w];
      realEdgesSize[v] += realEdgesSize[w];
      orgNodesInSubTree[v].insert(orgNodesInSubTree[w].begin(), orgNodesInSubTree[w].end());
      
      //      const NodePair& cutPair = spqr.getCutPair(e);
      orgPosNodesInSubTree[v].insert(orgPosNodesInSubTree[w].begin(), orgPosNodesInSubTree[w].end());
      //      if (score[cutPair.first] > 0)
      //      {
      //        orgPosNodesInSubTree[v].insert(cutPair.first);
      //      }
      //      if (score[cutPair.second] > 0)
      //      {
      //        orgPosNodesInSubTree[v].insert(cutPair.second);
      //      }
    }
  }

}
  
template<typename GR, typename WGHT, typename NLBL, typename EWGHT>
inline bool EnumSolverUnrooted<GR, WGHT, NLBL, EWGHT>::solveBlock(MwcsPreGraphType& mwcsGraph,
                                                                  BlockCutTreeType& bcTree,
                                                                  const BcTreeBoolBlockNodeMap& bcTreeNeg,
                                                                  BcTreeBlockNode b,
                                                                  int blockIndex,
                                                                  int nBlocks)
{
  const Graph& g = mwcsGraph.getGraph();
  const BcTree& T = bcTree.getBlockCutTree();
  
  typename BcTree::Edge e(typename BcTree::IncEdgeIt(T, b));
  BcTreeCutNode c = e != lemon::INVALID ? T.redNode(e) : lemon::INVALID;
  Node orgC = c != lemon::INVALID ? bcTree.getArticulationPoint(c) : lemon::INVALID;
  
  const NodeSet& nodesB = bcTree.getRealNodes(b);
//  mwcsGraph.print(std::cout);
  
  if (g_verbosity >= VERBOSE_ESSENTIAL)
  {

    
    std::cout << std::endl;
    std::cout << "// Considering block " << blockIndex + 1 << "/" << nBlocks
              << ": contains " << nodesB.size() << " nodes and "
              << bcTree.getRealEdges(b).size() << " edges" << std::endl;
//    std::cout << "Nodes in block:" << std::endl;
//    printNodeSet(mwcsGraph, nodesB);
//    if (orgC != lemon::INVALID)
//    {
//      std::cout << "Cut node: " << mwcsGraph.getLabel(orgC) << std::endl;
//    }
  }

  if (bcTreeNeg[b])
  {
    if (g_verbosity >= VERBOSE_ESSENTIAL)
    {
      std::cout << "// Removed block, as it consists of only negatively weighted nodes"
      << std::endl;
    }
    
    // don't remove cut node if it connects to other blocks
    NodeSet nodesToRemove = nodesB;
    if (c != lemon::INVALID && bcTree.getDegree(c) > 1)
    {
      nodesToRemove.erase(orgC);
    }
    else if (g_verbosity >= VERBOSE_ESSENTIAL)
    {
      std::cout << "// Removed corresponding cut vertex, as there are no other blocks"
      << std::endl;
    }
    
    // remove block from graph
    mwcsGraph.remove(nodesToRemove);
  }
  else
  {
    // create a new graph containing only block b
    Graph subG;
    DoubleNodeMap weightSubG(subG);
    LabelNodeMap labelSubG(subG);
    NodeMap mapToG(subG);
    NodeMap mapToSubG(g);
    MwcsPreGraphType mwcsSubGraph;
    
    BoolNodeMap sameBlock(g, false);
    for (NodeSetIt nodeIt = nodesB.begin(); nodeIt != nodesB.end(); ++nodeIt)
    {
      sameBlock[*nodeIt] = true;
    }
    
    while (processBlock(mwcsGraph, bcTree, b, sameBlock));
    
    printNodeSet(mwcsGraph, nodesB);
    
    initLocalGraph(g,
                   mwcsGraph.getScores(),
                   mwcsGraph.getLabels(),
                   sameBlock,
                   subG,
                   weightSubG,
                   labelSubG,
                   mapToG,
                   mapToSubG,
                   mwcsSubGraph);
    
    // solve the unrooted formulation first
    NodeSet subSolutionSet;
    double subSolutionScore;
    if (!solveBlockUnrooted(mwcsSubGraph, subSolutionSet, subSolutionScore))
    {
      return false;
    }
   
    // let's map the solution back to our node space
    NodeSet solutionSet, solutionComplementSet;
    map(mwcsSubGraph, mapToG, subSolutionSet, solutionSet);
    
    std::set_difference(nodesB.begin(), nodesB.end(),
                        solutionSet.begin(), solutionSet.end(),
                        std::inserter(solutionComplementSet, solutionComplementSet.begin()));
    
    // now check if the cut node is in the solution
    if (orgC == lemon::INVALID || solutionSet.find(orgC) != solutionSet.end())
    {
//      std::cout << std::endl << "Sol set" << std::endl;
//      printNodeSet(mwcsGraph, solutionSet);
//      std::cout << std::endl << "Complement set" << std::endl;
//      printNodeSet(mwcsGraph, solutionComplementSet);

      // no need to solve the rooted formulation!
      // just merge the solution nodes
      
      // cut node should be kept!!!! update merge function!
      if (orgC == lemon::INVALID)
      {
        mwcsGraph.merge(solutionSet);
      }
      else
      {
        mwcsGraph.merge(orgC, solutionSet);
      }
      
      // and remove the other nodes
      
      // but don't the remove cut node if it connects to other blocks
      if (c != lemon::INVALID && bcTree.getDegree(c) > 1)
      {
        solutionComplementSet.erase(orgC);
      }
      
      mwcsGraph.remove(solutionComplementSet);
    }
    else if (orgC != lemon::INVALID) // orgC is not in the solutionSet
    {
      // extract into an isolated node
      mwcsGraph.extract(solutionSet);
      
      // recreate subgraph
      initLocalGraph(g,
                     mwcsGraph.getScores(),
                     mwcsGraph.getLabels(),
                     sameBlock,
                     subG,
                     weightSubG,
                     labelSubG,
                     mapToG,
                     mapToSubG,
                     mwcsSubGraph);
      
      // solve rooted formulation
      subSolutionSet.clear();
      subSolutionScore = 0;
      if (!solveBlockRooted(mwcsSubGraph, mapToSubG[orgC], subSolutionSet, subSolutionScore))
      {
        return false;
      }
      
      solutionSet.clear();
      map(mwcsSubGraph, mapToG, subSolutionSet, solutionSet);
      
      assert(solutionSet.find(orgC) != solutionSet.end());
      
      solutionComplementSet.clear();
      std::set_difference(nodesB.begin(), nodesB.end(),
                          solutionSet.begin(), solutionSet.end(),
                          std::inserter(solutionComplementSet, solutionComplementSet.begin()));
      
      // if rooted solution is negative then discard
      // otherwise merge the solution (collapse into orgC)
      // and remove the nodes that are not in the solution
      if (subSolutionScore > 0)
      {
        mwcsGraph.merge(orgC, solutionSet);
        
        mwcsGraph.remove(solutionComplementSet);
      }
      else
      {
        // but don't the remove cut node if it connects to other blocks
        if (c != lemon::INVALID && bcTree.getDegree(c) > 1)
        {
          solutionComplementSet.erase(orgC);
        }
        
        mwcsGraph.remove(solutionComplementSet);
      }
    }
  }
  
  return true;
}
  
template<typename GR, typename WGHT, typename NLBL, typename EWGHT>
inline bool EnumSolverUnrooted<GR, WGHT, NLBL, EWGHT>::solveBlockUnrooted(MwcsPreGraphType& mwcsGraph,
                                                                          NodeSet& solutionSet,
                                                                          double& solutionScore)
{
  const Graph& g = mwcsGraph.getGraph();
  BoolNodeMap solutionMap(g, false);
  
  if (_preprocess)
  {
    // preprocess the graph
    mwcsGraph.preprocess(NodeSet());
  }
  
  _pImpl->init(mwcsGraph);
  return _pImpl->solve(solutionScore, solutionMap, solutionSet);
}
  
template<typename GR, typename WGHT, typename NLBL, typename EWGHT>
inline bool EnumSolverUnrooted<GR, WGHT, NLBL, EWGHT>::solveBlockRooted(MwcsPreGraphType& mwcsGraph,
                                                                        Node rootNode,
                                                                        NodeSet& solutionSet,
                                                                        double& solutionScore)
{
  const Graph& g = mwcsGraph.getGraph();
  BoolNodeMap solutionMap(g, false);
  
  NodeSet rootNodes = mwcsGraph.getPreNodes(rootNode);
//  std::cout << "Root nodes: " << std::endl;
//  printNodeSet(mwcsGraph, rootNodes);
  
  if (_preprocess)
  {
    // preprocess the graph
//    mwcsGraph.print(std::cout);
    mwcsGraph.preprocess(rootNodes);
  }
  
  _pRootedImpl->init(mwcsGraph, rootNodes);
  return _pRootedImpl->solve(solutionScore, solutionMap, solutionSet);
}

} // namespace mwcs
} // namespace nina

#endif // ENUMSOLVERUNROOTED_H
