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
  
  bool processBlock1(MwcsPreGraphType& mwcsGraph,
                     BlockCutTreeType& bcTree,
                     BcTreeBlockNode b,
                     BcTreeCutNode c,
                     Node orgC,
                     BoolNodeMap& sameBlock);

  bool processBlock2(MwcsPreGraphType& mwcsGraph,
                     BlockCutTreeType& bcTree,
                     BcTreeBlockNode b,
                     BcTreeCutNode c,
                     Node orgC,
                     BoolNodeMap& sameBlock);
  
  bool solveBlock(MwcsPreGraphType& mwcsGraph,
                  BlockCutTreeType& bcTree,
                  const BcTreeBoolBlockNodeMap& bcTreeNeg,
                  BcTreeBlockNode b,
                  int blockIndex,
                  int nBlocks);
  
  
  bool solveTriComp(MwcsPreGraphType& mwcsGraph,
                    const NodePair& cutPair,
                    const NodeSet& nodesTriComp,
                    BlockCutTreeType& bcTree,
                    BcTreeBlockNode b);
  
  bool solveUnrooted(MwcsPreGraphType& mwcsGraph,
                     NodeSet& solutionSet,
                     double& solutionScore);

  bool solveRooted(MwcsPreGraphType& mwcsGraph,
                   const NodeSet& rootNodes,
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
  
  for (int blockDegree = nBlocks > 1 ? 1 : 0; blockDegree >= 0; --blockDegree)
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
//  mwcsGraph.print(std::cout);

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
  
  assert(dijkstra.reached(cutPair.second));
  
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
inline bool EnumSolverUnrooted<GR, WGHT, NLBL, EWGHT>::processBlock2(MwcsPreGraphType& mwcsGraph,
                                                                     BlockCutTreeType& bcTree,
                                                                     BcTreeBlockNode b,
                                                                     BcTreeCutNode c,
                                                                     Node orgC,
                                                                     BoolNodeMap& sameBlock)
{
  bool result = false;
  
  const WeightNodeMap& score = mwcsGraph.getScores();
  const Graph& g = mwcsGraph.getGraph();
  
  SubGraph subG(g, sameBlock);
  
  if (lemon::countEdges(subG) <= 2)
  {
    return false;
  }
  
  SpqrType spqr(subG);
  bool spqrRes = spqr.run();
  assert(spqrRes);
  
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
  
  // Next, look for subtrees with more than 4 nodes
  typedef std::pair<int, SpqrTreeNode> IntSpqrTreeNodePair;
  typedef std::vector<IntSpqrTreeNodePair> IntSpqrTreeNodePairVector;
  typedef typename IntSpqrTreeNodePairVector::const_iterator IntSpqrTreeNodePairVectorIt;
  
  IntSpqrTreeNodePairVector triComponents;
  for (RootedSpqrTreeOutArcIt a(rootedT, rootNode); a != lemon::INVALID; ++a)
  {
    const SpqrTreeNode child = T.v(a);
    if (child != rootNode && orgNodesInSubTree[child].size() > 3)
    {
      assert(!negativeSubTree[child]);
      
//      const NodePair& cutPair = spqr.getCutPair(a);
//      std::cout << T.id(child) << " " << g.id(cutPair.first) << " " << g.id(cutPair.second)
//                << " " << orgNodesInSubTree[child].size() << std::endl;
      triComponents.push_back(std::make_pair(static_cast<int>(orgNodesInSubTree[child].size()), child));
    }
  }
  
  std::sort(triComponents.begin(), triComponents.end());
  int triCompIdx = 0;
  for (IntSpqrTreeNodePairVectorIt triCompIt = triComponents.begin();
       triCompIt != triComponents.end(); ++triCompIt, ++triCompIdx)
  {
    if (g_verbosity >= VERBOSE_ESSENTIAL)
    {
      std::cout << std::endl;
      std::cout << "// Considering triconnected component " << triCompIdx + 1 << "/" << triComponents.size()
                << ": contains " << orgNodesInSubTree[triCompIt->second].size() << " nodes and "
                << triCompIt->first << " edges" << std::endl;
    }
    
    const NodePair& cutPair = spqr.getCutPair(RootedSpqrTreeInArcIt(rootedT, triCompIt->second));
    result |= solveTriComp(mwcsGraph, cutPair, orgNodesInSubTree[triCompIt->second], bcTree, b);
  }
  
  NodeSet cutNodeSet;
  if (orgC != lemon::INVALID)
  {
    cutNodeSet.insert(orgC);
  }
  
  mwcsGraph.preprocess(cutNodeSet);
  bcTree.recomputeRealNodesAndEdges();
  
//  std::cout << T.id(triComponents.begin()->second) << " " << triComponents.begin()->first << std::endl;
  
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
//
//    std::ofstream nnOut("nodes.txt");
//    mwcsGraph.printNodeList(nnOut, false);
//    nnOut.close();
//    
//    std::ofstream eeOut("edges.txt");
//    mwcsGraph.printEdgeList(eeOut, false);
//    eeOut.close();
//  }
  
  return result;
}
  
template<typename GR, typename WGHT, typename NLBL, typename EWGHT>
inline bool EnumSolverUnrooted<GR, WGHT, NLBL, EWGHT>::solveTriComp(MwcsPreGraphType& mwcsGraph,
                                                                    const NodePair& cutPair,
                                                                    const NodeSet& nodesTriComp,
                                                                    BlockCutTreeType& bcTree,
                                                                    BcTreeBlockNode b)
{
  const Graph& g = mwcsGraph.getGraph();
  const DoubleNodeMap& score = mwcsGraph.getScores();
  
  // create a new graph induced by nodesTriComp
  Graph subG;
  DoubleNodeMap weightSubG(subG);
  LabelNodeMap labelSubG(subG);
  NodeMap mapToG(subG);
  NodeMap mapToSubG(g);
  MwcsPreGraphType mwcsSubGraph;
  
  BoolNodeMap sameTriComp(g, false);
  for (NodeSetIt nodeIt = nodesTriComp.begin(); nodeIt != nodesTriComp.end(); ++nodeIt)
  {
    sameTriComp[*nodeIt] = true;
  }
  
  initLocalGraph(g,
                 mwcsGraph.getScores(),
                 mwcsGraph.getLabels(),
                 sameTriComp,
                 subG,
                 weightSubG,
                 labelSubG,
                 mapToG,
                 mapToSubG,
                 mwcsSubGraph);
  
  // start by solving the unrooted formulation
//  printNodeSet(mwcsGraph, nodesTriComp);
//  std::cout << mwcsGraph.getLabel(cutPair.first) << " -- " << mwcsGraph.getLabel(cutPair.second) << std::endl;
//  mwcsSubGraph.print(std::cout);
//  assert(lemon::connected(subG));
  
  NodeSet subSolutionSet;
  double solutionScore;
  if (!solveUnrooted(mwcsSubGraph, subSolutionSet, solutionScore))
  {
    abort();
  }
  
  assert(nodesTriComp.find(cutPair.first) != nodesTriComp.end());
  assert(nodesTriComp.find(cutPair.second) != nodesTriComp.end());
  
  NodeSet V4;
  map(mwcsSubGraph, mapToG, subSolutionSet, V4);
  
//  printNodeSet(mwcsSubGraph, subSolutionSet);
//  printNodeSet(mwcsGraph, V4);
  
  // V1 is rooted at cutPair.first
  // V2 is rooted at cutPair.second
  // V3 is rooted at cutPair.first and cutPair.second
  NodeSet V1, V2, V3;
  
  if (V4.find(cutPair.first) != V4.end())
  {
    V1 = V4;
    V1.erase(cutPair.first);
  }
  else
  {
    // solve rooted at cutPair.first
    mwcsSubGraph.clear();
    
    subSolutionSet.clear();
    if (!solveRooted(mwcsSubGraph,
                     mwcsSubGraph.getPreNodes(mapToSubG[cutPair.first]),
                     subSolutionSet,
                     solutionScore))
    {
//      std::cout << mwcsGraph.getLabel(cutPair.first) << std::endl;
      abort();
    }
    map(mwcsSubGraph, mapToG, subSolutionSet, V1);
    V1.erase(cutPair.first);
    V1.erase(cutPair.second);
  }
  
  if (V4.find(cutPair.second) != V4.end())
  {
    V2 = V4;
    V2.erase(cutPair.second);
  }
  else
  {
    // solve rooted at cutPair.second
    mwcsSubGraph.clear();
    
    subSolutionSet.clear();
    if (!solveRooted(mwcsSubGraph,
                     mwcsSubGraph.getPreNodes(mapToSubG[cutPair.second]),
                     subSolutionSet,
                     solutionScore))
    {
//      mwcsSubGraph.print(std::cout);
//      std::cout << mwcsGraph.getLabel(cutPair.second) << std::endl;
      abort();
    }
    map(mwcsSubGraph, mapToG, subSolutionSet, V2);
    V2.erase(cutPair.first);
    V2.erase(cutPair.second);
  }
  
  if (V4.find(cutPair.first) != V4.end() && V4.find(cutPair.second) != V4.end())
  {
    V3 = V4;
    V3.erase(cutPair.first);
    V3.erase(cutPair.second);
  }
  else if (V1.find(cutPair.first) != V1.end() && V1.find(cutPair.second) != V1.end())
  {
    V3 = V1;
    V3.erase(cutPair.first);
    V3.erase(cutPair.second);
  }
  else if (V2.find(cutPair.first) != V2.end() && V2.find(cutPair.second) != V2.end())
  {
    V3 = V2;
    V3.erase(cutPair.first);
    V3.erase(cutPair.second);
  }
  else
  {
    // solve rooted at cutPair.first and cutPair.second
    mwcsSubGraph.clear();
    
    NodeSet rootNodes = mwcsSubGraph.getPreNodes(mapToSubG[cutPair.first]);
    const NodeSet& tmp = mwcsSubGraph.getPreNodes(mapToSubG[cutPair.second]);
    rootNodes.insert(tmp.begin(), tmp.end());

    assert(rootNodes.size() == 2);
    
//    mwcsSubGraph.print(std::cout);
    
    subSolutionSet.clear();
    if (!solveRooted(mwcsSubGraph,
                     rootNodes,
                     subSolutionSet,
                     solutionScore))
    {
//      mwcsSubGraph.print(std::cout);
      abort();
    }
    
    map(mwcsSubGraph, mapToG, subSolutionSet, V3);
    V3.erase(cutPair.first);
    V3.erase(cutPair.second);
  }
  
//  std::cout << "V1" << std::endl;
//  printNodeSet(mwcsGraph, V1);
//  std::cout << std::endl << "V2" << std::endl;
//  printNodeSet(mwcsGraph, V2);
//  std::cout << std::endl << "V3" << std::endl;
//  printNodeSet(mwcsGraph, V3);
//  std::cout << std::endl << "V4" << std::endl;
//  printNodeSet(mwcsGraph, V4);

  // introduce gadget
  mwcsGraph.extract(V4);
  
  NodeSet V1_minus_V2;
  std::set_difference(V1.begin(), V1.end(),
                      V2.begin(), V2.end(),
                      std::inserter(V1_minus_V2, V1_minus_V2.begin()));
  
  NodeSet V2_minus_V1;
  std::set_difference(V2.begin(), V2.end(),
                      V1.begin(), V1.end(),
                      std::inserter(V2_minus_V1, V2_minus_V1.begin()));
  
  NodeSet V1_cup_V2 = V1;
  V1_cup_V2.insert(V2.begin(), V2.end());
  
  NodeSet V1_cup_V2_cup_V3 = V1_cup_V2;
  V1_cup_V2_cup_V3.insert(V3.begin(), V3.end());
  
  NodeSet comp_V1_cup_V2_cup_V3;
  std::set_difference(nodesTriComp.begin(), nodesTriComp.end(),
                      V1_cup_V2_cup_V3.begin(), V1_cup_V2_cup_V3.end(),
                      std::inserter(comp_V1_cup_V2_cup_V3, comp_V1_cup_V2_cup_V3.begin()));
  comp_V1_cup_V2_cup_V3.erase(cutPair.first);
  comp_V1_cup_V2_cup_V3.erase(cutPair.second);
  
  NodeSet V3_minus_V1_cup_V2;
  std::set_difference(V3.begin(), V3.end(),
                      V1_cup_V2.begin(), V1_cup_V2.end(),
                      std::inserter(V3_minus_V1_cup_V2, V3_minus_V1_cup_V2.begin()));
  
  NodeSet V1_cap_V2;
  std::set_intersection(V1.begin(), V1.end(),
                        V2.begin(), V2.end(),
                        std::inserter(V1_cap_V2, V1_cap_V2.begin()));
  
  assert(V1.find(cutPair.first) == V1.end());
  assert(V2.find(cutPair.first) == V2.end());
  assert(V3.find(cutPair.first) == V3.end());
  assert(V1.find(cutPair.second) == V1.end());
  assert(V2.find(cutPair.second) == V2.end());
  assert(V3.find(cutPair.second) == V3.end());
  
  int res = 0;
  NodeSet gadget;
  
  Node nV1_minus_V2 = lemon::INVALID;
  if (!V1_minus_V2.empty())
  {
    nV1_minus_V2 = mwcsGraph.merge(V1_minus_V2);
    gadget.insert(nV1_minus_V2);
    ++res;
    
    // let's disconnect nV1_minus_V2
    while (IncEdgeIt(g, nV1_minus_V2) != lemon::INVALID)
      mwcsGraph.remove(IncEdgeIt(g, nV1_minus_V2));
  }
  
  Node nV2_minus_V1 = lemon::INVALID;
  if (!V2_minus_V1.empty())
  {
    nV2_minus_V1 = mwcsGraph.merge(V2_minus_V1);
    gadget.insert(nV2_minus_V1);
    ++res;
    
    // let's disconnect nV2_minus_V1
    while (IncEdgeIt(g, nV2_minus_V1) != lemon::INVALID)
      mwcsGraph.remove(IncEdgeIt(g, nV2_minus_V1));
  }
  
  Node nV1_cap_V2 = lemon::INVALID;
  if (!V1_cap_V2.empty())
  {
    nV1_cap_V2 = mwcsGraph.merge(V1_cap_V2);
    gadget.insert(nV1_cap_V2);
    ++res;
    
    // let's disconnect nV1_cap_V2
    while (IncEdgeIt(g, nV1_cap_V2) != lemon::INVALID)
      mwcsGraph.remove(IncEdgeIt(g, nV1_cap_V2));
  }
  
  Node nV3_minus_V1_cup_V2 = lemon::INVALID;
  if (!V3_minus_V1_cup_V2.empty())
  {
    nV3_minus_V1_cup_V2 = mwcsGraph.merge(V3_minus_V1_cup_V2);
    gadget.insert(nV3_minus_V1_cup_V2);
    ++res;
    
    // let's disconnect nV3_minus_V1_cup_V2
    while (IncEdgeIt(g, nV3_minus_V1_cup_V2) != lemon::INVALID)
      mwcsGraph.remove(IncEdgeIt(g, nV3_minus_V1_cup_V2));
  }
  
  if (!comp_V1_cup_V2_cup_V3.empty())
  {
    mwcsGraph.remove(comp_V1_cup_V2_cup_V3);
  }
  
  // let's connect up the gadget
  if (nV1_cap_V2 == lemon::INVALID)
  {
    if (nV1_minus_V2 != lemon::INVALID && nV3_minus_V1_cup_V2 != lemon::INVALID)
    {
      assert(score[nV1_minus_V2] >= 0);
      gadget.erase(nV1_minus_V2);
      nV1_minus_V2 = mwcsGraph.merge(nV1_minus_V2, cutPair.first);
    }
    if (nV2_minus_V1 != lemon::INVALID && nV3_minus_V1_cup_V2 != lemon::INVALID)
    {
      assert(score[nV2_minus_V1] >= 0);
      gadget.erase(nV2_minus_V1);
      nV2_minus_V1 = mwcsGraph.merge(nV2_minus_V1, cutPair.second);
    }
  }
  else
  {
    gadget.insert(nV1_cap_V2);
    if (nV1_minus_V2 != lemon::INVALID)
    {
      Edge e = mwcsGraph.addEdge(nV1_minus_V2, nV1_cap_V2);
      bcTree.assignEdgeToBlock(e, b);
    }
    else
    {
      Edge e = mwcsGraph.addEdge(cutPair.first, nV1_cap_V2);
      bcTree.assignEdgeToBlock(e, b);
    }
    if (nV2_minus_V1 != lemon::INVALID)
    {
      Edge e = mwcsGraph.addEdge(nV2_minus_V1, nV1_cap_V2);
      bcTree.assignEdgeToBlock(e, b);
    }
    else
    {
      Edge e = mwcsGraph.addEdge(cutPair.second, nV1_cap_V2);
      bcTree.assignEdgeToBlock(e, b);
    }
  }
  if (nV1_minus_V2 != lemon::INVALID && nV1_minus_V2 != cutPair.first)
  {
    gadget.insert(nV1_minus_V2);
    Edge e = mwcsGraph.addEdge(cutPair.first, nV1_minus_V2);
    bcTree.assignEdgeToBlock(e, b);
  }
  if (nV2_minus_V1 != lemon::INVALID && nV2_minus_V1 != cutPair.second)
  {
    Edge e = mwcsGraph.addEdge(cutPair.second, nV2_minus_V1);
    bcTree.assignEdgeToBlock(e, b);
  }
  if (nV3_minus_V1_cup_V2 != lemon::INVALID)
  {
    gadget.insert(nV3_minus_V1_cup_V2);
    Edge e1 = mwcsGraph.addEdge(cutPair.first, nV3_minus_V1_cup_V2);
    Edge e2 = mwcsGraph.addEdge(cutPair.second, nV3_minus_V1_cup_V2);
    bcTree.assignEdgeToBlock(e1, b);
    bcTree.assignEdgeToBlock(e2, b);
  }
  
  if (nV3_minus_V1_cup_V2 == lemon::INVALID && nV1_cap_V2 == lemon::INVALID)
  {
    if (nV1_minus_V2 != lemon::INVALID)
    {
      assert(nV2_minus_V1 == lemon::INVALID);
      Edge e = mwcsGraph.addEdge(cutPair.second, nV1_minus_V2);
      bcTree.assignEdgeToBlock(e, b);
    }
    if (nV2_minus_V1 != lemon::INVALID)
    {
      assert(nV1_minus_V2 == lemon::INVALID);
      Edge e = mwcsGraph.addEdge(cutPair.first, nV2_minus_V1);
      bcTree.assignEdgeToBlock(e, b);
    }
  }
  
  gadget.insert(cutPair.first);
  gadget.insert(cutPair.second);
  
  if (g_verbosity >= VERBOSE_ESSENTIAL)
  {
    std::cout << "// Replaced tricomponent of "
              << nodesTriComp.size()
              << " nodes by a gadget of "
              << gadget.size() << " nodes" << std::endl;
  }
  
  // print the gadget
  subG.clear();
  lemon::mapFill(g, sameTriComp, false);
  for (NodeSetIt nodeIt = gadget.begin(); nodeIt != gadget.end(); ++nodeIt)
  {
    sameTriComp[*nodeIt] = true;
  }
  
  mwcsSubGraph.clear();
  initLocalGraph(g,
                 mwcsGraph.getScores(),
                 mwcsGraph.getLabels(),
                 sameTriComp,
                 subG,
                 weightSubG,
                 labelSubG,
                 mapToG,
                 mapToSubG,
                 mwcsSubGraph);
  
//  std::cout << "gadget:" << std::endl;
//  mwcsSubGraph.print(std::cout);
  
  assert(lemon::connected(SubGraph(g, sameTriComp)));
  return gadget.size() < nodesTriComp.size();
}

template<typename GR, typename WGHT, typename NLBL, typename EWGHT>
inline bool EnumSolverUnrooted<GR, WGHT, NLBL, EWGHT>::processBlock1(MwcsPreGraphType& mwcsGraph,
                                                                     BlockCutTreeType& bcTree,
                                                                     BcTreeBlockNode b,
                                                                     BcTreeCutNode c,
                                                                     Node orgC,
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
  bool spqrRes = spqr.run();
  assert(spqrRes);
  
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
  
  // first let's replace the negative tri components by single nodes
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
      bcTree.removeFromBlockNode(b, nodesToRemoveFromB);
      
      mwcsGraph.remove(nodesToRemove);
      
      if (g_verbosity >= VERBOSE_ESSENTIAL)
      {
        if (singleNode != lemon::INVALID)
        {
          std::cout << "// Replacing negative triconnected component of "
                    << orgNodesInSubTree[*it].size() << " nodes"
                    << " and " << realEdgesSize[*it] << " edges with a single node"
                    << std::endl;
        }
        else
        {
          std::cout << "// Removed negative triconnected component of "
                    << orgNodesInSubTree[*it].size() << " nodes"
                    << " and " << realEdgesSize[*it] << " edges"
                    << std::endl;
        }
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
  }
  
  NodeSet cutNodeSet;
  if (orgC != lemon::INVALID)
  {
    cutNodeSet.insert(orgC);
  }
  
  mwcsGraph.preprocess(cutNodeSet);
  bcTree.recomputeRealNodesAndEdges();
  
  return result;
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
      
      orgNodesInSubTree[v].insert(orgNodesInSubTree[w].begin(),
                                  orgNodesInSubTree[w].end());
      
      orgPosNodesInSubTree[v].insert(orgPosNodesInSubTree[w].begin(),
                                     orgPosNodesInSubTree[w].end());
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
  
  if (g_verbosity >= VERBOSE_ESSENTIAL)
  {
    std::cout << std::endl;
    std::cout << "// Considering block " << blockIndex + 1 << "/" << nBlocks
              << ": contains " << nodesB.size() << " nodes and "
              << bcTree.getRealEdges(b).size() << " edges" << std::endl;
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
    
    // process this block
    while (processBlock1(mwcsGraph, bcTree, b, c, orgC, sameBlock)
           || processBlock2(mwcsGraph, bcTree, b, c, orgC, sameBlock));
    
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
    if (!solveUnrooted(mwcsSubGraph, subSolutionSet, subSolutionScore))
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
      // no need to solve the rooted formulation!
      // just merge the solution nodes
      
      // cut node should be kept!!!!
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
      if (!solveRooted(mwcsSubGraph,
                       mwcsSubGraph.getPreNodes(mapToSubG[orgC]),
                       subSolutionSet,
                       subSolutionScore))
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
inline bool EnumSolverUnrooted<GR, WGHT, NLBL, EWGHT>::solveUnrooted(MwcsPreGraphType& mwcsGraph,
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
  
  // are we done?
  if (mwcsGraph.getNodeCount() == 0)
  {
    solutionSet.clear();
    solutionScore = 0;
    return true;
  }
  else if (mwcsGraph.getNodeCount() == 1)
  {
    solutionSet.clear();
    Node v = NodeIt(mwcsGraph.getGraph());
    solutionSet.insert(v);
    solutionScore = mwcsGraph.getScore(v);
    return true;
  }
  
  _pImpl->init(mwcsGraph);
  return _pImpl->solve(solutionScore, solutionMap, solutionSet);
}
  
template<typename GR, typename WGHT, typename NLBL, typename EWGHT>
inline bool EnumSolverUnrooted<GR, WGHT, NLBL, EWGHT>::solveRooted(MwcsPreGraphType& mwcsGraph,
                                                                   const NodeSet& rootNodes,
                                                                   NodeSet& solutionSet,
                                                                   double& solutionScore)
{
  const Graph& g = mwcsGraph.getGraph();
  BoolNodeMap solutionMap(g, false);
  
  if (_preprocess)
  {
    // preprocess the graph
    mwcsGraph.preprocess(rootNodes);
  }
  
  // are we done?
  if (mwcsGraph.getNodeCount() == static_cast<int>(rootNodes.size()))
  {
    solutionSet.clear();
    solutionScore = 0;
    for (NodeSetIt rootIt = rootNodes.begin(); rootIt != rootNodes.end(); ++rootIt)
    {
      solutionSet.insert(*rootIt);
      solutionScore += mwcsGraph.getScore(*rootIt);
    }
    return true;
  }
  
  _pRootedImpl->init(mwcsGraph, rootNodes);
  return _pRootedImpl->solve(solutionScore, solutionMap, solutionSet);
}

} // namespace mwcs
} // namespace nina

#endif // ENUMSOLVERUNROOTED_H
