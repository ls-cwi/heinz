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
#include <assert.h>
#include <ostream>

#include "mwcs.h"
#include "mwcsgraph.h"
#include "mwcspreprocessedgraph.h"

#include "solver/solverunrooted.h"
#include "solver/impl/solverrootedimpl.h"

//#include "preprocessing/negdeg01.h"
//#include "preprocessing/posedge.h"
//#include "preprocessing/negedge.h"
//#include "preprocessing/negcircuit.h"
//#include "preprocessing/negdiamond.h"
//#include "preprocessing/negmirroredhubs.h"
//#include "preprocessing/posdeg01.h"
//#include "preprocessing/posdiamond.h"

#include "blockcuttree.h"

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

//  // preprocessing rules
//  typedef NegDeg01<Graph, WeightNodeMap> NegDeg01Type;
//  typedef PosEdge<Graph, WeightNodeMap> PosEdgeType;
//  typedef NegEdge<Graph, WeightNodeMap> NegEdgeType;
//  typedef NegCircuit<Graph, WeightNodeMap> NegCircuitType;
//  typedef NegDiamond<Graph, WeightNodeMap> NegDiamondType;
//  typedef NegMirroredHubs<Graph, WeightNodeMap> NegMirroredHubsType;
//  typedef PosDeg01<Graph, WeightNodeMap> PosDeg01Type;
//  typedef PosDiamond<Graph, WeightNodeMap> PosDiamondType;
  
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
  typedef typename SubGraph::NodeIt SubNodeIt;
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
  
private:
  SolverRootedImplType* _pRootedImpl;
  bool _preprocess;
  
  bool solveComponent(MwcsPreGraphType& mwcsGraph,
                      NodeSet& solutionSet,
                      double& solutionScore);
  
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
  
  void processBlock(MwcsPreGraphType& mwcsGraph,
                    SubGraph& subG);
  
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
inline void EnumSolverUnrooted<GR, WGHT, NLBL, EWGHT>::processBlock(MwcsPreGraphType& mwcsGraph,
                                                                    SubGraph& subG)
{
  
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
  
  mwcsGraph.print(std::cout);
  
  NodeSet nodesB;
  bcTree.getRealNodes(b, nodesB);
  
  if (g_verbosity >= VERBOSE_ESSENTIAL)
  {
    std::cout << std::endl;
    std::cout << "// Considering block " << blockIndex + 1 << "/" << nBlocks
              << ": contains " << nodesB.size() << " nodes and "
              << bcTree.getRealEdges(b).size() << " edges" << std::endl;
    std::cout << "Nodes in block:" << std::endl;
    printNodeSet(mwcsGraph, nodesB);
    if (orgC != lemon::INVALID)
    {
      std::cout << "Cut node: " << mwcsGraph.getLabel(orgC) << std::endl;
    }
  }

  if (bcTreeNeg[b])
  {
    if (g_verbosity >= VERBOSE_ESSENTIAL)
    {
      std::cout << "// Removed block, as it consists of only negatively weighted nodes"
      << std::endl;
    }
    
    // don't remove cut node if it connects to other blocks
    if (c != lemon::INVALID && bcTree.getDegree(c) > 1)
    {
      nodesB.erase(orgC);
    }
    else if (g_verbosity >= VERBOSE_ESSENTIAL)
    {
      std::cout << "// Removed corresponding cut vertex, as there are no other blocks"
      << std::endl;
    }
    
    // remove block from graph
    mwcsGraph.remove(nodesB);
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
      std::cout << std::endl << "Sol set" << std::endl;
      printNodeSet(mwcsGraph, solutionSet);
      std::cout << std::endl << "Complement set" << std::endl;
      printNodeSet(mwcsGraph, solutionComplementSet);

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
