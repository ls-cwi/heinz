/*
 * negbicomponent.h
 *
 *  Created on: 14-mar-2014
 *      Author: M. El-Kebir
 */

#ifndef NEGBICOMPONENT_H
#define NEGBICOMPONENT_H

#include <lemon/core.h>
#include <lemon/connectivity.h>
#include <string>
#include <vector>
#include <set>
#include "unrootedrule.h"
#include "solver/blocktree.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename WGHT = typename GR::template NodeMap<double> >
class NegBiComponent : public UnrootedRule<GR, WGHT>
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

  NegBiComponent();
  virtual ~NegBiComponent() {}
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

  virtual std::string name() const { return "NegBiComponent"; }
  
private:
  typedef BlockTree<Graph> BlockTreeType;
};

template<typename GR, typename WGHT>
inline NegBiComponent<GR, WGHT>::NegBiComponent()
  : Parent()
{
}

template<typename GR, typename WGHT>
inline int NegBiComponent<GR, WGHT>::apply(Graph& g,
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
  BlockTreeType blockTree(g);
  blockTree.run();
  
  typedef typename BlockTreeType::Tree Tree;
  typedef typename Tree::template NodeMap<int> IntTreeNodeMap;
  typedef typename Tree::template NodeMap<bool> BoolTreeNodeMap;

  const Tree& T = blockTree.getBlockTree();
  
  // determine block tree node degrees and negativity
  IntTreeNodeMap treeDeg(T, 0);
  BoolTreeNodeMap treeNeg(T, true);
  
  for (typename Tree::NodeIt v(T); v != lemon::INVALID; ++v)
  {
    for (typename Tree::IncEdgeIt e(T, v); e != lemon::INVALID; ++e)
    {
      ++treeDeg[v];
    }
    
    const typename BlockTreeType::EdgeVector& realEdges = blockTree.getRealEdges(v);
    for (typename BlockTreeType::EdgeVectorIt edgeIt = realEdges.begin();
         edgeIt != realEdges.end(); ++edgeIt)
    {
      treeNeg[v] = treeNeg[v] && (score[g.u(*edgeIt)] <= 0 && score[g.v(*edgeIt)] <= 0);
    }
  }
  
  for (typename Tree::NodeIt v(T); v != lemon::INVALID; ++v)
  {
    if (treeNeg[v])
    {
      std::cerr << T.id(v) << " with deg " << treeDeg[v] << " is neg" << std::endl;
    }
  }
  
  std::cerr << "#nodes in block tree: " << blockTree.getNumBlockTreeNodes() << std::endl;
  blockTree.printNodes(std::cerr);
  blockTree.printEdges(std::cerr);
  return 0;
  
  typedef std::set<int> IntSet;
  typedef typename Graph::template NodeMap<IntSet> IntSetMap;
  typedef std::vector<NodeSet> NodeSetVector;
  typedef std::vector<bool> BoolVector;
  typedef std::vector<int> IntVector;
  
  IntEdgeMap biCompMap(g, -1);
  IntSetMap biCompNodeMap(g, IntSet());
  
  int nBiComp = lemon::biNodeConnectedComponents(g, biCompMap);
  
  NodeSetVector biComponents(nBiComp, NodeSet());
  BoolVector negative(nBiComp, true);
  IntVector articulationCount(nBiComp, 0);
  
  // identify negative terminal biconnected components
  for (EdgeIt e(g); e != lemon::INVALID; ++e)
  {
    Node u = g.u(e);
    Node v = g.v(e);
    
    int biCompIdx = biCompMap[e];
    NodeSet& biComponent = biComponents[biCompIdx];
    
    if (biComponent.find(u) == biComponent.end())
    {
      biComponent.insert(u);
      biCompNodeMap[u].insert(biCompIdx);
      
      negative[biCompIdx] = negative[biCompIdx] && (score[u] <= 0);
    }
    if (biComponent.find(v) == biComponent.end())
    {
      biComponent.insert(v);
      biCompNodeMap[v].insert(biCompIdx);
      
      negative[biCompIdx] = negative[biCompIdx] && (score[v] <= 0);
    }
  }
  
  // count #articulation nodes per component
  for (int biCompIdx = 0; biCompIdx < nBiComp; ++biCompIdx)
  {
    const NodeSet& biComponent = biComponents[biCompIdx];
    for (NodeSetIt nodeIt = biComponent.begin(); nodeIt != biComponent.end(); ++nodeIt)
    {
      assert(biCompNodeMap[*nodeIt].size() >= 1);
      if (biCompNodeMap[*nodeIt].size() > 1)
      {
        ++articulationCount[biCompIdx];
      }
    }
  }

//  // let's print them for now
//  std::cerr << "[";
//  for (int biCompIdx = 0; biCompIdx < nBiComp; ++biCompIdx)
//  {
//    std::cerr << " " << biComponents[biCompIdx].size();
//    if (negative[biCompIdx])
//      std::cerr << "*";
//  }
//  std::cerr << " ]" << std::endl;
//  
//  std::cerr << "[";
//  for (int biCompIdx = 0; biCompIdx < nBiComp; ++biCompIdx)
//  {
//    std::cerr << " " << articulationCount[biCompIdx];
//  }
//  std::cerr << " ]" << std::endl << std::endl;
  
  // let's remove negative biconnected components with only one articulation node
  // note that we keep the articulation node!
  int res = 0;
  for (int biCompIdx = 0; biCompIdx < nBiComp; ++biCompIdx)
  {
    if (negative[biCompIdx] && articulationCount[biCompIdx] == 1)
    {
      const NodeSet& biComponent = biComponents[biCompIdx];
      for (NodeSetIt nodeIt = biComponent.begin(); nodeIt != biComponent.end(); ++nodeIt)
      {
        if (biCompNodeMap[*nodeIt].size() == 1)
        {
          remove(g, mapToPre, preOrigNodes, neighbors,
                 nNodes, nArcs, nEdges,
                 degree, degreeVector, *nodeIt);
          ++res;
        }
      }
    }
  }
  
  return res;
}

} // namespace mwcs
} // namespace nina

#endif // NEGBICOMPONENT_H
