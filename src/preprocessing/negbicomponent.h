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
#include "rule.h"
#include "solver/blockcuttree.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename WGHT = typename GR::template NodeMap<double> >
class NegBiComponent : public Rule<GR, WGHT>
{
public:
  typedef GR Graph;
  typedef WGHT WeightNodeMap;
  typedef Rule<GR, WGHT> Parent;
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
  virtual int apply(Graph& g, const NodeSet& rootNodes,
                    const ArcLookUpType& arcLookUp,
                    LabelNodeMap& label,
                    WeightNodeMap& score,
                    NodeSetMap& mapToPre,
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
  typedef BlockCutTree<Graph> BlockCutTreeType;
};

template<typename GR, typename WGHT>
inline NegBiComponent<GR, WGHT>::NegBiComponent()
  : Parent()
{
}

template<typename GR, typename WGHT>
inline int NegBiComponent<GR, WGHT>::apply(Graph& g,
                                           const NodeSet& rootNodes,
                                           const ArcLookUpType& arcLookUp,
                                           LabelNodeMap& label,
                                           WeightNodeMap& score,
                                           NodeSetMap& mapToPre,
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
  
  typedef typename BlockCutTreeType::Tree Tree;
  typedef typename Tree::template NodeMap<int> IntTreeNodeMap;
  typedef typename Tree::template BlueNodeMap<bool> BoolTreeNodeMap;

  const Tree& T = blockCutTree.getBlockCutTree();
  
  // determine block negativity
  BoolTreeNodeMap treeNeg(T, true);
  
  for (typename Tree::BlueNodeIt v(T); v != lemon::INVALID; ++v)
  {
    const typename BlockCutTreeType::EdgeVector& realEdges = blockCutTree.getRealEdges(v);
    for (typename BlockCutTreeType::EdgeVectorIt edgeIt = realEdges.begin();
         edgeIt != realEdges.end(); ++edgeIt)
    {
      treeNeg[v] = treeNeg[v] && (score[g.u(*edgeIt)] <= 0 && score[g.v(*edgeIt)] <= 0);
    }
  }

  // let's remove negative biconnected components with only one articulation node
  // note that we keep the articulation node!
  int res = 0;
  for (typename Tree::BlueNodeIt v(T); v != lemon::INVALID; ++v)
  {
    if (treeNeg[v] && blockCutTree.getDegree(v) == 1)
    {
      Node articulationNode = blockCutTree.getArticulationPoint(T.redNode(typename Tree::IncEdgeIt(T, v)));
      const typename BlockCutTreeType::EdgeVector& realEdges = blockCutTree.getRealEdges(v);
      
      // determine nodes to remove
      NodeSet nodesToRemove;
      for (typename BlockCutTreeType::EdgeVectorIt edgeIt = realEdges.begin(); edgeIt != realEdges.end(); ++edgeIt)
      {
        if (g.u(*edgeIt) != articulationNode)
        {
          nodesToRemove.insert(g.u(*edgeIt));
        }
        if (g.v(*edgeIt) != articulationNode)
        {
          nodesToRemove.insert(g.v(*edgeIt));
        }
      }
      
      for (NodeSetIt nodeIt = nodesToRemove.begin(); nodeIt != nodesToRemove.end(); ++nodeIt)
      {
        remove(g, mapToPre, preOrigNodes, neighbors,
               nNodes, nArcs, nEdges,
               degree, degreeVector, *nodeIt);
        ++res;
      }
    }
  }
  
  return res;
}

} // namespace mwcs
} // namespace nina

#endif // NEGBICOMPONENT_H
