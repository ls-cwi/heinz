/*
 * posdiamond.h
 *
 *  Created on: 10-mar-2014
 *      Author: M. El-Kebir
 */

#ifndef POSDIAMOND_H
#define POSDIAMOND_H

#include <lemon/core.h>
#include <string>
#include <vector>
#include <set>
#include "unrootedrule.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename WGHT = typename GR::template NodeMap<double> >
class PosDiamond : public UnrootedRule<GR, WGHT>
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

  PosDiamond();
  virtual ~PosDiamond() {}
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

  virtual std::string name() const { return "PosDiamond"; }
};

template<typename GR, typename WGHT>
inline PosDiamond<GR, WGHT>::PosDiamond()
  : Parent()
{
}

template<typename GR, typename WGHT>
inline int PosDiamond<GR, WGHT>::apply(Graph& g,
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
  typedef std::pair<double, Node> WeightNodePair;
  typedef std::set<WeightNodePair> WeightNodePairSet;
  typedef typename WeightNodePairSet::const_iterator WeightNodePairSetIt;
  
  typedef std::pair<Node, Node> NodePair;
  typedef std::map<NodePair, WeightNodePairSet> NodePairMap;
  typedef typename NodePairMap::const_iterator NodePairMapIt;

  if (degreeVector.size() <= 2)
  {
    // nothing to remove, there are no degree 2 nodes
    return 0;
  }
  
  NodePairMap posMap;
  const NodeSet& nodes = degreeVector[2];
  for (NodeSetIt nodeIt = nodes.begin(); nodeIt != nodes.end(); ++nodeIt)
  {
    Node v = *nodeIt;
    
    assert(degree[v] == 2);
    Edge e1 = IncEdgeIt(g, v);
    Edge e2 = ++IncEdgeIt(g, v);
    
    Node u = g.oppositeNode(v, e1);
    Node w = g.oppositeNode(v, e2);
    if (score[v] >= 0 && score[u] <= 0 && score[w] <= 0)
    {
      if (u < w)
      {
        posMap[std::make_pair(u, w)].insert(std::make_pair(score[v], v));
      }
      else
      {
        posMap[std::make_pair(w, u)].insert(std::make_pair(score[v], v));
      }
    }
  }
  
  int res = 0;
  for (NodePairMapIt it = posMap.begin(); it != posMap.end(); ++it)
  {
    const WeightNodePairSet& set = it->second;
    Node u = it->first.first;
    Node w = it->first.first;
    
    int set_size = static_cast<int>(set.size());
    if (degree[u] == set_size && degree[w] != set_size)
    {
      remove(g, mapToPre, preOrigNodes, neighbors,
             nNodes, nArcs, nEdges,
             degree, degreeVector, u);
      ++res;
    }
    else if (degree[u] != set_size && degree[w] == set_size)
    {
      remove(g, mapToPre, preOrigNodes, neighbors,
             nNodes, nArcs, nEdges,
             degree, degreeVector, w);
      ++res;
    }
    else if (degree[u] == set_size && degree[w] == set_size)
    {
      if (score[u] < score[w])
      {
        remove(g, mapToPre, preOrigNodes, neighbors,
               nNodes, nArcs, nEdges,
               degree, degreeVector, u);
        ++res;
      }
      else
      {
        remove(g, mapToPre, preOrigNodes, neighbors,
               nNodes, nArcs, nEdges,
               degree, degreeVector, w);
        ++res;
      }
    }
  }

  return res;
}

} // namespace mwcs
} // namespace nina

#endif // POSDIAMOND_H
