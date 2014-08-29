/*
 * negdeg01.h
 *
 *  Created on: 14-jan-2013
 *      Author: M. El-Kebir
 */

#ifndef NEGDEG01_H
#define NEGDEG01_H

#include <lemon/core.h>
#include <string>
#include <vector>
#include <set>
#include "rule.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename WGHT = typename GR::template NodeMap<double> >
class NegDeg01 : public Rule<GR, WGHT>
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

  NegDeg01();
  virtual ~NegDeg01() {}
  virtual int apply(Graph& g,
                    const NodeSet& rootNodes,
                    const ArcLookUpType& arcLookUp,
                    LabelNodeMap& label,
                    WeightNodeMap& score,
                    IntNodeMap& comp,
                    NodeSetMap& mapToPre,
                    NodeSetMap& preOrigNodes,
                    NodeSetMap& neighbors,
                    int& nNodes,
                    int& nArcs,
                    int& nEdges,
                    int& nComponents,
                    DegreeNodeMap& degree,
                    DegreeNodeSetVector& degreeVector,
                    double& LB);

  virtual std::string name() const { return "NegDeg01"; }

protected:
  int apply(Graph& g,
            const NodeSet& rootNodes,
            const ArcLookUpType& arcLookUp,
            LabelNodeMap& label,
            WeightNodeMap& score,
            IntNodeMap& comp,
            NodeSetMap& mapToPre,
            NodeSetMap& preOrigNodes,
            NodeSetMap& neighbors,
            int& nNodes,
            int& nArcs,
            int& nEdges,
            int& nComponents,
            DegreeNodeMap& degree,
            DegreeNodeSetVector& degreeVector,
            int d);
};

template<typename GR, typename WGHT>
inline NegDeg01<GR, WGHT>::NegDeg01()
  : Parent()
{
}

template<typename GR, typename WGHT>
inline int NegDeg01<GR, WGHT>::apply(Graph& g,
                                     const NodeSet& rootNodes,
                                     const ArcLookUpType& arcLookUp,
                                     LabelNodeMap& label,
                                     WeightNodeMap& score,
                                     IntNodeMap& comp,
                                     NodeSetMap& mapToPre,
                                     NodeSetMap& preOrigNodes,
                                     NodeSetMap& neighbors,
                                     int& nNodes,
                                     int& nArcs,
                                     int& nEdges,
                                     int& nComponents,
                                     DegreeNodeMap& degree,
                                     DegreeNodeSetVector& degreeVector,
                                     double& LB)
{
  return apply(g, rootNodes, arcLookUp,
               label, score, comp,
               mapToPre, preOrigNodes, neighbors,
               nNodes, nArcs, nEdges, nComponents,
               degree, degreeVector, 0)
      + apply(g, rootNodes, arcLookUp,
              label, score, comp,
              mapToPre, preOrigNodes, neighbors,
              nNodes, nArcs, nEdges, nComponents,
              degree, degreeVector, 1);
}

template<typename GR, typename WGHT>
inline int NegDeg01<GR, WGHT>::apply(Graph& g,
                                     const NodeSet& rootNodes,
                                     const ArcLookUpType& arcLookUp,
                                     LabelNodeMap& label,
                                     WeightNodeMap& score,
                                     IntNodeMap& comp,
                                     NodeSetMap& mapToPre,
                                     NodeSetMap& preOrigNodes,
                                     NodeSetMap& neighbors,
                                     int& nNodes,
                                     int& nArcs,
                                     int& nEdges,
                                     int& nComponents,
                                     DegreeNodeMap& degree,
                                     DegreeNodeSetVector& degreeVector,
                                     int d)
{
  if (static_cast<int>(degreeVector.size()) <= d)
    return 0;

  const NodeSet& nodes = degreeVector[d];

  for (NodeSetIt nodeIt = nodes.begin(); nodeIt != nodes.end();)
  {
    NodeSetIt nextNodeIt = nodeIt;
    nextNodeIt++;
    
    Node v = *nodeIt;
    // remove if negative and not the root node
    if (score[v] < 0 && rootNodes.find(v) == rootNodes.end())
    {
      remove(g, comp, mapToPre, preOrigNodes, neighbors,
             nNodes, nArcs, nEdges, nComponents,
             degree, degreeVector, v);
      return 1;
    }

    nodeIt = nextNodeIt;
  }

  return 0;
}

} // namespace mwcs
} // namespace nina

#endif // NEGDEG01_H
