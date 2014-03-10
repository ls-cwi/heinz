/*
 * posdeg01.h
 *
 *  Created on: 8-mar-2014
 *      Author: M. El-Kebir
 */

#ifndef POSDEG01_H
#define POSDEG01_H

#include <lemon/core.h>
#include <string>
#include <vector>
#include <set>
#include "mwcspreprocessrule.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename WGHT = typename GR::template NodeMap<double> >
class PosDeg01 : public MwcsPreprocessRule<GR, WGHT>
{
public:
  typedef GR Graph;
  typedef WGHT WeightNodeMap;
  typedef MwcsPreprocessRule<GR, WGHT> Parent;
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

  PosDeg01();
  virtual ~PosDeg01() {}
  virtual int apply(Graph& g,
                    const ArcLookUpType& arcLookUp,
                    LabelNodeMap& label,
                    WeightNodeMap& score,
                    NodeMap& mapToPre,
                    NodeSetMap& preOrigNodes,
                    int& nNodes,
                    int& nArcs,
                    int& nEdges,
                    DegreeNodeMap& degree,
                    DegreeNodeSetVector& degreeVector,
                    double& LB);

  virtual std::string name() const { return "Unrooted - PosDeg01"; }
};

template<typename GR, typename WGHT>
inline PosDeg01<GR, WGHT>::PosDeg01()
  : Parent()
{
}

template<typename GR, typename WGHT>
inline int PosDeg01<GR, WGHT>::apply(Graph& g,
                                     const ArcLookUpType& arcLookUp,
                                     LabelNodeMap& label,
                                     WeightNodeMap& score,
                                     NodeMap& mapToPre,
                                     NodeSetMap& preOrigNodes,
                                     int& nNodes,
                                     int& nArcs,
                                     int& nEdges,
                                     DegreeNodeMap& degree,
                                     DegreeNodeSetVector& degreeVector,
                                     double& LB)
{
  // positive deg 0 nodes smaller than LB are to be removed
  const NodeSet& nodes0 = degreeVector[0];
  for (NodeSetIt nodeIt = nodes0.begin(); nodeIt != nodes0.end(); ++nodeIt)
  {
    if (0 <= score[*nodeIt] && score[*nodeIt] < LB)
    {
      remove(g, mapToPre, preOrigNodes,
             nNodes, nArcs, nEdges,
             degree, degreeVector, *nodeIt);
      return 1;
    }
  }
  
  const NodeSet& nodes1 = degreeVector[1];
  for (NodeSetIt nodeIt = nodes1.begin(); nodeIt != nodes1.end(); ++nodeIt)
  {
    Node v = *nodeIt;
    if (0 <= score[v] && score[v] < LB)
    {
      Node u = g.oppositeNode(v, IncEdgeIt(g, v));
      
      merge(g, arcLookUp, label, score,
            mapToPre, preOrigNodes,
            nNodes, nArcs, nEdges,
            degree, degreeVector, u, v, LB);
      
      return 1;
    }
  }
              
  return 0;
}

} // namespace mwcs
} // namespace nina

#endif // POSDEG01_H
