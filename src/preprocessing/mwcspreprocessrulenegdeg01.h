/*
 * mwcspreprocessrulenegdeg01.h
 *
 *  Created on: 14-jan-2013
 *      Author: M. El-Kebir
 */

#ifndef MWCSPREPROCESSRULENEGDEG01_H
#define MWCSPREPROCESSRULENEGDEG01_H

#include <lemon/core.h>
#include <string>
#include <vector>
#include <set>
#include "mwcspreprocessrule.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename WGHT = typename GR::template NodeMap<double> >
class MwcsPreprocessRuleNegDeg01 : public MwcsPreprocessRule<GR, WGHT>
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

  MwcsPreprocessRuleNegDeg01();
  virtual ~MwcsPreprocessRuleNegDeg01() {}
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

  virtual std::string name() const { return "NegDeg01"; }

protected:
  int apply(Graph& g,
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
            int d);
};

template<typename GR, typename WGHT>
inline MwcsPreprocessRuleNegDeg01<GR, WGHT>::MwcsPreprocessRuleNegDeg01()
  : Parent()
{
}

template<typename GR, typename WGHT>
inline int MwcsPreprocessRuleNegDeg01<GR, WGHT>::apply(Graph& g,
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
  return apply(g, arcLookUp, label, score, mapToPre, preOrigNodes, nNodes, nArcs, nEdges, degree, degreeVector, 0)
      + apply(g, arcLookUp, label, score, mapToPre, preOrigNodes, nNodes, nArcs, nEdges, degree, degreeVector, 1);
}

template<typename GR, typename WGHT>
inline int MwcsPreprocessRuleNegDeg01<GR, WGHT>::apply(Graph& g,
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
                                                       int d)
{
  if (static_cast<int>(degreeVector.size()) <= d)
    return 0;

  const NodeSet& nodes = degreeVector[d];

  for (NodeSetIt nodeIt = nodes.begin(); nodeIt != nodes.end();)
  {
    NodeSetIt nextNodeIt = nodeIt;
    nextNodeIt++;

    // remove if negative
    if (score[*nodeIt] < 0)
    {
      remove(g, mapToPre, preOrigNodes,
             nNodes, nArcs, nEdges,
             degree, degreeVector, *nodeIt);
      return 1;
    }

    nodeIt = nextNodeIt;
  }

  return 0;
}

} // namespace mwcs
} // namespace nina

#endif // MWCSPREPROCESSRULENEGDEG01_H
