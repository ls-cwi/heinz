/*
 * mwcspreprocessrulenegehub.h
 *
 *  Created on: 13-may-2013
 *      Author: M. El-Kebir
 */

#ifndef MWCSPREPROCESSRULENEGHUB_H
#define MWCSPREPROCESSRULENEGHUB_H

#include <lemon/core.h>
#include <lemon/adaptors.h>
#include <lemon/connectivity.h>
#include <string>
#include <vector>
#include <set>
#include "mwcspreprocessrule.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename WGHT = typename GR::template NodeMap<double> >
class MwcsPreprocessRuleNegHub : public MwcsPreprocessRule<GR, WGHT>
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

  MwcsPreprocessRuleNegHub();
  virtual ~MwcsPreprocessRuleNegHub() {}
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

  virtual std::string name() const { return "NegHub"; }

private:
  typedef lemon::FilterNodes<const Graph> Subgraph;

  void determinePosNeighbors(const Graph& g,
                             const WeightNodeMap& score,
                             const LabelNodeMap& label,
                             Node negHub,
                             NodeSet& posNeighbors);
};

template<typename GR, typename WGHT>
inline void MwcsPreprocessRuleNegHub<GR, WGHT>::determinePosNeighbors(const Graph& g,
                                                                      const WeightNodeMap& score,
                                                                      const LabelNodeMap& label,
                                                                      Node negHub,
                                                                      NodeSet& posNeighbors)
{
  posNeighbors.clear();

  BoolNodeMap filter(g, true);
  for (IncEdgeIt e(g, negHub); e != lemon::INVALID; ++e)
  {
    Node posNeighbor = g.oppositeNode(negHub, e);
    if (score[posNeighbor] >= 0)
    {
      posNeighbors.insert(posNeighbor);
    }
  }

  if (posNeighbors.size() < 2)
  {
    posNeighbors.clear();
    return;
  }

  IntNodeMap comp(g);
  filter[negHub] = false;
  Subgraph subG(g, filter);
  int nComp = lemon::connectedComponents(subG, comp);

  if (nComp < 2)
  {
    posNeighbors.clear();
  }
  //else
  //{
  //  std::cout << "// NegHub: " << label[negHub] << "\t" << score[negHub] << std::endl;
  //  for (NodeSetIt posNodeIt = posNeighbors.begin();
  //       posNodeIt != posNeighbors.end(); posNodeIt++)
  //  {
  //    std::cout << "// " << label[*posNodeIt] << "\t" << score[*posNodeIt] << std::endl;
  //  }
  //}
}

template<typename GR, typename WGHT>
inline MwcsPreprocessRuleNegHub<GR, WGHT>::MwcsPreprocessRuleNegHub()
{
}

template<typename GR, typename WGHT>
inline int MwcsPreprocessRuleNegHub<GR, WGHT>::apply(Graph& g,
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
  int res = 0;

  NodeSet bestPosNeighbors;
  Node bestNegHub = lemon::INVALID;
  double bestGain = 0;

  for (NodeIt u(g); u != lemon::INVALID; ++u)
  {
    double score_u = score[u];

    if (score_u < 0)
    {
      NodeSet posNeighbors;
      determinePosNeighbors(g, score, label, u, posNeighbors);

      double maxScorePosNeighbor = 0;
      double scorePosNeighbors = 0;
      for (NodeSetIt nodeIt = posNeighbors.begin();
           nodeIt != posNeighbors.end(); nodeIt++)
      {
        Node v = *nodeIt;
        double score_v = score[v];
        assert(score_v >= 0);

        scorePosNeighbors += score_v;
        if (score_v > maxScorePosNeighbor)
          maxScorePosNeighbor = score_v;
      }

      double gain = scorePosNeighbors + score_u;

      if (gain > maxScorePosNeighbor && gain > bestGain)
      {
        bestGain = gain;
        bestPosNeighbors = posNeighbors;
        bestNegHub = u;
      }
    }
  }

  if (bestNegHub != lemon::INVALID)
  {
    std::cout << "// NegHub: " << label[bestNegHub] << "\t" << score[bestNegHub] << std::endl;
  }

  for (NodeSetIt posNodeIt = bestPosNeighbors.begin();
       posNodeIt != bestPosNeighbors.end(); posNodeIt++)
  {
    std::cout << "// " << label[*posNodeIt] << "\t" << score[*posNodeIt] << std::endl;

    bestNegHub = merge(g, arcLookUp, label, score,
                       mapToPre, preOrigNodes,
                       nNodes, nArcs, nEdges,
                       degree, degreeVector, bestNegHub, *posNodeIt, LB);
    res++;
  }

  return res;
}

} // namespace mwcs
} // namespace nina


#endif // MWCSPREPROCESSRULENEGHUB_H
