/*
 * negcircuit.h
 *
 *  Created on: 7-mar-2014
 *      Author: M. El-Kebir
 */

#ifndef NEGCIRCUIT_H
#define NEGCIRCUIT_H

#include <lemon/core.h>
#include <string>
#include <vector>
#include <set>
#include "unrootedrule.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename WGHT = typename GR::template NodeMap<double> >
class NegCircuit : public UnrootedRule<GR, WGHT>
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

  NegCircuit();
  virtual ~NegCircuit() {}
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

  virtual std::string name() const { return "NegCircuit"; }
};

template<typename GR, typename WGHT>
inline NegCircuit<GR, WGHT>::NegCircuit()
  : Parent()
{
}

template<typename GR, typename WGHT>
inline int NegCircuit<GR, WGHT>::apply(Graph& g,
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
  const NodeSet& nodes = degreeVector[2];
  
  for (NodeSetIt nodeIt = nodes.begin(); nodeIt != nodes.end(); ++nodeIt)
  {
    Node v = *nodeIt;
    if (score[v] <= 0)
    {
      assert(degree[v] == 2);
      Edge e1 = IncEdgeIt(g, v);
      Edge e2 = ++IncEdgeIt(g, v);

      Node u = g.oppositeNode(v, e1);
      Node w = g.oppositeNode(v, e2);
      if (arcLookUp(u, w) != lemon::INVALID)
      {
        remove(g, mapToPre, preOrigNodes, neighbors,
               nNodes, nArcs, nEdges,
               degree, degreeVector, v);
        return 1;
      }
    }
  }
              
  return 0;
}

} // namespace mwcs
} // namespace nina

#endif // NEGCIRCUIT_H
