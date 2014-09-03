/*
 * posedge.h
 *
 *  Created on: 14-jan-2013
 *      Author: M. El-Kebir
 */

#ifndef POSEDGE_H
#define POSEDGE_H

#include <lemon/core.h>
#include <string>
#include <vector>
#include <set>
#include "rule.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename WGHT = typename GR::template NodeMap<double> >
class PosEdge : public Rule<GR, WGHT>
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

  PosEdge();
  virtual ~PosEdge() {}
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

  virtual std::string name() const { return "PosEdge"; }
};

template<typename GR, typename WGHT>
inline PosEdge<GR, WGHT>::PosEdge()
  : Parent()
{
}

template<typename GR, typename WGHT>
inline int PosEdge<GR, WGHT>::apply(Graph& g,
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
  for (EdgeIt e(g); e != lemon::INVALID; ++e)
  {
    Node u = g.u(e);
    Node v = g.v(e);

    if (score[u] >= 0 && score[v] >= 0
        && rootNodes.find(u) == rootNodes.end()
        && rootNodes.find(v) == rootNodes.end())
    {
      merge(g, arcLookUp, label, score,
            mapToPre, preOrigNodes, neighbors,
            nNodes, nArcs, nEdges,
            degree, degreeVector, u, v, LB);
      return 1;
    }
  }
  return 0;
}

} // namespace mwcs
} // namespace nina

#endif // POSEDGE_H
