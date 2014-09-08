/*
 * negedge.h
 *
 *  Created on: 14-jan-2013
 *      Author: M. El-Kebir
 */

#ifndef NEGEDGE_H
#define NEGEDGE_H

#include <lemon/core.h>
#include <string>
#include <vector>
#include <set>
#include "rule.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename WGHT = typename GR::template NodeMap<double> >
class NegEdge : public Rule<GR, WGHT>
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

  TEMPLATE_GRAPH_TYPEDEFS(Graph);

  using Parent::remove;
  using Parent::merge;

  NegEdge();
  virtual ~NegEdge() {}
  virtual int apply(Graph& g,
                    const NodeSet& rootNodes,
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

  virtual std::string name() const { return "NegEdge"; }
};

template<typename GR, typename WGHT>
inline NegEdge<GR, WGHT>::NegEdge()
  : Parent()
{
}

template<typename GR, typename WGHT>
inline int NegEdge<GR, WGHT>::apply(Graph& g,
                                    const NodeSet& rootNodes,
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
  int res = 0;

  for (EdgeIt e(g); e != lemon::INVALID; ++e)
  {
    Node u = g.u(e);
    Node v = g.v(e);

    if (score[u] <= 0 && score[v] <= 0 && degree[u] == 2 && degree[v] == 2)
    {
      // don't merge if both u and v are root nodes
      // otherwise ensure that root node is kept
      if (rootNodes.find(u) == rootNodes.end() && rootNodes.find(v) == rootNodes.end())
      {
        res++;
        merge(g, label, score,
              mapToPre, preOrigNodes, neighbors,
              nNodes, nArcs, nEdges,
              degree, degreeVector, u, v, LB);
      }
//      else if (rootNodes.find(u) != rootNodes.end() && rootNodes.find(v) == rootNodes.end())
//      {
//        res++;
//        merge(g, arcLookUp, label, score,
//              mapToPre, preOrigNodes, neighbors,
//              nNodes, nArcs, nEdges,
//              degree, degreeVector, v, u, LB);
//      }
    }
  }

  return res;
}

} // namespace mwcs
} // namespace nina

#endif // NEGEDGE_H
