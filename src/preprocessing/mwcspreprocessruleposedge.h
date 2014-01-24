/*
 * mwcspreprocessruleposedge.h
 *
 *  Created on: 14-jan-2013
 *      Author: M. El-Kebir
 */

#ifndef MWCSPREPROCESSRULEPOSEDGE_H
#define MWCSPREPROCESSRULEPOSEDGE_H

#include <lemon/core.h>
#include <string>
#include <vector>
#include <set>
#include "mwcspreprocessrule.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename WGHT = typename GR::template NodeMap<double> >
class MwcsPreprocessRulePosEdge : public MwcsPreprocessRule<GR, WGHT>
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

  MwcsPreprocessRulePosEdge();
  virtual ~MwcsPreprocessRulePosEdge() {}
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
                    DegreeNodeSetVector& degreeVector);

  virtual std::string name() const { return "PosEdge"; }
};

template<typename GR, typename WGHT>
inline MwcsPreprocessRulePosEdge<GR, WGHT>::MwcsPreprocessRulePosEdge()
  : Parent()
{
}

template<typename GR, typename WGHT>
inline int MwcsPreprocessRulePosEdge<GR, WGHT>::apply(Graph& g,
                                                      const ArcLookUpType& arcLookUp,
                                                      LabelNodeMap& label,
                                                      WeightNodeMap& score,
                                                      NodeMap& mapToPre,
                                                      NodeSetMap& preOrigNodes,
                                                      int& nNodes,
                                                      int& nArcs,
                                                      int& nEdges,
                                                      DegreeNodeMap& degree,
                                                      DegreeNodeSetVector& degreeVector)
{
  for (EdgeIt e(g); e != lemon::INVALID; ++e)
  {
    Node u = g.u(e);
    Node v = g.v(e);

    if (score[u] >= 0 && score[v] >= 0)
    {
      merge(g, arcLookUp, label, score,
            mapToPre, preOrigNodes,
            nNodes, nArcs, nEdges,
            degree, degreeVector, u, v);
      return 1;
    }
  }
  return 0;
}

} // namespace mwcs
} // namespace nina

#endif // MWCSPREPROCESSRULEPOSEDGE_H
