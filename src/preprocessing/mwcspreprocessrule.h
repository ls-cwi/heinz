/*
 * mwcspreprocessrule.h
 *
 *  Created on: 12-jan-2013
 *      Author: M. El-Kebir
 */

#ifndef MWCSPREPROCESSRULE_H
#define MWCSPREPROCESSRULE_H

#include "mwcspreprocessrulebase.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename WGHT = typename GR::template NodeMap<double> >
class MwcsPreprocessRule : public MwcsPreprocessRuleBase<GR, WGHT>
{
public:
  typedef GR Graph;
  typedef WGHT WeightNodeMap;
  typedef MwcsPreprocessRuleBase<GR, WGHT> Parent;

  TEMPLATE_GRAPH_TYPEDEFS(Graph);

  typedef typename Parent::LabelNodeMap LabelNodeMap;
  typedef typename Parent::DegreeNodeMap DegreeNodeMap;
  typedef typename Parent::NodeMap NodeMap;
  typedef typename Parent::NodeSet NodeSet;
  typedef typename Parent::NodeSetIt NodeSetIt;
  typedef typename Parent::NodeSetMap NodeSetMap;
  typedef typename Parent::DegreeNodeSetVector DegreeNodeSetVector;
  typedef typename Parent::ArcLookUpType ArcLookUpType;

  using Parent::merge;
  using Parent::remove;

public:
  MwcsPreprocessRule();
  virtual ~MwcsPreprocessRule() {}
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
                    double& LB) = 0;
  virtual std::string name() const = 0;
};

template<typename GR, typename WGHT>
inline MwcsPreprocessRule<GR, WGHT>::MwcsPreprocessRule()
  : Parent()
{
}

} // namespace mwcs
} // namespace nina

#endif // MWCSPREPROCESSRULE_H
