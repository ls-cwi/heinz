/*
 * mwcspreprocessrootrule.h
 *
 *  Created on: 21-jan-2013
 *      Author: M. El-Kebir
 */

#ifndef MWCSPREPROCESSROOTRULE_H
#define MWCSPREPROCESSROOTRULE_H

#include "mwcspreprocessrulebase.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename WGHT = typename GR::template NodeMap<double> >
class MwcsPreprocessRootRule : public MwcsPreprocessRuleBase<GR, WGHT>
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
  MwcsPreprocessRootRule();
  virtual ~MwcsPreprocessRootRule() {}
  virtual int apply(Graph& g,
                    Node& root,
                    const ArcLookUpType& arcLookUp,
                    LabelNodeMap& label,
                    WeightNodeMap& score,
                    NodeMap& mapToPre,
                    NodeSetMap& preOrigNodes,
                    int& nNodes,
                    int& nArcs,
                    int& nEdges,
                    DegreeNodeMap& degree,
                    DegreeNodeSetVector& degreeVector) = 0;
  virtual std::string name() const = 0;
};

template<typename GR, typename WGHT>
inline MwcsPreprocessRootRule<GR, WGHT>::MwcsPreprocessRootRule()
  : Parent()
{
}

} // namespace mwcs
} // namespace nina

#endif // MWCSPREPROCESSROOTRULE_H
