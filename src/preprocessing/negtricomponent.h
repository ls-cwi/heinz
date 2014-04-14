/*
 * negtricomponent.h
 *
 *  Created on: 30-mar-2014
 *      Author: M. El-Kebir
 */

#ifndef NEGTRICOMPONENT_H
#define NEGTRICOMPONENT_H

#include <lemon/core.h>
#include <lemon/connectivity.h>
#include <string>
#include <vector>
#include <set>
#include "unrootedrule.h"
#include "solver/spqrtree.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename WGHT = typename GR::template NodeMap<double> >
class NegTriComponent : public UnrootedRule<GR, WGHT>
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
  
  typedef typename nina::SpqrTree<Graph> SpqrType;

  TEMPLATE_GRAPH_TYPEDEFS(Graph);

  using Parent::remove;
  using Parent::merge;

  NegTriComponent();
  virtual ~NegTriComponent() {}
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

  virtual std::string name() const { return "NegTriComponent"; }
};

template<typename GR, typename WGHT>
inline NegTriComponent<GR, WGHT>::NegTriComponent()
  : Parent()
{
}

template<typename GR, typename WGHT>
inline int NegTriComponent<GR, WGHT>::apply(Graph& g,
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
  // todo do this on the biconnected components
  SpqrType spqr(g);
  spqr.run();
  
  return 0;
}

} // namespace mwcs
} // namespace nina

#endif // NEGTRICOMPONENT_H
