/*
 * solverrootedimpl.h
 *
 *  Created on: 25-aug-2014
 *      Author: M. El-Kebir, G. W. Klau
 */

#ifndef SOLVERROOTEDIMPL_H
#define SOLVERROOTEDIMPL_H

#include <set>
#include "solverimpl.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class SolverRootedImpl : public SolverImpl<GR, NWGHT, NLBL, EWGHT>
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;
  
  typedef SolverImpl<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> Parent;

  typedef typename Parent::MwcsGraphType MwcsGraphType;
  typedef typename Parent::NodeSet NodeSet;
  typedef typename Parent::NodeSetIt NodeSetIt;
  typedef typename Parent::NodeSetNonConstIt NodeSetNonConstIt;
  
  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  
  using Parent::_pMwcsGraph;

protected:
  NodeSet _rootNodes;
  
public:
  SolverRootedImpl()
    : Parent()
    , _rootNodes()
  {
  }
  
  virtual ~SolverRootedImpl()
  {
  }
  
  virtual void init(const MwcsGraphType& mwcsGraph,
                    const NodeSet& rootNodes)
  {
    _pMwcsGraph = &mwcsGraph;
    _rootNodes = rootNodes;
  }  
};

} // namespace mwcs
} // namespace nina

#endif // SOLVERROOTEDIMPL_H
