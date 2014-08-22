/*
 * solverrooted.h
 *
 *  Created on: 26-jul-2014
 *      Author: M. El-Kebir
 */

#ifndef SOLVERROOTED_H
#define SOLVERROOTED_H

#include "solver.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class SolverRooted : public virtual Solver<GR, NWGHT, NLBL, EWGHT>
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;
  
  typedef Solver<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> Parent;
  typedef typename Parent::MwcsGraphType MwcsGraphType;
  typedef typename Parent::NodeSet NodeSet;
  typedef typename Parent::NodeSetIt NodeSetIt;
  typedef typename Parent::NodeVector NodeVector;
  typedef typename Parent::NodeVectorIt NodeVectorIt;
  
public:
  SolverRooted(const MwcsGraphType& mwcsGraph,
               const NodeSet& rootNodes)
    : Parent(mwcsGraph)
    , _rootNodes(rootNodes)
  {
    assert(rootNodes.size() > 0);
  }
  
  virtual ~SolverRooted()
  {
  }
  
protected:
  NodeSet _rootNodes;
};
  
} // namespace mwcs
} // namespace nina

#endif // SOLVERROOTED_H
