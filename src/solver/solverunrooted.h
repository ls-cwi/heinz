/*
 * solverunrooted.h
 *
 *  Created on: 7-aug-2012
 *     Authors: C.I. Bucur, M. El-Kebir
 */

#ifndef SOLVERUNROOTED_H
#define SOLVERUNROOTED_H

#include "solver.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class SolverUnrooted : public virtual Solver<GR, NWGHT, NLBL, EWGHT>
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
  SolverUnrooted(const MwcsGraphType& mwcsGraph)
    : Parent(mwcsGraph)
  {
  }

  virtual ~SolverUnrooted()
  {
  }
};


} // namespace mwcs
} // namespace nina

#endif // MWCSSOLVERUNROOTED_H
