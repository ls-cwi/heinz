/*
 * solverrooted.h
 *
 *  Created on: 26-jul-2014
 *      Author: M. El-Kebir, G. W. Klau
 */

#ifndef SOLVERROOTED_H
#define SOLVERROOTED_H

#include "solver.h"
#include "impl/solverrootedimpl.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class SolverRooted : public Solver<GR, NWGHT, NLBL, EWGHT>
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;

  typedef Solver<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> Parent;
  typedef SolverRootedImpl<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> SolverRootedImplType;
  typedef typename Parent::MwcsGraphType MwcsGraphType;
  typedef typename Parent::NodeSet NodeSet;
  typedef typename Parent::NodeSetIt NodeSetIt;
  typedef typename Parent::NodeVector NodeVector;
  typedef typename Parent::NodeVectorIt NodeVectorIt;

  TEMPLATE_GRAPH_TYPEDEFS(Graph);

  using Parent::_score;
  using Parent::_scoreUB;
  using Parent::_pSolutionMap;
  using Parent::_solutionSet;

public:
  SolverRooted(SolverRootedImplType* pImpl)
    : _pImpl(pImpl)
  {
  }

  virtual ~SolverRooted()
  {
    delete _pImpl;
  }

  virtual bool solve(const MwcsGraphType& mwcsGraph,
                     const NodeSet& rootNodes)
  {
    delete _pSolutionMap;
    _pSolutionMap = new BoolNodeMap(mwcsGraph.getGraph(), false);

    _pImpl->init(mwcsGraph, rootNodes);
    return _pImpl->solve(_score, _scoreUB, *_pSolutionMap, _solutionSet);
  }

protected:
  SolverRootedImplType* _pImpl;
};

} // namespace mwcs
} // namespace nina

#endif // SOLVERROOTED_H
