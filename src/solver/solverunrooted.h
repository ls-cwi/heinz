/*
 * solverunrooted.h
 *
 *  Created on: 7-aug-2012
 *     Authors: C.I. Bucur, M. El-Kebir, G. W. Klau
 */

#ifndef SOLVERUNROOTED_H
#define SOLVERUNROOTED_H

#include "solver.h"
#include "impl/solverunrootedimpl.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class SolverUnrooted : public Solver<GR, NWGHT, NLBL, EWGHT>
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;

  typedef Solver<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> Parent;
  typedef SolverUnrootedImpl<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> SolverUnrootedImplType;
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
  SolverUnrooted(SolverUnrootedImplType* pImpl)
    : _pImpl(pImpl)
  {
  }

  virtual ~SolverUnrooted()
  {
    delete _pImpl;
  }

  virtual bool solve(const MwcsGraphType& mwcsGraph)
  {
    delete _pSolutionMap;
    _pSolutionMap = new BoolNodeMap(mwcsGraph.getGraph(), false);

    _pImpl->init(mwcsGraph);
    return _pImpl->solve(_score, _scoreUB, *_pSolutionMap, _solutionSet);
  }

protected:
  SolverUnrootedImplType* _pImpl;
};


} // namespace mwcs
} // namespace nina

#endif // MWCSSOLVERUNROOTED_H
