/*
 * solverimpl.h
 *
 *  Created on: 25-aug-2014
 *      Author: M. El-Kebir
 */

#ifndef SOLVERIMPL_H
#define SOLVERIMPL_H

#include <set>

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class SolverImpl
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;
  
  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  
  typedef std::set<Node> NodeSet;
  typedef typename NodeSet::const_iterator NodeSetIt;
  typedef typename NodeSet::iterator NodeSetNonConstIt;
  
  typedef MwcsGraph<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> MwcsGraphType;
  
protected:
  const MwcsGraphType* _pMwcsGraph;
  
public:
  SolverImpl()
    : _pMwcsGraph(NULL)
  {
  }
  
  virtual ~SolverImpl()
  {
  }
  
  virtual bool solve(double& score, BoolNodeMap& solutionMap, NodeSet& solutionSet) = 0;  
};

} // namespace mwcs
} // namespace nina

#endif // SOLVERIMPL_H
