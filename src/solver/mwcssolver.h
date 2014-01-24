/*
 * mwcssolver.h
 *
 *  Created on: 7-aug-2012
 *     Authors: C.I. Bucur, M. El-Kebir
 */

#ifndef MWCSSOLVER_H
#define MWCSSOLVER_H

#include <set>
#include <lemon/maps.h>
#include "mwcsgraph.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class MwcsSolver
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;

  typedef MwcsSolver<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> MwcsSolverType;
  typedef MwcsGraph<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> MwcsGraphType;

  TEMPLATE_GRAPH_TYPEDEFS(Graph);

  typedef typename std::set<Node> NodeSet;
  typedef typename NodeSet::const_iterator NodeSetIt;

  typedef std::vector<Node> NodeVector;
  typedef typename NodeVector::const_iterator NodeVectorIt;
  typedef typename std::vector<NodeSet> ModuleVector;
  typedef typename ModuleVector::const_iterator ModuleVectorIt;

public:
  MwcsSolver(const MwcsGraphType& mwcsGraph);

  virtual ~MwcsSolver();
  virtual void init(Node root = lemon::INVALID) = 0;
  virtual bool solve() = 0;
  virtual void setLowerBound(double LB) = 0;
  void printSolution(std::ostream& out, bool moduleOnly) const;

protected:
  const MwcsGraphType& _mwcsGraph;

  Node _root;
  double _score;

  BoolNodeMap _solutionMap;
  NodeSet _solutionSet;

public:
  const MwcsGraphType& getMwcsGraph() const
  {
    return _mwcsGraph;
  }

  double getSolutionWeight() const
  {
    return _score;
  }

  const BoolNodeMap& getSolutionNodeMap() const
  {
    return _solutionMap;
  }

  const NodeSet& getSolutionModule() const
  {
    return _solutionSet;
  }

  bool isNodeInSolution(Node n) const
  {
    assert(n != lemon::INVALID);
    return _solutionMap[n];
  }
};

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline MwcsSolver<GR, NWGHT, NLBL, EWGHT>::MwcsSolver(const MwcsGraphType& mwcsGraph)
  : _mwcsGraph(mwcsGraph)
  , _root(lemon::INVALID)
  , _score(0)
  , _solutionMap(mwcsGraph.getGraph(), false)
  , _solutionSet()
{
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline MwcsSolver<GR, NWGHT, NLBL, EWGHT>::~MwcsSolver()
{
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsSolver<GR, NWGHT, NLBL, EWGHT>::printSolution(std::ostream& out, bool moduleOnly) const
{
  const Graph& g = _mwcsGraph.getGraph();

  for (NodeIt n(g); n != lemon::INVALID; ++n)
  {
    if (isNodeInSolution(n))
    {
      out << "// " << _mwcsGraph.getLabel(n) << " "
          << _mwcsGraph.getScore(n) << std::endl;
    }
    else if (!moduleOnly)
    {
      out << "// " << _mwcsGraph.getLabel(n) << " 0" << std::endl;
    }
  }
}

} // namespace mwcs
} // namespace nina

#endif // MWCSSOLVER_H
