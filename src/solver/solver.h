/*
 * solver.h
 *
 *  Created on: 21-aug-2014
 *      Author: M. El-Kebir
 */

#ifndef SOLVER_H
#define SOLVER_H

#include <set>
#include <vector>
#include <ostream>

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class Solver
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;
  
  typedef MwcsGraph<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> MwcsGraphType;
  
  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  
  typedef std::set<Node> NodeSet;
  typedef typename NodeSet::const_iterator NodeSetIt;
  typedef typename NodeSet::iterator NodeSetNonConstIt;
  typedef std::vector<Node> NodeVector;
  typedef typename NodeVector::const_iterator NodeVectorIt;
  typedef typename NodeVector::iterator NodeVectorNonConstIt;
  
public:
  Solver(const MwcsGraphType& mwcsGraph)
    : _mwcsGraph(mwcsGraph)
    , _score(0)
    , _solutionMap(mwcsGraph.getGraph(), false)
    , _solutionSet()
  {
  }
  
  virtual ~Solver()
  {
  }
  
protected:
  const MwcsGraphType& _mwcsGraph;
  
  double _score;
  
  BoolNodeMap _solutionMap;
  NodeSet _solutionSet;
  
public:
  virtual void init() = 0;
  
  virtual bool solve() = 0;

  void printSolution(std::ostream& out, bool moduleOnly) const
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
  
} // namespace mwcs
} // namespace nina

#endif // SOLVER_H