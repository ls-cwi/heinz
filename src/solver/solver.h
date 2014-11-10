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
#include <limits>

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
  Solver()
    : _score(0)
    , _scoreUB(std::numeric_limits<double>::max())
    , _pSolutionMap(NULL)
    , _solutionSet()
  {
  }
  
  virtual ~Solver()
  {
    delete _pSolutionMap;
  }
  
protected:
  double _score;
  double _scoreUB;
  BoolNodeMap* _pSolutionMap;
  NodeSet _solutionSet;
  
public:
  void printSolution(const MwcsGraphType& mwcsGraph,
                     std::ostream& out,
                     bool moduleOnly) const
  {
    const Graph& g = mwcsGraph.getGraph();
    
    for (NodeIt n(g); n != lemon::INVALID; ++n)
    {
      if (isNodeInSolution(n))
      {
        out << "// " << mwcsGraph.getLabel(n) << " "
            << mwcsGraph.getScore(n) << std::endl;
      }
      else if (!moduleOnly)
      {
        out << "// " << mwcsGraph.getLabel(n) << " 0" << std::endl;
      }
    }
  }

  double getSolutionWeight() const
  {
    return _score;
  }
  
  double getSolutionWeightUB() const
  {
    return _scoreUB;
  }
  
  const BoolNodeMap& getSolutionNodeMap() const
  {
    return *_pSolutionMap;
  }
  
  const NodeSet& getSolutionModule() const
  {
    return _solutionSet;
  }
  
  bool isNodeInSolution(Node n) const
  {
    assert(n != lemon::INVALID);
    return (*_pSolutionMap)[n];
  }
};
  
} // namespace mwcs
} // namespace nina

#endif // SOLVER_H