/*
 * cplexsolver.cpp
 *
 *  Created on: 10-aug-2012
 *     Authors: M. El-Kebir
 */

#ifndef CPLEXSOLVER_H
#define CPLEXSOLVER_H

#include "solver.h"
#include "cplex_cut/backoff.h"
#include "analysis.h"

// ILOG stuff
#include <ilconcert/iloalg.h>
#include <ilcplex/ilocplex.h>

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class CplexSolver : public Solver<GR, NWGHT, NLBL, EWGHT>
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;

  typedef Solver<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> Parent;
  typedef typename Parent::MwcsGraphType MwcsGraphType;
  typedef MwcsAnalyze<Graph> MwcsAnalyzeType;
  
  typedef typename Parent::NodeSet NodeSet;
  typedef typename Parent::NodeSetIt NodeSetIt;
  typedef typename Parent::NodeSetNonConstIt NodeSetNonConstIt;
  typedef typename Parent::NodeVector NodeVector;
  typedef typename Parent::NodeVectorIt NodeVectorIt;
  typedef typename Parent::NodeVectorNonConstIt NodeVectorNonConstIt;

  TEMPLATE_GRAPH_TYPEDEFS(Graph);

  typedef std::vector<Node> InvNodeIntMap;
  typedef std::vector<Arc> InvArcIntMap;

  using Parent::_mwcsGraph;
  using Parent::_score;
  using Parent::_solutionMap;
  using Parent::_solutionSet;
  using Parent::printSolution;
  using Parent::getSolutionWeight;
  using Parent::getSolutionNodeMap;
  using Parent::getSolutionModule;
  using Parent::isNodeInSolution;
  using Parent::init;
  
  struct Options
  {
    Options(const BackOff& backOff,
            int maxNumberOfCuts = -1,
            int timeLimit = -1,
            int multiThreading = 1)
      : _backOff(backOff)
      , _maxNumberOfCuts(maxNumberOfCuts)
      , _timeLimit(timeLimit)
      , _multiThreading(multiThreading)
    {
    }
    
    const BackOff& _backOff;
    int _maxNumberOfCuts;
    int _timeLimit;
    int _multiThreading;
  };

public:
  CplexSolver(const MwcsGraphType& mwcsGraph,
              const Options& options,
              const MwcsAnalyzeType& analysis)
    : Parent(mwcsGraph)
    , _options(options)
    , _analysis(analysis)
    , _n(mwcsGraph.getNodeCount())
    , _m(mwcsGraph.getArcCount())
    , _pNode(NULL)
    , _invNode()
    , _env()
    , _model(_env)
    , _cplex(_model)
    , _x()
  {
  }
  
  virtual ~CplexSolver()
  {
    _env.end();
  }

  void exportModel(const std::string& filename)
  {
    _cplex.exportModel(filename.c_str());
  }

  virtual bool solve();
  virtual void init();
  virtual void printVariables(std::ostream& out);

protected:
  const Options& _options;
  const MwcsAnalyzeType& _analysis;
  int _n;
  int _m;
  IntNodeMap* _pNode;
  InvNodeIntMap _invNode;
  IloEnv _env;
  IloModel _model;
  IloCplex _cplex;
  IloBoolVarArray _x;

  virtual void initVariables();
  virtual void initConstraints();
  virtual void clean();

  virtual bool solveCplex() = 0;

private:
  struct NodesDegComp
  {
  private:
    const IntNodeMap& _deg;

  public:
    NodesDegComp(const IntNodeMap& deg)
      : _deg(deg)
    {
    }

    bool operator ()(Node u, Node v)
    {
      return _deg[u] > _deg[v];
    }
  };
};

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
void CplexSolver<GR, NWGHT, NLBL, EWGHT>::clean()
{
  _cplex.end();
  _env.end();
  _env = IloEnv();
  _model = IloModel(_env);
  _cplex = IloCplex(_model);
  delete _pNode;
  _pNode = NULL;
}


template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void CplexSolver<GR, NWGHT, NLBL, EWGHT>::init()
{
  initVariables();
  initConstraints();
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void CplexSolver<GR, NWGHT, NLBL, EWGHT>::printVariables(std::ostream& out)
{
  for (int id_v = 0; id_v < _n; id_v++)
  {
    out << "Node '" << _mwcsGraph.getLabel(_invNode[id_v])
        << "', x_" << id_v << " = " << _cplex.getValue(_x[id_v])
        << ", w = " << _mwcsGraph.getScore(_invNode[id_v])
        << std::endl;
  }
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void CplexSolver<GR, NWGHT, NLBL, EWGHT>::initVariables()
{
  const Graph& g = _mwcsGraph.getGraph();
  _n = _mwcsGraph.getNodeCount();
  _m = _mwcsGraph.getArcCount();

  delete _pNode;
  _pNode = new IntNodeMap(g);

  _x = IloBoolVarArray(_env, _n);

  // let's sort the nodes on degree
  IntNodeMap deg(g);
  _invNode.clear();
  _invNode.reserve(_n);
  for (NodeIt v(g); v != lemon::INVALID; ++v)
  {
    int d = 0;
    for (IncEdgeIt e(g, v); e != lemon::INVALID; ++e) ++d;
    deg[v] = d;
    _invNode.push_back(v);
  }

  NodesDegComp comp(deg);
  std::sort(_invNode.begin(), _invNode.end(), comp);

  char buf[1024];

  int i = 0;
  for (NodeVectorIt it = _invNode.begin(); it != _invNode.end(); ++it, ++i)
  {
    Node v = *it;

    // x_i = 0 if node i is not in the subgraph
    // x_i = 1 if node i is the subgraph
    snprintf(buf, 1024, "x_%s", _mwcsGraph.getLabel(v).c_str());
//    snprintf(buf, 1024, "x_%d", g.id(v));
    _x[i].setName(buf);

    (*_pNode)[v] = i;
  }
}
  
template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void CplexSolver<GR, NWGHT, NLBL, EWGHT>::initConstraints()
{
  const Graph& g = _mwcsGraph.getGraph();
  const WeightNodeMap& weight = _mwcsGraph.getScores();
  
  // objective function
  IloExpr expr(_env);
  for (int i = 0; i < _n ; i++)
  {
    expr += _x[i] * _mwcsGraph.getScore(_invNode[i]);
  }
  _model.add(IloObjective(_env, expr, IloObjective::Maximize));

  // add equality constraints
  int nAnalyzeConstraints = 0;
  for (NodeIt i(g); i != lemon::INVALID; ++i)
  {
    if (_mwcsGraph.getScore(i) > 0)
    {
      for (NodeIt j(g); j != lemon::INVALID; ++j)
      {
        if (_mwcsGraph.getScore(j) > 0 && _analysis.ok(i, j))
        {
          _model.add(_x[(*_pNode)[i]] <= _x[(*_pNode)[j]]);
          nAnalyzeConstraints++;
        }
      }
    }
  }

  if (g_verbosity >= VERBOSE_ESSENTIAL)
  {
    std::cout << "// Added " << nAnalyzeConstraints << " analyze constraints" << std::endl;
  }
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline bool CplexSolver<GR, NWGHT, NLBL, EWGHT>::solve()
{
  // shut up cplex
  if (g_verbosity < VERBOSE_NON_ESSENTIAL)
  {
    _cplex.setOut(_env.getNullStream());
    _cplex.setWarning(_env.getNullStream());
    _cplex.setError(_env.getNullStream());
  }
  else
  {
    _cplex.setOut(std::cerr);
    _cplex.setWarning(std::cerr);
    _cplex.setError(std::cerr);
  }
  
  if (_options._timeLimit > 0)
  {
    _cplex.setParam(IloCplex::TiLim, _options._timeLimit);
  }
  
  if (_options._multiThreading > 1)
  {
    _cplex.setParam(IloCplex::ParallelMode, -1);
    _cplex.setParam(IloCplex::Threads, _options._multiThreading);
  }

  bool optimal = solveCplex();
  if (!optimal)
  {
    if (_cplex.getStatus() == IloAlgorithm::Infeasible)
    {
      clean();
      return false;
    }
    else
    {
      if (g_verbosity >= VERBOSE_ESSENTIAL)
      {
        std::cerr << _cplex.getStatus() << std::endl;
        std::cerr << "Optimization problems. CPLEX status code " << _cplex.getStatus();
      }

      clean();
      return false;
    }
  }

  lemon::Tolerance<double> tol(1e-5);

  // solution
  _solutionSet.clear();
  for (int i = 0; i < _n ; i++)
  {
    Node node = _invNode[i];
    _solutionMap[node] = _cplex.getValue(_x[i]);

    if (tol.nonZero(_cplex.getValue(_x[i])))
    {
      _solutionSet.insert(node);
    }
  }

  _score = _cplex.getObjValue();
  clean();

  return true;
}

} // namespace mwcs
} // namespace nina

#endif // CPLEXSOLVER_H
