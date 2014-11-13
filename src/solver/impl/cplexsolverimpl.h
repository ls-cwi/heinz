/*
 * cplexsolverimpl.h
 *
 *  Created on: 10-aug-2012
 *     Authors: M. El-Kebir
 */

#ifndef CPLEXSOLVERIMPL_H
#define CPLEXSOLVERIMPL_H

#include "cplex_cut/backoff.h"
#include "analysis.h"

#include <set>
#include <vector>

// ILOG stuff
#include <ilconcert/iloalg.h>
#include <ilcplex/ilocplex.h>

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class CplexSolverImpl
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;

  typedef MwcsGraph<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> MwcsGraphType;
  typedef MwcsAnalyze<Graph> MwcsAnalyzeType;
  
  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  
  typedef std::set<Node> NodeSet;
  typedef typename NodeSet::const_iterator NodeSetIt;
  typedef std::vector<Node> NodeVector;
  typedef typename NodeVector::const_iterator NodeVectorIt;
  typedef std::vector<Node> InvNodeIntMap;
  typedef std::vector<Arc> InvArcIntMap;
  
  struct Options
  {
    Options(const BackOff& backOff,
            bool analysis,
            int maxNumberOfCuts,
            int timeLimit,
            int multiThreading,
            int memoryLimit)
      : _backOff(backOff)
      , _analysis(analysis)
      , _maxNumberOfCuts(maxNumberOfCuts)
      , _timeLimit(timeLimit)
      , _multiThreading(multiThreading)
      , _memoryLimit(memoryLimit)
    {
    }
    
    const BackOff& _backOff;
    bool _analysis;
    int _maxNumberOfCuts;
    int _timeLimit;
    int _multiThreading;
    int _memoryLimit;
  };

protected:
  CplexSolverImpl(const Options& options)
    : _options(options)
    , _pAnalysis(NULL)
    , _n(0)
    , _m(0)
    , _pNode(NULL)
    , _invNode()
    , _env()
    , _model(_env)
    , _cplex(_model)
    , _x()
  {
  }
  
  virtual ~CplexSolverImpl()
  {
    _env.end();
  }

  void exportModel(const std::string& filename)
  {
    _cplex.exportModel(filename.c_str());
  }

  virtual void printVariables(const MwcsGraphType& mwcsGraph, std::ostream& out)
  {
    for (int id_v = 0; id_v < _n; id_v++)
    {
      out << "Node '" << mwcsGraph.getLabel(_invNode[id_v])
      << "', x_" << id_v << " = " << _cplex.getValue(_x[id_v])
      << ", w = " << mwcsGraph.getScore(_invNode[id_v])
      << std::endl;
    }
  }

protected:
  const Options& _options;
  MwcsAnalyzeType* _pAnalysis;

  int _n;
  int _m;
  IntNodeMap* _pNode;
  InvNodeIntMap _invNode;
  IloEnv _env;
  IloModel _model;
  IloCplex _cplex;
  IloBoolVarArray _x;

  virtual void initVariables(const MwcsGraphType& mwcsGraph);
  virtual void initConstraints(const MwcsGraphType& mwcsGraph);
  virtual void clean();

  virtual bool solveCplex(const MwcsGraphType& mwcsGraph,
                          double& score,
                          double& scoreUB,
                          BoolNodeMap& solutionMap,
                          NodeSet& solutionSet);
  
  virtual bool solveModel() = 0;

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
  
  struct NodesWeightComp
  {
  private:
    const DoubleNodeMap& _weight;
    
  public:
    NodesWeightComp(const DoubleNodeMap& weight)
    : _weight(weight)
    {
    }
    
    bool operator ()(Node u, Node v)
    {
      return _weight[u] > _weight[v];
    }
  };
};

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
void CplexSolverImpl<GR, NWGHT, NLBL, EWGHT>::clean()
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
inline void CplexSolverImpl<GR, NWGHT, NLBL, EWGHT>::initVariables(const MwcsGraphType& mwcsGraph)
{
  const Graph& g = mwcsGraph.getGraph();
  _n = mwcsGraph.getNodeCount();
  _m = mwcsGraph.getArcCount();

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
//  NodesWeightComp comp(mwcsGraph.getScores());
  std::sort(_invNode.begin(), _invNode.end(), comp);

  char buf[1024];

  int i = 0;
  for (NodeVectorIt it = _invNode.begin(); it != _invNode.end(); ++it, ++i)
  {
    Node v = *it;

    // x_i = 0 if node i is not in the subgraph
    // x_i = 1 if node i is the subgraph
    snprintf(buf, 1024, "x_%s", mwcsGraph.getLabel(v).c_str());
//    snprintf(buf, 1024, "x_%d", g.id(v));
    _x[i].setName(buf);

    (*_pNode)[v] = i;
  }
}
  
template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void CplexSolverImpl<GR, NWGHT, NLBL, EWGHT>::initConstraints(const MwcsGraphType& mwcsGraph)
{
  const Graph& g = mwcsGraph.getGraph();
  const WeightNodeMap& weight = mwcsGraph.getScores();
  
  // objective function
  IloExpr expr(_env);
  for (int i = 0; i < _n ; i++)
  {
    expr += _x[i] * weight[_invNode[i]];
  }
  _model.add(IloObjective(_env, expr, IloObjective::Maximize));

  // add equality constraints
  if (_pAnalysis)
  {
    int nAnalyzeConstraints = 0;
    for (NodeIt i(g); i != lemon::INVALID; ++i)
    {
      if (weight[i] > 0)
      {
        for (NodeIt j(g); j != lemon::INVALID; ++j)
        {
          if (weight[j] > 0 && _pAnalysis->ok(i, j))
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
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline bool CplexSolverImpl<GR, NWGHT, NLBL, EWGHT>::solveCplex(const MwcsGraphType& mwcsGraph,
                                                                double& score,
                                                                double& scoreUB,
                                                                BoolNodeMap& solutionMap,
                                                                NodeSet& solutionSet)
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
    int limit = _options._timeLimit - g_timer.realTime();
    limit = std::max(1, limit);
    _cplex.setParam(IloCplex::TiLim, limit);
  }
  
  if (_options._memoryLimit > 0)
  {
    _cplex.setParam(IloCplex::TreLim, _options._memoryLimit);
  }
  
  if (_options._multiThreading > 1)
  {
    _cplex.setParam(IloCplex::ParallelMode, -1);
    _cplex.setParam(IloCplex::Threads, _options._multiThreading);
  }

  bool optimal = solveModel();
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
  solutionSet.clear();
  for (int i = 0; i < _n ; i++)
  {
    Node node = _invNode[i];
    solutionMap[node] = _cplex.getValue(_x[i]);

    if (tol.nonZero(_cplex.getValue(_x[i])))
    {
      solutionSet.insert(node);
    }
  }

  score = _cplex.getObjValue();
  scoreUB = _cplex.getBestObjValue();
  clean();

  return true;
}

} // namespace mwcs
} // namespace nina

#endif // CPLEXSOLVER_H
