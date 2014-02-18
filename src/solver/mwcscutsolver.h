/*
 * mwcscutsolver.h
 *
 *  Created on: 10-dec-2012
 *      Author: M. El-Kebir
 */

#ifndef MWCSCUTSOLVER_H
#define MWCSCUTSOLVER_H

#include <vector>
#include <set>
#include <assert.h>

#include <ilcplex/ilocplex.h>
#include <ilcplex/ilocplexi.h>
#include <ilconcert/ilothread.h>

#include <lemon/hao_orlin.h>
#include <lemon/gomory_hu.h>
#include <lemon/bfs.h>
#include <lemon/kruskal.h>
#include <lemon/adaptors.h>
#include <lemon/preflow.h>

#include "cplex_heuristic/mwcscutsolverheuristic.h"
#include "cplex_cut/nodecutrootedbkcallback.h"
#include "cplex_cut/nodecutunrootedbkcallback.h"
#include "mwcscplexsolver.h"
#include "mwcsanalyze.h"
#include "parser/identityparser.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class MwcsCutSolver : public MwcsCplexSolver<GR, NWGHT, NLBL, EWGHT>
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;

  typedef MwcsCplexSolver<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> Parent;
  typedef typename Parent::MwcsGraphType MwcsGraphType;
  typedef MwcsCutSolverHeuristic<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> MwcsCutSolverHeuristicType;

  TEMPLATE_GRAPH_TYPEDEFS(Graph);

  using Parent::_mwcsGraph;
  using Parent::_root;
  using Parent::_score;
  using Parent::_solutionMap;
  using Parent::_solutionSet;
  using Parent::_n;
  using Parent::_m;
  using Parent::_pNode;
  using Parent::_invNode;
  using Parent::_env;
  using Parent::_model;
  using Parent::_cplex;
  using Parent::_x;
  using Parent::_y;
  using Parent::_LB;
  using Parent::printSolution;
  using Parent::getSolutionWeight;
  using Parent::getSolutionNodeMap;
  using Parent::getSolutionModule;
  using Parent::isNodeInSolution;
  using Parent::init;
  using Parent::printVariables;
  using Parent::initConstraints;
  using Parent::clean;
  using Parent::solveCplex;
  using Parent::exportModel;

public:
  MwcsCutSolver(const MwcsGraphType& mwcsGraph,
                int maxNumberOfCuts = -1,
                int timeLimit = -1,
                int multiThreading = 1);

  bool solveCplex();

protected:
  virtual void initConstraints();

private:
  typedef std::vector<bool> BoolVector;
  typedef std::vector<Node> NodeVector;
  typedef std::set<Node> NodeSet;
  typedef typename NodeSet::const_iterator NodeSetIt;
  typedef MwcsAnalyze<Graph> MwcsAnalyzeType;

  int _maxNumberOfCuts;
  int _timeLimit;
  int _multiThreading;
};



template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline MwcsCutSolver<GR, NWGHT, NLBL, EWGHT>::MwcsCutSolver(const MwcsGraphType& mwcsGraph,
                                                            int maxNumberOfCuts,
                                                            int timeLimit,
                                                            int multiThreading)
  : Parent(mwcsGraph)
  , _maxNumberOfCuts(maxNumberOfCuts)
  , _timeLimit(timeLimit)
  , _multiThreading(multiThreading)
{
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline bool MwcsCutSolver<GR, NWGHT, NLBL, EWGHT>::solveCplex()
{
  const Graph& g = _mwcsGraph.getGraph();
  const WeightNodeMap& weight = _mwcsGraph.getScores();

  IloFastMutex* pMutex = NULL;
  if (_multiThreading > 1)
  {
    pMutex = new IloFastMutex();
  }

  IloCplex::LazyConstraintCallbackI* pLazyCut = NULL;
  IloCplex::UserCutCallbackI* pUserCut = NULL;
  if (_root != lemon::INVALID)
  {
    _cplex.setParam( IloCplex::HeurFreq      , -1 );
    _cplex.setParam( IloCplex::Cliques       , -1 );
    _cplex.setParam( IloCplex::Covers        , -1 );
    _cplex.setParam( IloCplex::FlowCovers    , -1 );
    _cplex.setParam( IloCplex::GUBCovers     , -1 );
    _cplex.setParam( IloCplex::FracCuts      , -1 );
    _cplex.setParam( IloCplex::MIRCuts       , -1 );
    _cplex.setParam( IloCplex::FlowPaths     , -1 );
    _cplex.setParam( IloCplex::ImplBd        , -1 );
    _cplex.setParam( IloCplex::DisjCuts      , -1 );
    _cplex.setParam( IloCplex::ZeroHalfCuts  , -1 );
    _cplex.setParam( IloCplex::MCFCuts       , -1 );
    _cplex.setParam( IloCplex::AggFill       ,  0 );
    _cplex.setParam( IloCplex::PreInd        ,  0 );
    _cplex.setParam( IloCplex::RelaxPreInd   ,  0 );
    _cplex.setParam( IloCplex::PreslvNd      , -1 );
    _cplex.setParam( IloCplex::RepeatPresolve,  0 );

    pLazyCut = new (_env) NodeCutRootedBkLazyCallback<GR, NWGHT, NLBL, EWGHT>(_env, _x, g, weight, _root, *_pNode,
                                                                              _n, _m, _maxNumberOfCuts, _mwcsGraph.getComponentMap(),
                                                                              pMutex);
    pUserCut = new (_env) NodeCutRootedBkUserCallback<GR, NWGHT, NLBL, EWGHT>(_env, _x, g, weight, _root, *_pNode,
                                                                              _n, _m, _maxNumberOfCuts, _mwcsGraph.getComponentMap(),
                                                                              pMutex);
  }
  else
  {
    _cplex.setParam( IloCplex::HeurFreq      , -1 );
    _cplex.setParam( IloCplex::Cliques       , -1 );
    _cplex.setParam( IloCplex::Covers        , -1 );
    _cplex.setParam( IloCplex::FlowCovers    , -1 );
    _cplex.setParam( IloCplex::GUBCovers     , -1 );
    _cplex.setParam( IloCplex::FracCuts      , -1 );
    _cplex.setParam( IloCplex::MIRCuts       , -1 );
    _cplex.setParam( IloCplex::FlowPaths     , -1 );
    _cplex.setParam( IloCplex::ImplBd        , -1 );
    _cplex.setParam( IloCplex::DisjCuts      , -1 );
    _cplex.setParam( IloCplex::ZeroHalfCuts  , -1 );
    _cplex.setParam( IloCplex::MCFCuts       , -1 );
    _cplex.setParam( IloCplex::AggFill       ,  0 );
    _cplex.setParam( IloCplex::PreInd        ,  0 );
    _cplex.setParam( IloCplex::RelaxPreInd   ,  0 );
    _cplex.setParam( IloCplex::PreslvNd      , -1 );
    _cplex.setParam( IloCplex::RepeatPresolve,  0 );

    pLazyCut = new (_env) NodeCutUnrootedBkLazyCallback<GR, NWGHT, NLBL, EWGHT>(_env, _x, _y, g, weight, *_pNode,
                                                                                _n, _m, _maxNumberOfCuts, _mwcsGraph.getComponentMap(),
                                                                                pMutex);
    //pUserCut = new (_env) NodeCutUnrootedBkUserCallback<GR, NWGHT, NLBL, EWGHT>(_env, _x, _y, g, weight, *_pNode,
    //                                                                            _n, _m, _maxNumberOfCuts, _mwcsGraph.getComponentMap(),
    //                                                                            pMutex);

  }

  if (_timeLimit > 0)
  {
    _cplex.setParam(IloCplex::TiLim, _timeLimit);
  }

  if (_multiThreading > 1)
  {
    _cplex.setParam(IloCplex::ParallelMode, -1);
    _cplex.setParam(IloCplex::Threads, _multiThreading);
  }

  IloCplex::Callback cb(pLazyCut);
  _cplex.use(cb);

  IloCplex::Callback cb2(new (_env) MwcsCutSolverHeuristicType(_env, _x, _y,
                                                               g, weight, _root,
                                                               *_pNode,
                                                               _n, _m, pMutex));
  _cplex.use(cb2);

  IloCplex::Callback cb3(pUserCut);
  _cplex.use(cb3);

  //exportModel("/tmp/model.lp");
  bool res = _cplex.solve();
  cb.end();
  cb2.end();
  cb3.end();

  _cplex.setParam(IloCplex::MIPInterval, 2);
  _cplex.setParam(IloCplex::MIPDisplay, 5);

  if (res)
  {
    std::cerr << "[" << _cplex.getObjValue() << ", "
              << _cplex.getBestObjValue() << "]" << std::endl;
  }
  else
  {
    std::cerr << "[0, 0]" << std::endl;
  }

  //printVariables(std::cerr);
  delete pMutex;
  return res;
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsCutSolver<GR, NWGHT, NLBL, EWGHT>::initConstraints()
{
  Parent::initConstraints();

  const Graph& g = _mwcsGraph.getGraph();
  const WeightNodeMap& weight = _mwcsGraph.getScores();

  // objective must be positive
  IloExpr expr(_env);
  for (int i = 0; i < _n ; i++)
  {
    expr += _x[i] * weight[_invNode[i]];
  }
  _model.add(expr >= _LB);

  // if you pick a node then it must be the root node
  // or at least one of its direct neighbors must be part
  // of the solution as well
  for (NodeIt i(g); i != lemon::INVALID; ++i)
  {
    if (i == _root)
      continue;

    expr.clear();
    for (IncEdgeIt e(g, i); e != lemon::INVALID; ++e)
    {
      Node j = g.oppositeNode(i, e);
      expr += _x[(*_pNode)[j]];

      // if i is negative then its positive neighbors must be in
      if (weight[i] < 0 && weight[j] > 0)
        _model.add(_x[(*_pNode)[i]] <= _x[(*_pNode)[j]]);
    }

    if (_root == lemon::INVALID)
      expr += _y[(*_pNode)[i]];

    _model.add(_x[(*_pNode)[i]] <= expr);
  }

  if (_root == lemon::INVALID)
  {
    // if you pick a negative node, then at least one of its direct neighbors
    // must be part of the solution as well
    for (NodeIt i(g); i != lemon::INVALID; ++i)
    {
      if (_mwcsGraph.getScore(i) <= 0)
      {
        expr.clear();
        for (IncEdgeIt e(g, i); e != lemon::INVALID; ++e)
        {
          Node j = g.oppositeNode(i, e);
          expr += _x[(*_pNode)[j]];
        }
        _model.add(_x[(*_pNode)[i]] <= expr);
      }
      else
      {
        // symmetry breaking
        int id_i = (*_pNode)[i];

        for (int id_j = 0; id_j < id_i; ++id_j)
        {
          Node j = _invNode[id_j];
          if (_mwcsGraph.getScore(j) < 0) continue;
          _model.add(_y[id_i] <= 1 - _x[id_j]);
        }
      }
    }
  }

  if (_root != lemon::INVALID)
  {
    // nodes i that are not in the same component as the root get x_i = 0
    const int rootComp = _mwcsGraph.getComponent(_root);
    for (NodeIt i(g); i != lemon::INVALID; ++i)
    {
      if (_mwcsGraph.getComponent(i) != rootComp)
      {
        _model.add(_x[(*_pNode)[i]] == 0);
      }
    }
  }
  else
  {
    // nodes in different components can never occur together
    //for (NodeIt i(g); i != lemon::INVALID; ++i)
    //{
    //  for (NodeIt j = i; j != lemon::INVALID; ++j)
    //  {
    //    if (_mwcsGraph.getComponent(i) != _mwcsGraph.getComponent(j))
    //    {
    //      _model.add(_x[(*_pNode)[i]] <= 1 - _x[(*_pNode)[j]]);
    //      _model.add(_x[(*_pNode)[j]] <= 1 - _x[(*_pNode)[i]]);
    //    }
    //  }
    //}
  }

  // add equality constraints
  int nAnalyzeConstraints = 0;
  MwcsAnalyzeType analyze(_mwcsGraph);
  analyze.analyze();
  for (NodeIt i(g); i != lemon::INVALID; ++i)
  {
    if (_mwcsGraph.getScore(i) > 0)
    {
      for (NodeIt j(g); j != lemon::INVALID; ++j)
      {
        if (_mwcsGraph.getScore(j) > 0 && analyze.ok(i, j))
        {
          _model.add(_x[(*_pNode)[i]] <= _x[(*_pNode)[j]]);
          nAnalyzeConstraints++;
        }
      }
    }

    // set priorities
    //if (_mwcsGraph.getScore(i) > 0)
    //  _cplex.setPriority(_x[(*_pNode)[i]], 1);
  }

  if (g_verbosity >= VERBOSE_ESSENTIAL)
  {
    std::cout << "// Added " << nAnalyzeConstraints << " analyze constraints" << std::endl;
  }

  //for (int i = 0; i < analyze.getEqClassesCount(); i++)
  //{
  //  const NodeSet& nodes = analyze.getEqClasses()[i];
  //  for (NodeSetIt nodeIt1 = nodes.begin(); nodeIt1 != nodes.end(); nodeIt1++)
  //  {
  //    //std::cout << g.id(*nodeIt1) << "\t" << _mwcsGraph.getLabel(*nodeIt1) << "\t" << _mwcsGraph.getScore(*nodeIt1) << std::endl;
  //    for (NodeSetIt nodeIt2 = nodeIt1; nodeIt2 != nodes.end(); nodeIt2++)
  //    {
  //      if (nodeIt1 == nodeIt2) continue;
  //      _model.add(_x[(*_pNode)[*nodeIt1]] == _x[(*_pNode)[*nodeIt2]]);
  //    }
  //  }
  //}
}

} // namespace mwcs
} // namespace nina

#endif // MWCSCUTSOLVER_H
