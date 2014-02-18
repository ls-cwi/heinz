/*
 * mwcssizecutsolver.h
 *
 *  Created on: 02-jul-2013
 *      Author: M. El-Kebir
 */

#ifndef MWCSSIZECUTSOLVER_H
#define MWCSSIZECUTSOLVER_H

#include <vector>
#include <set>
#include <assert.h>

#include <lemon/hao_orlin.h>
#include <lemon/gomory_hu.h>
#include <lemon/bfs.h>
#include <lemon/kruskal.h>
#include <lemon/adaptors.h>
#include <lemon/preflow.h>

#include <ilcplex/ilocplex.h>
#include <ilcplex/ilocplexi.h>
#include <ilconcert/ilothread.h>

#include "mwcscplexsolver.h"
#include "mwcsanalyze.h"
#include "cplex_cut/nodecutrootedbkcallback.h"
#include "cplex_cut/nodecutunrootedbkcallback.h"
#include "parser/identityparser.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class MwcsSizeCutSolver : public MwcsCplexSolver<GR, NWGHT, NLBL, EWGHT>
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;

  typedef MwcsCplexSolver<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> Parent;
  typedef typename Parent::MwcsGraphType MwcsGraphType;

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
  MwcsSizeCutSolver(const MwcsGraphType& mwcsGraph,
                    int moduleSize,
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

  int _moduleSize;
  int _maxNumberOfCuts;
  int _timeLimit;
  int _multiThreading;
};



template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline MwcsSizeCutSolver<GR, NWGHT, NLBL, EWGHT>::MwcsSizeCutSolver(const MwcsGraphType& mwcsGraph,
                                                                    int moduleSize,
                                                                    int maxNumberOfCuts,
                                                                    int timeLimit,
                                                                    int multiThreading)
  : Parent(mwcsGraph)
  , _moduleSize(moduleSize)
  , _maxNumberOfCuts(maxNumberOfCuts)
  , _timeLimit(timeLimit)
  , _multiThreading(multiThreading)
{
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline bool MwcsSizeCutSolver<GR, NWGHT, NLBL, EWGHT>::solveCplex()
{
  const Graph& g = _mwcsGraph.getGraph();
  const WeightNodeMap& weight = _mwcsGraph.getScores();

  IloFastMutex* pMutex = NULL;
  if (_multiThreading > 1)
  {
    pMutex = new IloFastMutex();
  }


  IloCplex::LazyConstraintCallbackI* pCut = NULL;
  if (_root != lemon::INVALID)
  {
    //_cplex.setParam(IloCplex::Cliques, -1);
    //_cplex.setParam(IloCplex::Covers, -1);
    //_cplex.setParam(IloCplex::FlowCovers, -1);
    //_cplex.setParam(IloCplex::GUBCovers, -1);
    //_cplex.setParam(IloCplex::FracCuts, -1);
    //_cplex.setParam(IloCplex::MIRCuts, -1);
    //_cplex.setParam(IloCplex::FlowPaths, -1);
    //_cplex.setParam(IloCplex::ImplBd, -1);
    //_cplex.setParam(IloCplex::DisjCuts, -1);
    //_cplex.setParam(IloCplex::ZeroHalfCuts, -1);
    //_cplex.setParam(IloCplex::MCFCuts, -1);

    pCut = new (_env) NodeCutRootedBkLazyCallback<GR, NWGHT, NLBL, EWGHT>(_env, _x, g, weight, _root, *_pNode,
                                                                          _n, _m, _maxNumberOfCuts, _mwcsGraph.getComponentMap(),
                                                                          pMutex);
  }
  else
  {
    //_cplex.setParam(IloCplex::Cliques, -1);
    //_cplex.setParam(IloCplex::Covers, -1);
    //_cplex.setParam(IloCplex::FlowCovers, -1);
    //_cplex.setParam(IloCplex::GUBCovers, -1);
    //_cplex.setParam(IloCplex::FracCuts, -1);
    //_cplex.setParam(IloCplex::MIRCuts, -1);
    //_cplex.setParam(IloCplex::FlowPaths, -1);
    //_cplex.setParam(IloCplex::ImplBd, -1);
    //_cplex.setParam(IloCplex::DisjCuts, -1);
    //_cplex.setParam(IloCplex::ZeroHalfCuts, -1);
    //_cplex.setParam(IloCplex::MCFCuts, -1);
//
    //_cplex.setParam(IloCplex::AggFill, 0);
    //_cplex.setParam(IloCplex::PreInd, 0);
    //_cplex.setParam(IloCplex::RelaxPreInd, 0);
    //_cplex.setParam(IloCplex::PreslvNd, -1);
    //_cplex.setParam(IloCplex::RepeatPresolve, 0);

    pCut = new (_env) NodeCutUnrootedBkLazyCallback<GR, NWGHT, NLBL, EWGHT>(_env, _x, _y, g, weight, *_pNode,
                                                                            _n, _m, _maxNumberOfCuts, _mwcsGraph.getComponentMap(),
                                                                            pMutex);
  }

  IloCplex::Callback cb(pCut);
  _cplex.use(cb);

  if (_timeLimit > 0)
  {
    _cplex.setParam(IloCplex::TiLim, _timeLimit);
  }

  if (_multiThreading > 1)
  {
    _cplex.setParam(IloCplex::ParallelMode, -1);
    _cplex.setParam(IloCplex::Threads, _multiThreading);
  }

  //exportModel("/tmp/model.lp");
  bool res = _cplex.solve();
  cb.end();

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
inline void MwcsSizeCutSolver<GR, NWGHT, NLBL, EWGHT>::initConstraints()
{
  //Parent::initConstraints();

  const Graph& g = _mwcsGraph.getGraph();

  if (_root == lemon::INVALID)
  {
    IloExpr expr(_env);
    // there is at most one root node
    // \sum_{i \in V} y_i <= 1
    for (int i = 0; i < _n; i++)
    {
      expr += _y[i];
    }
    IloConstraint c1;
    _model.add(c1 = (expr == 1));
    c1.setName("one_root");

    // the root node has to be one of the selected nodes
    // in the final graph
    // y_i <= x_i for all nodes i in V
    for (int i = 0; i < _n; i++)
    {
      IloConstraint c2;
      _model.add(c2 = (_y[i] <= _x[i]));
      c2.setName("root_in_x");
    }
  }
  else
  {
    // x_r = 1
    int r = (*_pNode)[_root];
    IloConstraint c2;
    _model.add(c2 = (_x[r] == 1));
    c2.setName("root");
  }

  //lemon::Bfs<Graph> bfs(g);

  // module size
  IloExpr expr(_env);
  for (NodeIt i(g); i != lemon::INVALID; ++i)
  {
    expr += _x[(*_pNode)[i]];
  }
  _model.add(expr == _moduleSize);

  // if you pick a node then it must be the root node
  // or at least one of its direct neighbors must be part
  // of the solution as well
  //for (NodeIt i(g); i != lemon::INVALID; ++i)
  //{
  //  if (i == _root)
  //    continue;
//
  //  expr.clear();
  //  for (IncEdgeIt e(g, i); e != lemon::INVALID; ++e)
  //  {
  //    Node j = g.oppositeNode(i, e);
  //    expr += _x[(*_pNode)[j]];
  //  }
//
  //  if (_root == lemon::INVALID)
  //    expr += _y[(*_pNode)[i]];
//
  //  _model.add(_x[(*_pNode)[i]] <= expr);
  //}

  //if (_root == lemon::INVALID)
  //{
  //  // symmetry breaking
  //  for (NodeIt i(g); i != lemon::INVALID; ++i)
  //  {
  //    for (NodeIt j(g); j != i; ++j)
  //    {
  //      if (j == i) continue;
  //      _model.add(_y[(*_pNode)[i]] <= 1 - _x[(*_pNode)[j]]);
  //    }
  //  }
  //}

  //lemon::Bfs<Graph> bfs(g);
  //for (NodeIt i(g); i != lemon::INVALID; ++i)
  //{
  //  // do a bfs until level > k
  //  bfs.init();
  //  bfs.addSource(i);
  //  while (!bfs.emptyQueue())
  //  {
  //    Node node = bfs.processNextNode();
  //    int d = bfs.dist(node);
  //    if (d > _moduleSize)
  //    {
  //      break;
  //    }
  //  }
//
  //  for (NodeIt j(g); j != lemon::INVALID; ++j)
  //  {
  //    if (!bfs.reached(j))
  //    {
  //      _model.add(_x[(*_pNode)[i]] <= 1 - _x[(*_pNode)[j]]);
  //      _model.add(_x[(*_pNode)[j]] <= 1 - _x[(*_pNode)[i]]);
  //    }
  //  }
  //}
}

} // namespace mwcs
} // namespace nina

#endif // MWCSSIZECUTSOLVER_H
