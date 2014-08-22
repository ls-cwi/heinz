/*
 * cutsolverrooted.h
 *
 *  Created on: 22-aug-2014
 *      Author: M. El-Kebir
 */

#ifdef CUTSOLVERUNROOTED_H
#define CUTSOLVERUNROOTED_H

#include "solverrooted.h"
#include "cplexsolver.h"

#include <ilconcert/ilothread.h>

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class CutSolverRooted : public SolverRooted<GR, NWGHT, NLBL, EWGHT>,
                        public CplexSolver<GR, NWGHT, NLBL, EWGHT>
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;
  
  typedef Solver<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> GrandParent;
  typedef SolverRooted<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> Parent1;
  typedef CplexSolver<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> Parent2;

  typedef typename Parent1::MwcsGraphType MwcsGraphType;
  typedef typename Parent1::NodeSet NodeSet;
  typedef typename Parent1::NodeSetIt NodeSetIt;
  typedef typename Parent1::NodeVector NodeVector;
  typedef typename Parent1::NodeVectorIt NodeVectorIt;
  
  typedef typename Parent2::MwcsAnalyzeType MwcsAnalyzeType;
  typedef typename Parent2::InvNodeIntMap InvNodeIntMap;
  typedef typename Parent2::InvArcIntMap InvArcIntMap;
  typedef typename Parent2::Options Options;
  
  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  
  using Parent1::_mwcsGraph;
  using Parent1::_score;
  using Parent1::_solutionMap;
  using Parent1::_solutionSet;
  using Parent1::_rootNodes;
  using Parent1::printSolution;
  using Parent1::getSolutionWeight;
  using Parent1::getSolutionNodeMap;
  using Parent1::getSolutionModule;
  using Parent1::isNodeInSolution;

  using Parent2::_options;
  using Parent2::_n;
  using Parent2::_m;
  using Parent2::_pNode;
  using Parent2::_invNode;
  using Parent2::_env;
  using Parent2::_model;
  using Parent2::_cplex;
  using Parent2::_x;
  using Parent2::init;
  using Parent2::initVariables;
  using Parent2::initConstraints;
  using Parent2::clean;
  using Parent2::solveCplex;
  
public:
  CutSolverRooted(const MwcsGraphType& mwcsGraph,
                  const NodeSet& rootNodes,
                  const Options& options,
                  const MwcsAnalyzeType& analysis)
    : GrandParent(mwcsGraph)
    , Parent1(mwcsGraph, rootNodes)
    , Parent2(mwcsGraph, options, analysis)
  {
  }
  
  virtual ~CutSolverUnrooted()
  {
  }
  
  bool solveCplex();
  
protected:
  virtual void initConstraints();
};

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void CutSolverRooted<GR, NWGHT, NLBL, EWGHT>::initConstraints()
{
  Parent2::initConstraints();
  
  char buf[1024];
  
  int i = 0;
  for (NodeSetIt nodeIt = _rootNodes.begin(); nodeIt != _rootNodes.end(); ++nodeIt, ++i)
  {
    // x_r = 1
    int r = (*_pNode)[*nodeIt];
    IloConstraint c;
    _model.add(c = (_x[r] == 1));
    
    snprintf(buf, 1024, "root_%d", i);
    c.setName(buf);
  }
  
  // if you pick a non-root node then at least one
  // of its direct neighbors must be part of the solution as well
  for (NodeIt i(g); i != lemon::INVALID; ++i)
  {
    if (i == _rootNodes.find(i) != _rootNodes.end())
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
    
    _model.add(_x[(*_pNode)[i]] <= expr);
  }
  
  for (NodeSetIt rootIt = _rootNodes.begin(); rootIt != _rootNodes.end(); ++rootIt)
  {
    // nodes i that are not in the same component as the root get x_i = 0
    const int rootComp = _mwcsGraph.getComponent(*rootIt);
    for (NodeIt i(g); i != lemon::INVALID; ++i)
    {
      if (_mwcsGraph.getComponent(i) != rootComp)
      {
        _model.add(_x[(*_pNode)[i]] == 0);
      }
    }
  }
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
  IloCplex::HeuristicCallbackI* pHeuristic = NULL;
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
//  _cplex.setParam( IloCplex::ZeroHalfCuts  , -1 );
  _cplex.setParam( IloCplex::MCFCuts       , -1 );
//  _cplex.setParam( IloCplex::AggFill       ,  0 );
//  _cplex.setParam( IloCplex::PreInd        ,  0 );
//  _cplex.setParam( IloCplex::RelaxPreInd   ,  0 );
//  _cplex.setParam( IloCplex::PreslvNd      , -1 );
//  _cplex.setParam( IloCplex::RepeatPresolve,  0 );
  _cplex.setParam( IloCplex::MIPEmphasis, IloCplex::MIPEmphasisBestBound );

  pLazyCut = new (_env) NodeCutRootedLazyConstraint<GR, NWGHT, NLBL, EWGHT>(_env, _x, g, weight, _root, *_pNode,
                                                                            _n, _options._maxNumberOfCuts, pMutex);
  pUserCut = new (_env) NodeCutRootedUserCut<GR, NWGHT, NLBL, EWGHT>(_env, _x, g, weight, _root, *_pNode,
                                                                     _n, _options._maxNumberOfCuts, pMutex,
                                                                     _options._backOff);
    
  pHeuristic = new (_env) HeuristicRootedType(_env, _x, //_z,
                                              g, weight, _root,
                                              *_pNode, //*_pEdge,
                                              _n, _m, pMutex);

  _cplex.setParam(IloCplex::MIPInterval, 1);

  IloCplex::Callback cb(pLazyCut);
  _cplex.use(cb);

  IloCplex::Callback cb2(pHeuristic);
  _cplex.use(cb2);

  IloCplex::Callback cb3(pUserCut);
  _cplex.use(cb3);
  
  // determine degrees
//  IntNodeMap deg(g, 0);
//  IntNodeMap posDeg(g, 0);
//  for (NodeIt i(g); i != lemon::INVALID; ++i)
//  {
//    for (IncEdgeIt e(g, i); e != lemon::INVALID; ++e)
//    {
//      Node j = g.oppositeNode(i, e);
//      ++deg[i];
//      if (weight[j] >= 0) ++posDeg[i];
//    }
//  }
//  
//  IloCplex::BranchCallbackI* pBranch = new Branch<GR, NWGHT, NLBL, EWGHT>(_env, _x, g, weight, *_pNode, _n, deg, posDeg);
//  IloCplex::Callback cb4(pBranch);
////  _cplex.use(cb4);
  
  bool res = _cplex.solve();
  cb.end();
  cb2.end();
  cb3.end();
  cb4.end();

  if (res)
  {
    std::cerr << "[" << _cplex.getObjValue() << ", "
              << _cplex.getBestObjValue() << "]" << std::endl;
    // print overview
    std::cerr << "# Cover cuts: " << _cplex.getNcuts(IloCplex::CutCover) << std::endl;
    std::cerr << "# GUB cover cuts: " << _cplex.getNcuts(IloCplex::CutGubCover) << std::endl;
    std::cerr << "# Flow cover cuts: " << _cplex.getNcuts(IloCplex::CutFlowCover) << std::endl;
    std::cerr << "# Clique cuts: " << _cplex.getNcuts(IloCplex::CutClique) << std::endl;
    std::cerr << "# Fractional cuts: " << _cplex.getNcuts(IloCplex::CutFrac) << std::endl;
    std::cerr << "# MCF cuts: " << _cplex.getNcuts(IloCplex::CutMCF) << std::endl;
    std::cerr << "# MIR cuts: " << _cplex.getNcuts(IloCplex::CutMir) << std::endl;
    std::cerr << "# Flow path cuts: " << _cplex.getNcuts(IloCplex::CutFlowPath) << std::endl;
    std::cerr << "# Implied bound cuts: " << _cplex.getNcuts(IloCplex::CutImplBd) << std::endl;
    std::cerr << "# Zero-half cuts: " << _cplex.getNcuts(IloCplex::CutZeroHalf) << std::endl;
    std::cerr << "# Local cover cuts: " << _cplex.getNcuts(IloCplex::CutLocalCover) << std::endl;
    std::cerr << "# Tighten cuts: " << _cplex.getNcuts(IloCplex::CutTighten) << std::endl;
    std::cerr << "# Obj disj cuts: " << _cplex.getNcuts(IloCplex::CutObjDisj) << std::endl;
    std::cerr << "# Lift-and-project cuts: " << _cplex.getNcuts(IloCplex::CutLiftProj) << std::endl;
    std::cerr << "# User cuts: " << _cplex.getNcuts(IloCplex::CutUser) << std::endl;
    std::cerr << "# Table cuts: " << _cplex.getNcuts(IloCplex::CutTable) << std::endl;
    std::cerr << "# Soln pool cuts: " << _cplex.getNcuts(IloCplex::CutSolnPool) << std::endl;
  }
  else
  {
    std::cerr << "[0, 0]" << std::endl;
  }

  //printVariables(std::cerr);
  delete pMutex;
  return res;
}

} // namespace mwcs
} // namespace nina

#endif // CUTSOLVERROOTED_H
