/*
 * cutsolverunrootedimpl.h
 *
 *  Created on: 22-aug-2014
 *      Author: M. El-Kebir
 */

#ifndef CUTSOLVERUNROOTEDIMPL_H
#define CUTSOLVERUNROOTEDIMPL_H

#include "solverunrootedimpl.h"
#include "cplexsolverimpl.h"
#include "cplex_cut/nodecutunrooted.h"
#include "cplex_heuristic/heuristicunrooted.h"
#include "cplex_incumbent/incumbent.h"
#include "cplex_incumbent/pcstincumbent.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class CutSolverUnrootedImpl : public SolverUnrootedImpl<GR, NWGHT, NLBL, EWGHT>,
                              public CplexSolverImpl<GR, NWGHT, NLBL, EWGHT>
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;
  
  typedef SolverUnrootedImpl<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> Parent1;
  typedef CplexSolverImpl<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> Parent2;

  typedef typename Parent1::MwcsGraphType MwcsGraphType;
  typedef typename Parent1::NodeSet NodeSet;
  typedef typename Parent1::NodeSetIt NodeSetIt;
  
  typedef typename Parent2::NodeVector NodeVector;
  typedef typename Parent2::NodeVectorIt NodeVectorIt;
  typedef typename Parent2::MwcsAnalyzeType MwcsAnalyzeType;
  typedef typename Parent2::InvNodeIntMap InvNodeIntMap;
  typedef typename Parent2::InvArcIntMap InvArcIntMap;
  typedef typename Parent2::Options Options;
  
  typedef NodeCutUnrootedLazyConstraint<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> NodeCutUnrootedLazyConstraintType;
  typedef NodeCutUnrootedUserCut<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> NodeCutUnrootedUserCutType;
  typedef HeuristicUnrooted<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> HeuristicUnrootedType;
  typedef PcstIncumbent<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap>  PcstIncumbentType;
  
  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  
  using Parent1::_pMwcsGraph;
  using Parent2::_options;
  using Parent2::_n;
  using Parent2::_m;
  using Parent2::_pNode;
  using Parent2::_invNode;
  using Parent2::_env;
  using Parent2::_model;
  using Parent2::_cplex;
  using Parent2::_x;
  using Parent2::initVariables;
  using Parent2::initConstraints;
  using Parent2::clean;
  using Parent2::solveCplex;
  
public:
  CutSolverUnrootedImpl(const Options& options)
    : Parent1()
    , Parent2(options)
    , _y()
  {
  }
  
  virtual ~CutSolverUnrootedImpl()
  {
  }
  
  virtual void printVariables(const MwcsGraphType& mwcsGraph,
                              std::ostream& out)
  {
    for (int id_v = 0; id_v < _n; id_v++)
    {
      out << "Node '" << mwcsGraph.getLabel(_invNode[id_v])
      << "', x_" << id_v << " = " << _cplex.getValue(_x[id_v])
      << ", y = " << _cplex.getValue(_y[id_v])
      << ", w = " << mwcsGraph.getScore(_invNode[id_v])
      << std::endl;
    }
  }
  
  void init(const MwcsGraphType& mwcsGraph)
  {
    Parent1::init(mwcsGraph);
    initVariables(mwcsGraph);
    initConstraints(mwcsGraph);
  }
  
  bool solve(double& score, double& scoreUB, BoolNodeMap& solutionMap, NodeSet& solutionSet)
  {
    return Parent2::solveCplex(*_pMwcsGraph, score, scoreUB, solutionMap, solutionSet);
  }

protected:
  IloBoolVarArray _y;
  
  virtual void initVariables(const MwcsGraphType& mwcsGraph);
  virtual void initConstraints(const MwcsGraphType& mwcsGraph);
  
  bool solveModel();
};
  
template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void CutSolverUnrootedImpl<GR, NWGHT, NLBL, EWGHT>::initVariables(const MwcsGraphType& mwcsGraph)
{
  Parent2::initVariables(mwcsGraph);
  
  _n = mwcsGraph.getNodeCount();
  _y = IloBoolVarArray(_env, _n);
  
  char buf[1024];
  int i = 0;
  for (NodeVectorIt it = _invNode.begin(); it != _invNode.end(); ++it, ++i)
  {
    Node v = *it;

    // y_i = 0 if node i is not the root node
    // y_i = 1 if node i is picked as the root node
    snprintf(buf, 1024, "y_%s", mwcsGraph.getLabel(v).c_str());
    //snprintf(buf, 1024, "y_%d", g.id(v));
    _y[i].setName(buf);
  }
}
  
template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void CutSolverUnrootedImpl<GR, NWGHT, NLBL, EWGHT>::initConstraints(const MwcsGraphType& mwcsGraph)
{
  Parent2::initConstraints(mwcsGraph);
  
  const Graph& g = mwcsGraph.getGraph();
  const WeightNodeMap& weight = mwcsGraph.getScores();
  
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
  
  // TODO maybe equality up there, but need to double check whether there are positive nodes
  
  // the root node has to be one of the selected nodes
  // in the final graph
  // y_i <= x_i for all nodes i in V
  for (int i = 0; i < _n; i++)
  {
    IloConstraint c2;
    _model.add(c2 = (_y[i] <= _x[i]));
    c2.setName("root_in_x");
  }
  
  // root node has to be positive
  expr.clear();
  for (int i = 0; i < _n; i++)
  {
    double weight_i = weight[_invNode[i]];
    if (weight_i < 0)
    {
      expr += _y[i];
    }
  }
  _model.add(expr == 0);
  
  // objective must be positive
  expr.clear();
  for (int i = 0; i < _n ; i++)
  {
    expr += _x[i] * weight[_invNode[i]];
  }
  _model.add(expr >= 0);
  
  // if you pick a node then it must be the root node
  // or at least one of its direct neighbors must be part
  // of the solution as well
  for (NodeIt i(g); i != lemon::INVALID; ++i)
  {
    expr.clear();
    for (IncEdgeIt e(g, i); e != lemon::INVALID; ++e)
    {
      Node j = g.oppositeNode(i, e);
      expr += _x[(*_pNode)[j]];
      
      // if i is negative then its positive neighbors must be in
      if (weight[i] < 0 && weight[j] > 0)
        _model.add(_x[(*_pNode)[i]] <= _x[(*_pNode)[j]]);
    }
    
    expr += _y[(*_pNode)[i]];
    
    _model.add(_x[(*_pNode)[i]] <= expr);
  }
  
  // if you pick a negative node, then at least two of its direct neighbors
  // must be part of the solution as well
  // if you get in, you have to get out as well
  // BIG FAT WARNING: not true for xHeinz!!!
  if (_n < 5000 || !_options._pcst)
  {
    if (g_verbosity >= VERBOSE_DEBUG)
    {
      std::cout << std::endl;
    }
    int idx = 1;
    for (NodeIt i(g); i != lemon::INVALID; ++i, ++idx)
    {
      if (g_verbosity >= VERBOSE_DEBUG)
      {
        std::cout << "\rAdding symmetry breaking constraints " << idx << "/" << _n << std::flush;
      }
      if (weight[i] <= 0)
      {
        expr.clear();
        for (IncEdgeIt e(g, i); e != lemon::INVALID; ++e)
        {
          Node j = g.oppositeNode(i, e);
          expr += _x[(*_pNode)[j]];
        }
        _model.add(2 * _x[(*_pNode)[i]] <= expr);
      }
      else
      {
        // symmetry breaking
        // sum_{j > i, w_j > 0} y_j <= 1 - x_i
        int id_i = (*_pNode)[i];
        expr.clear();
        for (int id_j = id_i + 1; id_j < _n; ++id_j)
        {
          Node j = _invNode[id_j];
          if (weight[j] < 0) continue;
          expr += _y[id_j];
          _model.add(_y[id_j] <= 1 - _x[id_i]);
        }
        _model.add(expr <= 1 - _x[id_i]);
      }
    }
    if (g_verbosity >= VERBOSE_DEBUG)
    {
      std::cout << std::endl;
    }
  }
}
  
template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline bool CutSolverUnrootedImpl<GR, NWGHT, NLBL, EWGHT>::solveModel()
{
  const Graph& g = _pMwcsGraph->getGraph();
  const WeightNodeMap& weight = _pMwcsGraph->getScores();

  IloFastMutex* pMutex = NULL;
  if (_options._multiThreading > 1)
  {
    pMutex = new IloFastMutex();
  }

  IloCplex::LazyConstraintCallbackI* pLazyCut = NULL;
  IloCplex::UserCutCallbackI* pUserCut = NULL;
  IloCplex::HeuristicCallbackI* pHeuristic = NULL;
  IloCplex::IncumbentCallbackI* pIncumbent = NULL;

  _cplex.setParam( IloCplex::HeurFreq      , -1 );
  _cplex.setParam( IloCplex::Cliques       , -1 );
//  _cplex.setParam( IloCplex::Covers        , -1 );
//  _cplex.setParam( IloCplex::FlowCovers    , -1 );
//  _cplex.setParam( IloCplex::GUBCovers     , -1 );
//  _cplex.setParam( IloCplex::FracCuts      , -1 );
//  _cplex.setParam( IloCplex::MIRCuts       , -1 );
//  _cplex.setParam( IloCplex::FlowPaths     , -1 );
//  _cplex.setParam( IloCplex::ImplBd        , -1 );
//  _cplex.setParam( IloCplex::DisjCuts      , -1 );
//  _cplex.setParam( IloCplex::ZeroHalfCuts  , -1 );
//    _cplex.setParam( IloCplex::LiftProjCuts  , -1 );
  _cplex.setParam( IloCplex::MCFCuts       , -1 );
  _cplex.setParam( IloCplex::MIPEmphasis, IloCplex::MIPEmphasisBestBound );
//    _cplex.setParam( IloCplex::AggFill       ,  0 );
//    _cplex.setParam( IloCplex::PreInd        ,  0 );
//    _cplex.setParam( IloCplex::RelaxPreInd   ,  0 );
//    _cplex.setParam( IloCplex::PreslvNd      , -1 );
//    _cplex.setParam( IloCplex::RepeatPresolve,  0 );

  pLazyCut = new (_env) NodeCutUnrootedLazyConstraint<GR, NWGHT, NLBL, EWGHT>(_env, _x, _y, g, weight, *_pNode,
                                                                              _n, _options._maxNumberOfCuts, pMutex);
  pUserCut = new (_env) NodeCutUnrootedUserCut<GR, NWGHT, NLBL, EWGHT>(_env, _x, _y, g, weight, *_pNode,
                                                                       _n, _options._maxNumberOfCuts, pMutex,
                                                                       _options._backOff);

  pHeuristic = new (_env) HeuristicUnrootedType(_env, _x, _y, //_z,
                                                g, weight,
                                                *_pNode, //*_pEdge,
                                                _n, _m, pMutex);
  
  if (g_pOut)
  {
    if (_options._pcst)
      pIncumbent = new (_env) PcstIncumbentType(_env, _pMwcsGraph->getTotalNodeProfitPCST(), pMutex);
    else
      pIncumbent = new (_env) Incumbent(_env, pMutex);
  }

  _cplex.setParam(IloCplex::MIPInterval, 1);
  
  IloCplex::Callback cb(pLazyCut);
  _cplex.use(cb);

  IloCplex::Callback cb2(pHeuristic);
  _cplex.use(cb2);

  IloCplex::Callback cb3(pUserCut);
  _cplex.use(cb3);
  
  IloCplex::Callback cb4(pIncumbent);
  if (pIncumbent)
    _cplex.use(cb4);
  
//  // determine degrees
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
  
  //exportModel("/tmp/model.lp");
  bool res = _cplex.solve();
  cb.end();
  cb2.end();
  cb3.end();
  if (g_pOut)
  {
    cb4.end();
  }
  
  if (res)
  {
    if (g_verbosity > VERBOSE_NONE)
    {
      std::cerr << "[" << _cplex.getObjValue() << ", "
                << _cplex.getBestObjValue() << "]" << std::endl;
    }
    // print overview
//    std::cerr << "# Cover cuts: " << _cplex.getNcuts(IloCplex::CutCover) << std::endl;
//    std::cerr << "# GUB cover cuts: " << _cplex.getNcuts(IloCplex::CutGubCover) << std::endl;
//    std::cerr << "# Flow cover cuts: " << _cplex.getNcuts(IloCplex::CutFlowCover) << std::endl;
//    std::cerr << "# Clique cuts: " << _cplex.getNcuts(IloCplex::CutClique) << std::endl;
//    std::cerr << "# Fractional cuts: " << _cplex.getNcuts(IloCplex::CutFrac) << std::endl;
//    std::cerr << "# MCF cuts: " << _cplex.getNcuts(IloCplex::CutMCF) << std::endl;
//    std::cerr << "# MIR cuts: " << _cplex.getNcuts(IloCplex::CutMir) << std::endl;
//    std::cerr << "# Flow path cuts: " << _cplex.getNcuts(IloCplex::CutFlowPath) << std::endl;
//    std::cerr << "# Implied bound cuts: " << _cplex.getNcuts(IloCplex::CutImplBd) << std::endl;
//    std::cerr << "# Zero-half cuts: " << _cplex.getNcuts(IloCplex::CutZeroHalf) << std::endl;
//    std::cerr << "# Local cover cuts: " << _cplex.getNcuts(IloCplex::CutLocalCover) << std::endl;
//    std::cerr << "# Tighten cuts: " << _cplex.getNcuts(IloCplex::CutTighten) << std::endl;
//    std::cerr << "# Obj disj cuts: " << _cplex.getNcuts(IloCplex::CutObjDisj) << std::endl;
//    std::cerr << "# Lift-and-project cuts: " << _cplex.getNcuts(IloCplex::CutLiftProj) << std::endl;
//    std::cerr << "# User cuts: " << _cplex.getNcuts(IloCplex::CutUser) << std::endl;
//    std::cerr << "# Table cuts: " << _cplex.getNcuts(IloCplex::CutTable) << std::endl;
//    std::cerr << "# Soln pool cuts: " << _cplex.getNcuts(IloCplex::CutSolnPool) << std::endl;
  }
  else
  {
    if (g_verbosity > VERBOSE_NONE)
    {
      std::cerr << "[0, 0]" << std::endl;
    }
  }

  //printVariables(std::cerr);
  delete pMutex;
  return res;
}


} // namespace mwcs
} // namespace nina

#endif // CUTSOLVERUNROOTEDIMPL_H
