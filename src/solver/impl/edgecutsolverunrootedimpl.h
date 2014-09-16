/*
 * edgecutsolverunrootedimpl.h
 *
 *  Created on: 10-sep-2014
 *      Author: M. El-Kebir
 */

#ifndef EDGECUTSOLVERUNROOTEDIMPL_H
#define EDGECUTSOLVERUNROOTEDIMPL_H

#include "solverunrootedimpl.h"
#include "edgecutsolverimpl.h"
#include "cplex_cut/edgecut.h"
#include "cplex_cut/edgecutlazy.h"
#include "cplex_cut/edgecutuser.h"
#include "cplex_cut/edgecutunrooted.h"

#include <ilconcert/ilothread.h>

namespace nina {
namespace mwcs {
  
template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class EdgeCutSolverUnrootedImpl : public SolverUnrootedImpl<GR, NWGHT, NLBL, EWGHT>,
                                  public EdgeCutSolverImpl<GR, NWGHT, NLBL, EWGHT>
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;
  
  typedef SolverUnrootedImpl<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> Parent1;
  typedef EdgeCutSolverImpl<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> Parent2;
  
  typedef typename Parent1::MwcsGraphType MwcsGraphType;
  typedef typename Parent1::NodeSet NodeSet;
  typedef typename Parent1::NodeSetIt NodeSetIt;
  
  typedef typename Parent2::NodeVector NodeVector;
  typedef typename Parent2::NodeVectorIt NodeVectorIt;
  typedef typename Parent2::MwcsAnalyzeType MwcsAnalyzeType;
  typedef typename Parent2::InvNodeIntMap InvNodeIntMap;
  typedef typename Parent2::InvArcIntMap InvArcIntMap;
  typedef typename Parent2::Options Options;
  
//  typedef HeuristicUnrooted<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> HeuristicUnrootedType;
  
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
  using Parent2::_z;
  using Parent2::_d;
  using Parent2::_pToD;
  using Parent2::_toG;
  using Parent2::_weight;
  using Parent2::_arcCost;
  using Parent2::_label;
  using Parent2::_diNode;
  using Parent2::_diArc;
  using Parent2::_invDiNode;
  using Parent2::_invDiArc;

  using Parent2::initVariables;
  using Parent2::initConstraints;
  using Parent2::clean;
  using Parent2::solveCplex;
  using Parent2::printVariables;
  using Parent2::exportModel;
  
public:
  EdgeCutSolverUnrootedImpl(const Options& options)
    : Parent1()
    , Parent2(options)
    , _root(lemon::INVALID)
  {
  }
  
  virtual ~EdgeCutSolverUnrootedImpl()
  {
  }
  
  void init(const MwcsGraphType& mwcsGraph)
  {
    Parent1::init(mwcsGraph);
    initGraph(mwcsGraph);
    initVariables(mwcsGraph);
    initConstraints(mwcsGraph);
  }
  
  bool solve(double& score, BoolNodeMap& solutionMap, NodeSet& solutionSet)
  {
    return Parent2::solveCplex(*_pMwcsGraph, score, solutionMap, solutionSet);
  }
  
protected:
  typedef typename Parent2::DiGraph DiGraph;
  typedef typename Parent2::DiNode DiNode;
  typedef typename Parent2::DiNodeIt DiNodeIt;
  typedef typename Parent2::DiArc DiArc;
  typedef typename Parent2::DiArcIt DiArcIt;
  typedef typename Parent2::DiInArcIt DiInArcIt;
  typedef typename Parent2::DiOutArcIt DiOutArcIt;
  typedef typename Parent2::DiNodeSet DiNodeSet;
  typedef typename Parent2::DiNodeVector DiNodeVector;
  typedef typename Parent2::DiNodeVectorIt DiNodeVectorIt;
  typedef typename Parent2::DiArcVector DiArcVector;
  typedef typename Parent2::DiArcVectorIt DiArcVectorIt;
  typedef typename Parent2::InvDiNodeIntMap InvDiNodeIntMap;
  typedef typename Parent2::InvDiArcIntMap InvDiArcIntMap;
  typedef typename Parent2::IntDiNodeMap IntDiNodeMap;
  typedef typename Parent2::IntDiArcMap IntDiArcMap;
  typedef typename Parent2::DoubleDiNodeMap DoubleDiNodeMap;
  typedef typename Parent2::LabelDiNodeMap LabelDiNodeMap;
  typedef typename Parent2::NodeDiNodeMap NodeDiNodeMap;
  typedef typename Parent2::DiNodeNodeMap DiNodeNodeMap;
  typedef EdgeCutUnrootedLazyConstraint<DiGraph, DoubleDiNodeMap> EdgeCutUnrootedLazyConstraintType;
  typedef EdgeCutUnrootedUserCut<DiGraph, DoubleDiNodeMap> EdgeCutUnrootedUserCutType;
  
  bool solveModel();
  void initGraph(const MwcsGraphType& mwcsGraph);
  virtual void initConstraints(const MwcsGraphType& mwcsGraph);
  
  DiNode _root;
};

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void EdgeCutSolverUnrootedImpl<GR, NWGHT, NLBL, EWGHT>::initConstraints(const MwcsGraphType& mwcsGraph)
{
  Parent2::initConstraints(mwcsGraph);

  IloExpr expr(_env);
  for (DiNodeIt v(_d); v != lemon::INVALID; ++v)
  {
    if (v == _root)
      continue;
    
    expr.clear();
    for (DiInArcIt a(_d, v); a != lemon::INVALID; ++a)
    {
      expr += _z[_diArc[a]];
    }
    
    IloConstraint c;
    _model.add(c = (expr == _x[_diNode[v]]));
    c.setName("one_in_arc");
  }
  
  expr.clear();
  for (DiOutArcIt a(_d, _root); a != lemon::INVALID; ++a)
  {
    int id_i = _diNode[_d.target(a)];
    expr += _z[_diArc[a]];
    
    // symmetry breaking
    for (int id_j = 0; id_j < id_i; ++id_j)
    {
      DiNode j = _invDiNode[id_j];
      if (_weight[j] < 0 || j == _root) continue;
      _model.add(_z[_diArc[a]] <= 1 - _x[id_j]);
    }
  }
  
  IloConstraint c2;
  _model.add(c2 = (expr == 1));
  c2.setName("one_out_arc_from_root");
  
  IloConstraint c3;
  _model.add(c3 = (_x[_diNode[_root]] == 1));
  c3.setName("root_in_sol");
}
  
template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void EdgeCutSolverUnrootedImpl<GR, NWGHT, NLBL, EWGHT>::initGraph(const MwcsGraphType& mwcsGraph)
{
  const Graph& g = mwcsGraph.getGraph();
  const DoubleNodeMap& score = mwcsGraph.getScores();
  const LabelNodeMap& label = mwcsGraph.getLabels();
  
  // let's construct the graph d
  delete _pToD;
  _pToD = new DiNodeNodeMap(g);
  
  _d.clear();
  _d.reserveNode(mwcsGraph.getNodeCount() + 1);
  _d.reserveArc(mwcsGraph.getArcCount() + mwcsGraph.getNodeCount());
  
  // determine _arcCost
  _arcCost = lemon::mapMinValue(g, score);
  
  // introduce nodes
  _root = _d.addNode();
  _toG[_root] = lemon::INVALID;
  _weight[_root] = 0;
  _label[_root] = "UBER_ROOT";
  
  for (NodeIt v(g); v != lemon::INVALID; ++v)
  {
    DiNode di_v = _d.addNode();
    _toG[di_v] = v;
    _pToD->set(v, di_v);
    
    _weight[di_v] = score[v] - _arcCost;
    _label[di_v] = label[v];
  }
  
  // introduce arcs
  for (ArcIt a(g); a != lemon::INVALID; ++a)
  {
    Node u = g.source(a);
    Node v = g.target(a);
    _d.addArc((*_pToD)[u], (*_pToD)[v]);
  }
  
  // introduce arcs from root to positively weighted nodes
  for (NodeIt v(g); v != lemon::INVALID; ++v)
  {
    if (score[v] > 0)
    {
      _d.addArc(_root, (*_pToD)[v]);
    }
  }
  
  _n = lemon::countNodes(_d);
  _m = lemon::countArcs(_d);
}
  
template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline bool EdgeCutSolverUnrootedImpl<GR, NWGHT, NLBL, EWGHT>::solveModel()
{
  IloFastMutex* pMutex = NULL;
  if (_options._multiThreading > 1)
  {
    
    pMutex = new IloFastMutex();
  }

  IloCplex::LazyConstraintCallbackI* pLazyCut = NULL;
  IloCplex::UserCutCallbackI* pUserCut = NULL;
  IloCplex::HeuristicCallbackI* pHeuristic = NULL;

  _cplex.setParam( IloCplex::HeurFreq      , -1 );
  _cplex.setParam( IloCplex::MCFCuts       , -1 );
  _cplex.setParam( IloCplex::MIPEmphasis, IloCplex::MIPEmphasisBestBound );

  pLazyCut = new (_env) EdgeCutUnrootedLazyConstraintType(_env, _x, _z, _d, _root, _weight, _diNode, _diArc,
                                                          _n, _m, _options._maxNumberOfCuts, pMutex);
  pUserCut = new (_env) EdgeCutUnrootedUserCutType(_env, _x, _z, _d, _root, _weight, _diNode, _diArc,
                                                   _n, _m, _options._maxNumberOfCuts, pMutex,
                                                   _options._backOff);
//
//  pHeuristic = new (_env) HeuristicUnrootedType(_env, _x, _y, //_z,
//                                                g, weight,
//                                                *_pNode, //*_pEdge,
//                                                _n, _m, pMutex);

  _cplex.setParam(IloCplex::MIPInterval, 1);
  
  IloCplex::Callback cb(pLazyCut);
  _cplex.use(cb);
//
//  IloCplex::Callback cb2(pHeuristic);
//  _cplex.use(cb2);
//
  IloCplex::Callback cb3(pUserCut);
//  _cplex.use(cb3);
  
  bool res = _cplex.solve();
  cb.end();
//  cb2.end();
  cb3.end();

  if (res)
  {
    std::cerr << "[" << _cplex.getObjValue() << ", "
              << _cplex.getBestObjValue() << "]" << std::endl;
  }
  else
  {
    std::cerr << "[0, 0]" << std::endl;
  }
  
  delete pMutex;
  return res;
}

} // namespace mwcs
} // namespace nina

#endif // EDGECUTSOLVERUNROOTEDIMPL_H
