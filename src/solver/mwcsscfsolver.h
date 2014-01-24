/*
 * mwcsscfsolver.h
 *
 *  Created on: 7-aug-2012
 *     Authors: C.I. Bucur and Mohammed El-Kebir
 */

#ifndef MWCSSCFSOLVER_H
#define MWCSSCFSOLVER_H

#include "mwcsflowsolver.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class MwcsSCFSolver : public MwcsFlowSolver<GR, NWGHT, NLBL, EWGHT>
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;

  typedef MwcsFlowSolver<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> Parent;
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
  using Parent::_pArc;
  using Parent::_invArc;
  using Parent::_env;
  using Parent::_model;
  using Parent::_cplex;
  using Parent::_x;
  using Parent::_y;
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

public:
  MwcsSCFSolver(const MwcsGraphType& mwcsGraph);
  MwcsSCFSolver(const Graph& g,
                const WeightNodeMap& weight,
                const LabelNodeMap* pLabel);
  void printVariables(std::ostream& out);

protected:
  virtual void initVariables();
  virtual void initConstraints();
  void initVariablesNoRoot();
  void initConstraintsNoRoot();
  void initConstraintsRoot();

private:
  IloIntVarArray _f;
  IloIntVarArray _h;
};

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline MwcsSCFSolver<GR, NWGHT, NLBL, EWGHT>::MwcsSCFSolver(const MwcsGraphType& mwcsGraph)
  : Parent(mwcsGraph)
  , _f()
  , _h()
{
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline MwcsSCFSolver<GR, NWGHT, NLBL, EWGHT>::MwcsSCFSolver(const Graph& g,
                                        const WeightNodeMap& weight,
                                        const LabelNodeMap* pLabel)
  : Parent(g, weight, pLabel)
  , _f()
  , _h()
{
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsSCFSolver<GR, NWGHT, NLBL, EWGHT>::printVariables(std::ostream& out)
{
  const Graph& g = _mwcsGraph.getGraph();

  Parent::printVariables(out);

  if (_root == lemon::INVALID)
  {
    for (int i = 0; i < _n; i++)
    {
      out << "Node h_" << i << " = " << _cplex.getValue(_h[i]) << std::endl;
    }
  }

  for (ArcIt a(g); a != lemon::INVALID; ++a)
  {
    out << "f_" << (*_pArc)[a] << " = " << _cplex.getValue(_f[(*_pArc)[a]])
        << ", s = " << (*_pNode)[g.source(a)]
        << ", t = " << (*_pNode)[g.target(a)] << std::endl;
  }
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsSCFSolver<GR, NWGHT, NLBL, EWGHT>::initVariables()
{
  Parent::initVariables();

  char buf[1024];

  // flow variables for every arc
  _f = IloIntVarArray(_env, _m);
  for (int i = 0; i < _m; i++)
  {
    _f[i] = IloIntVar(_env);
    sprintf(buf, "f_%d", i);
    _f[i].setName(buf);
  }

  if (_root == lemon::INVALID)
    initVariablesNoRoot();
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsSCFSolver<GR, NWGHT, NLBL, EWGHT>::initVariablesNoRoot()
{
  char buf[1024];

  // flow variables for the root node
  // the root node has an incoming arc with flow attached
  // such that the flow conservation holds for this arc as well
  _h = IloIntVarArray(_env, _n);
  for (int i = 0; i < _n; i++)
  {
    _h[i] = IloIntVar(_env);
    sprintf(buf, "h_%d", i);
    _h[i].setName(buf);
  }
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsSCFSolver<GR, NWGHT, NLBL, EWGHT>::initConstraintsNoRoot()
{
  const Graph& g = _mwcsGraph.getGraph();

  char buf[1024];
  IloExpr expr(_env);

  // maximum incoming flow for the root node
  // 0 <= h_i <= |V| * y_i, for all nodes i in V
  for (int i = 0; i < _n; i++)
  {
    IloConstraint c3a, c3b;

    _model.add(c3a = (0 <= _h[i]));
    sprintf(buf, "root_capacity_a_%d", i);
    c3a.setName(buf);

    _model.add(c3b = (_h[i] <= _n * _y[i]));
    sprintf(buf, "root_capacity_b_%d", i);
    c3b.setName(buf);
  }

  // flow conservation
  // incoming flow = outgoing flow, except for the root node
  // incoming flow of the root node is artificially created using h
  // h_j + \sum_{j \in \delta^-(j) f_{ij} = \sum_{k \in \delta^+(j) f_{jk} for all nodes j
  for (int j = 0; j < _n; j++)
  {
    expr.clear();

    expr += _h[j];

    for (InArcIt a(g, _invNode[j]); a != lemon::INVALID; ++a)
      expr += _f[(*_pArc)[a]];
    for (OutArcIt  a(g, _invNode[j]); a != lemon::INVALID; ++a)
      expr -= _f[(*_pArc)[a]];

    IloConstraint c4;
    _model.add(c4 = (expr == _x[j]));
    sprintf(buf, "flow_cons_%d", j);
    c4.setName(buf);
  }

  // extra constraint for performance reasons
  // \sum_{i \in V} h_i = \sum_{i \in V} x_i
  /*expr.clear();
  for (int i = 0; i < n; i++)
  {
    expr += _h[i] - _x[i];
  }
  _model.add(expr == 0);*/
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsSCFSolver<GR, NWGHT, NLBL, EWGHT>::initConstraintsRoot()
{
  const Graph& g = _mwcsGraph.getGraph();
  const int r = _root == lemon::INVALID ? -1 : (*_pNode)[_root];

  char buf[1024];
  IloExpr expr(_env);

  // flow conservation
  // incoming flow = outgoing flow, except for the root node
  // incoming flow of the root node is artificially created using h
  // \sum_{j \in \delta^-(j) f_{ij} = \sum_{k \in \delta^+(j) f_{jk} + x_j
  // for all nodes j
  for (int j = 0; j < _n; j++)
  {
    if (j == r)
    {
      continue;
    }
    else
    {
      expr.clear();

      for (InArcIt a(g, _invNode[j]); a != lemon::INVALID; ++a)
        expr += _f[(*_pArc)[a]];
      for (OutArcIt  a(g, _invNode[j]); a != lemon::INVALID; ++a)
        expr -= _f[(*_pArc)[a]];

      IloConstraint c4;
      _model.add(c4 = (expr == _x[j]));
      sprintf(buf, "flow_cons_%d", j);
      c4.setName(buf);
    }
  }

  // flow conservation for the root node
  expr.clear();

  for (int j = 0; j < _n; j++)
    expr += _x[j];
  for (OutArcIt  a(g, _root); a != lemon::INVALID; ++a)
    expr -= _f[(*_pArc)[a]];

  IloConstraint c5;
  _model.add(c5 = (expr == _x[r]));
  sprintf(buf, "flow_cons_root_%d", r);
  c5.setName(buf);
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsSCFSolver<GR, NWGHT, NLBL, EWGHT>::initConstraints()
{
  const Graph& g = _mwcsGraph.getGraph();
  Parent::initConstraints();

  char buf[1024];

  if (_root == lemon::INVALID)
  {
    initConstraintsNoRoot();
  }
  else
  {
    initConstraintsRoot();
  }

  // flow should not exceed capacity
  // 0 <= f_{ij} <= |V| * x_j, for all arcs (i,j) in A
  for (int ij = 0; ij < _m; ij++)
  {
    IloConstraint c5a, c5b;

    _model.add(c5a = (0 <= _f[ij]));
    sprintf(buf, "capacity_a_%d", ij);
    c5a.setName(buf);

    _model.add(c5b = (_f[ij] <= _n * _x[(*_pNode)[g.target(_invArc[ij])]]));
    sprintf(buf, "capacity_b_%d", ij);
    c5b.setName(buf);
  }
}

} // namespace mwcs
} // namespace nina

#endif // MWCSSCFSOLVER_H
