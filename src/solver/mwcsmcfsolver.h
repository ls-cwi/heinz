/*
 * mwcsscfsolver.h
 *
 *  Created on: 9-aug-2012
 *     Authors: C.I. Bucur, M. El-Kebir
 */

#ifndef MWCSMCFSOLVER_H
#define MWCSMCFSOLVER_H

#include "mwcsflowsolver.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class MwcsMCFSolver : public MwcsFlowSolver<GR, NWGHT, NLBL, EWGHT>
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;

  typedef MwcsFlowSolver<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> Parent;
  typedef typename Parent::MwcsGraphType MwcsGraphType;

  TEMPLATE_GRAPH_TYPEDEFS(Graph);

  typedef IloArray<IloBoolVarArray> IloBoolVarMatrix;

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
  MwcsMCFSolver(const MwcsGraphType& mwcsGraph);
  void printVariables(std::ostream& out);

protected:
  virtual void initVariables();
  virtual void initConstraints();
  void initVariablesNoRoot();
  void initConstraintsNoRoot();
  void initConstraintsRoot();

private:
  IloBoolVarMatrix _f;
  IloBoolVarMatrix _h;
};

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline MwcsMCFSolver<GR, NWGHT, NLBL, EWGHT>::MwcsMCFSolver(const MwcsGraphType& mwcsGraph)
  : Parent(mwcsGraph)
  , _f()
  , _h()
{
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsMCFSolver<GR, NWGHT, NLBL, EWGHT>::printVariables(std::ostream &out)
{
  const Graph& g = _mwcsGraph.getGraph();
  Parent::printVariables(out);

  // get solution
  for (int j = 0; j < _n; j++)
  {
    out << "Node x_" << j << " = " << _cplex.getValue(_x[j]);

    if (_root == lemon::INVALID)
      out << ", y = " << _cplex.getValue(_y[j])
          << ", w = " << _mwcsGraph.getScore(_invNode[j]);

    out << std::endl;
  }


  if (_root == lemon::INVALID)
  {
    for (int j = 0; j < _n; j++)
    {
      for (int k = 0; k < _n; k++)
        out << "h_" << j << "_" << k
            << " = " << _cplex.getValue(_h[j][k]) << std::endl;
    }
  }

  for (int ij = 0; ij < _m; ij++)
  {
    for (int k = 0; k < _n; k++)
      out << "f_" << ij << "_" << k << " = " << _cplex.getValue(_f[ij][k])
          << ", s = " << (*_pNode)[g.source(_invArc[ij])]
          << ", t = " << (*_pNode)[g.target(_invArc[ij])] << std::endl;
  }
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsMCFSolver<GR, NWGHT, NLBL, EWGHT>::initVariables()
{
  Parent::initVariables();

  char buf[1024];

  // flow variables f_{ijk} for arc (i,j) and commodity k
  _f = IloBoolVarMatrix(_env, _m);
  for (int ij = 0; ij < _m; ij++)
  {
    _f[ij] = IloBoolVarArray(_env, _n);
    for (int k = 0; k < _n; k++)
    {
      _f[ij][k] = IloBoolVar(_env);
      sprintf(buf, "f_%d_%d", ij, k);
      _f[ij][k].setName(buf);
    }
  }

  if (_root == lemon::INVALID)
    initVariablesNoRoot();
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsMCFSolver<GR, NWGHT, NLBL, EWGHT>::initVariablesNoRoot()
{
  char buf[1024];

  // artificial flow variables h_{jk} for node j and commodity k
  _h = IloBoolVarMatrix(_env, _n);
  for (int j = 0; j < _n; j++)
  {
    _h[j] = IloBoolVarArray(_env, _n);
    for (int k = 0; k < _n; k++)
    {
      _h[j][k] = IloBoolVar(_env);
      sprintf(buf, "h_%d_%d", j, k);
      _h[j][k].setName(buf);
    }
  }
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsMCFSolver<GR, NWGHT, NLBL, EWGHT>::initConstraintsRoot()
{
  const Graph& g = _mwcsGraph.getGraph();
  const int r = _root == lemon::INVALID ? -1 : (*_pNode)[_root];

  char buf[1024];

  // \sum_{j \in delta^+(r)} f_{rjk} = x_k for all k in V, k != r
  for (int k = 0; k < _n; k++)
  {
    if (k == r)
      continue;

    IloExpr expr(_env);

    for (OutArcIt a(g, _root); a != lemon::INVALID; ++a)
      expr += _f[(*_pArc)[a]][k];

    IloConstraint c1;
    _model.add(c1 = (expr == _x[k]));
    sprintf(buf, "root_flow_%d", k);
    c1.setName(buf);
  }


  // flow conservation
  for (int j = 0; j < _n; j++)
  {
    if (j == r)
      continue;

    for (int k = 0; k < _n; k++)
    {
      if (k == r || j == k)
        continue;

      IloExpr expr(_env);

      for (InArcIt a(g, _invNode[j]); a != lemon::INVALID; ++a)
        expr += _f[(*_pArc)[a]][k];

      for (OutArcIt a(g, _invNode[j]); a != lemon::INVALID; ++a)
        expr -= _f[(*_pArc)[a]][k];

      IloConstraint c3;
      _model.add(c3 = (expr == 0));
      sprintf(buf, "flow_cons_jk_%d_%d", j, k);
      c3.setName(buf);
    }
  }

  // flow conservation
  for (int k = 0; k < _n; k++)
  {
    if (k == r)
      continue;

    IloExpr expr(_env);

    for (InArcIt a(g, _invNode[k]); a != lemon::INVALID; ++a)
      expr += _f[(*_pArc)[a]][k];

    for (OutArcIt a(g, _invNode[k]); a != lemon::INVALID; ++a)
      expr -= _f[(*_pArc)[a]][k];

    expr -= _x[k];

    IloConstraint c4;
    _model.add(c4 = (expr == 0));
    sprintf(buf, "flow_cons_k_%d", k);
    c4.setName(buf);
  }
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsMCFSolver<GR, NWGHT, NLBL, EWGHT>::initConstraintsNoRoot()
{
  const Graph& g = _mwcsGraph.getGraph();
  char buf[1024];

  // maximum incoming flow for one good
  // h_{jk} <= y_j, for all nodes j and all commodities k in V
  for (int j = 0; j < _n; j++)
  {
    for (int k = 0; k < _n; k++)
    {
      IloConstraint c3;
      _model.add(c3 = (_h[j][k] <= _y[j]));
      sprintf(buf, "root_capacity_h_%d_%d", j, k);
      c3.setName(buf);
    }
  }

  IloExpr expr(_env);

  // flow conservation
  // incoming flow = outgoing flow, except for the root node
  // incoming flow of the root node is artificially created using h
  // h_j + \sum_{j \in \delta^-(j) f_{ij} = \sum_{k \in \delta^+(j) f_{jk} for all nodes j
  for (int j = 0; j < _n; j++)
  {
    for (int k = 0; k < _n; k++)
    {
      if (j != k)
      {
        expr.clear();
        expr += _h[j][k];

        for (InArcIt a(g, _invNode[j]); a != lemon::INVALID; ++a)
          expr += _f[(*_pArc)[a]][k];
        for (OutArcIt a(g, _invNode[j]); a != lemon::INVALID; ++a)
          expr -= _f[(*_pArc)[a]][k];

        IloConstraint c5;
        _model.add(c5 = (expr == 0));
        sprintf(buf, "flow_cons_%d_%d", j, k);
        c5.setName(buf);
      }
      else
      {
        expr.clear();
        expr += _h[k][k];

        for (InArcIt a(g, _invNode[k]); a != lemon::INVALID; ++a)
          expr += _f[(*_pArc)[a]][k];

        expr -= _x[k];

        for (OutArcIt a(g, _invNode[k]); a != lemon::INVALID; ++a)
          expr -= _f[(*_pArc)[a]][k];

        IloConstraint c6;
        _model.add(c6 = (expr == 0));
        sprintf(buf, "flow_cons_%d_%d", k, k);
        c6.setName(buf);
      }
    }
  }
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsMCFSolver<GR, NWGHT, NLBL, EWGHT>::initConstraints()
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
  // 0 <= f_{ijk} <= |V| * x_j, for all arcs (i,j) in A
  for (int ij = 0; ij < _m; ij++)
  {
    for (int k = 0; k < _n; k++)
    {
      int i = (*_pNode)[g.source(_invArc[ij])];

      //flow should not go through a node that is not in the solution
      IloConstraint c4;
      _model.add(c4 = (_f[ij][k] <= _x[i]));
      sprintf(buf, "capacity_%d_%d_%d", ij, k, i);
      c4.setName(buf);

      /*model.add(c4 = (_f[ij][k] <= _x[j]));
      sprintf(buf, "capacity_%d_%d_%d", ij, k, j);
      c4.setName(buf);*/

      // prevent the circulation of something that is not used
      /*_model.add(c4 = (_f[ij][k] <= _x[k]));
      sprintf(buf, "capacity_%d_%d_%d", ij, k, k);
      c4.setName(buf);*/
    }
  }
}

} // namespace mwcs
} // namespace nina

#endif // MWCSMCFSOLVER_H
