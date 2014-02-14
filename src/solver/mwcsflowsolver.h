/*
 * mwcsflowsolver.cpp
 *
 *  Created on: 10-aug-2012
 *     Authors: M. El-Kebir
 */

#ifndef MWCSFLOWSOLVER_H
#define MWCSFLOWSOLVER_H

#include "mwcssolver.h"

// ILOG stuff
#include <ilconcert/iloalg.h>
#include <ilcplex/ilocplex.h>

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class MwcsFlowSolver : public MwcsSolver<GR, NWGHT, NLBL, EWGHT>
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;

  typedef MwcsSolver<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> Parent;
  typedef typename Parent::MwcsGraphType MwcsGraphType;

  TEMPLATE_GRAPH_TYPEDEFS(Graph);

  typedef std::vector<Node> InvNodeIntMap;
  typedef std::vector<Arc> InvArcIntMap;
  typedef typename Parent::NodeVector NodeVector;
  typedef typename Parent::NodeVectorIt NodeVectorIt;

  using Parent::_mwcsGraph;
  using Parent::_root;
  using Parent::_score;
  using Parent::_solutionMap;
  using Parent::_solutionSet;
  using Parent::printSolution;
  using Parent::getSolutionWeight;
  using Parent::getSolutionNodeMap;
  using Parent::getSolutionModule;
  using Parent::isNodeInSolution;
  using Parent::init;

public:
  MwcsFlowSolver(const MwcsGraphType& mwcsGraph);
  virtual ~MwcsFlowSolver();

  void exportModel(const std::string& filename)
  {
    _cplex.exportModel(filename.c_str());
  }

  virtual bool solve();
  virtual void init(Node root);
  virtual void printVariables(std::ostream& out);

protected:
  int _n;
  int _m;
  IntNodeMap* _pNode;
  InvNodeIntMap _invNode;
  IntArcMap* _pArc;
  InvArcIntMap _invArc;
  IloEnv _env;
  IloModel _model;
  IloCplex _cplex;
  IloBoolVarArray _x;
  IloBoolVarArray _y;
  double _LB;

  virtual void initVariables();
  virtual void initConstraints();
  virtual void clean();
  virtual void setLowerBound(double LB)
  {
    _LB = LB;
  }

  virtual bool solveCplex()
  {
    return _cplex.solve();
  }

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
inline MwcsFlowSolver<GR, NWGHT, NLBL, EWGHT>::MwcsFlowSolver(const MwcsGraphType& mwcsGraph)
  : Parent(mwcsGraph)
  , _n(mwcsGraph.getNodeCount())
  , _m(mwcsGraph.getArcCount())
  , _pNode(NULL)
  , _invNode()
  , _pArc(NULL)
  , _invArc()
  , _env()
  , _model(_env)
  , _cplex(_model)
  , _x()
  , _y()
  , _LB(0)
{
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline MwcsFlowSolver<GR, NWGHT, NLBL, EWGHT>::~MwcsFlowSolver()
{
  _env.end();
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
void MwcsFlowSolver<GR, NWGHT, NLBL, EWGHT>::clean()
{
  _cplex.end();
  _env.end();
  _env = IloEnv();
  _model = IloModel(_env);
  _cplex = IloCplex(_model);
  delete _pArc;
  delete _pNode;
  _pArc = NULL;
  _pNode = NULL;
}


template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsFlowSolver<GR, NWGHT, NLBL, EWGHT>::init(Node root)
{
  _root = root;

  initVariables();
  initConstraints();
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsFlowSolver<GR, NWGHT, NLBL, EWGHT>::printVariables(std::ostream& out)
{
  for (int id_v = 0; id_v < _n; id_v++)
  {
    out << "Node '" << _mwcsGraph.getLabel(_invNode[id_v])
        << "', x_" << id_v << " = " << _cplex.getValue(_x[id_v]);

    if (_root == lemon::INVALID)
      out << ", y = " << _cplex.getValue(_y[id_v])
          << ", w = " << _mwcsGraph.getScore(_invNode[id_v]);

    out << std::endl;
  }
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsFlowSolver<GR, NWGHT, NLBL, EWGHT>::initVariables()
{
  const Graph& g = _mwcsGraph.getGraph();
  _n = _mwcsGraph.getNodeCount();
  _m = _mwcsGraph.getArcCount();

  delete _pNode;
  delete _pArc;
  _pNode = new IntNodeMap(g);
  _pArc = new IntArcMap(g);

  _x = IloBoolVarArray(_env, _n);
  if (_root == lemon::INVALID)
    _y = IloBoolVarArray(_env, _n);


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
    _x[i].setName(buf);

    if (_root == lemon::INVALID)
    {
      // y_i = 0 if node i is not the root node
      // y_i = 1 if node i is picked as the root node
      snprintf(buf, 1024, "y_%s", _mwcsGraph.getLabel(v).c_str());
      _y[i].setName(buf);
    }

    (*_pNode)[v] = i;
  }

  _invArc.clear();
  _invArc.reserve(_m);
  for (ArcIt a(g); a != lemon::INVALID; ++a)
  {
    (*_pArc)[a] = _invArc.size();
    _invArc.push_back(a);
  }
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsFlowSolver<GR, NWGHT, NLBL, EWGHT>::initConstraints()
{
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
    for (int i = 0; i < _n; i++)
    {
      double weight = _mwcsGraph.getScore(_invNode[i]);
      if (weight < 0)
      {
        _model.add(_y[i] == 0);
      }
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
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline bool MwcsFlowSolver<GR, NWGHT, NLBL, EWGHT>::solve()
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

  // objective function
  IloExpr expr(_env);
  for (int i = 0; i < _n ; i++)
  {
    expr += _x[i] * _mwcsGraph.getScore(_invNode[i]);
  }
  _model.add(IloObjective(_env, expr, IloObjective::Maximize));

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

  lemon::Tolerance<double> tol(1e-6);

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

#endif // MWCSFLOWSOLVER_H
