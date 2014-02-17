/*
 * mwcscutsolverheuristic.h
 *
 *  Created on: 2-jul-2013
 *      Author: M. El-Kebir
 */

#ifndef MWCSCUTSOLVERHEURISTIC_H
#define MWCSCUTSOLVERHEURISTIC_H

#include <ilcplex/ilocplex.h>
#include <ilcplex/ilocplex.h>
#include <ilconcert/ilothread.h>
#include <lemon/adaptors.h>
#include <lemon/bfs.h>
#include "../mwcstreesolver.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class MwcsCutSolverHeuristic : public IloCplex::HeuristicCallbackI
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;
  TEMPLATE_GRAPH_TYPEDEFS(Graph);

  MwcsCutSolverHeuristic(IloEnv env,
                         IloBoolVarArray x,
                         IloBoolVarArray y,
                         const Graph& g,
                         const WeightNodeMap& weight,
                         Node root,
                         const IntNodeMap& nodeMap,
                         int n,
                         int m,
                         IloFastMutex* pMutex);
  MwcsCutSolverHeuristic(const MwcsCutSolverHeuristic& other);
  ~MwcsCutSolverHeuristic();

protected:
  void main();
  IloCplex::CallbackI* duplicateCallback() const
  {
    return (new (getEnv()) MwcsCutSolverHeuristic(*this));
  }

private:
  typedef lemon::FilterEdges<const Graph, const BoolEdgeMap> SubGraphType;
  typedef MwcsGraph<const SubGraphType, const WeightNodeMap, LabelNodeMap, DoubleEdgeMap> MwcsSubGraphType;
  typedef MwcsTreeSolver<const SubGraphType, const WeightNodeMap, LabelNodeMap, DoubleEdgeMap> MwcsSubTreeSolver;

  IloBoolVarArray _x;
  IloBoolVarArray _y;
  const Graph& _g;
  const WeightNodeMap& _weight;
  const Node _root;
  const IntNodeMap& _nodeMap;
  const int _n;
  const int _m;
  DoubleEdgeMap* _pEdgeCost;
  BoolEdgeMap* _pEdgeFilterMap;
  const SubGraphType* _pSubG;
  MwcsSubGraphType* _pMwcsSubGraph;
  MwcsSubTreeSolver* _pMwcsSubTreeSolver;
  IloFastMutex* _pMutex;

  void lock()
  {
    if (_pMutex)
      _pMutex->lock();
  }

  void unlock()
  {
    if (_pMutex)
      _pMutex->unlock();
  }
};

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline MwcsCutSolverHeuristic<GR, NWGHT, NLBL, EWGHT>::MwcsCutSolverHeuristic(
      IloEnv env,
      IloBoolVarArray x,
      IloBoolVarArray y,
      const Graph& g,
      const WeightNodeMap& weight,
      Node root,
      const IntNodeMap& nodeMap,
      int n,
      int m,
      IloFastMutex* pMutex)
  : IloCplex::HeuristicCallbackI(env)
  , _x(x)
  , _y(y)
  , _g(g)
  , _weight(weight)
  , _root(root)
  , _nodeMap(nodeMap)
  , _n(n)
  , _m(m)
  , _pEdgeCost(NULL)
  , _pEdgeFilterMap(NULL)
  , _pSubG(NULL)
  , _pMwcsSubGraph(NULL)
  , _pMwcsSubTreeSolver(NULL)
  , _pMutex(pMutex)
{
  lock();
  _pEdgeCost = new DoubleEdgeMap(_g);
  _pEdgeFilterMap = new BoolEdgeMap(_g, false);
  _pSubG = new SubGraphType(_g, *_pEdgeFilterMap);
  _pMwcsSubGraph = new MwcsSubGraphType();
  _pMwcsSubGraph->init(_pSubG, NULL, &_weight, NULL);
  _pMwcsSubTreeSolver = new MwcsSubTreeSolver(*_pMwcsSubGraph);
  unlock();
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline MwcsCutSolverHeuristic<GR, NWGHT, NLBL, EWGHT>::MwcsCutSolverHeuristic(
    const MwcsCutSolverHeuristic<GR, NWGHT, NLBL, EWGHT>& other)
  : IloCplex::HeuristicCallbackI(other._env)
  , _x(other._x)
  , _y(other._y)
  , _g(other._g)
  , _weight(other._weight)
  , _root(other._root)
  , _nodeMap(other._nodeMap)
  , _n(other._n)
  , _m(other._m)
  , _pEdgeCost(NULL)
  , _pEdgeFilterMap(NULL)
  , _pSubG(NULL)
  , _pMwcsSubGraph(NULL)
  , _pMwcsSubTreeSolver(NULL)
  , _pMutex(other._pMutex)
{
  lock();
  _pEdgeCost = new DoubleEdgeMap(_g);
  _pEdgeFilterMap = new BoolEdgeMap(_g, false);
  _pSubG = new SubGraphType(_g, *_pEdgeFilterMap);
  _pMwcsSubGraph = new MwcsSubGraphType();
  _pMwcsSubGraph->init(_pSubG, NULL, &_weight, NULL);
  _pMwcsSubTreeSolver = new MwcsSubTreeSolver(*_pMwcsSubGraph);
  unlock();
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline MwcsCutSolverHeuristic<GR, NWGHT, NLBL, EWGHT>::~MwcsCutSolverHeuristic()
{
  lock();
  delete _pMwcsSubTreeSolver;
  delete _pMwcsSubGraph;
  delete _pEdgeCost;
  delete _pEdgeFilterMap;
  delete _pSubG;
  unlock();
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsCutSolverHeuristic<GR, NWGHT, NLBL, EWGHT>::main()
{
  // 1. compute edge weights
  IloNumArray x_values(getEnv(), _n);
  getValues(x_values, _x);

  for (EdgeIt e(_g); e != lemon::INVALID; ++e)
  {
    Node u = _g.u(e);
    Node v = _g.v(e);

    double x_u = x_values[_nodeMap[u]];
    double x_v = x_values[_nodeMap[v]];

    _pEdgeCost->set(e, 2 - (x_u + x_v));
  }

  x_values.end();

  // 2. compute minimum spanning tree
  lock();
  lemon::kruskal(_g, *_pEdgeCost, *_pEdgeFilterMap);
  unlock();

  // 3. solve dp on the resulting spanning tree
  _pMwcsSubGraph->resetCounts();
  //_pMwcsSubGraph->init(_pSubG, NULL, &_weight, NULL);

  IloBoolVarArray solutionVar(_env, 0);
  solutionVar.add(_x);
  if (_root == lemon::INVALID)
    solutionVar.add(_y);
  IloNumArray solution(_env, solutionVar.getSize());

  double solutionWeight = hasIncumbent() ? getIncumbentObjValue() : -1;
  bool foundSolution = false;
  if (_root == lemon::INVALID)
  {
    for (NodeIt i(_g); i != lemon::INVALID; ++i)
    {
      if (_weight[i] > 0)
      {
        _pMwcsSubTreeSolver->init(i);
        _pMwcsSubTreeSolver->solve();

        if (_pMwcsSubTreeSolver->getSolutionWeight() > solutionWeight)
        {
          foundSolution = true;
          int foundRoot = _n;
          for (NodeIt n(_g); n != lemon::INVALID; ++n)
          {
            bool inSolution = _pMwcsSubTreeSolver->isNodeInSolution(n);
            solution[_nodeMap[n]] = (inSolution ? 1 : 0);

            int id_n = _nodeMap[n];
            solution[_n + id_n] = 0;
            if (inSolution && _weight[n] >= 0 && id_n < foundRoot)
            {
              foundRoot = id_n;
            }
          }
          assert(foundRoot != -1);
          solution[_n + foundRoot] = 1;

          solutionWeight = _pMwcsSubTreeSolver->getSolutionWeight();
        }
      }
    }
  }
  else
  {
    _pMwcsSubTreeSolver->init(_root);
    _pMwcsSubTreeSolver->solve();
    if (_pMwcsSubTreeSolver->getSolutionWeight() > solutionWeight)
    {
      solutionWeight = _pMwcsSubTreeSolver->getSolutionWeight();
      foundSolution = true;
      for (NodeIt n(_g); n != lemon::INVALID; ++n)
      {
        solution[_nodeMap[n]] = _pMwcsSubTreeSolver->isNodeInSolution(n);
      }
    }
  }

  // 4. pass the solution on to cplex
  if (foundSolution)
  {
    setSolution(solutionVar, solution, solutionWeight);
  }

  solutionVar.end();
  solution.end();
}

} // namespace mwcs
} // namespace nina

#endif // MWCSCUTSOLVERHEURISTIC_H
