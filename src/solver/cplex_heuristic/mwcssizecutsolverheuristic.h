/*
 * mwcssizecutsolverheuristic.h
 *
 *  Created on: 3-jul-2013
 *      Author: M. El-Kebir
 */

#ifndef MWCSSIZECUTSOLVERHEURISTIC_H
#define MWCSSIZECUTSOLVERHEURISTIC_H

#include <lemon/adaptors.h>
#include <lemon/bfs.h>
#include "mwcssizetreememsolver.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class MwcsSizeCutSolverHeuristic : public IloCplex::HeuristicCallbackI
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;
  TEMPLATE_GRAPH_TYPEDEFS(Graph);

  MwcsSizeCutSolverHeuristic(IloEnv env,
                             IloBoolVarArray x,
                             IloBoolVarArray y,
                             const Graph& g,
                             const WeightNodeMap& weight,
                             Node root,
                             const IntNodeMap& nodeMap,
                             const IntArcMap& arcMap,
                             int n,
                             int m,
                             int k);
  ~MwcsSizeCutSolverHeuristic()
  {
  }

protected:
  virtual void main();
  virtual IloCplex::CallbackI* duplicateCallback() const
  {
    return (new (getEnv()) MwcsSizeCutSolverHeuristic(*this));
  }

private:
  IloBoolVarArray _x;
  IloBoolVarArray _y;
  const Graph& _g;
  const WeightNodeMap& _weight;
  const Node _root;
  const IntNodeMap& _nodeMap;
  const IntArcMap& _arcMap;
  const int _n;
  const int _m;
  const int _k;
};

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline MwcsSizeCutSolverHeuristic<GR, NWGHT, NLBL, EWGHT>::MwcsSizeCutSolverHeuristic(
      IloEnv env,
      IloBoolVarArray x,
      IloBoolVarArray y,
      const Graph& g,
      const WeightNodeMap& weight,
      Node root,
      const IntNodeMap& nodeMap,
      const IntArcMap& arcMap,
      int n,
      int m,
      int k)
  : IloCplex::HeuristicCallbackI(env)
  , _x(x)
  , _y(y)
  , _g(g)
  , _weight(weight)
  , _root(root)
  , _nodeMap(nodeMap)
  , _arcMap(arcMap)
  , _n(n)
  , _m(m)
  , _k(k)
{
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsSizeCutSolverHeuristic<GR, NWGHT, NLBL, EWGHT>::main()
{
  // 1. compute edge weights and determine best node
  DoubleEdgeMap edgeCost(_g);

  Node root = _root;
  float root_val = root != lemon::INVALID ? 1 : 0;

  for (EdgeIt e(_g); e != lemon::INVALID; ++e)
  {
    Node u = _g.u(e);
    Node v = _g.v(e);

    double x_u = getValue(_x[_nodeMap[u]]);
    double x_v = getValue(_x[_nodeMap[v]]);
    double y_u = getValue(_y[_nodeMap[u]]);
    double y_v = getValue(_y[_nodeMap[v]]);

    if (y_u > root_val)
    {
      root_val = y_u;
      root = u;
    }
    if (y_v > root_val)
    {
      root_val = y_v;
      root = v;
    }

    edgeCost[e] = 2 - (x_u + x_v);
  }

  // 2. compute minimum spanning tree
  BoolEdgeMap edgeFilterMap(_g, false);
  lemon::kruskal(_g, edgeCost, edgeFilterMap);

  // 3. solve dp on the resulting spanning tree
  typedef lemon::FilterEdges<const Graph, const BoolEdgeMap> SubGraphType;
  const SubGraphType subG(_g, edgeFilterMap);

  typedef MwcsGraph<const SubGraphType, const WeightNodeMap, LabelNodeMap, DoubleEdgeMap> MwcsSubGraphType;
  MwcsSubGraphType mwcsSubGraph;
  mwcsSubGraph.init(&subG, NULL, &_weight, NULL);

  MwcsSizeTreeMemSolver<const SubGraphType, const WeightNodeMap, LabelNodeMap, DoubleEdgeMap> solver(mwcsSubGraph, _k);
  IloBoolVarArray solutionVar(_env, 0);
  solutionVar.add(_x);
  if (_root == lemon::INVALID)
    solutionVar.add(_y);
  IloNumArray solution(_env, solutionVar.getSize());

  double solutionWeight = hasIncumbent() ? getIncumbentObjValue() : -1;
  bool foundSolution = false;
  if (_root == lemon::INVALID)
  {
    solver.init(root);
    bool res = solver.solve();
    if (res && solver.getSolutionWeight() > solutionWeight)
    {
      foundSolution = true;
      bool foundRoot = false;
      for (NodeIt n(_g); n != lemon::INVALID; ++n)
      {
        bool inSolution = solver.isNodeInSolution(n);
        solution[_nodeMap[n]] = (inSolution ? 1 : 0);
        if (inSolution && !foundRoot)
        {
          solution[_n + _nodeMap[n]] = 1;
          foundRoot = true;
        }
        else
        {
          solution[_n + _nodeMap[n]] = 0;
        }
      }

      solutionWeight = solver.getSolutionWeight();
      //std::cout << "Found solution with weight: " << solutionWeight << std::endl;
    }
  }
  else
  {
    solver.init(_root);
    solver.solve();
    if (solver.getSolutionWeight() > solutionWeight)
    {
      solutionWeight = solver.getSolutionWeight();
      foundSolution = true;
      for (NodeIt n(_g); n != lemon::INVALID; ++n)
      {
        solution[_nodeMap[n]] = solver.isNodeInSolution(n);
      }
    }
  }

  // 4. pass the solution on to cplex
  if (foundSolution)
  {
    setSolution(solutionVar, solution, solutionWeight);
  }
}

} // namespace mwcs
} // namespace nina

#endif // MWCSSIZECUTSOLVERHEURISTIC_H
