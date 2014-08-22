/*
 * treeheurisiticsolver.h
 *
 *  Created on: 11-jan-2013
 *      Author: M. El-Kebir
 */

#ifndef TREEHEURISTICSOLVER_H
#define TREEHEURISTICSOLVER_H

#include <set>
#include <vector>
#include <map>
#include <limits>
#include <assert.h>
#include <lemon/bfs.h>
#include <lemon/connectivity.h>
#include <lemon/adaptors.h>
#include <lemon/kruskal.h>
#include <lemon/random.h>
#include "mwcsgraph.h"
#include "solver.h"
#include "treesolverrooted.h"
#include "analysis.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class TreeHeuristicSolver : public virtual Solver<GR, NWGHT, NLBL, EWGHT>
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;

  typedef Solver<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> Parent;
  typedef typename Parent::MwcsGraphType MwcsGraphType;
  typedef MwcsAnalyze<Graph> MwcsAnalyzeType;
  
  typedef typename Parent::NodeSet NodeSet;
  typedef typename Parent::NodeSetIt NodeSetIt;
  typedef typename Parent::NodeSetNonConstIt NodeSetNonConstIt;
  typedef typename Parent::NodeVector NodeVector;
  typedef typename Parent::NodeVectorIt NodeVectorIt;
  typedef typename Parent::NodeVectorNonConstIt NodeVectorNonConstIt;
  
  TEMPLATE_GRAPH_TYPEDEFS(Graph);

  typedef enum {
    EDGE_COST_FIXED,
    EDGE_COST_UNIFORM_RANDOM,
    EDGE_COST_RANDOM,
    EDGE_COST_MIN_MAX,
  } EdgeHeuristic;

protected:
  typedef lemon::FilterEdges<const Graph, const BoolEdgeMap> SubGraphType;
  typedef MwcsGraph<const SubGraphType, const WeightNodeMap, LabelNodeMap, DoubleEdgeMap> MwcsSubGraphType;
  typedef TreeSolverRooted<const SubGraphType, const WeightNodeMap, LabelNodeMap, DoubleEdgeMap> MwcsSubTreeSolverType;

  BoolEdgeMap _filterMap;
  SubGraphType _subG;
  MwcsSubGraphType* _pMwcsSubGraph;
  MwcsSubTreeSolverType* _pMwcsSubSolver;
  DoubleEdgeMap* _pEdgeCost;
  
private:
  bool _localEdgeCost;

public:
  TreeHeuristicSolver(const MwcsGraphType& mwcsGraph);
  TreeHeuristicSolver(const MwcsGraphType& mwcsGraph,
                      EdgeHeuristic heuristic,
                      MwcsAnalyzeType* pAnalyze);
  virtual ~TreeHeuristicSolver()
  {
    delete _pMwcsSubSolver;
    delete _pMwcsSubGraph;
    if (_localEdgeCost)
      delete _pEdgeCost;
  }

  void computeEdgeWeights(const MwcsGraphType& mwcsGraph,
                          EdgeHeuristic heuristic,
                          MwcsAnalyzeType* pAnalyze);

  void setEdgeCosts(DoubleEdgeMap& cost)
  {
    delete _pEdgeCost;
    _localEdgeCost = false;
    _pEdgeCost = &cost;
  }
};

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline TreeHeuristicSolver<GR, NWGHT, NLBL, EWGHT>::TreeHeuristicSolver(const MwcsGraphType& mwcsGraph)
  : Parent(mwcsGraph)
  , _filterMap(mwcsGraph.getGraph())
  , _subG(mwcsGraph.getGraph(), _filterMap)
  , _pMwcsSubGraph(NULL)
  , _pMwcsSubSolver(NULL)
  , _pEdgeCost(NULL)
  , _localEdgeCost(true)
{
}
  
template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline TreeHeuristicSolver<GR, NWGHT, NLBL, EWGHT>::TreeHeuristicSolver(const MwcsGraphType& mwcsGraph,
                                                                        EdgeHeuristic heuristic,
                                                                        MwcsAnalyzeType* pAnalyze)
  : Parent(mwcsGraph)
  , _filterMap(mwcsGraph.getGraph())
  , _subG(mwcsGraph.getGraph(), _filterMap)
  , _pMwcsSubGraph(NULL)
  , _pMwcsSubSolver(NULL)
  , _pEdgeCost(NULL)
  , _localEdgeCost(true)
{
  computeEdgeWeights(mwcsGraph, heuristic, pAnalyze);
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void TreeHeuristicSolver<GR, NWGHT, NLBL, EWGHT>::computeEdgeWeights(const MwcsGraphType& mwcsGraph,
                                                                            EdgeHeuristic heuristic,
                                                                            MwcsAnalyzeType* pAnalyze)
{
  const Graph& g = mwcsGraph.getGraph();

  delete _pEdgeCost;
  _pEdgeCost = new DoubleEdgeMap(g);

  switch (heuristic)
  {
    case EDGE_COST_MIN_MAX:
      {
        const DoubleNodeMap& score = mwcsGraph.getScores();
        double minScore = lemon::mapMinValue(g, score);
        double maxScore = lemon::mapMaxValue(g, score);

        DoubleNodeMap nodeProb(g);
        for (NodeIt n(g); n != lemon::INVALID; ++n)
        {
          double score_n = score[n];
          // let's do min-max normalization
          double norm_score_n = (score_n - minScore) / (maxScore - minScore);
          nodeProb[n] = norm_score_n;
        }

        for (EdgeIt e(g); e != lemon::INVALID; ++e)
        {
          Node u = g.u(e);
          Node v = g.v(e);

          double w_u = nodeProb[u];
          double w_v = nodeProb[v];

          _pEdgeCost->set(e, 2 - (w_u + w_v));
        }
      }
      break;
    case EDGE_COST_FIXED:
      for (EdgeIt e(g); e != lemon::INVALID; ++e)
      {
        Node u = g.u(e);
        Node v = g.v(e);

        double w_u = mwcsGraph.getScore(u);
        double w_v = mwcsGraph.getScore(v);

        _pEdgeCost->set(e, - (w_u + w_v));
      }
      break;
    case EDGE_COST_RANDOM:
      {
        DoubleNodeMap nodeProb(g);
        for (NodeIt n(g); n != lemon::INVALID; ++n)
        {
          nodeProb[n] = lemon::rnd.real();
        }

        for (EdgeIt e(g); e != lemon::INVALID; ++e)
        {
          Node u = g.u(e);
          Node v = g.v(e);

          double w_u = nodeProb[u];
          double w_v = nodeProb[v];

          _pEdgeCost->set(e, 2 - (w_u + w_v));
        }
      }
      break;
    case EDGE_COST_UNIFORM_RANDOM:
      {
        const DoubleNodeMap& score = mwcsGraph.getScores();
        double minScore = lemon::mapMinValue(g, score);
        double maxScore = lemon::mapMaxValue(g, score);

        DoubleNodeMap nodeProb(g);
        for (NodeIt n(g); n != lemon::INVALID; ++n)
        {
          double score_n = score[n];
          // let's do min-max normalization
          double norm_score_n = (score_n - minScore) / (maxScore - minScore);
          nodeProb[n] = lemon::rnd(norm_score_n);
        }

        for (EdgeIt e(g); e != lemon::INVALID; ++e)
        {
          Node u = g.u(e);
          Node v = g.v(e);

          double w_u = nodeProb[u];
          double w_v = nodeProb[v];

          _pEdgeCost->set(e, 2 - (w_u + w_v));
        }
      }
      break;
  }

  if (pAnalyze)
  {
    NodeVector rouletteWheel = pAnalyze->getRouletteWheel();
    int nDraws = pAnalyze->getNumberOfBeneficialNegHubs() / 10;
    if (nDraws == 0 && pAnalyze->getNumberOfBeneficialNegHubs() > 0)
      nDraws = 1;

    for (int i = 0; i < nDraws; i++)
    {
      int spot = lemon::rnd.integer(static_cast<int>(rouletteWheel.size()));
      Node negHub = rouletteWheel[spot];

      // remove everything from rouletteWheel corresponding to negHub;
      NodeVectorNonConstIt beginIt = std::find(rouletteWheel.begin(), rouletteWheel.end(), negHub);
      NodeVectorNonConstIt endIt = std::find(rouletteWheel.rbegin(), rouletteWheel.rend(), negHub).base();
      rouletteWheel.erase(beginIt, endIt);

      for (IncEdgeIt e(g, negHub); e != lemon::INVALID; ++e)
      {
        Node posNeighbor = g.oppositeNode(negHub, e);
        if (mwcsGraph.getScore(posNeighbor) >= 0)
        {
          double c = lemon::rnd(1e-3);
          _pEdgeCost->set(e, c);
        }
      }
    }

    const NodeSet& negHubs = pAnalyze->getBeneficialNegHubs();
    for (NodeSetIt negHubIt = negHubs.begin(); negHubIt != negHubs.end(); negHubIt++)
    {
      if (lemon::rnd(1) < 0.1)
      {
        Node negHub = *negHubIt;
        for (IncEdgeIt e(g, negHub); e != lemon::INVALID; ++e)
        {
          Node posNeighbor = g.oppositeNode(negHub, e);
          if (mwcsGraph.getScore(posNeighbor) >= 0)
          {
            double c = lemon::rnd(1e-3);
            _pEdgeCost->set(e, c);
          }
        }
      }
    }
  }
}

} // namespace mwcs
} // namespace nina

#endif // TREEHEURISTICSOLVERROOTED_H
