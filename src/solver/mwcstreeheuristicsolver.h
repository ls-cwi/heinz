/*
 * mwcstreeheurisiticsolver.h
 *
 *  Created on: 11-jan-2013
 *      Author: M. El-Kebir
 */

#ifndef MWCSTREEHEURISTICSOLVER_H
#define MWCSTREEHEURISTICSOLVER_H

#include "mwcssolver.h"
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
#include "mwcstreesolver.h"
#include "analysis.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class MwcsTreeHeuristicSolver : public MwcsSolver<GR, NWGHT, NLBL, EWGHT>
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;

  typedef MwcsSolver<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> Parent;
  typedef typename Parent::MwcsGraphType MwcsGraphType;
  typedef MwcsAnalyze<Graph> MwcsAnalyzeType;

  TEMPLATE_GRAPH_TYPEDEFS(Graph);

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

  typedef enum {
    EDGE_COST_FIXED,
    EDGE_COST_UNIFORM_RANDOM,
    EDGE_COST_RANDOM,
    EDGE_COST_MIN_MAX,
  } EdgeHeuristic;

private:
  typedef lemon::FilterEdges<const Graph, const BoolEdgeMap> SubGraphType;
  typedef MwcsGraph<const SubGraphType, const WeightNodeMap, LabelNodeMap, DoubleEdgeMap> MwcsSubGraphType;
  typedef MwcsTreeSolver<const SubGraphType, const WeightNodeMap, LabelNodeMap, DoubleEdgeMap> MwcsSubTreeSolverType;
  typedef typename MwcsAnalyzeType::PathType PathType;
  typedef std::set<Node> NodeSet;
  typedef typename NodeSet::const_iterator NodeSetIt;
  typedef std::vector<Node> NodeVector;
  typedef typename NodeVector::const_iterator NodeVectorIt;
  typedef typename NodeVector::iterator NodeVectorNonConstIt;

  BoolEdgeMap _filterMap;
  SubGraphType _subG;
  MwcsSubGraphType* _pMwcsSubGraph;
  MwcsSubTreeSolverType* _pMwcsSubSolver;
  DoubleEdgeMap* _pEdgeCost;
  bool _localEdgeCost;

public:
  MwcsTreeHeuristicSolver(const MwcsGraphType& mwcsGraph);
  virtual ~MwcsTreeHeuristicSolver()
  {
    delete _pMwcsSubSolver;
    delete _pMwcsSubGraph;
    if (_localEdgeCost)
      delete _pEdgeCost;
  }

  virtual void init(Node root);
  virtual bool solve();
  virtual void setLowerBound(double LB) {}

  //void setEdgeHeuristic(EdgeHeuristic heuristic)
  //{
  //  _edgeHeuristic = heuristic;
  //}
//
  //EdgeHeuristic getEdgeHeuristic() const
  //{
  //  return _edgeHeuristic;
  //}

  void computeEdgeWeights(EdgeHeuristic heuristic,
                          MwcsAnalyzeType* pAnalyze = NULL)
  {
    const Graph& g = _mwcsGraph.getGraph();

    delete _pEdgeCost;
    _pEdgeCost = new DoubleEdgeMap(g);

    switch (heuristic)
    {
      case EDGE_COST_MIN_MAX:
        {
          const DoubleNodeMap& score = _mwcsGraph.getScores();
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

          double w_u = _mwcsGraph.getScore(u);
          double w_v = _mwcsGraph.getScore(v);

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
          const DoubleNodeMap& score = _mwcsGraph.getScores();
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
    //  for (NodeIt v(g); v != lemon::INVALID; ++v)
    //  {
    //    for (NodeIt w(g); w != lemon::INVALID; ++w)
    //    {
    //      if (v == w) continue;
    //      if (pAnalyze->ok(v, w) && pAnalyze->ok(w, v))
    //      {
    //        PathType p = pAnalyze->path(v, w);
    //        for (typename PathType::ArcIt a(p); a != lemon::INVALID; ++a)
    //        {
    //          Arc g_a = a;
    //          Edge e = g_a;
    //          _pEdgeCost->set(e, -100);
    //        }
    //      }
    //    }
    //  }

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
          if (_mwcsGraph.getScore(posNeighbor) >= 0)
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
            if (_mwcsGraph.getScore(posNeighbor) >= 0)
            {
              double c = lemon::rnd(1e-3);
              _pEdgeCost->set(e, c);
            }
          }
        }
      }

      //for (NodeIt u(g); u != lemon::INVALID; ++u)
      //{
      //  double score_u = _mwcsGraph.getScore(u);
//
      //  if (score_u < 0)
      //  {
      //    NodeSet posNeighbors;
      //    for (IncEdgeIt e(g, u); e != lemon::INVALID; ++e)
      //    {
      //      Node posNeighbor = g.oppositeNode(u, e);
      //      if (_mwcsGraph.getScore(posNeighbor) >= 0)
      //      {
      //        posNeighbors.insert(posNeighbor);
      //      }
      //    }
//
      //    double maxScorePosNeighbor = 0;
      //    double scorePosNeighbors = 0;
      //    for (NodeSetIt nodeIt = posNeighbors.begin();
      //         nodeIt != posNeighbors.end(); nodeIt++)
      //    {
      //      Node v = *nodeIt;
      //      double score_v = _mwcsGraph.getScore(v);
      //      assert(score_v >= 0);
//
      //      scorePosNeighbors += score_v;
      //      if (score_v > maxScorePosNeighbor)
      //        maxScorePosNeighbor = score_v;
      //    }
//
      //    double gain = scorePosNeighbors + score_u;
//
      //    if (gain > maxScorePosNeighbor && lemon::rnd(1) < 0.1)
      //    {
      //      for (IncEdgeIt e(g, u); e != lemon::INVALID; ++e)
      //      {
      //        Node posNeighbor = g.oppositeNode(u, e);
      //        if (_mwcsGraph.getScore(posNeighbor) >= 0)
      //        {
      //          double c = lemon::rnd(1e-3);
      //          _pEdgeCost->set(e, c);
      //        }
      //      }
      //    }
      //  }
      //}


    }
  }

  void setEdgeCosts(DoubleEdgeMap& cost)
  {
    delete _pEdgeCost;
    _localEdgeCost = false;
    _pEdgeCost = &cost;
  }
};

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline MwcsTreeHeuristicSolver<GR, NWGHT, NLBL, EWGHT>::MwcsTreeHeuristicSolver(const MwcsGraphType& mwcsGraph)
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
inline void MwcsTreeHeuristicSolver<GR, NWGHT, NLBL, EWGHT>::init(Node root)
{
  _root = root;

  // construct mwcssubgraph instance
  delete _pMwcsSubGraph;
  _pMwcsSubGraph = new MwcsSubGraphType();
  _pMwcsSubGraph->init(&_subG, NULL, &_mwcsGraph.getScores(), NULL);

  // construct mwcssubsolver instance
  delete _pMwcsSubSolver;
  _pMwcsSubSolver = new MwcsSubTreeSolverType(*_pMwcsSubGraph);
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline bool MwcsTreeHeuristicSolver<GR, NWGHT, NLBL, EWGHT>::solve()
{
  const Graph& g = _mwcsGraph.getGraph();

  if (!_pEdgeCost)
    computeEdgeWeights(EDGE_COST_FIXED);

  // compute minimum spanning tree
  lemon::kruskal(g, *_pEdgeCost, _filterMap);

  bool res = false;

  if (_root == lemon::INVALID)
  {
    for (NodeIt v(g); v != lemon::INVALID; ++v)
    {
      if (_mwcsGraph.getScore(v) > 0)
      {
        _pMwcsSubSolver->init(v);
        res = _pMwcsSubSolver->solve();
        if (_score < _pMwcsSubSolver->getSolutionWeight())
        {
          _score = _pMwcsSubSolver->getSolutionWeight();
          lemon::mapCopy(g, _pMwcsSubSolver->getSolutionNodeMap(), _solutionMap);
          _solutionSet = _pMwcsSubSolver->getSolutionModule();
        }
      }
    }
  }
  else
  {
    _pMwcsSubSolver->init(_root);
    res = _pMwcsSubSolver->solve();
    _score = _pMwcsSubSolver->getSolutionWeight();
    lemon::mapCopy(g, _pMwcsSubSolver->getSolutionNodeMap(), _solutionMap);
    _solutionSet = _pMwcsSubSolver->getSolutionModule();
  }

  return res;
}

} // namespace mwcs
} // namespace nina

#endif // MWCSTREEHEURISTICSOLVER_H
