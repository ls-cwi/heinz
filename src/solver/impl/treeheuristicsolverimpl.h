/*
 * treeheurisiticsolverimpl.h
 *
 *  Created on: 11-jan-2013
 *      Author: M. El-Kebir
 */

#ifndef TREEHEURISTICSOLVERIMPL_H
#define TREEHEURISTICSOLVERIMPL_H

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
#include "../solver.h"
#include "treesolverimpl.h"
#include "analysis.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class TreeHeuristicSolverImpl
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;
  
  typedef MwcsGraph<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> MwcsGraphType;
  typedef MwcsAnalyze<Graph> MwcsAnalyzeType;
  
  TEMPLATE_GRAPH_TYPEDEFS(Graph);

  typedef enum {
    EDGE_COST_FIXED,
    EDGE_COST_UNIFORM_RANDOM,
    EDGE_COST_RANDOM,
    EDGE_COST_MIN_MAX,
  } EdgeHeuristic;
  
  struct Options
  {
    Options(EdgeHeuristic heuristic,
            bool analysis,
            int nRepetitions)
      : _heuristic(heuristic)
      , _analysis(analysis)
      , _nRepetitions(nRepetitions)
    {
    }
    
    EdgeHeuristic _heuristic;
    bool _analysis;
    int _nRepetitions;
  };

protected:
  typedef lemon::FilterEdges<const Graph, const BoolEdgeMap> SubGraphType;
  typedef MwcsGraph<const SubGraphType, const WeightNodeMap, LabelNodeMap, DoubleEdgeMap> MwcsSubGraphType;
  typedef TreeSolverImpl<const SubGraphType, const WeightNodeMap, LabelNodeMap, DoubleEdgeMap> MwcsSubTreeSolverType;
  typedef typename MwcsSubTreeSolverType::BoolNodeMap SubBoolNodeMap;
  
  typedef std::set<Node> NodeSet;
  typedef typename NodeSet::const_iterator NodeSetIt;
  
  typedef std::vector<Node> NodeVector;
  typedef typename NodeVector::const_iterator NodeVectorIt;
  typedef typename NodeVector::iterator NodeVectorNonConstIt;

  Options _options;
  BoolEdgeMap* _pFilterMap;
  SubGraphType* _pSubG;
  MwcsSubGraphType* _pMwcsSubGraph;
  DoubleEdgeMap* _pEdgeCost;
  
private:
  MwcsAnalyzeType* _pAnalysis;

protected:
  TreeHeuristicSolverImpl(Options options)
    : _options(options)
    , _pFilterMap(NULL)
    , _pSubG(NULL)
    , _pMwcsSubGraph(NULL)
    , _pEdgeCost(NULL)
    , _pAnalysis(NULL)
  {
  }

  virtual ~TreeHeuristicSolverImpl()
  {
    delete _pFilterMap;
    delete _pMwcsSubGraph;
    
    delete _pEdgeCost;
    
    delete _pSubG;
    
    delete _pAnalysis;
  }

  void init(const MwcsGraphType& mwcsGraph);
  
  bool solveMonteCarlo(const MwcsGraphType& mwcsGraph,
                       MwcsSubTreeSolverType& subTreeSolver,
                       double& score,
                       BoolNodeMap& solutionMap,
                       NodeSet& solutionSet);
  
  void computeEdgeCosts(const MwcsGraphType& mwcsGraph);
  
  virtual void initSubTreeSolver() = 0;
  
  virtual bool solveSubTree(double& score, SubBoolNodeMap& solutionMap, NodeSet& solutionSet) = 0;
  
private:
  void analyzeNegHubs(const MwcsGraphType& mwcsGraph);
};

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void TreeHeuristicSolverImpl<GR, NWGHT, NLBL, EWGHT>::init(const MwcsGraphType& mwcsGraph)
{
  const Graph& g = mwcsGraph.getGraph();
  
  delete _pFilterMap;
  _pFilterMap = new BoolEdgeMap(g);
  
  delete _pSubG;
  _pSubG = new SubGraphType(g, *_pFilterMap);
  
  if (_options._analysis)
  {
    _pAnalysis = new MwcsAnalyzeType(mwcsGraph);
    _pAnalysis->analyzeNegHubs();
    
    if (g_verbosity >= VERBOSE_ESSENTIAL)
    {
      std::cout << "// Number of beneficial negative hubs: "
                << _pAnalysis->getNumberOfBeneficialNegHubs()
                << std::endl;
    }
  }
  
  delete _pEdgeCost;
  _pEdgeCost = new DoubleEdgeMap(g);
  
  delete _pMwcsSubGraph;
  _pMwcsSubGraph = new MwcsSubGraphType();
  _pMwcsSubGraph->init(_pSubG, NULL, &mwcsGraph.getScores(), NULL);
}
  
template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void TreeHeuristicSolverImpl<GR, NWGHT, NLBL, EWGHT>::computeEdgeCosts(const MwcsGraphType& mwcsGraph)
{
  const Graph& g = mwcsGraph.getGraph();
  const DoubleNodeMap& score = mwcsGraph.getScores();
  
  switch (_options._heuristic)
  {
    case EDGE_COST_MIN_MAX:
      {
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
  
  if (_options._analysis)
  {
    assert(_pAnalysis);
    analyzeNegHubs(mwcsGraph);
  }
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void TreeHeuristicSolverImpl<GR, NWGHT, NLBL, EWGHT>::analyzeNegHubs(const MwcsGraphType& mwcsGraph)
{
  const Graph& g = mwcsGraph.getGraph();
  
  NodeVector rouletteWheel = _pAnalysis->getRouletteWheel();
  int nDraws = _pAnalysis->getNumberOfBeneficialNegHubs() / 10;
  if (nDraws == 0 && _pAnalysis->getNumberOfBeneficialNegHubs() > 0)
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
  
  const NodeSet& negHubs = _pAnalysis->getBeneficialNegHubs();
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
  
template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline bool TreeHeuristicSolverImpl<GR, NWGHT, NLBL, EWGHT>::solveMonteCarlo(const MwcsGraphType& mwcsGraph,
                                                                             MwcsSubTreeSolverType& subTreeSolver,
                                                                             double& score,
                                                                             BoolNodeMap& solutionMap,
                                                                             NodeSet& solutionSet)
{
  typedef typename MwcsSubTreeSolverType::BoolNodeMap SubBoolNodeMap;
  
  const Graph& g = mwcsGraph.getGraph();
  
  double newScore;
  double newScoreUB;
  SubBoolNodeMap newSolutionMap(_pMwcsSubGraph->getGraph());
  NodeSet newSolutionSet;

  for (int i = 0; i < _options._nRepetitions; ++i)
  {
    if (g_verbosity >= VERBOSE_NON_ESSENTIAL)
    {
      std::cerr << "\rIteration " << i << ": " << std::flush;
    }
    
    computeEdgeCosts(mwcsGraph);

    // compute minimum spanning tree
    lemon::kruskal(g, *_pEdgeCost, *_pFilterMap);

    initSubTreeSolver();
    if (!subTreeSolver.solve(newScore, newScoreUB, newSolutionMap, newSolutionSet))
    {
      return false;
    }
    
    if (score < newScore)
    {
      score = newScore;
      
      lemon::mapCopy(_pMwcsSubGraph->getGraph(), newSolutionMap, solutionMap);
      solutionSet = newSolutionSet;
    }
    
    if (g_verbosity >= VERBOSE_NON_ESSENTIAL)
    {
      std::cerr << score << std::flush;
    }
  }
  
  if (g_verbosity >= VERBOSE_NON_ESSENTIAL)
  {
    std::cerr << std::endl;
  }
  
  return true;
}

  
} // namespace mwcs
} // namespace nina

#endif // TREEHEURISTICSOLVERIMPL_H
