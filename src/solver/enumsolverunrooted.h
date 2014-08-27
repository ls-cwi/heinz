/*
 * enumsolverunrooted.h
 *
 *  Created on: 30-jan-2013
 *      Author: M. El-Kebir
 */

#ifndef ENUMSOLVERUNROOTED_H
#define ENUMSOLVERUNROOTED_H

#include <set>
#include <assert.h>
#include <ostream>
#include "mwcs.h"
#include "mwcsgraph.h"
#include "mwcspreprocessedgraph.h"

#include "solver/solver.h"
#include "solver/solverrooted.h"
#include "solver/solverunrooted.h"
#include "solver/cutsolverrooted.h"
#include "solver/cutsolverunrooted.h"
#include "solver/cplex_cut/backoff.h"

#include "preprocessing/negdeg01.h"
#include "preprocessing/posedge.h"
#include "preprocessing/negedge.h"
#include "preprocessing/rootedposdeg01.h"
#include "preprocessing/negcircuit.h"
#include "preprocessing/negdiamond.h"
#include "preprocessing/negmirroredhubs.h"
#include "preprocessing/posdeg01.h"
#include "preprocessing/posdiamond.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class EnumSolverUnrooted : public SolverUnrooted<GR, NWGHT, NLBL, EWGHT>
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;
  
  typedef SolverUnrooted<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> Parent;
  typedef typename Parent::MwcsGraphType MwcsGraphType;
  typedef typename Parent::NodeSet NodeSet;
  typedef typename Parent::NodeSetIt NodeSetIt;
  
  typedef MwcsPreprocessedGraph<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> MwcsPreGraphType;
  typedef Solver<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> SolverType;
  typedef CutSolverRooted<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> CutSolverRootedType;
  typedef CutSolverUnrooted<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> CutSolverUnrootedType;
  typedef typename CutSolverRooted::Options Options;
  typedef typename CutSolverRooted::MwcsAnalyzeType MwcsAnalyzeType;

  TEMPLATE_GRAPH_TYPEDEFS(Graph);

  typedef typename std::vector<NodeSet> NodeSetVector;
  typedef typename NodeSetVector::const_iterator NodeSetVectorIt;

  typedef lemon::FilterNodes<Graph, BoolNodeMap> SubGraph;
  typedef typename SubGraph::NodeIt SubNodeIt;

  // pre processing
  typedef NegDeg01<Graph, WeightNodeMap> NegDeg01Type;
  typedef PosEdge<Graph, WeightNodeMap> PosEdgeType;
  typedef NegEdge<Graph, WeightNodeMap> NegEdgeType;
  typedef RootedPosDeg01<Graph, WeightNodeMap> RootedPosDeg01Type;
  typedef NegCircuit<Graph, WeightNodeMap> NegCircuitType;
  typedef NegDiamond<Graph, WeightNodeMap> NegDiamondType;
  typedef NegMirroredHubs<Graph, WeightNodeMap> NegMirroredHubsType;
  typedef PosDeg01<Graph, WeightNodeMap> PosDeg01Type;
  typedef PosDiamond<Graph, WeightNodeMap> PosDiamondType;

public:
  EnumSolverUnrooted(MwcsGraphType& mwcsGraph,
                     const Options& options,
                     const MwcsAnalyzeType& analysis)
    : Parent(mwcsGraph)
    , _options(options)
    , _analysis(analysis)
  {
  }
  
  virtual ~EnumSolverUnrooted()
  {
  }
  
  virtual void init();
  virtual void solve();
  
  void printOutput(std::ostream& out) const;

  const NodeSetVector& getModules() const
  {
    return _modules;
  }

  const NodeSet& getModule(int idx) const
  {
    assert(0 <= idx && idx < static_cast<int>(_modules.size()));
    return _modules[idx];
  }

  int getModuleIndex(Node node) const
  {
    assert(node != lemon::INVALID);
    return _moduleIdx[node];
  }

  double getModuleWeight(Node node) const
  {
    assert(node != lemon::INVALID);
    return _moduleWeight[node];
  }

protected:
  const Options& _options;
  const MwcsAnalyzeType& _analysis;

  ModuleVector _modules;
  IntNodeMap _moduleIdx;
  WeightNodeMap _moduleWeight;

  typedef typename Graph::template NodeMap<Node> NodeMap;

  void processModule(const Module& module, double moduleWeight)
  {
    int idx = static_cast<int>(_modules.size());

    _modules.push_back(NodeSet());

    for (ModuleIt nodeIt = module.begin(); nodeIt != module.end(); nodeIt++)
    {
      NodeSet orgNodes = _mwcsGraph.getOrgNodes(*nodeIt);
      _modules.back().insert(orgNodes.begin(), orgNodes.end());

      assert(orgNodes.size() != 0);
      assert(_moduleWeight[*orgNodes.begin()] < moduleWeight);

      for (NodeSetIt orgNodeIt = orgNodes.begin(); orgNodeIt != orgNodes.end(); orgNodeIt++)
      {
        _moduleIdx[*orgNodeIt] = idx;
        _moduleWeight[*orgNodeIt] = moduleWeight;
      }
    }
  }

  MwcsSolverType* createSolver(MwcsGraphType* pMwcsSubGraph,
                               MwcsSolverEnum solver)
  {
    MwcsSolverType* pResult = NULL;

    switch (solver)
    {
      case MwcsSolverCutNodeSeparatorBk:
        pResult = new MwcsCutSolverType(*pMwcsSubGraph,
                                        _backOff,
                                        _maxNumberOfCuts,
                                        _timeLimit,
                                        _multiThreading);
        break;
      case MwcsSolverTreeDP:
        pResult = new MwcsTreeSolverType(*pMwcsSubGraph);
        break;
      case MwcsSizeSolverTreeDP:
        pResult = new MwcsSizeTreeSolverType(*pMwcsSubGraph, _moduleSize);
        break;
      case MwcsSizeSolverCutNodeSeparatorBk:
        pResult = new MwcsSizeCutSolverType(*pMwcsSubGraph,
                                            _moduleSize,
                                            _maxNumberOfCuts,
                                            _timeLimit,
                                            _multiThreading);
        break;
    }

    return pResult;
  }

  MwcsGraphType* createMwcsGraph(bool preprocess)
  {
    if (preprocess)
    {
      MwcsPreprocessedGraphType* pPreprocessedMwcs = new MwcsPreprocessedGraphType();
      pPreprocessedMwcs->addPreprocessRule(1, new NegDeg01Type());
      pPreprocessedMwcs->addPreprocessRule(1, new PosEdgeType());
      pPreprocessedMwcs->addPreprocessRule(1, new NegEdgeType());
      pPreprocessedMwcs->addPreprocessRootRule(1, new RootedPosDeg01Type());
      pPreprocessedMwcs->addPreprocessRule(1, new NegCircuitType());
      pPreprocessedMwcs->addPreprocessRule(1, new NegDiamondType());
      pPreprocessedMwcs->addPreprocessRule(1, new PosDeg01Type());
      pPreprocessedMwcs->addPreprocessRule(1, new PosDiamondType());
      pPreprocessedMwcs->addPreprocessRule(2, new NegMirroredHubsType());

      return pPreprocessedMwcs;
    }
    else
    {
      return new MwcsGraphType();
    }
  }

  Module mapModule(const Module& module, const NodeMap& mapToG) const
  {
    Module result;
    for (ModuleIt nodeIt = module.begin(); nodeIt != module.end(); nodeIt++)
    {
      result.insert(mapToG[*nodeIt]);
    }
    return result;
  }

  virtual bool solveMWCS(MwcsGraphType* pMwcsSubGraph,
                         const NodeMap& mapToG,
                         MwcsSolverEnum solver,
                         NodeSet& pickedNodes)
  {
    bool result;

    MwcsSolverType* pSolver = createSolver(pMwcsSubGraph, solver);
    pSolver->init();

    if (pSolver->solve() && (pSolver->getSolutionWeight() > 0 || _moduleSize > 0))
    {
      Module mappedModule =
          mapModule(pMwcsSubGraph->getOrgNodes(pSolver->getSolutionModule()), mapToG);

      pickedNodes.insert(mappedModule.begin(), mappedModule.end());

      processModule(mappedModule, pSolver->getSolutionWeight());

      if (g_verbosity >= VERBOSE_ESSENTIAL)
      {
        std::cout << "// Solution with weight " << pSolver->getSolutionWeight()
                  << " and " << mappedModule.size() << " nodes found" << std::endl;
      }

      result = true;
    }
    else
    {
      if (g_verbosity >= VERBOSE_ESSENTIAL)
      {
        std::cout << "// No feasible solution found" << std::endl;
      }

      // add the empty module
      Module mappedModule;
      processModule(mappedModule, 0);

      result = false;
    }

    delete pSolver;
    return result;
  }
};

template<typename GR, typename WGHT>
inline MwcsEnumerate<GR, WGHT>::MwcsEnumerate(MwcsGraphType& mwcsGraph)
  : _mwcsGraph(mwcsGraph)
  , _modules()
  , _moduleIdx(_mwcsGraph.getOrgGraph(), -1)
  , _moduleWeight(_mwcsGraph.getOrgGraph(), 0)
  , _timeLimit(-1)
  , _multiThreading(1)
  , _moduleSize(-1)
  , _backOff(1)
{
}

template<typename GR, typename WGHT>
inline void MwcsEnumerate<GR, WGHT>::printOutput(std::ostream& out) const
{
  for (NodeIt node(_mwcsGraph.getOrgGraph()); node != lemon::INVALID; ++node)
  {
    out << _mwcsGraph.getOrgLabel(node) << "\t"
        << _moduleIdx[node] << "\t"
        << _moduleWeight[node] << std::endl;
  }
}

template<typename GR, typename WGHT>
inline void MwcsEnumerate<GR, WGHT>::enumerate(MwcsSolverEnum solver, bool preprocess)
{
  Graph subG;
  DoubleNodeMap weightSubG(subG);
  LabelNodeMap labelSubG(subG);
  NodeMap mapToG(subG);
  
  MwcsGraphType* pMwcsSubGraph = createMwcsGraph(preprocess);
  
  // contains the set of picked nodes (not necessarily in the original graph)
  NodeSet pickedNodes;
  
  // 1. construct the subgraph
  Graph& g = _mwcsGraph.getGraph();
  
  // 1a. determine allowed nodes
  BoolNodeMap allowedNodes(g, true);
  for (NodeSetIt nodeIt = pickedNodes.begin(); nodeIt != pickedNodes.end(); nodeIt++)
  {
    allowedNodes[*nodeIt] = false;
  }
  
  // 1b. construct subgraph
  SubGraph subTmpG(g, allowedNodes);
  
  // 1c. determine connected components in subG
  IntNodeMap comp(g, -1);
  int nComponents = lemon::connectedComponents(subTmpG, comp);
  
  for (int compIdx = 0; compIdx < nComponents; compIdx++)
  {
    // 2a. determine nodes in same component
    int nNodesComp = 0;
    BoolNodeMap allowedNodesSameComp(g, false);
    for (SubNodeIt node(subTmpG); node != lemon::INVALID; ++node)
    {
      if (comp[node] == compIdx)
      {
        allowedNodesSameComp[node] = true;
        nNodesComp++;
      }
    }
    
    // 2b. create and preprocess subgraph
    SubGraph subTmpSameCompG(g, allowedNodesSameComp);
    lemon::graphCopy(subTmpSameCompG, subG)
    .nodeMap(_mwcsGraph.getScores(), weightSubG)
    .nodeMap(_mwcsGraph.getLabels(), labelSubG)
    .nodeCrossRef(mapToG)
    .run();
    
    if (g_verbosity >= VERBOSE_ESSENTIAL)
    {
      std::cout << std::endl;
      std::cout << "// Considering component " << compIdx + 1 << "/" << nComponents
      << ": contains " << nNodesComp << " nodes" << std::endl;
    }
    pMwcsSubGraph->init(&subG, &labelSubG, &weightSubG, NULL);
    
    // 3. solve
    solveMWCS(pMwcsSubGraph, mapToG, solver, pickedNodes);
  }
  
  delete pMwcsSubGraph;
}

} // namespace mwcs
} // namespace nina

#endif // ENUMSOLVERUNROOTED_H
