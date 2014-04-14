/*
 * mwcsenumerate.h
 *
 *  Created on: 30-jan-2013
 *      Author: M. El-Kebir
 */

#ifndef MWCSENUMERATE_H
#define MWCSENUMERATE_H

#include <set>
#include <assert.h>
#include <ostream>
#include "mwcs.h"
#include "mwcsgraph.h"
#include "mwcspreprocessedgraph.h"
#include "solver/cplex_cut/backoff.h"
#include "solver/mwcssolver.h"
#include "solver/mwcscutsolver.h"
#include "solver/mwcssizecutsolver.h"
#include "solver/mwcstreesolver.h"
#include "solver/mwcssizetreesolver.h"
//#include "solver/mwcssizetreememsolver.h"
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
         typename WGHT = typename GR::template NodeMap<double> >
class MwcsEnumerate
{
public:
  typedef GR Graph;
  typedef WGHT WeightNodeMap;
  typedef MwcsGraph<Graph, WeightNodeMap> MwcsGraphType;
  typedef MwcsPreprocessedGraph<Graph, WeightNodeMap> MwcsPreprocessedGraphType;
  typedef MwcsSolver<Graph, WeightNodeMap> MwcsSolverType;
  typedef typename MwcsGraphType::LabelNodeMap LabelNodeMap;
  typedef typename MwcsGraphType::WeightEdgeMap WeightEdgeMap;

  TEMPLATE_GRAPH_TYPEDEFS(Graph);

  typedef typename std::set<Node> Module;
  typedef typename Module::const_iterator ModuleIt;
  typedef typename std::vector<Module> ModuleVector;
  typedef typename ModuleVector::const_iterator ModuleVectorIt;

  typedef typename std::set<Node> NodeSet;
  typedef typename NodeSet::const_iterator NodeSetIt;
  typedef typename std::vector<NodeSet> NodeSetVector;
  typedef typename NodeSetVector::const_iterator NodeSetVectorIt;

  typedef lemon::FilterNodes<Graph, BoolNodeMap> SubGraph;
  typedef typename SubGraph::NodeIt SubNodeIt;

  typedef MwcsCutSolver<Graph, WeightNodeMap> MwcsCutSolverType;
  typedef MwcsTreeSolver<Graph, WeightNodeMap> MwcsTreeSolverType;
  typedef MwcsSizeCutSolver<Graph, WeightNodeMap> MwcsSizeCutSolverType;
  typedef MwcsSizeTreeSolver<Graph, WeightNodeMap> MwcsSizeTreeSolverType;
  //typedef MwcsSizeTreeMemSolver<Graph, WeightNodeMap> MwcsSizeTreeMemSolverType;

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
  MwcsEnumerate(MwcsGraphType& mwcsGraph);
  virtual ~MwcsEnumerate() {}
  virtual void enumerate(MwcsSolverEnum solver, bool preprocess);
  void printOutput(std::ostream& out) const;

  const ModuleVector& getModules() const
  {
    return _modules;
  }

  const Module& getModule(int idx) const
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

  void setTimeLimit(int timeLimit)
  {
    _timeLimit = timeLimit;
  }

  void setMultiThreading(int multiThreading)
  {
    _multiThreading = multiThreading;
  }

  void setModuleSize(int moduleSize)
  {
    _moduleSize = moduleSize;
    if (moduleSize > 0)
      lemon::mapFill(_mwcsGraph.getGraph(),
                     _moduleWeight,
                     -std::numeric_limits<double>::max());
  }
  
  void setMaxNumberOfCuts(int maxNumberOfCuts)
  {
    _maxNumberOfCuts = maxNumberOfCuts;
  }
  
  void setBackOff(BackOff backOff)
  {
    _backOff = backOff;
  }

protected:
  MwcsGraphType& _mwcsGraph;

  ModuleVector _modules;
  IntNodeMap _moduleIdx;
  WeightNodeMap _moduleWeight;
  int _timeLimit;
  int _multiThreading;
  int _moduleSize;
  int _maxNumberOfCuts;
  BackOff _backOff;

  typedef typename Graph::template NodeMap<Node> NodeMap;

protected:
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

#endif // MWCSENUMERATE_H
