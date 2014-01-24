/*
 * mwcsenumeratecomp.h
 *
 *  Created on: 07-may-2013
 *      Author: M. El-Kebir
 */

#ifndef MWCSENUMERATECOMP_H
#define MWCSENUMERATECOMP_H

#include <set>
#include <assert.h>
#include "solver/mwcssolver.h"
#include "mwcsgraph.h"
#include "mwcsenumerate.h"
#include "mwcsanalyze.h"
#include "mwcspreprocessedgraph.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename WGHT = typename GR::template NodeMap<double> >
class MwcsEnumerateComp : public MwcsEnumerate<GR, WGHT>
{
public:
  typedef GR Graph;
  typedef WGHT WeightNodeMap;
  typedef MwcsEnumerate<Graph, WeightNodeMap> Parent;
  typedef typename Parent::MwcsGraphType MwcsGraphType;
  typedef typename Parent::MwcsPreprocessedGraphType MwcsPreProcessedGraphType;
  typedef typename Parent::MwcsSolverType MwcsSolverType;
  typedef typename Parent::LabelNodeMap LabelNodeMap;
  typedef typename Parent::WeightEdgeMap WeightEdgeMap;

  TEMPLATE_GRAPH_TYPEDEFS(Graph);

  typedef typename Parent::Module Module;
  typedef typename Parent::ModuleIt ModuleIt;
  typedef typename Parent::ModuleVector ModuleVector;
  typedef typename Parent::ModuleVectorIt ModuleVectorIt;
  typedef typename ModuleVector::iterator ModuleVectorNonConstIt;

  typedef typename Parent::NodeSet NodeSet;
  typedef typename Parent::NodeSetIt NodeSetIt;
  typedef typename Parent::NodeSetVector NodeSetVector;
  typedef typename Parent::NodeSetVectorIt NodeSetVectorIt;

  typedef typename Parent::SubGraph SubGraph;
  typedef typename Parent::SubNodeIt SubNodeIt;
  typedef typename Parent::NodeMap NodeMap;

  typedef MwcsAnalyze<Graph, WeightNodeMap> MwcsAnalyzeType;

  using Parent::createSolver;
  using Parent::createMwcsGraph;
  using Parent::mapModule;
  using Parent::processModule;
  using Parent::solveMWCS;
  using Parent::_mwcsGraph;
  using Parent::_modules;
  using Parent::_moduleIdx;
  using Parent::_moduleWeight;
  using Parent::_moduleSize;

public:
  MwcsEnumerateComp(MwcsGraphType& mwcsGraph)
    : Parent(mwcsGraph)
  {
  }

  virtual ~MwcsEnumerateComp() {}
  virtual void enumerate(MwcsSolverEnum solver, bool preprocess);
};

template<typename GR, typename WGHT>
inline void MwcsEnumerateComp<GR, WGHT>::enumerate(MwcsSolverEnum solver, bool preprocess)
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

#endif // MWCSENUMERATECOMP_H
