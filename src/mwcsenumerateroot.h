/*
 * mwcsenumerate.h
 *
 *  Created on: 30-jan-2013
 *      Author: M. El-Kebir
 */

#ifndef MWCSENUMERATEROOT_H
#define MWCSENUMERATEROOT_H

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
class MwcsEnumerateRoot : public MwcsEnumerate<GR, WGHT>
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
  using Parent::mapModule;
  using Parent::processModule;
  using Parent::_mwcsGraph;
  using Parent::_modules;
  using Parent::_moduleIdx;
  using Parent::_moduleWeight;
  using Parent::_moduleSize;

public:
  MwcsEnumerateRoot(MwcsGraphType& mwcsGraph)
    : Parent(mwcsGraph)
  {
  }

  virtual ~MwcsEnumerateRoot() {}

protected:
  virtual bool solveMWCS(MwcsGraphType* pMwcsSubGraph,
                         const NodeMap& mapToG,
                         MwcsSolverEnum solver,
                         NodeSet& pickedNodes)
  {
    MwcsAnalyzeType mwcsAnalyze(*pMwcsSubGraph);
    mwcsAnalyze.analyze(true);

    const NodeSetVector& eqClasses = mwcsAnalyze.getEqClasses();

    MwcsSolverType* pSolver = createSolver(pMwcsSubGraph, solver);

    Module mappedModule;

    double LB;
    if (_moduleSize > 0)
      LB = -std::numeric_limits<double>::max();
    else
      LB = 0;

    for (NodeSetVectorIt eqClassIt = eqClasses.begin();
         eqClassIt != eqClasses.end(); eqClassIt++)
    {
      // ToDo: something more involved here might pay off
      Node root = *eqClassIt->begin();
      Node mappedRoot = pMwcsSubGraph->init(root);

      const Graph& g = pMwcsSubGraph->getGraph();

      // first initialize
      pSolver->setLowerBound(LB);
      pSolver->init(mappedRoot);

      bool solved = pSolver->solve();
      double moduleWeight = pSolver->getSolutionWeight();

      if (solved && moduleWeight > LB)
      {
        LB = moduleWeight;
        mappedModule =
          mapModule(pMwcsSubGraph->getOrgNodes(pSolver->getSolutionModule()), mapToG);
        if (g_verbosity >= VERBOSE_ESSENTIAL)
        {
          std::cout << "// Intermediate solution with weight " << LB
                    << " and " << mappedModule.size() << " nodes found" << std::endl;
        }
      }

      pMwcsSubGraph->deinit();
    }

    bool result;
    if (LB > 0 || (LB != -std::numeric_limits<double>::max() && _moduleSize > 0))
    {
      pickedNodes.insert(mappedModule.begin(), mappedModule.end());
      processModule(mappedModule, LB);

      if (g_verbosity >= VERBOSE_ESSENTIAL)
      {
        std::cout << "// Solution with weight " << LB
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

      result = false;
    }

    delete pSolver;
    return result;
  }
};

} // namespace mwcs
} // namespace nina

#endif // MWCSENUMERATEROOT_H
