/*
 * mwcsgraphparser.h
 *
 *  Created on: 30-jan-2013
 *      Author: M. El-Kebir
 */

#ifndef MWCSGRAPHPARSER_H
#define MWCSGRAPHPARSER_H

#include <map>
#include <assert.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <math.h>
#include <lemon/core.h>
#include <lemon/lgf_writer.h>
#include <lemon/connectivity.h>
#include "mwcsgraph.h"
#include "parser/parser.h"
#include "verbose.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class MwcsGraphParser : public MwcsGraph<GR, NWGHT, NLBL, EWGHT>
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;
  typedef MwcsGraph<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> Parent;

  TEMPLATE_GRAPH_TYPEDEFS(Graph);

public:
  typedef typename Parent::ParserType ParserType;
  typedef typename Parent::InvLabelNodeMap InvLabelNodeMap;
  typedef typename Parent::InvLabelNodeMapIt InvLabelNodeMapIt;

  typedef typename Parent::NodeSet NodeSet;
  typedef typename Parent::NodeSetIt NodeSetIt;
  typedef typename Parent::NodeSetVector NodeSetVector;
  typedef typename Parent::NodeSetVectorIt NodeSetVectorIt;

  using Parent::getGraph;
  using Parent::getScores;
  using Parent::getOrgArcCount;
  using Parent::getOrgComponent;
  using Parent::getOrgComponentCount;
  using Parent::getOrgComponentMap;
  using Parent::getOrgEdgeCount;
  using Parent::getOrgGraph;
  using Parent::getOrgLabel;
  using Parent::getOrgLabels;
  using Parent::getOrgNodeByLabel;
  using Parent::getOrgNodeCount;
  using Parent::getOrgScore;
  using Parent::getOrgScores;
  using Parent::getOrgPValue;

public:
  MwcsGraphParser()
    : Parent()
  {
  }

  virtual ~MwcsGraphParser()
  {
  }

  virtual void computeScores(double tau)
  {
    const double log_tau = log(tau);

    if (g_verbosity >= VERBOSE_ESSENTIAL)
    {
      std::cerr << "// Offset: " << tau << std::endl;
    }

    const Graph& g = getOrgGraph();
    WeightNodeMap& score = getOrgScores();

    int nPos = 0;
    int nNeg = 0;
    for (NodeIt v(g); v != lemon::INVALID; ++v)
    {
      double x = getOrgPValue(v);
      double score_x = -log(x) + log_tau;

      if (score_x > 0)
        nPos++;
      else
        nNeg++;

      score.set(v, score_x);
    }

    if (g_verbosity >= VERBOSE_ESSENTIAL)
    {
      std::cerr << "// Number of positive nodes: " << nPos << std::endl;
      std::cerr << "// Number of negative nodes: " << nNeg << std::endl;
      std::cerr << "// Fraction of positive nodes: "
                << (double)nPos / (double)getOrgNodeCount() << std::endl;
    }
  }

  virtual void computeScores(double lambda, double a, double FDR)
  {
    double tau = pow(
          ((lambda + (1 - lambda) * a) - FDR * lambda) / (FDR * (1 - lambda)),
          1 / (a - 1));
    double a_log_tau = (a-1) * log(tau);

    if (g_verbosity >= VERBOSE_ESSENTIAL)
    {
      std::cerr << "// FDR offset: " << -a_log_tau << std::endl;
    }

    const Graph& g = getOrgGraph();
    WeightNodeMap& score = getOrgScores();

    int nPos = 0;
    int nNeg = 0;
    for (NodeIt v(g); v != lemon::INVALID; ++v)
    {
      double x = getOrgPValue(v);
      double score_x = (a - 1) * log(x) - a_log_tau;

      if (score_x > 0)
        nPos++;
      else
        nNeg++;

      score.set(v, score_x);
    }

    if (g_verbosity >= VERBOSE_ESSENTIAL)
    {
      std::cerr << "// Number of positive nodes: " << nPos << std::endl;
      std::cerr << "// Number of negative nodes: " << nNeg << std::endl;
      std::cerr << "// Fraction of positive nodes: "
                << (double)nPos / (double)getOrgNodeCount() << std::endl;
    }
  }

protected:
  virtual void initParserMembers(Graph*& pG,
                                 LabelNodeMap*& pLabel,
                                 WeightNodeMap*& pPVal,
                                 WeightNodeMap*& pScore)
  {
    pG = new Graph();
    pLabel = new LabelNodeMap(*pG);
    pPVal = new WeightNodeMap(*pG);
    pScore = new WeightNodeMap(*pG);
  }
};

} // namespace mwcs
} // namespace nina

#endif // MWCSGRAPHPARSER_H
