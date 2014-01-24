/*
 * mainR.cpp
 *
 *  Created on: 26-oct-2012
 *     Authors: M. El-Kebir
 */


#include <Rcpp.h>

// ILOG stuff
#include <ilconcert/iloalg.h>
#include <ilcplex/ilocplex.h>

#include "io/input/lgfparser.h"
#include "io/input/identityparser.h"
#include "mwcsparser.h"
#include "mwcsgraph.h"
#include "mwcsscfsolver.h"
#include "mwcsmcfsolver.h"
#include "utils.h"

// make -f Makefile.R mwcsR.so
typedef IdentityParser<Graph> ParserType;
typedef LgfParser<Graph> LgfParserType;
typedef MwcsParser<Graph> MwcsParserType;

RcppExport SEXP solveMWCS(SEXP nodesR, SEXP nodeWeightsR,
                          SEXP sourcesR, SEXP targetsR, SEXP edgeWeightsR)
{
  Rcpp::CharacterVector xnodes(nodesR);
  Rcpp::CharacterVector xsources(sourcesR);
  Rcpp::CharacterVector xtargets(targetsR);
  Rcpp::NumericVector xnodeWeights(nodeWeightsR);
  Rcpp::NumericVector xedgeWeights(edgeWeightsR);

  Graph g;
  g.reserveNode(static_cast<int>(xnodes.size()));
  g.reserveEdge(static_cast<int>(xsources.size()));

  g_verbosity = VERBOSE_ESSENTIAL;

  // the LEMON node/edge maps
  ParserType::IdNodeMap label(g);
  ParserType::WeightNodeMap weight(g);
  ParserType::WeightEdgeMap weightEdge(g);

  // let's construct the graph, starting with the nodes
  ParserType::InvIdNodeMap invNodeMap;
  for (int i = 0; i < xnodes.size(); i++)
  {
    Node n = g.addNode();
    std::string nodeLabel = (std::string) xnodes[i];

    invNodeMap[nodeLabel] = n;
    label[n] = nodeLabel;
    weight[n] = xnodeWeights[i];
  }

  // and now the edges
  std::vector<Edge> edgeMap;
  for (int i = 0; i < xsources.size(); i++)
  {
    Edge e = g.addEdge(invNodeMap[(std::string)xsources[i]],
                       invNodeMap[(std::string)xtargets[i]]);

    edgeMap.push_back(e);
    weightEdge[e] = xedgeWeights[i];
  }

  // now solve the instance
  ParserType parser(g, &weight, &weightEdge, &label);

  MwcsGraph mwcsGraph;
  mwcsGraph.init(&parser);

  MwcsSCFSolver solver(mwcsGraph);
  solver.init();
  solver.solve();

  // construct the solution
  std::vector<bool> solution;
  for (int i = 0; i < xnodes.size(); i++)
  {
    std::string nodeLabel = (std::string) xnodes[i];
    Node n = parser.getNodeRefMap()[invNodeMap[nodeLabel]];
    solution.push_back(solver.isNodeInSolution(n));
  }

  return Rcpp::wrap(solution);
}
