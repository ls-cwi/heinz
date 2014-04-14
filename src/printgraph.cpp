/*
 *  printgraph.cpp
 *
 *   Created on: 18-jan-2013
 *       Author: M. El-Kebir
 */

#include <iostream>
#include <fstream>
#include <lemon/arg_parser.h>
#include "parser/mwcsparser.h"
#include "mwcsgraph.h"
#include "mwcsgraphparser.h"
#include "mwcspreprocessedgraph.h"
#include "preprocessing/negdeg01.h"
#include "preprocessing/posedge.h"
#include "preprocessing/negedge.h"
#include "preprocessing/rootedposdeg01.h"
#include "preprocessing/negcircuit.h"
#include "preprocessing/negdiamond.h"
#include "preprocessing/negmirroredhubs.h"
#include "preprocessing/posdeg01.h"
#include "preprocessing/posdiamond.h"
#include "preprocessing/negbicomponent.h"
#include "preprocessing/negtricomponent.h"
#include "utils.h"

using namespace nina;
using namespace nina::mwcs;

typedef Parser<Graph> ParserType;
typedef MwcsParser<Graph> MwcsParserType;
typedef MwcsGraphParser<Graph> MwcsGraphType;
typedef MwcsPreprocessedGraph<Graph> MwcsPreprocessedGraphType;
typedef NegDeg01<Graph> NegDeg01Type;
typedef PosEdge<Graph> PosEdgeType;
typedef NegEdge<Graph> NegEdgeType;
typedef RootedPosDeg01<Graph> RootedPosDeg01Type;
typedef NegCircuit<Graph> NegCircuitType;
typedef NegDiamond<Graph> NegDiamondType;
typedef NegMirroredHubs<Graph> NegMirroredHubsType;
typedef PosDeg01<Graph> PosDeg01Type;
typedef PosDiamond<Graph> PosDiamondType;
typedef NegBiComponent<Graph> NegBiComponentType;
typedef NegTriComponent<Graph> NegTriComponentType;

void pairs(const MwcsGraphType& mwcsGraph)
{
  const Graph& g = mwcsGraph.getGraph();
  const MwcsGraphType::ArcLookUpType& arcLookUp = mwcsGraph.getArcLookUp();
  BoolNodeMap present(g, true);
  IntNodeMap comp(g, -1);
  lemon::FilterNodes<const Graph> subG(g, present);
  
  for (NodeIt v(g); v != lemon::INVALID; ++v)
  {
    present[v] = false;
    for (NodeIt w = v; w != lemon::INVALID; ++w)
    {
      if (w == v) continue;
      if (arcLookUp(w, v) != lemon::INVALID) continue;
      
      present[w] = false;
      
      // determine connected components...
      int nComp = lemon::connectedComponents(subG, comp);
      if (nComp > 1)
      {
        // determine smallest component size
        std::vector<int> compSizes(nComp, 0);
        for (NodeIt u(g); u != lemon::INVALID; ++u)
        {
          if (present[u]) ++compSizes[comp[u]];
        }
        
        int maxComp = 0;
        for (size_t i = 0; i < compSizes.size(); ++i)
        {
          if (maxComp < compSizes[i]) maxComp = compSizes[i];
        }
        
        if (mwcsGraph.getNodeCount() - 2 - maxComp > 1)
        {
          std::cerr << mwcsGraph.getNodeCount() - 2 - maxComp << "\t" << mwcsGraph.getLabel(v) << " -- " << mwcsGraph.getLabel(w) << std::endl;
        }
      }
      
      present[w] = true;
      
    }
    present[v] = true;
  }
}

int main (int argc, char** argv)
{
  // parse command line arguments
  bool preprocess = false;
  double lambda = 0;
  double a = 0;
  double fdr = 0;
  int verbosityLevel = 2;
  
  std::string nodeListOutput;
  std::string edgeListOutput;

  lemon::ArgParser ap(argc, argv);

  ap
    .refOption("p", "Enable preprocessing", preprocess, false)
    .refOption("lambda", "Specifies lambda", lambda, false)
    .refOption("a", "Specifies a", a, false)
    .refOption("FDR", "Specifies fdr", fdr, false)
    .refOption("v", "Specifies the verbosity level:\n"
               "     0 - No output\n"
               "     1 - Only necessary output\n"
               "     2 - More verbose output (default)\n"
               "     3 - Debug output", verbosityLevel, false)
    .refOption("no", "Node list output filename", nodeListOutput, false)
    .refOption("eo", "Edge list output filename", edgeListOutput, false);
  ap.parse();
  
  g_verbosity = static_cast<VerbosityLevel>(verbosityLevel);

  bool pval = ap.given("FDR") && ap.given("lambda") && ap.given("a");
  if (pval)
  {
    // check if ok
    if (!(0 <= fdr && fdr <= 1))
    {
      std::cerr << "Value of FDR should be in the range [0,1]" << std::endl;
      return 1;
    }
    if (!(0 <= lambda && lambda <= 1))
    {
      std::cerr << "Value of lambda should be in the range [0,1]" << std::endl;
      return 1;
    }
    if (!(0 <= a && a <= 1))
    {
      std::cerr << "Value of a should be in the range [0,1]" << std::endl;
      return 1;
    }
  }

  // Construct parser
  ParserType* pParser = NULL;
  if (ap.files().size() < 2)
  {
    std::cerr << "Please provide two input files: first for the nodes "
              << "and second for the edges" << std::endl;
    return 1;
  }
  pParser = new MwcsParserType(ap.files()[0], ap.files()[1]);

  // Parse the input graph file and preprocess
  MwcsGraphType* pMwcs;
  if (preprocess)
  {
    MwcsPreprocessedGraphType* pPreprocessedMwcs = new MwcsPreprocessedGraphType();
    pMwcs = pPreprocessedMwcs;
    pPreprocessedMwcs->addPreprocessRule(1, new NegDeg01Type());
    pPreprocessedMwcs->addPreprocessRule(1, new PosEdgeType());
    pPreprocessedMwcs->addPreprocessRule(1, new NegEdgeType());
    pPreprocessedMwcs->addPreprocessRootRule(1, new RootedPosDeg01Type());
    pPreprocessedMwcs->addPreprocessRule(1, new NegCircuitType());
    pPreprocessedMwcs->addPreprocessRule(1, new NegDiamondType());
    pPreprocessedMwcs->addPreprocessRule(1, new PosDeg01Type());
    pPreprocessedMwcs->addPreprocessRule(1, new PosDiamondType());
    pPreprocessedMwcs->addPreprocessRule(2, new NegMirroredHubsType());
    pPreprocessedMwcs->addPreprocessRule(3, new NegBiComponentType());
    pPreprocessedMwcs->addPreprocessRule(4, new NegTriComponentType());
  }
  else
  {
    pMwcs = new MwcsGraphType();
  }

  if (!pMwcs->init(pParser, pval))
  {
    delete pParser;
    return 1;
  }

  if (pval)
  {
    pMwcs->computeScores(lambda, a, fdr);
  }

  // Now let's print the graph
  //pairs(*pMwcs);
//  BoolNodeMap cut(pMwcs->getGraph());
//  std::cerr << "#articulation nodes: " << lemon::biNodeConnectedCutNodes(pMwcs->getGraph(), cut) << std::endl;
//  for (NodeIt v(pMwcs->getGraph()); v != lemon::INVALID; ++v)
//  {
//    if (cut[v])
//    {
//      std::cerr << pMwcs->getLabel(v) << std::endl;
//    }
//  }
//  std::cerr << "Bi-edge connected components: " << lemon::biEdgeConnected(pMwcs->getGraph()) << std::endl;
  pMwcs->print(std::cout);
  
  if (!edgeListOutput.empty())
  {
    std::ofstream outFile(edgeListOutput.c_str());
    pMwcs->printEdgeList(outFile);
    outFile.flush();
    outFile.close();
  }
  
  if (!nodeListOutput.empty())
  {
    std::ofstream outFile(nodeListOutput.c_str());
    pMwcs->printNodeList(outFile);
    outFile.flush();
    outFile.close();
  }

  delete pMwcs;
  delete pParser;

  return 0;
}
