/*
 *  mwcs-mc.cpp
 *
 *   Created on: 29-apr-2013
 *       Author: M. El-Kebir
 */

#include <iostream>
#include <lemon/arg_parser.h>
#include <lemon/time_measure.h>

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
#include "solver/solver.h"
#include "solver/solverrooted.h"
#include "solver/solverunrooted.h"
#include "solver/impl/treeheuristicsolverrootedimpl.h"
#include "solver/impl/treeheuristicsolverunrootedimpl.h"
#include "mwcs.h"
#include "utils.h"
#include "config.h"

using namespace nina::mwcs;
using namespace nina;

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

typedef Solver<Graph> SolverType;
typedef SolverRooted<Graph> SolverRootedType;
typedef SolverUnrooted<Graph> SolverUnrootedType;
typedef TreeHeuristicSolverImpl<Graph> TreeHeuristicSolverImplType;
typedef TreeHeuristicSolverUnrootedImpl<Graph> TreeHeuristicSolverUnrootedImplType;
typedef TreeHeuristicSolverRootedImpl<Graph> TreeHeuristicSolverRootedImplType;
typedef typename SolverType::NodeSet NodeSet;
typedef typename TreeHeuristicSolverImplType::Options Options;

int main (int argc, char** argv)
{
  // parse command line arguments
  int verbosityLevel = 2;
  bool noPreprocess = false;
  bool noAnalyze = false;
  std::string root;
  std::string outputFile;
  double lambda = 0;
  double a = 0;
  double fdr = 0;
  std::string nodeFile;
  std::string edgeFile;
  int n = 100;
  int h = 3;

  lemon::ArgParser ap(argc, argv);

  ap
    .boolOption("version", "Show version number")
    .refOption("p", "Disable preprocessing", noPreprocess, false)
    .refOption("e", "Edge list file", edgeFile, true)
    .refOption("n", "Node file", nodeFile, true)
    .refOption("v", "Specifies the verbosity level:\n"
                    "     0 - No output\n"
                    "     1 - Only necessary output\n"
                    "     2 - More verbose output (default)\n"
                    "     3 - Debug output", verbosityLevel, false)
    .synonym("-verbosity", "v")
    .refOption("r", "Specifies the root node (optional)", root, false)
    .refOption("lambda", "Specifies lambda", lambda, false)
    .refOption("a", "Specifies a", a, false)
    .refOption("FDR", "Specifies fdr", fdr, false)
    .refOption("o", "Output file", outputFile, false)
    .refOption("h", "Specifies heuristic\n"
                    "     0 - fixed_edge\n"
                    "     1 - random_edge\n"
                    "     2 - uniform_edge\n"
                    "     3 - min_max_edge (default)", h, false)
    .refOption("m", "Number of Monte Carlo iterations", n, false)
    .refOption("z", "Disable negative hubs sampling", noAnalyze, false);
  ap.parse();

  if (ap.given("version"))
  {
    std::cout << "Version number: " << HEINZ_VERSION << std::endl;
    return 0;
  }

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
  
  Options options(static_cast<TreeHeuristicSolverImplType::EdgeHeuristic>(h),
                  !noAnalyze, n);

  g_verbosity = static_cast<VerbosityLevel>(verbosityLevel);

  // Construct parser
  ParserType* pParser = new MwcsParserType(nodeFile, edgeFile);

  // Parse the input graph file and preprocess
  MwcsGraphType* pMwcs;
  if (!noPreprocess)
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

  // compute scores
  if (pval)
  {
    pMwcs->computeScores(lambda, a, fdr);
  }

  // Solve
  lemon::Timer t;

  const Node rootNode = pMwcs->getNodeByLabel(root);
  SolverType* pSolver = NULL;
  if (rootNode != lemon::INVALID)
  {
    NodeSet rootNodes;
    rootNodes.insert(rootNode);
    
    SolverRootedType* pSolverRooted = new SolverRootedType(new TreeHeuristicSolverRootedImplType(options));
    pSolverRooted->solve(*pMwcs, rootNodes);
    pSolver = pSolverRooted;
  }
  else
  {
    SolverUnrootedType* pSolverUnrooted = new SolverUnrootedType(new TreeHeuristicSolverUnrootedImplType(options));
    pSolverUnrooted->solve(*pMwcs);
    pSolver = pSolverUnrooted;
  }

//  for (int i = 0; i < n; i++)
//  {
//    std::cerr << "\rIteration " << i << ": " << std::flush;
//    pTreeSolver->computeEdgeCosts(*pMwcs);
//    pTreeSolver->solve();
//    std::cerr << pTreeSolver->getSolutionWeight() << std::flush;
//  }
//  std::cerr << std::endl;

  if (outputFile != "-" && !outputFile.empty())
  {
    std::ofstream outFile(outputFile.c_str());
    pMwcs->printHeinz(pSolver->getSolutionModule(), outFile);
    pMwcs->printModule(pSolver->getSolutionModule(), std::cout, true);
  }
  else if (outputFile == "-")
  {
    pMwcs->printHeinz(pSolver->getSolutionModule(), std::cout);
  }
  else
  {
    pMwcs->printModule(pSolver->getSolutionModule(), std::cout);
  }

  std::cerr << "Score: " << pSolver->getSolutionWeight() << std::endl;
  std::cerr << "Time: " << t.realTime() << "s" << std::endl;

  delete pSolver;
  delete pParser;
  delete pMwcs;

  return 0;
}
