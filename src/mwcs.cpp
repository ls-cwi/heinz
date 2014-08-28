/*
 *  mwcs.cpp
 *
 *   Created on: 27-jul-2012
 *      Authors: C.I. Bucur and M. El-Kebir
 */

// TODO: determine which preprocessing rules are invalid for rooted formulation!

#include <iostream>
#include <lemon/arg_parser.h>
#include <lemon/time_measure.h>

// ILOG stuff
#include <ilconcert/iloalg.h>
#include <ilcplex/ilocplex.h>

#include "parser/mwcsparser.h"
#include "parser/stpparser.h"
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

#include "solver/solver.h"
#include "solver/solverrooted.h"
#include "solver/solverunrooted.h"
#include "solver/impl/cplexsolverimpl.h"
#include "solver/impl/cutsolverrootedimpl.h"
#include "solver/impl/cutsolverunrootedimpl.h"
#include "solver/impl/cplex_cut/backoff.h"

//#include "solver/enumsolverunrooted.h"
#include "mwcs.h"
#include "utils.h"
#include "config.h"

using namespace nina::mwcs;
using namespace nina;

typedef Parser<Graph> ParserType;
typedef MwcsParser<Graph> MwcsParserType;
typedef StpParser<Graph> StpParserType;
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

typedef Solver<Graph> SolverType;
typedef SolverRooted<Graph> SolverRootedType;
typedef SolverUnrooted<Graph> SolverUnrootedType;
typedef CplexSolverImpl<Graph> CplexSolverImplType;
typedef typename CplexSolverImplType::Options Options;
typedef CutSolverRootedImpl<Graph> CutSolverRootedImplType;
typedef CutSolverUnrootedImpl<Graph> CutSolverUnrootedImplType;
typedef typename SolverType::NodeSet NodeSet;

BackOff createBackOff(int function, int period)
{
  BackOff::Function f = static_cast<BackOff::Function>(function);
  switch (f)
  {
    case nina::mwcs::BackOff::ConstantWaiting:
      return BackOff(period);
    default:
      return BackOff(f);
  }
}

int main(int argc, char** argv)
{
  // parse command line arguments
  int formulation = 5;
  int verbosityLevel = 2;
  int maxNumberOfCuts = 3;
  int timeLimit = -1;
  bool noPreprocess = false;
  int multiThreading = 1;
  int backOffFunction = 1;
  int backOffPeriod = 1;
  std::string root;
  std::string outputFile;
  double lambda = 0;
  double a = 0;
  double fdr = 0;
  std::string stpNodeFile;
  std::string nodeFile;
  std::string edgeFile;

  lemon::ArgParser ap(argc, argv);

  ap
    .boolOption("version", "Show version number")
    .refOption("t", "Time limit (in seconds, default: -1)", timeLimit, false)
    .refOption("e", "Edge list file", edgeFile, true)
    .refOption("n", "Node file", nodeFile, true)
    .refOption("period", "Back-off period (default: 1)", backOffPeriod, false)
    .refOption("b", "Back-off function:\n"
                        "     0 - Constant waiting (period: 1, override with '-period')\n"
                        "     1 - Linear waiting (default)\n"
                        "     2 - Quadratic waiting\n"
                        "     3 - Exponential waiting\n"
                        "     4 - Infinite waiting", backOffFunction, false)
    .refOption("p", "Disable preprocessing", noPreprocess, false)
    .refOption("s", "STP node file", stpNodeFile, false)
    .refOption("v", "Specifies the verbosity level:\n"
                    "     0 - No output\n"
                    "     1 - Only necessary output\n"
                    "     2 - More verbose output (default)\n"
                    "     3 - Debug output", verbosityLevel, false)
    .refOption("o", "Output file", outputFile, false)
    .refOption("m", "Specifies number of threads (default: 1)", multiThreading, false)
    .synonym("-verbosity", "v")
    .refOption("r", "Specifies the root node (optional)", root, false)
    .refOption("lambda", "Specifies lambda", lambda, false)
    .refOption("a", "Specifies a", a, false)
    .refOption("FDR", "Specifies fdr", fdr, false)
    .refOption("maxCuts", "Specifies the number of cut iterations per node in the B&B tree (default: 3)\n",
               maxNumberOfCuts, false);
  ap.parse();

  if (ap.given("version"))
  {
    std::cout << "Version number: " << HEINZ_VERSION << std::endl;
    return 0;
  }

  bool pval = ap.given("FDR");
  if (pval)
  {
    // check if ok
    if (!(0 <= fdr && fdr <= 1))
    {
      std::cerr << "Value of FDR should be in the range [0,1]" << std::endl;
      return 1;
    }
    if (ap.given("lambda") && !(0 <= lambda && lambda <= 1))
    {
      std::cerr << "Value of lambda should be in the range [0,1]" << std::endl;
      return 1;
    }
    if (ap.given("a") && !(0 <= a && a <= 1))
    {
      std::cerr << "Value of a should be in the range [0,1]" << std::endl;
      return 1;
    }
  }

  g_verbosity = static_cast<VerbosityLevel>(verbosityLevel);

  // Construct parser
  ParserType* pParser = NULL;
  if (!stpNodeFile.empty())
  {
    pParser = new StpParserType(stpNodeFile);
  }
  else
  {
    pParser = new MwcsParserType(nodeFile, edgeFile);
  }

  // Parse the input graph file and preprocess
  MwcsGraphType* pMwcs;
  if (!noPreprocess)
  {
    MwcsPreprocessedGraphType* pPreprocessedMwcs = new MwcsPreprocessedGraphType();
    pMwcs = pPreprocessedMwcs;
    pPreprocessedMwcs->addPreprocessRule(1, new NegDeg01Type());
    pPreprocessedMwcs->addPreprocessRule(1, new PosEdgeType());
    pPreprocessedMwcs->addPreprocessRule(1, new NegEdgeType());
    pPreprocessedMwcs->addPreprocessRule(1, new NegCircuitType());
    pPreprocessedMwcs->addPreprocessRule(1, new NegDiamondType());
    pPreprocessedMwcs->addPreprocessRule(1, new PosDiamondType());
    pPreprocessedMwcs->addPreprocessRule(2, new NegMirroredHubsType());
//    pPreprocessedMwcs->addPreprocessRule(3, new NegBiComponentType());
//    pPreprocessedMwcs->addPreprocessRule(4, new NegTriComponentType());
    if (root.empty())
    {
      pPreprocessedMwcs->addPreprocessRule(1, new PosDeg01Type());
    }
    else
    {
      pPreprocessedMwcs->addPreprocessRootRule(1, new RootedPosDeg01Type());
    }
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
    if (ap.given("a") && ap.given("lambda"))
      pMwcs->computeScores(lambda, a, fdr);
    else
      pMwcs->computeScores(fdr);
  }

  // Solve
  lemon::Timer t;
  const Node rootNode = pMwcs->getNodeByLabel(root);
  SolverType* pSolver = NULL;
  
  Options options(createBackOff(backOffFunction, backOffPeriod),
                  true,
                  maxNumberOfCuts,
                  timeLimit,
                  multiThreading);
  
  if (rootNode != lemon::INVALID)
  {
    NodeSet rootNodes;
    rootNodes.insert(rootNode);
    
    SolverRootedType* pSolverRooted = new SolverRootedType(new CutSolverRootedImplType(options));
    pSolverRooted->solve(*pMwcs, rootNodes);
    pSolver = pSolverRooted;
  }
  else
  {
    SolverUnrootedType* pSolverUnrooted = new SolverUnrootedType(new CutSolverUnrootedImplType(options));
    pSolverUnrooted->solve(*pMwcs);
    pSolver = pSolverUnrooted;
  }
  
  if (rootNode == lemon::INVALID && !root.empty())
  {
    std::cerr << "No node with label '" << root
              << "' present. Defaulting to unrooted formulation." << std::endl;
  }
  
  if (outputFile != "-" && !outputFile.empty())
  {
    std::ofstream outFile(outputFile.c_str());
    pMwcs->printHeinz(pSolver->getSolutionModule(), outFile);
    pMwcs->printModule(pSolver->getSolutionModule(), std::cout, false);
  }
  else if (outputFile == "-")
  {
    pMwcs->printHeinz(pSolver->getSolutionModule(), std::cout);
  }
  else
  {
    pMwcs->printModule(pSolver->getSolutionModule(), std::cout, false);
  }

  delete pSolver;
  
//  if (rootNode == lemon::INVALID)
//  {
//    MwcsEnumerateType mwcsEnumerate(*pMwcs);
//    mwcsEnumerate.setTimeLimit(timeLimit);
//    mwcsEnumerate.setMultiThreading(multiThreading);
//    mwcsEnumerate.setMaxNumberOfCuts(maxNumberOfCuts);
//    mwcsEnumerate.setBackOff(createBackOff(backOffFunction, backOffPeriod));
//    mwcsEnumerate.enumerate(static_cast<MwcsSolverEnum>(formulation), !noPreprocess);
//    
//    const Graph& g = pMwcs->getOrgGraph();
//    double maxScore = 0;
//    int maxIdx = -1;
//    for (NodeIt v(g); v != lemon::INVALID; ++v)
//    {
//      double score = mwcsEnumerate.getModuleWeight(v);
//      if (score >= maxScore)
//      {
//        maxScore = score;
//        maxIdx = mwcsEnumerate.getModuleIndex(v);
//      }
//    }
//    
//    if (maxIdx != -1)
//    {
//      if (outputFile != "-" && !outputFile.empty())
//      {
//        std::ofstream outFile(outputFile.c_str());
//        pMwcs->printHeinzOrg(mwcsEnumerate.getModule(maxIdx), outFile);
//        pMwcs->printModule(mwcsEnumerate.getModule(maxIdx), std::cout, true);
//      }
//      else if (outputFile == "-")
//      {
//        pMwcs->printHeinzOrg(mwcsEnumerate.getModule(maxIdx), std::cout);
//      }
//      else
//      {
//        pMwcs->printModule(mwcsEnumerate.getModule(maxIdx), std::cout, true);
//      }
//    }
//  }
//  else
//  {
//    // Initialize solver
//    const Node rootNode = pMwcs->getNodeByLabel(root);
//    if (rootNode == lemon::INVALID && !root.empty())
//    {
//      std::cerr << "No node with label '" << root
//                << "' present. Defaulting to unrooted formulation." << std::endl;
//    }
//
//    MwcsSolverType* pSolver = new MwcsCutSolverType(*pMwcs,
//                                                    createBackOff(backOffFunction, backOffPeriod),
//                                                    maxNumberOfCuts, timeLimit, multiThreading);
//    pSolver->init(rootNode);
//    pSolver->solve();
//    if (outputFile != "-" && !outputFile.empty())
//    {
//      std::ofstream outFile(outputFile.c_str());
//      pMwcs->printHeinz(pSolver->getSolutionModule(), outFile);
//      pMwcs->printModule(pSolver->getSolutionModule(), std::cout, false);
//    }
//    else if (outputFile == "-")
//    {
//      pMwcs->printHeinz(pSolver->getSolutionModule(), std::cout);
//    }
//    else
//    {
//      pMwcs->printModule(pSolver->getSolutionModule(), std::cout, false);
//    }
//
//    delete pSolver;
//  }

  std::cerr << "Time: " << t.realTime() << "s" << std::endl;

  delete pParser;
  delete pMwcs;

  return 0;
}
