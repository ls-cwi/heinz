/*
 *  mwcs.cpp
 *
 *   Created on: 27-jul-2012
 *      Authors: C.I. Bucur, M. El-Kebir, G. W. Klau
 */

#include <iostream>
#include <lemon/arg_parser.h>
#include <lemon/time_measure.h>

// ILOG stuff
#include <ilconcert/iloalg.h>
#include <ilcplex/ilocplex.h>

#include "parser/mwcsparser.h"
#include "parser/stpparser.h"
#include "parser/stppcstparser.h"

#include "mwcsgraph.h"
#include "mwcsgraphparser.h"
#include "mwcspreprocessedgraph.h"

#include "solver/solver.h"
#include "solver/solverrooted.h"
#include "solver/solverunrooted.h"
#include "solver/enumsolverunrooted.h"
#include "solver/impl/cplexsolverimpl.h"
#include "solver/impl/cutsolverrootedimpl.h"
#include "solver/impl/cutsolverunrootedimpl.h"
#include "solver/impl/cplex_cut/backoff.h"

#include "mwcs.h"
#include "utils.h"
#include "config.h"

using namespace nina::mwcs;
using namespace nina;

typedef Parser<Graph> ParserType;
typedef MwcsParser<Graph> MwcsParserType;
typedef StpParser<Graph> StpParserType;
typedef StpPcstParser<Graph> StpPcstParserType;

typedef MwcsGraphParser<Graph> MwcsGraphType;
typedef MwcsPreprocessedGraph<Graph> MwcsPreprocessedGraphType;

typedef Solver<Graph> SolverType;
typedef SolverRooted<Graph> SolverRootedType;
typedef SolverUnrooted<Graph> SolverUnrootedType;
typedef EnumSolverUnrooted<Graph> EnumSolverUnrootedType;
typedef CplexSolverImpl<Graph> CplexSolverImplType;
typedef CplexSolverImplType::Options Options;
typedef CutSolverRootedImpl<Graph> CutSolverRootedImplType;
typedef CutSolverUnrootedImpl<Graph> CutSolverUnrootedImplType;
typedef SolverType::NodeSet NodeSet;

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
  int verbosityLevel = 2;
  int maxNumberOfCuts = 3;
  int timeLimit = -1;
  int memoryLimit = -1;
  int k = 1;
  bool noPreprocess = false;
  bool noEnum = false;
  int multiThreading = 1;
  int backOffFunction = 1;
  int backOffPeriod = 1;
  std::string root;
  std::string outputFile;
  double lambda = 0;
  double a = 0;
  double fdr = 0;
  std::string stpFile;
  std::string stpPcstFile;
  std::string nodeFile;
  std::string edgeFile;

  lemon::ArgParser ap(argc, argv);

  ap
    .boolOption("version", "Show version number")
    .refOption("t", "Time limit (in seconds, default: -1)", timeLimit, false)
    .refOption("ml", "Memory limit (in MB, default: -1)", memoryLimit, false)
    .refOption("e", "Edge list file", edgeFile, false)
    .refOption("n", "Node file", nodeFile, false)
    .refOption("k", "number of modules (1)", k, false)
    .refOption("period", "Back-off period (default: 1)", backOffPeriod, false)
    .refOption("b", "Back-off function:\n"
                        "     0 - Constant waiting (period: 1, override with '-period')\n"
                        "     1 - Linear waiting (default)\n"
                        "     2 - Quadratic waiting\n"
                        "     3 - Exponential waiting\n"
                        "     4 - Infinite waiting", backOffFunction, false)
    .refOption("no-pre", "Disable preprocessing", noPreprocess, false)
    .refOption("no-enum", "Disable enumerator", noEnum, false)
    .refOption("stp", "STP file", stpFile, false)
    .refOption("stp-pcst", "STP-PCST file", stpPcstFile, false)
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

  if (!(ap.given("n") && ap.given("e")) && !ap.given("stp") &&  !ap.given("stp-pcst"))
  {
    std::cerr << "Please specify either '-n' and '-e', or '-stp', or '-stp-pcst'" << std::endl;
    return 1;
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
  if (!stpFile.empty())
  {
    pParser = new StpParserType(stpFile);
  }
  else if (!stpPcstFile.empty())
  {
    pParser = new StpPcstParserType(stpPcstFile);
  }
  else
  {
    pParser = new MwcsParserType(nodeFile, edgeFile);
  }

  // Parse the input graph file and preprocess
  MwcsGraphType* pMwcs;
  MwcsPreprocessedGraphType* pPreprocessedMwcs = NULL;
  if (!noPreprocess)
  {
    pMwcs = pPreprocessedMwcs = new MwcsPreprocessedGraphType();
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
  const NodeSet rootNodeSet = pMwcs->getNodeByLabel(root);
  assert(rootNodeSet.size() == 0 || rootNodeSet.size() == 1);

  if (pPreprocessedMwcs && (noEnum || rootNodeSet.size() > 0))
  {
    pPreprocessedMwcs->preprocess(rootNodeSet);
  }

  SolverType* pSolver = NULL;

  Options options(createBackOff(backOffFunction, backOffPeriod),
                  true,
                  maxNumberOfCuts,
                  timeLimit,
                  multiThreading,
                  memoryLimit,
                  !stpPcstFile.empty(),
                  k);

  if (rootNodeSet.size() == 0 && !root.empty())
  {
    std::cerr << "No node with label '" << root
              << "' present. Defaulting to unrooted formulation." << std::endl;
  }

  if (rootNodeSet.size() == 1)
  {
    SolverRootedType* pSolverRooted = new SolverRootedType(new CutSolverRootedImplType(options));
    pSolverRooted->solve(*pMwcs, rootNodeSet);
    pSolver = pSolverRooted;
  }
  else if (noEnum)
  {
    SolverUnrootedType* pSolverUnrooted = new SolverUnrootedType(new CutSolverUnrootedImplType(options));
    pSolverUnrooted->solve(*pMwcs);
    pSolver = pSolverUnrooted;
  }
  else
  {
    SolverUnrootedType* pSolverUnrooted = new EnumSolverUnrootedType(new CutSolverUnrootedImplType(options),
                                                                     new CutSolverRootedImplType(options),
                                                                     !noPreprocess);
    pSolverUnrooted->solve(*pMwcs);
    pSolver = pSolverUnrooted;
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

  std::cerr << "Time: " << g_timer.realTime() << "s" << std::endl;

  delete pParser;
  delete pMwcs;

  return 0;
}
