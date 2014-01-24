/*
 *  montecarlo.cpp
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
#include "preprocessing/mwcspreprocessrulenegdeg01.h"
#include "preprocessing/mwcspreprocessruleposedge.h"
#include "preprocessing/mwcspreprocessrulenegedge.h"
#include "preprocessing/mwcspreprocessruleposdeg01.h"
#include "solver/mwcstreeheuristicsolver.h"
#include "mwcs.h"
#include "utils.h"
#include "config.h"

using namespace nina::mwcs;
using namespace nina;

typedef Parser<Graph> ParserType;
typedef MwcsParser<Graph> MwcsParserType;
typedef MwcsGraphParser<Graph> MwcsGraphType;
typedef MwcsPreprocessedGraph<Graph> MwcsPreprocessedGraphType;
typedef MwcsPreprocessRuleNegDeg01<Graph> MwcsPreprocessRuleNegDeg01Type;
typedef MwcsPreprocessRulePosEdge<Graph> MwcsPreprocessRulePosEdgeType;
typedef MwcsPreprocessRuleNegEdge<Graph> MwcsPreprocessRuleNegEdgeType;
typedef MwcsPreprocessRulePosDeg01<Graph> MwcsPreprocessRulePosDeg01Type;
typedef MwcsSolver<Graph> MwcsSolverType;
typedef MwcsTreeHeuristicSolver<Graph> MwcsTreeHeuristicSolverType;
typedef MwcsTreeHeuristicSolverType::MwcsAnalyzeType MwcsAnalyzeType;

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
    .refOption("h", "Specifies heuristic", h, false)
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

  g_verbosity = static_cast<VerbosityLevel>(verbosityLevel);

  // Construct parser
  ParserType* pParser = new MwcsParserType(nodeFile, edgeFile);

  // Parse the input graph file and preprocess
  MwcsGraphType* pMwcs;
  if (!noPreprocess)
  {
    MwcsPreprocessedGraphType* pPreprocessedMwcs = new MwcsPreprocessedGraphType();
    pMwcs = pPreprocessedMwcs;
    pPreprocessedMwcs->addPreprocessRule(new MwcsPreprocessRuleNegDeg01Type());
    pPreprocessedMwcs->addPreprocessRule(new MwcsPreprocessRulePosEdgeType());
    pPreprocessedMwcs->addPreprocessRule(new MwcsPreprocessRuleNegEdgeType());
    pPreprocessedMwcs->addPreprocessRootRule(new MwcsPreprocessRulePosDeg01Type());
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
  MwcsTreeHeuristicSolverType* pTreeSolver = new MwcsTreeHeuristicSolverType(*pMwcs);

  MwcsAnalyzeType* pAnalyze = NULL;
  if (!noAnalyze)
  {
    pAnalyze = new MwcsAnalyzeType(*pMwcs);
    pAnalyze->analyzeNegHubs();
    std::cout << "// Number of beneficial negative hubs: "
              << pAnalyze->getNumberOfBeneficialNegHubs() << std::endl;
  }

  const Node rootNode = pMwcs->getNodeByLabel(root);
  pTreeSolver->init(rootNode);
  for (int i = 0; i < n; i++)
  {
    std::cerr << "\rIteration " << i << ": " << std::flush;
    pTreeSolver->computeEdgeWeights(static_cast<MwcsTreeHeuristicSolverType::EdgeHeuristic>(h),
                                    pAnalyze);
    pTreeSolver->solve();
    std::cerr << pTreeSolver->getSolutionWeight() << std::flush;
  }
  std::cerr << std::endl;

  if (outputFile != "-" && !outputFile.empty())
  {
    std::ofstream outFile(outputFile.c_str());
    pMwcs->printHeinz(pTreeSolver->getSolutionModule(), outFile);
    pMwcs->printModule(pTreeSolver->getSolutionModule(), std::cout, true);
  }
  else if (outputFile == "-")
  {
    pMwcs->printHeinz(pTreeSolver->getSolutionModule(), std::cout);
  }
  else
  {
    pMwcs->printModule(pTreeSolver->getSolutionModule(), std::cout);
  }

  delete pTreeSolver;
  std::cerr << "Score: " << pTreeSolver->getSolutionWeight() << std::endl;
  std::cerr << "Time: " << t.realTime() << "s" << std::endl;

  delete pParser;
  delete pMwcs;

  return 0;
}
