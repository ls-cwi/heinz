/*
 *  mwcs.cpp
 *
 *   Created on: 27-jul-2012
 *      Authors: C.I. Bucur and M. El-Kebir
 */

#include <iostream>
#include <lemon/arg_parser.h>
#include <lemon/time_measure.h>

// ILOG stuff
#include <ilconcert/iloalg.h>
#include <ilcplex/ilocplex.h>

#include "parser/mwcsparser.h"
#include "parser/mwcsstpparser.h"
#include "mwcsgraph.h"
#include "mwcsgraphparser.h"
#include "mwcspreprocessedgraph.h"
#include "preprocessing/mwcspreprocessrulenegdeg01.h"
#include "preprocessing/mwcspreprocessruleposedge.h"
#include "preprocessing/mwcspreprocessrulenegedge.h"
#include "preprocessing/mwcspreprocessruleposdeg01.h"
#include "preprocessing/mwcspreprocessruleneghub.h"
#include "solver/mwcsscfsolver.h"
#include "solver/mwcsmcfsolver.h"
#include "solver/mwcscutsolver.h"
#include "solver/mwcstreesolver.h"
#include "solver/mwcstreeheuristicsolver.h"
#include "mwcsanalyze.h"
#include "mwcsenumerate.h"
#include "mwcsenumeratecomp.h"
#include "mwcsenumerateroot.h"
#include "mwcs.h"
#include "utils.h"
#include "config.h"

using namespace nina::mwcs;
using namespace nina;

typedef Parser<Graph> ParserType;
typedef MwcsParser<Graph> MwcsParserType;
typedef MwcsStpParser<Graph> MwcsStpParserType;
typedef MwcsGraphParser<Graph> MwcsGraphType;
typedef MwcsPreprocessedGraph<Graph> MwcsPreprocessedGraphType;
typedef MwcsPreprocessRuleNegDeg01<Graph> MwcsPreprocessRuleNegDeg01Type;
typedef MwcsPreprocessRulePosEdge<Graph> MwcsPreprocessRulePosEdgeType;
typedef MwcsPreprocessRuleNegEdge<Graph> MwcsPreprocessRuleNegEdgeType;
typedef MwcsPreprocessRulePosDeg01<Graph> MwcsPreprocessRulePosDeg01Type;
typedef MwcsPreprocessRuleNegHub<Graph> MwcsPreprocessRuleNegHubType;
typedef MwcsSolver<Graph> MwcsSolverType;
typedef MwcsSCFSolver<Graph> MwcsSCFSolverType;
typedef MwcsMCFSolver<Graph> MwcsMCFSolverType;
typedef MwcsCutSolver<Graph> MwcsCutSolverType;
typedef MwcsTreeSolver<Graph> MwcsTreeSolverType;
typedef MwcsTreeHeuristicSolver<Graph> MwcsTreeHeuristicSolverType;
typedef MwcsAnalyze<Graph> MwcsAnalyzeType;
typedef MwcsEnumerate<Graph> MwcsEnumerateType;
typedef MwcsEnumerateComp<Graph> MwcsEnumerateCompType;
typedef MwcsEnumerateRoot<Graph> MwcsEnumerateRootType;

int main(int argc, char** argv)
{
  // parse command line arguments
  int formulation = 5;
  int verbosityLevel = 2;
  int maxNumberOfCuts = -1;
  int enumerationMode = 3;
  int timeLimit = -1;
  bool preprocess = false;
  int multiThreading = 1;
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
    .refOption("enum", "Enumeration mode:\n"
                        "     0 - No enumeration\n"
                        "     1 - No root\n"
                        "     2 - Fix root\n"
                        "     3 - No root per component (default)",
               enumerationMode, false)
    .refOption("t", "Time limit (in seconds, default: -1)", timeLimit, false)
    .refOption("e", "Edge list file", edgeFile, false)
    .refOption("n", "Node file", nodeFile, false)
    .refOption("f", "Formulation of the problem:\n"
                        "     0 - Single Commodity Flow\n"
                        "     1 - Multi Commodity Flow\n"
                        "     2 - Cut formulation (Flow) \n"
                        "     3 - Cut formulation (Flow-min)\n"
                        "     4 - Cut formulation (Node-separator)\n"
                        "     5 - Cut formulation (Node-separator, BK, default)\n"
                        "     6 - Tree DP\n"
                        "     7 - Tree DP heuristic (fixed_edge)\n"
                        "     8 - Tree DP heuristic (random_edge)\n"
                        "     9 - Tree DP heuristic (uniform_edge)",
               formulation, false)
    .refOption("p", "Enable preprocessing", preprocess, false)
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
    .refOption("maxCuts", "Specifies the maximum number of cuts per step\n"
                          "     (only in conjuction with -f 2, optional, default: -1)", maxNumberOfCuts, false);
  ap.parse();

  if (ap.given("version"))
  {
    std::cout << "Version number: " << HEINZ_VERSION << std::endl;
    return 0;
  }

  bool pval = ap.given("FDR");// && ap.given("lambda") && ap.given("a");
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
    pParser = new MwcsStpParserType(stpNodeFile);
  }
  else
  {
    pParser = new MwcsParserType(nodeFile, edgeFile);
  }

  // Parse the input graph file and preprocess
  MwcsGraphType* pMwcs;
  if (preprocess)
  {
    MwcsPreprocessedGraphType* pPreprocessedMwcs = new MwcsPreprocessedGraphType();
    pMwcs = pPreprocessedMwcs;
    pPreprocessedMwcs->addPreprocessRule(new MwcsPreprocessRuleNegDeg01Type());
    pPreprocessedMwcs->addPreprocessRule(new MwcsPreprocessRulePosEdgeType());
    pPreprocessedMwcs->addPreprocessRule(new MwcsPreprocessRuleNegEdgeType());
    //pPreprocessedMwcs->addPreprocessRule(new MwcsPreprocessRuleNegHubType());
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
    if (ap.given("a") && ap.given("lambda"))
      pMwcs->computeScores(lambda, a, fdr);
    else
      pMwcs->computeScores(fdr);
  }

  if (pMwcs->getComponentCount() == 1 && enumerationMode == 1)
    enumerationMode = 0;
  if (!root.empty())
    enumerationMode = 0;

  // Solve
  lemon::Timer t;
  switch (enumerationMode)
  {
    case 1:
      {
        MwcsEnumerateType mwcsEnumerate(*pMwcs);
        mwcsEnumerate.enumerate(static_cast<MwcsSolverEnum>(formulation), preprocess);
        const MwcsEnumerateType::ModuleVector& modules = mwcsEnumerate.getModules();
        pMwcs->printModules(modules, std::cout, true);

        if (!outputFile.empty())
        {
          std::ofstream outFile(outputFile.c_str());
          if (outFile)
          {
            mwcsEnumerate.printOutput(outFile);
          }
          outFile.close();
        }
      }
      break;
    case 2:
      {
        MwcsEnumerateRootType mwcsEnumerate(*pMwcs);
        mwcsEnumerate.enumerate(static_cast<MwcsSolverEnum>(formulation), preprocess);
        const MwcsEnumerateType::ModuleVector& modules = mwcsEnumerate.getModules();
        pMwcs->printModules(modules, std::cout, true);

        if (!outputFile.empty())
        {
          std::ofstream outFile(outputFile.c_str());
          if (outFile)
          {
            mwcsEnumerate.printOutput(outFile);
          }
          outFile.close();
        }
      }
      break;
    case 3:
      {
        MwcsEnumerateCompType mwcsEnumerate(*pMwcs);
        mwcsEnumerate.setTimeLimit(timeLimit);
        mwcsEnumerate.setMultiThreading(multiThreading);
        mwcsEnumerate.enumerate(static_cast<MwcsSolverEnum>(formulation), preprocess);

        const Graph& g = pMwcs->getOrgGraph();
        double maxScore = 0;
        int maxIdx = -1;
        for (NodeIt v(g); v != lemon::INVALID; ++v)
        {
          double score = mwcsEnumerate.getModuleWeight(v);
          if (score >= maxScore)
          {
            maxScore = score;
            maxIdx = mwcsEnumerate.getModuleIndex(v);
          }
        }

        if (maxIdx != -1)
        {
          if (outputFile != "-" && !outputFile.empty())
          {
            std::ofstream outFile(outputFile.c_str());
            pMwcs->printHeinzOrg(mwcsEnumerate.getModule(maxIdx), outFile);
            pMwcs->printModule(mwcsEnumerate.getModule(maxIdx), std::cout, true);
          }
          else if (outputFile == "-")
          {
            pMwcs->printHeinzOrg(mwcsEnumerate.getModule(maxIdx), std::cout);
          }
          else
          {
            pMwcs->printModule(mwcsEnumerate.getModule(maxIdx), std::cout, true);
          }
        }
      }
      break;
    case 0:
      {
        // Initialize solver
        const Node rootNode = pMwcs->getNodeByLabel(root);
        if (rootNode == lemon::INVALID && !root.empty())
        {
          std::cerr << "No node with label '" << root
                    << "' present. Defaulting to unrooted formulation." << std::endl;
        }

        MwcsSolverType* pSolver = NULL;
        switch (formulation)
        {
          case 0:
            pSolver = new MwcsSCFSolverType(*pMwcs);
            break;
          case 1:
            pSolver = new MwcsMCFSolverType(*pMwcs);
            break;
          case 2:
            pSolver = new MwcsCutSolverType(*pMwcs, MwcsCutSolverType::MWCS_CUT_FLOW, maxNumberOfCuts, timeLimit, multiThreading);
            break;
          case 3:
            pSolver = new MwcsCutSolverType(*pMwcs, MwcsCutSolverType::MWCS_CUT_FLOW_MIN, maxNumberOfCuts, timeLimit, multiThreading);
            break;
          case 4:
            pSolver = new MwcsCutSolverType(*pMwcs, MwcsCutSolverType::MWCS_CUT_NODE_SEPARATOR, maxNumberOfCuts, timeLimit, multiThreading);
            break;
          case 5:
            pSolver = new MwcsCutSolverType(*pMwcs, MwcsCutSolverType::MWCS_CUT_NODE_SEPARATOR_BK, maxNumberOfCuts, timeLimit, multiThreading);
            break;
          case 6:
            if (rootNode == lemon::INVALID && !enumerationMode)
            {
              std::cerr << "Tree DP can only be used with a specified root node" << std::endl;
              return 1;
            }
            pSolver = new MwcsTreeSolverType(*pMwcs);
            break;
          case 7:
          case 8:
          case 9:
            {
              MwcsTreeHeuristicSolverType* pTreeSolver = new MwcsTreeHeuristicSolverType(*pMwcs);
              if (formulation == 7)
                pTreeSolver->computeEdgeWeights(MwcsTreeHeuristicSolverType::EDGE_COST_FIXED);
              else if (formulation == 8)
                pTreeSolver->computeEdgeWeights(MwcsTreeHeuristicSolverType::EDGE_COST_RANDOM);
              else if (formulation == 9)
                pTreeSolver->computeEdgeWeights(MwcsTreeHeuristicSolverType::EDGE_COST_UNIFORM_RANDOM);
              pSolver = pTreeSolver;
            }
            break;
          default:
          {
            std::cerr << "Wrong formulation specified" << std::endl;
            return 1;
          }
        }

        pSolver->init(rootNode);
        pSolver->solve();
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
      }
      break;
  }
  std::cerr << "Time: " << t.realTime() << "s" << std::endl;

  delete pParser;
  delete pMwcs;

  return 0;
}
