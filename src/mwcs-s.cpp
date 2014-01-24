/*
 *  mwcs-s.cpp
 *
 *   Created on: 2-jul-2013
 *       Author: M. El-Kebir
 */

#include <iostream>
#include <lemon/arg_parser.h>

// ILOG stuff
#include <ilconcert/iloalg.h>
#include <ilcplex/ilocplex.h>

#include "common/io/input/lgfparser.h"
#include "mwcsparser.h"
#include "mwcsgraph.h"
#include "mwcsgraphparser.h"
#include "mwcssizecutsolver.h"
#include "mwcssizetreesolver.h"
#include "mwcssizetreememsolver.h"
#include "mwcsenumerate.h"
#include "mwcsenumeratecomp.h"
#include "mwcsenumerateroot.h"
#include "utils.h"
#include <lemon/time_measure.h>

using namespace nina::mwcs;
using namespace nina;

typedef Parser<Graph> ParserType;
typedef LgfParser<Graph> LgfParserType;
typedef MwcsParser<Graph> MwcsParserType;
typedef MwcsGraphParser<Graph> MwcsGraphType;
typedef MwcsSolver<Graph> MwcsSolverType;
typedef MwcsSizeCutSolver<Graph> MwcsSizeCutSolverType;
typedef MwcsSizeTreeSolver<Graph> MwcsSizeTreeSolverType;
typedef MwcsSizeTreeMemSolver<Graph> MwcsSizeTreeMemSolverType;
typedef MwcsEnumerate<Graph> MwcsEnumerateType;
typedef MwcsEnumerateComp<Graph> MwcsEnumerateCompType;
typedef MwcsEnumerateRoot<Graph> MwcsEnumerateRootType;

int main (int argc, char** argv)
{
  // parse command line arguments
  int verbosityLevel = 2;
  int inputFormat = 1;
  int maxNumberOfCuts = -1;
  int enumerationMode = 3;
  int timeLimit = -1;
  int multiThreading = 1;
  std::string root;
  std::string outputFile;
  double lambda = 0;
  double a = 0;
  double fdr = 0;
  std::string nodeFile;
  std::string edgeFile;
  int moduleSize = -1;
  bool tree = false;

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
    .refOption("e", "Edge list file", edgeFile, true)
    .refOption("n", "Node file", nodeFile, true)
    .refOption("s", "Module size", moduleSize, true)
    .refOption("tree", "Tree", tree, false)
    .refOption("v", "Specifies the verbosity level:\n"
                    "     0 - No output\n"
                    "     1 - Only necessary output\n"
                    "     2 - More verbose output (default)\n"
                    "     3 - Debug output", verbosityLevel, false)
    .refOption("o", "Output file", outputFile, false)
    .refOption("m", "Specifies number of threads (default: 1)", multiThreading, false)
    .synonym("-verbosity", "v")
    .refOption("if", "Specifies the input file format:\n"
                     "     0 - LGF format\n"
                     "     1 - MWCS format (default)", inputFormat, false)
    .refOption("r", "Specifies the root node (optional)", root, false)
    .refOption("lambda", "Specifies lambda", lambda, false)
    .refOption("a", "Specifies a", a, false)
    .refOption("FDR", "Specifies fdr", fdr, false)
    .refOption("maxCuts", "Specifies the maximum number of cuts per step\n"
                          "     (only in conjuction with -f 2, optional, default: -1)", maxNumberOfCuts, false);
  ap.parse();

  if (ap.given("version"))
  {
    std::cout << "Version number: " << SVN_REV << std::endl;
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
  MwcsGraphType* pMwcs = new MwcsGraphType();

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

  if (pMwcs->getComponentCount() == 1)
    enumerationMode = 0;

  // Solve
  lemon::Timer t;
  switch (enumerationMode)
  {
    case 3:
      {
        MwcsEnumerateCompType mwcsEnumerate(*pMwcs);
        mwcsEnumerate.setTimeLimit(timeLimit);
        mwcsEnumerate.setMultiThreading(multiThreading);
        mwcsEnumerate.setModuleSize(moduleSize);
        if (tree)
          mwcsEnumerate.enumerate(MwcsSizeSolverTreeDP, false);
        else
          mwcsEnumerate.enumerate(MwcsSizeSolverCutNodeSeparatorBk, false);

        const Graph& g = pMwcs->getGraph();
        double maxScore = -std::numeric_limits<double>::max();
        double maxIdx = -1;
        for (NodeIt v(g); v != lemon::INVALID; ++v)
        {
          double score = mwcsEnumerate.getModuleWeight(v);
          if (score > maxScore)
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

        MwcsSolverType* pSolver = NULL;
        if (tree)
          pSolver = new MwcsSizeTreeMemSolverType(*pMwcs, moduleSize);
        else
          pSolver = new MwcsSizeCutSolverType(*pMwcs, moduleSize,
                                              MwcsSizeCutSolverType::MWCS_CUT_NODE_SEPARATOR_BK,
                                              maxNumberOfCuts, timeLimit,
                                              multiThreading);

        pSolver->init(rootNode);
        if (pSolver->solve())
        {
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
        }
        else
        {
          std::cerr << "No feasible solution" << std::endl;
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
