/*
 *  heinz_mwcs_mc.cpp
 *
 *   Created on: 8-nov-2014
 *       Author: M. El-Kebir
 */

#include <iostream>
#include <fstream>
#include <lemon/arg_parser.h>
#include <lemon/time_measure.h>

#include "utils.h"
#include "config.h"
#include "parser/stpparser.h"
#include "mwcspreprocessedgraph.h"

#include "solver/solver.h"
#include "solver/solverrooted.h"
#include "solver/solverunrooted.h"
#include "solver/impl/treeheuristicsolverrootedimpl.h"
#include "solver/impl/treeheuristicsolverunrootedimpl.h"

#define PROBLEM "MWCS"
#define METHOD "heinz-mwcs-mc"

using namespace nina::mwcs;

typedef StpParser<Graph> StpParserType;
typedef MwcsPreprocessedGraph<Graph> MwcsPreprocessedGraphType;
typedef Solver<Graph> SolverType;
typedef SolverRooted<Graph> SolverRootedType;
typedef SolverUnrooted<Graph> SolverUnrootedType;
typedef TreeHeuristicSolverImpl<Graph> TreeHeuristicSolverImplType;
typedef TreeHeuristicSolverUnrootedImpl<Graph> TreeHeuristicSolverUnrootedImplType;
typedef TreeHeuristicSolverRootedImpl<Graph> TreeHeuristicSolverRootedImplType;
typedef TreeHeuristicSolverImplType::Options Options;
typedef SolverType::NodeSet NodeSet;

void printUsage(std::ostream& out, const char* argv0)
{
  out << "Usage: " << argv0 << " filename time threads outputfile" << std::endl;
}

int main(int argc, char** argv)
{
  if (argc != 5)
  {
    printUsage(std::cerr, argv[0]);
    return 1;
  }
  
  const std::string input = argv[1];
  const int timelimit = atoi(argv[2]);
  const int threads = atoi(argv[3]);
  const std::string output = argv[4];
  
  if (timelimit <= 0)
  {
    std::cerr << "Invalid timelimit '" << timelimit << "'" << std::endl;
    return 1;
  }
  
  if (threads <= 0)
  {
    std::cerr << "Invalid thread count '" << threads << "'" << std::endl;
    return 1;
  }
  
  bool std_out_used = false;
  if (output != "-")
  {
    g_pOut = new std::ofstream(output.c_str());
    if (!g_pOut->good())
    {
      std::cerr << "Could not open file '" << output << "' for writing" << std::endl;
      delete g_pOut;
      return 1;
    }
  }
  else
  {
    std_out_used = true;
    g_pOut = &std::cout;
  }
  g_pOut->setf(std::ios_base::fixed);
  g_pOut->precision(9);

  g_verbosity = VERBOSE_NONE;
  g_verbosity = VERBOSE_NON_ESSENTIAL;
  StpParserType parser(input);
  
  MwcsPreprocessedGraphType instance;
  if (!instance.init(&parser, false))
  {
    return 1;
  }
  instance.preprocess(NodeSet());
  
  Options options(TreeHeuristicSolverImplType::EDGE_COST_RANDOM,
                  false,
                  -1,
                  timelimit - g_timer.realTime());

  printCommentSection(parser.getName(), PROBLEM, METHOD, HEINZ_VERSION);
  
  SolverUnrootedType solver(new TreeHeuristicSolverUnrootedImplType(options));
  
  *g_pOut << "SECTION Solutions" << std::endl;
  solver.solve(instance);
  *g_pOut << "End" << std::endl << std::endl;
  
  printRunSection(1, solver.getSolutionWeight(), -1);
  
  *g_pOut << "SECTION Finalsolution" << std::endl;
  instance.printMwcsDimacs(solver.getSolutionModule(), *g_pOut);
  *g_pOut << "End" << std::endl;
  
  if (!std_out_used)
    delete g_pOut;
  
  return 0;
}
