/*
 *  heinz_mwcs_dc.cpp
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
#include "solver/impl/cplexsolverimpl.h"
#include "solver/impl/cutsolverrootedimpl.h"
#include "solver/impl/cutsolverunrootedimpl.h"
#include "solver/impl/cplex_cut/backoff.h"

#define PROBLEM "MWCS"
#define METHOD "heinz-mwcs-no-dc"
#define MEMORY_LIMIT 30*1024 // 30 GB

using namespace nina::mwcs;

typedef StpParser<Graph> StpParserType;
typedef MwcsPreprocessedGraph<Graph> MwcsPreprocessedGraphType;
typedef Solver<Graph> SolverType;
typedef SolverUnrooted<Graph> SolverUnrootedType;
typedef CplexSolverImpl<Graph> CplexSolverImplType;
typedef CplexSolverImplType::Options Options;
typedef CutSolverUnrootedImpl<Graph> CutSolverUnrootedImplType;
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
  
  g_verbosity = VERBOSE_NONE;
  g_verbosity = VERBOSE_NON_ESSENTIAL;
  StpParserType parser(input);
  
  MwcsPreprocessedGraphType instance;
  if (!instance.init(&parser, false))
  {
    return 1;
  }
  instance.preprocess(NodeSet());
  
  Options options(BackOff(1), // linear waiting
                  true,
                  10,
                  timelimit,
                  threads,
                  MEMORY_LIMIT);

  printCommentSection(parser.getName(), PROBLEM, METHOD, HEINZ_VERSION);
  
  SolverUnrootedType solver(new CutSolverUnrootedImplType(options));
  
  *g_pOut << "SECTION Solutions" << std::endl;
  solver.solve(instance);
  *g_pOut << "End" << std::endl << std::endl;
  
  printRunSection(threads, solver.getSolutionWeight(), solver.getSolutionWeightUB());
  
  *g_pOut << "SECTION Finalsolution" << std::endl;
  instance.printMwcsDimacs(solver.getSolutionModule(), *g_pOut);
  *g_pOut << "End" << std::endl;
  
  if (!std_out_used)
    delete g_pOut;
  
  return 0;
}
