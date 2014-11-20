/*
 *  check_solution.cpp
 *
 *   Created on: 19-nov-2014
 *       Author: M. El-Kebir
 */

#include <iostream>
#include <lemon/arg_parser.h>

#include "parser/mwcsparser.h"
#include "parser/stpparser.h"
#include "parser/stppcstparser.h"
#include "parser/dimacsparser.h"

#include "mwcs.h"
#include "utils.h"
#include "config.h"

#include "mwcsgraph.h"
#include "mwcsgraphparser.h"

using namespace nina::mwcs;
using namespace nina;

typedef Parser<Graph> ParserType;
typedef MwcsParser<Graph> MwcsParserType;
typedef StpParser<Graph> StpParserType;
typedef StpPcstParser<Graph> StpPcstParserType;
typedef DimacsParser<Graph> DimacsParserType;
typedef MwcsGraphParser<Graph> MwcsGraphType;

int main(int argc, char** argv)
{
  // parse command line arguments
  int verbosityLevel = 2;
  std::string stpFile;
  std::string stpPcstFile;
  std::string nodeFile;
  std::string edgeFile;
  std::string dimacsFile;

  lemon::ArgParser ap(argc, argv);
  
  ap
    .boolOption("version", "Show version number")
    .refOption("e", "Edge list file", edgeFile, false)
    .refOption("n", "Node file", nodeFile, false)
    .refOption("stp", "STP file", stpFile, false)
    .refOption("stp-pcst", "STP-PCST file", stpPcstFile, false)
    .refOption("v", "Specifies the verbosity level:\n"
               "     0 - No output\n"
               "     1 - Only necessary output\n"
               "     2 - More verbose output (default)\n"
               "     3 - Debug output", verbosityLevel, false)
    .refOption("s", "DIMACS solution file", dimacsFile, true);
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
  
  g_verbosity = static_cast<VerbosityLevel>(verbosityLevel);
  
  // Construct parser
  ParserType* pParserInput = NULL;
  if (!stpFile.empty())
  {
    pParserInput = new StpParserType(stpFile);
  }
  else if (!stpPcstFile.empty())
  {
    pParserInput = new StpPcstParserType(stpPcstFile);
  }
  else
  {
    pParserInput = new MwcsParserType(nodeFile, edgeFile);
  }
  
  // Parse the input graph file
  MwcsGraphType mwcsInput;
  if (!mwcsInput.init(pParserInput, false))
  {
    delete pParserInput;
    return 1;
  }
  
  // Parse solution
  ParserType* pParserSolution = new DimacsParserType(dimacsFile);
  MwcsGraphType mwcsSolution;
  if (!mwcsSolution.init(pParserSolution, false))
  {
    delete pParserInput;
    delete pParserSolution;
    return 1;
  }
  
  // check if the solution graph is a subgraph of the input graph
  for (NodeIt v(mwcsSolution.getGraph()); v != lemon::INVALID; ++v)
  {
    if (mwcsInput.getNodeByLabel(mwcsSolution.getLabel(v)).size() == 0)
      return 1;
  }
  
  // check if the solution graph is connected
  if (!lemon::connected(mwcsSolution.getGraph()))
  {
    return 1;
  }
  
  delete pParserInput;
  delete pParserSolution;
  return 0;
}
