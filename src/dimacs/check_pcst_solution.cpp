/*
 *  check_pcst_solution.cpp
 *
 *   Created on: 20-nov-2014
 *       Author: M. El-Kebir
 */

#include <iostream>
#include <lemon/arg_parser.h>

#include "parser/mwcsparser.h"
#include "parser/stppcstparser.h"
#include "parser/dimacsparser.h"

#include "mwcs.h"
#include "utils.h"
#include "config.h"

#include "mwcsgraph.h"
#include "mwcsgraphparser.h"

#include <set>

using namespace nina::mwcs;
using namespace nina;

typedef Parser<Graph> ParserType;
typedef MwcsParser<Graph> MwcsParserType;
typedef StpPcstParser<Graph> StpPcstParserType;
typedef DimacsParser<Graph> DimacsParserType;
typedef MwcsGraphParser<Graph> MwcsGraphType;
typedef std::set<Node> NodeSet;
typedef NodeSet::const_iterator NodeSetIt;

int main(int argc, char** argv)
{
  // parse command line arguments
  int verbosityLevel = 0;
  std::string stpFile;
  std::string dimacsFile;

  lemon::ArgParser ap(argc, argv);
  
  ap
    .boolOption("version", "Show version number")
    .refOption("stp", "STP file", stpFile, true)
    .refOption("v", "Specifies the verbosity level:\n"
               "     0 - No output\n"
               "     1 - Only necessary output (default)\n"
               "     2 - More verbose output\n"
               "     3 - Debug output", verbosityLevel, false)
    .refOption("s", "DIMACS solution file", dimacsFile, true);
  ap.parse();
  
  if (ap.given("version"))
  {
    std::cout << "Version number: " << HEINZ_VERSION << std::endl;
    return 0;
  }
  
  g_verbosity = static_cast<VerbosityLevel>(verbosityLevel);
  
  // Construct parser
  StpPcstParserType* pParserInput = new StpPcstParserType(stpFile);
  
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
  const Graph& g = mwcsSolution.getGraph();
  for (NodeIt v(g); v != lemon::INVALID; ++v)
  {
    if (mwcsInput.getNodeByLabel(mwcsSolution.getLabel(v)).size() != 1)
    {
      std::cerr << "Node " << mwcsSolution.getLabel(v) << " missing in input file" << std::endl;
      return 1;
    }
  }

  const Graph& gg = mwcsInput.getGraph();
  for (EdgeIt e(g); e != lemon::INVALID; ++e)
  {
    Node u = g.u(e);
    Node v = g.v(e);
    
    const std::string& label_u = mwcsSolution.getLabel(u);
    const std::string& label_v = mwcsSolution.getLabel(v);

    if (mwcsInput.getNodeByLabel(label_u + "--" + label_v).size() == 0 &&
        mwcsInput.getNodeByLabel(label_v + "--" + label_u).size() == 0)
    {
      std::cerr << "Edge (" << mwcsSolution.getLabel(u) << "," << mwcsSolution.getLabel(v) << ")" << " missing in input file" << std::endl;
      return 1;
    }
  }
  
  // check if the solution graph is connected
  if (!lemon::connected(g))
  {
    std::cerr << "Solution graph is not connected" << std::endl;
    return 1;
  }
  
  // check if the solution graph contains the root nodes
  const NodeSet& roots = pParserInput->getRootNodes();
  for (NodeSetIt rootIt = roots.begin(); rootIt != roots.end(); ++rootIt)
  {
    Node root = *rootIt;
    if (mwcsSolution.getNodeByLabel(mwcsInput.getLabel(root)).size() == 0)
    {
      std::cerr << "Missing root node '" << mwcsInput.getLabel(root) << "'" << std::endl;
      return 1;
    }
  }
  
  delete pParserInput;
  delete pParserSolution;
  return 0;
}
