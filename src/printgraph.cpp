/*
 *  printgraph.cpp
 *
 *   Created on: 18-jan-2013
 *       Author: M. El-Kebir
 */

#include <iostream>
#include <lemon/arg_parser.h>

#include "parser/mwcsparser.h"
#include "parser/stpparser.h"
#include "parser/stppcstparser.h"

#include "mwcsgraph.h"
#include "mwcsgraphparser.h"
#include "mwcspreprocessedgraph.h"

#include "utils.h"
#include "config.h"

#include <set>

using namespace nina;
using namespace nina::mwcs;

typedef Parser<Graph> ParserType;
typedef MwcsParser<Graph> MwcsParserType;
typedef StpParser<Graph> StpParserType;
typedef StpPcstParser<Graph> StpPcstParserType;

typedef MwcsGraphParser<Graph> MwcsGraphType;
typedef MwcsPreprocessedGraph<Graph> MwcsPreprocessedGraphType;
typedef std::set<Node> NodeSet;

int main (int argc, char** argv)
{
  // parse command line arguments
  bool noPreprocess = false;
  double lambda = 0;
  double a = 0;
  double fdr = 0;
  int verbosityLevel = 2;
  std::string stpFile;
  std::string stpPcstFile;
  std::string nodeFile;
  std::string edgeFile;
  
  std::string nodeListOutput;
  std::string edgeListOutput;
  std::string stpOutput;

  lemon::ArgParser ap(argc, argv);

  ap
    .boolOption("version", "Show version number")
    .refOption("p", "Disable preprocessing", noPreprocess, false)
    .refOption("lambda", "Specifies lambda", lambda, false)
    .refOption("a", "Specifies a", a, false)
    .refOption("FDR", "Specifies fdr", fdr, false)
    .refOption("v", "Specifies the verbosity level:\n"
               "     0 - No output\n"
               "     1 - Only necessary output\n"
               "     2 - More verbose output (default)\n"
               "     3 - Debug output", verbosityLevel, false)
    .refOption("stp", "STP file", stpFile, false)
    .refOption("stp-pcst", "STP-PCST file", stpPcstFile, false)
    .refOption("e", "Edge list file", edgeFile, false)
    .refOption("n", "Node file", nodeFile, false)
    .refOption("no", "Node list output filename", nodeListOutput, false)
    .refOption("eo", "Edge list output filename", edgeListOutput, false)
    .refOption("so", "STP output filename", stpOutput, false);
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
  
  if (pPreprocessedMwcs)
  {
    pPreprocessedMwcs->preprocess(NodeSet());
  }

  // compute scores
  if (pval)
  {
    if (ap.given("a") && ap.given("lambda"))
      pMwcs->computeScores(lambda, a, fdr);
    else
      pMwcs->computeScores(fdr);
  }
  
  if (!edgeListOutput.empty())
  {
    std::ofstream outFile(edgeListOutput.c_str());
    pMwcs->printEdgeList(outFile);
    outFile.flush();
    outFile.close();
  }
  
  if (!nodeListOutput.empty())
  {
    std::ofstream outFile(nodeListOutput.c_str());
    pMwcs->printNodeList(outFile);
    outFile.flush();
    outFile.close();
  }
  
  if (!stpOutput.empty())
  {
    std::ofstream outFile(stpOutput.c_str());
    pMwcs->printSTP(stpOutput, "MEK", outFile);
    outFile.flush();
    outFile.close();
  }

  delete pMwcs;
  delete pParser;

  return 0;
}
