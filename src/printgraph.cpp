/*
 *  printgraph.cpp
 *
 *   Created on: 18-jan-2013
 *       Author: M. El-Kebir
 */

#include <iostream>
#include <lemon/arg_parser.h>
#include "parser/mwcsparser.h"
#include "mwcsgraph.h"
#include "mwcsgraphparser.h"
#include "mwcspreprocessedgraph.h"
#include "preprocessing/mwcspreprocessrulenegdeg01.h"
#include "preprocessing/mwcspreprocessruleposedge.h"
#include "preprocessing/mwcspreprocessrulenegedge.h"
#include "preprocessing/mwcspreprocessruleposdeg01.h"
#include "preprocessing/negcircuit.h"
#include "preprocessing/negdiamond.h"
#include "preprocessing/negmirroredhubs.h"
#include "preprocessing/posdeg01.h"
#include "preprocessing/posdiamond.h"
#include "utils.h"

using namespace nina;
using namespace nina::mwcs;

typedef Parser<Graph> ParserType;
typedef MwcsParser<Graph> MwcsParserType;
typedef MwcsGraphParser<Graph> MwcsGraphType;
typedef MwcsPreprocessedGraph<Graph> MwcsPreprocessedGraphType;
typedef MwcsPreprocessRuleNegDeg01<Graph> MwcsPreprocessRuleNegDeg01Type;
typedef MwcsPreprocessRulePosEdge<Graph> MwcsPreprocessRulePosEdgeType;
typedef MwcsPreprocessRuleNegEdge<Graph> MwcsPreprocessRuleNegEdgeType;
typedef MwcsPreprocessRulePosDeg01<Graph> MwcsPreprocessRulePosDeg01Type;
typedef NegCircuit<Graph> NegCircuitType;
typedef NegDiamond<Graph> NegDiamondType;
typedef NegMirroredHubs<Graph> NegMirroredHubsType;
typedef PosDeg01<Graph> PosDeg01Type;
typedef PosDiamond<Graph> PosDiamondType;

int main (int argc, char** argv)
{
  // parse command line arguments
  bool preprocess = false;
  double lambda = 0;
  double a = 0;
  double fdr = 0;

  lemon::ArgParser ap(argc, argv);

  ap
    .refOption("p", "Enable preprocessing", preprocess, false)
    .refOption("lambda", "Specifies lambda", lambda, false)
    .refOption("a", "Specifies a", a, false)
    .refOption("FDR", "Specifies fdr", fdr, false);
  ap.parse();

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


  // Construct parser
  ParserType* pParser = NULL;
  if (ap.files().size() < 2)
  {
    std::cerr << "Please provide two input files: first for the nodes "
              << "and second for the edges" << std::endl;
    return 1;
  }
  pParser = new MwcsParserType(ap.files()[0], ap.files()[1]);

  // Parse the input graph file and preprocess
  MwcsGraphType* pMwcs;
  if (preprocess)
  {
    MwcsPreprocessedGraphType* pPreprocessedMwcs = new MwcsPreprocessedGraphType();
    pMwcs = pPreprocessedMwcs;
    pPreprocessedMwcs->addPreprocessRule(new MwcsPreprocessRuleNegDeg01Type());
    pPreprocessedMwcs->addPreprocessRule(new MwcsPreprocessRulePosEdgeType());
    pPreprocessedMwcs->addPreprocessRule(new MwcsPreprocessRuleNegEdgeType());
    pPreprocessedMwcs->addPreprocessRootRule(new MwcsPreprocessRulePosDeg01Type());
    pPreprocessedMwcs->addPreprocessRule(new NegCircuitType());
    pPreprocessedMwcs->addPreprocessRule(new NegDiamondType());
    pPreprocessedMwcs->addPreprocessRule(new NegMirroredHubsType());
    pPreprocessedMwcs->addPreprocessRule(new PosDeg01Type());
    pPreprocessedMwcs->addPreprocessRule(new PosDiamondType());
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

  if (pval)
  {
    pMwcs->computeScores(lambda, a, fdr);
  }

  // Now let's print the graph
  pMwcs->print(std::cout);

  delete pMwcs;
  delete pParser;

  return 0;
}
