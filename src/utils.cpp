#include "utils.h"
#include <assert.h>

using namespace nina;
using namespace nina::mwcs;

VerbosityLevel nina::mwcs::g_verbosity = VERBOSE_ESSENTIAL;

lemon::Timer nina::mwcs::g_timer;

std::ostream* nina::mwcs::g_pOut = NULL;

void nina::mwcs::generateRandomGraph(Graph& g, Graph::NodeMap<int>& weight, int nNodes, int nEdges)
{
  for (int i = 0; i < nNodes; i++)
  {
    weight[g.addNode()] = rand() % 41 -20;
  }

  double p = ((double)nEdges) / ((nNodes - 1) * nNodes / 2);

  for (Graph::NodeIt v(g); v != lemon::INVALID; ++v)
    for (Graph::NodeIt w(g); w != lemon::INVALID; ++w)
    {
      if (((double)rand()) / RAND_MAX < p && v < w)
      {
        g.addEdge(v, w);
      }
    }
}

void nina::mwcs::generatePCSTGraph(Graph &graph, Graph::NodeMap<int> &weight)
{
  Node a = graph.addNode(); weight[a] = 5;
  Node b = graph.addNode(); weight[b] = -4;
  Node c = graph.addNode(); weight[c] = -5;
  Node d = graph.addNode(); weight[d] = 30;
  Node e = graph.addNode(); weight[e] = -5;
  Node f = graph.addNode(); weight[f] = -15;
  Node g = graph.addNode(); weight[g] = 10;
  Node h = graph.addNode(); weight[h] = 10;

  graph.addEdge(a, b);
  graph.addEdge(b, c);
  graph.addEdge(a, c);
  graph.addEdge(c, d);
  graph.addEdge(b, d);
  graph.addEdge(d, f);
  graph.addEdge(c, e);
  graph.addEdge(f, h);
  graph.addEdge(f, g);
  graph.addEdge(e, f);
}

void nina::mwcs::generateSimpleGraph(Graph &graph, Graph::NodeMap<int> &weight)
{
  Node a = graph.addNode(); weight[a] = 5;
  Node b = graph.addNode(); weight[b] = -1;
  Node c = graph.addNode(); weight[c] = 5;
  graph.addEdge(a, b);
  graph.addEdge(c, b);
}

void nina::mwcs::printCommentSection(const std::string& name,
                                     const std::string& problem,
                                     const std::string& method,
                                     const std::string& version)
{
  assert(g_pOut);
  *g_pOut << "SECTION Comment" << std::endl;
  *g_pOut << "Name " << name << std::endl;
  *g_pOut << "Problem \"" << problem << "\"" << std::endl;
  *g_pOut << "Program \"" << method << "\"" << std::endl;
  *g_pOut << "Version \"" << version << "\"" << std::endl;
  *g_pOut << "End" << std::endl;
  *g_pOut << std::endl;
}

void nina::mwcs::printRunSection(int threads, double primalObjValue, double dualObjValue)
{
  assert(g_pOut);
  *g_pOut << "SECTION Run" << std::endl;
  *g_pOut << "Threads " << threads << std::endl;
  *g_pOut << "Time " << g_timer.realTime() << std::endl;
  if (dualObjValue != -1)
  {
    *g_pOut << "Dual " << dualObjValue << std::endl;
  }
  *g_pOut << "Primal " << primalObjValue << std::endl;
  *g_pOut << "End" << std::endl;
  *g_pOut << std::endl;
}
