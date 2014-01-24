#include "utils.h"

using namespace nina;
using namespace nina::mwcs;

void generateRandomGraph(Graph& g, Graph::NodeMap<int>& weight, int nNodes, int nEdges)
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

void generatePCSTGraph(Graph &graph, Graph::NodeMap<int> &weight)
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

void generateSimpleGraph(Graph &graph, Graph::NodeMap<int> &weight)
{
  Node a = graph.addNode(); weight[a] = 5;
  Node b = graph.addNode(); weight[b] = -1;
  Node c = graph.addNode(); weight[c] = 5;
  graph.addEdge(a, b);
  graph.addEdge(c, b);
}
