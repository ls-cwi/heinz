#ifndef UTILS_H
#define UTILS_H

#include <lemon/list_graph.h>

namespace nina {
namespace mwcs {

typedef lemon::ListGraph Graph;
GRAPH_TYPEDEFS(Graph);

void generateRandomGraph(Graph& g, Graph::NodeMap<int> &weight, int nNodes, int nEdges);
void generatePCSTGraph(Graph &graph, Graph::NodeMap<int> &weight);
void generateSimpleGraph(Graph &graph, Graph::NodeMap<int> &weight);

} // namespace mwcs
} // namespace nina

#endif // UTILS_H
