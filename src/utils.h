#ifndef UTILS_H
#define UTILS_H

#include <lemon/list_graph.h>
#include <lemon/time_measure.h>
#include <ostream>

namespace nina {
namespace mwcs {

typedef lemon::ListGraph Graph;
GRAPH_TYPEDEFS(Graph);

void generateRandomGraph(Graph& g,
                         Graph::NodeMap<int> &weight,
                         int nNodes, int nEdges);
void generatePCSTGraph(Graph &graph,
                       Graph::NodeMap<int> &weight);
void generateSimpleGraph(Graph &graph,
                         Graph::NodeMap<int> &weight);
  
void printCommentSection(const std::string& name,
                         const std::string& problem,
                         const std::string& method,
                         const std::string& version);

void printRunSection(int threads, double primalObjValue, double dualObjValue);
  
typedef enum {
               VERBOSE_NONE = 0,
               VERBOSE_ESSENTIAL = 1,
               VERBOSE_NON_ESSENTIAL = 2,
               VERBOSE_DEBUG = 3
             } VerbosityLevel;

extern VerbosityLevel g_verbosity;
  
extern lemon::Timer g_timer;
  
extern std::ostream* g_pOut;
  
} // namespace mwcs
} // namespace nina

#endif // UTILS_H
