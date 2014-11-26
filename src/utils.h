#ifndef UTILS_H
#define UTILS_H

#include <lemon/list_graph.h>
#include <lemon/time_measure.h>
#include <ostream>
#include <set>

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

template<typename MWCSGR>
inline double reEvaluatePCST(const MWCSGR& mwcsGraph,
                             const std::set<Node>& solution)
{
  typedef std::set<Node> NodeSet;
  typedef NodeSet::const_iterator NodeSetIt;
  
  NodeSet solutionEdges;
  NodeSet solutionNodes;
  
  double cost = 0;
  for (NodeSetIt nodeIt = solution.begin(); nodeIt != solution.end(); ++nodeIt)
  {
    const NodeSet& orgNodes = mwcsGraph.getOrgNodes(*nodeIt);
    for (NodeSetIt orgNodeIt = orgNodes.begin(); orgNodeIt != orgNodes.end(); ++orgNodeIt)
    {
      Node orgNode = *orgNodeIt;
      const std::string& label = mwcsGraph.getOrgLabel(orgNode);
      int idU = -1, idV = -1;
      char c = '\0';
      if (sscanf(label.c_str(), "%d--%d%c", &idU, &idV, &c) == 2)
      {
        // edge
        solutionEdges.insert(orgNode);
        cost += -1 * mwcsGraph.getOrgScore(orgNode);
      }
      else if (sscanf(label.c_str(), "%d%c", &idU, &c) == 1)
      {
        // node
        solutionNodes.insert(*orgNodeIt);
      }
    }
  }
  
  const Graph& orgG = mwcsGraph.getOrgGraph();
  for (NodeIt v(orgG); v != lemon::INVALID; ++v)
  {
    const std::string& label = mwcsGraph.getOrgLabel(v);
    
    int idU = -1;
    char c = '\0';
    if (solutionEdges.find(v) == solutionEdges.end() &&
        solutionNodes.find(v) == solutionNodes.end() &&
        sscanf(label.c_str(), "%d%c", &idU, &c) == 1)
    {
      cost += mwcsGraph.getOrgScore(v);
    }
  }
  
  return cost;
}
  
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
