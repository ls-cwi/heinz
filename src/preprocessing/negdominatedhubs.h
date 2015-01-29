/*
 * negdominatedhubs.h
 *
 *  Created on: 24-jan-2015
 *      Author: G. Klau
 */

#ifndef NEGDOMINATEDHUBS_H
#define NEGDOMINATEDHUBS_H

#include "rule.h"

namespace nina {
  namespace mwcs {
    
    template<typename GR,
    typename WGHT = typename GR::template NodeMap<double> >
    class NegDominatedHubs : public Rule<GR, WGHT>
    {
    public:
      typedef GR Graph;
      typedef WGHT WeightNodeMap;
      typedef Rule<GR, WGHT> Parent;
      typedef typename Parent::NodeMap NodeMap;
      typedef typename Parent::NodeSet NodeSet;
      typedef typename Parent::NodeSetIt NodeSetIt;
      typedef typename Parent::NodeSetMap NodeSetMap;
      typedef typename Parent::DegreeNodeMap DegreeNodeMap;
      typedef typename Parent::DegreeNodeSetVector DegreeNodeSetVector;
      typedef typename Parent::LabelNodeMap LabelNodeMap;
      
      TEMPLATE_GRAPH_TYPEDEFS(Graph);
      
      using Parent::remove;
      using Parent::merge;
      
      NegDominatedHubs();
      virtual ~NegDominatedHubs() {}
      virtual int apply(Graph& g,
                        const NodeSet& rootNodes,
                        LabelNodeMap& label,
                        WeightNodeMap& score,
                        NodeSetMap& mapToPre,
                        NodeSetMap& preOrigNodes,
                        NodeSetMap& neighbors,
                        int& nNodes,
                        int& nArcs,
                        int& nEdges,
                        DegreeNodeMap& degree,
                        DegreeNodeSetVector& degreeVector,
                        double& LB);
      
      virtual std::string name() const { return "NegDominatedHubs"; }
    };
    
    template<typename GR, typename WGHT>
    inline NegDominatedHubs<GR, WGHT>::NegDominatedHubs()
    : Parent()
    {
    }
    
    template<typename GR, typename WGHT>
    inline int NegDominatedHubs<GR, WGHT>::apply(Graph& g,
                                                  const NodeSet& rootNodes,
                                                  LabelNodeMap& label,
                                                  WeightNodeMap& score,
                                                  NodeSetMap& mapToPre,
                                                  NodeSetMap& preOrigNodes,
                                                  NodeSetMap& neighbors,
                                                  int& nNodes,
                                                  int& nArcs,
                                                  int& nEdges,
                                                  DegreeNodeMap& degree,
                                                  DegreeNodeSetVector& degreeVector,
                                                  double& LB)
    {
      NodeSet negHubsToRemove;
      for (NodeIt u(g); u != lemon::INVALID; ++u)
      {
        if (score[u] > 0) continue;
        for (NodeIt v(g); v != lemon::INVALID; ++v)
        {
          if (u == v) continue;
          //if (score[v] >= 0) continue; // not really needed, but the whole thing was too slow. did not help...
          if (score[u] >= score[v]) continue;
          if (degree[u] > degree[v]) continue;
          
          // now check subset:
          const NodeSet& neighbors_u = neighbors[u], neighbors_v = neighbors[v];
          if (std::includes(neighbors_v.begin(), neighbors_v.end(),
                            neighbors_u.begin(), neighbors_u.end()))
            negHubsToRemove.insert(u);
        }
      }
      
      for (NodeSetIt nodeIt = rootNodes.begin(); nodeIt != rootNodes.end(); ++nodeIt)
      {
        // only remove if not a root node
        negHubsToRemove.erase(*nodeIt);
      }
      
      for (NodeSetIt nodeIt = negHubsToRemove.begin(); nodeIt != negHubsToRemove.end(); ++nodeIt)
      {
        Node v = *nodeIt;
        
        assert(rootNodes.find(v) == rootNodes.end());
        remove(g, mapToPre, preOrigNodes, neighbors,
               nNodes, nArcs, nEdges,
               degree, degreeVector, v);
      }
      
      return static_cast<int>(negHubsToRemove.size());
    }
    
  } // namespace mwcs
} // namespace nina

#endif // NEGDOMINATEDHUBS_H
