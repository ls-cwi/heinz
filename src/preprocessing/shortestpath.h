/*
 * shortestpath.h
 *
 *  Created on: 2-sep-2014
 *      Author: M. El-Kebir
 */

#ifndef SHORTESTPATH_H
#define SHORTESTPATH_H

#include <lemon/core.h>
#include <lemon/adaptors.h>
#include <lemon/dijkstra.h>
#include <string>
#include <vector>
#include <set>
#include <list>
#include "rule.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename WGHT = typename GR::template NodeMap<double> >
class ShortestPath : public Rule<GR, WGHT>
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
  typedef typename Parent::ArcLookUpType ArcLookUpType;
  
  TEMPLATE_GRAPH_TYPEDEFS(Graph);

  using Parent::remove;
  using Parent::merge;

  ShortestPath();
  virtual ~ShortestPath() {}
  virtual int apply(Graph& g,
                    const NodeSet& rootNodes,
                    const ArcLookUpType& arcLookUp,
                    LabelNodeMap& label,
                    WeightNodeMap& score,
                    IntNodeMap& comp,
                    NodeSetMap& mapToPre,
                    NodeSetMap& preOrigNodes,
                    NodeSetMap& neighbors,
                    int& nNodes,
                    int& nArcs,
                    int& nEdges,
                    int& nComponents,
                    DegreeNodeMap& degree,
                    DegreeNodeSetVector& degreeVector,
                    double& LB);

  virtual std::string name() const { return "ShortestPath"; }
  
private:
  typedef lemon::FilterNodes<const Graph, BoolNodeMap> SubGraph;
  typedef typename SubGraph::Node SubNode;
  typedef std::list<Node> NodeList;
  typedef typename NodeList::const_iterator NodeListIt;
  
  double shortCircuit(SubGraph& g,
                      const DoubleArcMap& arcCost,
                      const WeightNodeMap& score,
                      Node u,
                      Node v,
                      Node w,
                      NodeList& path);
};

template<typename GR, typename WGHT>
inline ShortestPath<GR, WGHT>::ShortestPath()
  : Parent()
{
}
  
template<typename GR, typename WGHT>
inline double ShortestPath<GR, WGHT>::shortCircuit(SubGraph& g,
                                                   const DoubleArcMap& arcCost,
                                                   const WeightNodeMap& score,
                                                   Node u,
                                                   Node v,
                                                   Node w,
                                                   NodeList& path)
{
  g.disable(v);
  
  double pathLength = 0;
  
  // and now let's do a Dijkstra
  lemon::Dijkstra<SubGraph, DoubleArcMap> dijkstra(g, arcCost);
  dijkstra.run(u);
  
  if (!dijkstra.reached(w))
  {
    pathLength = -std::numeric_limits<double>::max();
  }
  else
  {
    SubNode subNode = w;
    while (dijkstra.predArc(subNode) != lemon::INVALID)
    {
      subNode = dijkstra.predNode(subNode);
      if (subNode != u)
      {
        path.push_front(subNode);
        if (score[subNode] < 0)
          pathLength += score[subNode];
      }
    }
  }
  
  g.enable(v);
  
  return pathLength;
}
  
template<typename GR, typename WGHT>
inline int ShortestPath<GR, WGHT>::apply(Graph& g,
                                         const NodeSet& rootNodes,
                                         const ArcLookUpType& arcLookUp,
                                         LabelNodeMap& label,
                                         WeightNodeMap& score,
                                         IntNodeMap& comp,
                                         NodeSetMap& mapToPre,
                                         NodeSetMap& preOrigNodes,
                                         NodeSetMap& neighbors,
                                         int& nNodes,
                                         int& nArcs,
                                         int& nEdges,
                                         int& nComponents,
                                         DegreeNodeMap& degree,
                                         DegreeNodeSetVector& degreeVector,
                                         double& LB)
{
  int res = 0;
  
  if (degreeVector.size() <= 2)
  {
    // nothing to remove, there are no degree 2 nodes
    return 0;
  }
  
  // let's construct the arc costs
  DoubleArcMap arcCost(g);
  for (ArcIt a(g); a != lemon::INVALID; ++a)
  {
    arcCost[a] = score[g.target(a)] > 0 ? 0 : -score[g.target(a)];
    assert(arcCost[a] >= 0);
  }

  BoolNodeMap filter(g, true);
  SubGraph subG(g, filter);
  
  for (NodeIt v(g); v != lemon::INVALID;)
  {
    if (degree[v] == 2 && score[v] <= 0 && rootNodes.find(v) == rootNodes.end())
    {
      Edge e1 = IncEdgeIt(g, v);
      Edge e2 = ++IncEdgeIt(g, v);
      
      Node u = g.oppositeNode(v, e1);
      Node w = g.oppositeNode(v, e2);
      
      NodeList path;
      double pathLength = shortCircuit(subG, arcCost, score, u, v, w, path);
      
      if (pathLength >= score[v])
      {
//          std::cout << label[v] << " ( " << score[v] << " ) : "
//                    << label[u] << " (" << score[u] << ")";
//          for (NodeListIt nodeIt = path.begin(); nodeIt != path.end(); ++nodeIt)
//          {
//            std::cout << " -- " << label[*nodeIt] << " (" << score[*nodeIt] << ")";
//          }
//          std::cout << " -- " << label[w] << " (" << score[w] << ")" << std::endl;
        
        Node curV = v;
        ++v;
        remove(g, comp, mapToPre, preOrigNodes, neighbors,
               nNodes, nArcs, nEdges, nComponents,
               degree, degreeVector, curV);
        ++res;
      }
      else
      {
        ++v;
      }
    }
    else
    {
      ++v;
    }
  }
  
  return res;
}
  
} // namespace mwcs
} // namespace nina

#endif // SHORTESTPATH_H