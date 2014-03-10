/*
 * negdiamond.h
 *
 *  Created on: 7-mar-2014
 *      Author: M. El-Kebir
 */

#ifndef NEGDIAMOND_H
#define NEGDIAMOND_H

#include <lemon/core.h>
#include <string>
#include <vector>
#include <set>
#include "mwcspreprocessrule.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename WGHT = typename GR::template NodeMap<double> >
class NegDiamond : public MwcsPreprocessRule<GR, WGHT>
{
public:
  typedef GR Graph;
  typedef WGHT WeightNodeMap;
  typedef MwcsPreprocessRule<GR, WGHT> Parent;
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

  NegDiamond();
  virtual ~NegDiamond() {}
  virtual int apply(Graph& g,
                    const ArcLookUpType& arcLookUp,
                    LabelNodeMap& label,
                    WeightNodeMap& score,
                    NodeMap& mapToPre,
                    NodeSetMap& preOrigNodes,
                    int& nNodes,
                    int& nArcs,
                    int& nEdges,
                    DegreeNodeMap& degree,
                    DegreeNodeSetVector& degreeVector,
                    double& LB);

  virtual std::string name() const { return "NegDiamond"; }
};

template<typename GR, typename WGHT>
inline NegDiamond<GR, WGHT>::NegDiamond()
  : Parent()
{
}

template<typename GR, typename WGHT>
inline int NegDiamond<GR, WGHT>::apply(Graph& g,
                                       const ArcLookUpType& arcLookUp,
                                       LabelNodeMap& label,
                                       WeightNodeMap& score,
                                       NodeMap& mapToPre,
                                       NodeSetMap& preOrigNodes,
                                       int& nNodes,
                                       int& nArcs,
                                       int& nEdges,
                                       DegreeNodeMap& degree,
                                       DegreeNodeSetVector& degreeVector,
                                       double& LB)
{
  typedef std::pair<double, Node> WeightNodePair;
  typedef std::set<WeightNodePair> WeightNodePairSet;
  typedef typename WeightNodePairSet::const_iterator WeightNodePairSetIt;
  
  typedef std::pair<Node, Node> NodePair;
  typedef std::map<NodePair, WeightNodePairSet> NodePairMap;
  typedef typename NodePairMap::const_iterator NodePairMapIt;

  NodePairMap map;
  NodePairMap posMap;

  const NodeSet& nodes = degreeVector[2];
  for (NodeSetIt nodeIt = nodes.begin(); nodeIt != nodes.end(); ++nodeIt)
  {
    Node v = *nodeIt;
    
    assert(degree[v] == 2);
    Edge e1 = IncEdgeIt(g, v);
    Edge e2 = ++IncEdgeIt(g, v);
    
    Node u = g.oppositeNode(v, e1);
    Node w = g.oppositeNode(v, e2);
    if (score[v] <= 0)
    {
      if (u < w)
      {
        map[std::make_pair(u, w)].insert(std::make_pair(score[v], v));
      }
      else
      {
        map[std::make_pair(w, u)].insert(std::make_pair(score[v], v));
      }
    }
    else
    {
      if (u < w)
      {
        posMap[std::make_pair(u, w)].insert(std::make_pair(score[v], v));
      }
      else
      {
        posMap[std::make_pair(w, u)].insert(std::make_pair(score[v], v));
      }
    }
  }
  
  int res = 0;
  for (NodePairMapIt it = map.begin(); it != map.end(); ++it)
  {
    const WeightNodePairSet& set = it->second;
//    if (set.size() > 1)
    {
      // if there is a positive center node in the diamond remove all the negative center nodes
      WeightNodePairSetIt it_end = posMap.find(it->first) != posMap.end() ? set.end() : --set.end();
      for (WeightNodePairSetIt it = set.begin(); it != it_end; ++it)
      {
        remove(g, mapToPre, preOrigNodes,
               nNodes, nArcs, nEdges,
               degree, degreeVector, it->second);
//        std::cerr << it->first << "\t" << label[it->second] << std::endl;
        ++res;
      }
//      std::cerr << "(" << label[it->first.first] << "," << label[it->first.second] << "): " << it->second.size() << " " << i << std::endl;
    }
  }

  return res;
}

} // namespace mwcs
} // namespace nina

#endif // NEGDIAMOND_H
