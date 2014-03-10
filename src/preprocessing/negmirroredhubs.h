/*
 * negmirroredhubs.h
 *
 *  Created on: 8-mar-2014
 *      Author: M. El-Kebir
 */

#ifndef NEGMIRROREDHUBS_H
#define NEGMIRROREDHUBS_H

namespace nina {
namespace mwcs {

template<typename GR,
         typename WGHT = typename GR::template NodeMap<double> >
class NegMirroredHubs : public MwcsPreprocessRule<GR, WGHT>
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

  NegMirroredHubs();
  virtual ~NegMirroredHubs() {}
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

  virtual std::string name() const { return "NegMirroredHubs"; }
};

template<typename GR, typename WGHT>
inline NegMirroredHubs<GR, WGHT>::NegMirroredHubs()
  : Parent()
{
}

template<typename GR, typename WGHT>
inline int NegMirroredHubs<GR, WGHT>::apply(Graph& g,
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

  // TODO maintain neighborlist!
  NodeSet negHubsToRemove;
  for (size_t d = 3; d < degreeVector.size(); ++d)
  {
    const NodeSet& nodes = degreeVector[d];
    for (NodeSetIt nodeIt1 = nodes.begin(); nodeIt1 != nodes.end(); ++nodeIt1)
    {
      Node u = *nodeIt1;
      if (score[u] > 0) continue;
      
      NodeSet neighbors_u;
      for (IncEdgeIt e(g, u); e != lemon::INVALID; ++e)
      {
        neighbors_u.insert(g.oppositeNode(u, e));
      }

      for (NodeSetIt nodeIt2 = nodeIt1; nodeIt2 != nodes.end(); ++nodeIt2)
      {
        Node v = *nodeIt2;
        if (u == v) continue;
        if (score[v] > 0) continue;
//        if (neighbors_u.find(v) != neighbors_u.end()) continue;
        
        bool sameSets = true;
        for (IncEdgeIt e(g, v); e != lemon::INVALID && sameSets; ++e)
        {
          if (neighbors_u.find(g.oppositeNode(v, e)) == neighbors_u.end())
          {
            sameSets = false;
          }
        }
        
        if (sameSets)
        {
          // either u or v needs to go
          if (score[u] < score[v])
            negHubsToRemove.insert(u);
          else
            negHubsToRemove.insert(v);
        }
      }
    }
  }
  
  for (NodeSetIt nodeIt = negHubsToRemove.begin();
       nodeIt != negHubsToRemove.end(); ++nodeIt)
  {
    remove(g, mapToPre, preOrigNodes,
           nNodes, nArcs, nEdges,
           degree, degreeVector, *nodeIt);
//    std::cerr << "Removing neg hub " << label[*nodeIt] << " with deg " << degree[*nodeIt] << " and weight " << score[*nodeIt] << std::endl;
  }
  
//  int res = 0;
//  for (NodePairMapIt it = map.begin(); it != map.end(); ++it)
//  {
//    const WeightNodePairSet& set = it->second;
//    if (set.size() > 1)
//    {
//      for (WeightNodePairSetIt it = set.begin(), it_end = --set.end(); it != it_end; ++it)
//      {
//        remove(g, mapToPre, preOrigNodes,
//               nNodes, nArcs, nEdges,
//               degree, degreeVector, it->second);
////        std::cerr << it->first << "\t" << label[it->second] << std::endl;
//        ++res;
//      }
////      std::cerr << "(" << label[it->first.first] << "," << label[it->first.second] << "): " << it->second.size() << " " << i << std::endl;
//    }
//  }
  
//  std::cerr << "Done: " << negHubsToRemove.size() << std::endl;

  return static_cast<int>(negHubsToRemove.size());
}

} // namespace mwcs
} // namespace nina

#endif // NEGMIRROREDHUBS_H
