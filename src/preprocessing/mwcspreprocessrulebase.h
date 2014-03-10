/*
 * mwcspreprocessrulebase.h
 *
 *  Created on: 21-jan-2013
 *      Author: M. El-Kebir
 */

#ifndef MWCSPREPROCESSRULEBASE_H
#define MWCSPREPROCESSRULEBASE_H

#include <lemon/core.h>
#include <string>
#include <vector>
#include <set>

namespace nina {
namespace mwcs {

template<typename GR,
         typename WGHT = typename GR::template NodeMap<double> >
class MwcsPreprocessRuleBase
{
public:
  typedef GR Graph;
  typedef WGHT WeightNodeMap;

  TEMPLATE_GRAPH_TYPEDEFS(Graph);

  typedef typename Graph::template NodeMap<std::string> LabelNodeMap;
  typedef typename Graph::template NodeMap<int> DegreeNodeMap;
  typedef typename Graph::template NodeMap<Node> NodeMap;
  typedef std::set<Node> NodeSet;
  typedef typename NodeSet::iterator NodeSetIt;
  typedef typename Graph::template NodeMap<NodeSet> NodeSetMap;
  typedef typename std::vector<NodeSet> DegreeNodeSetVector;
  typedef lemon::DynArcLookUp<Graph> ArcLookUpType;

protected:
  void remove(Graph& g,
              NodeMap& mapToPre,
              NodeSetMap& preOrigNodes,
              int& nNodes,
              int& nArcs,
              int& nEdges,
              DegreeNodeMap& degree,
              DegreeNodeSetVector& degreeVector,
              Node node);
  Node merge(Graph& g,
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
             Node node1,
             Node node2,
             double& LB);

public:
  MwcsPreprocessRuleBase() {}
  virtual ~MwcsPreprocessRuleBase() {}
  virtual std::string name() const = 0;
};


template<typename GR, typename WGHT>
inline void MwcsPreprocessRuleBase<GR, WGHT>::remove(Graph& g,
                                                     NodeMap& mapToPre,
                                                     NodeSetMap& preOrigNodes,
                                                     int& nNodes,
                                                     int& nArcs,
                                                     int& nEdges,
                                                     DegreeNodeMap& degree,
                                                     DegreeNodeSetVector& degreeVector,
                                                     Node node)
{
  // decrease the degrees of adjacent nodes
  for (IncEdgeIt e(g, node); e != lemon::INVALID; ++e)
  {
    Node adjNode = g.oppositeNode(node, e);
    int d = degree[adjNode]--;
    degreeVector[d].erase(adjNode);
    degreeVector[d-1].insert(adjNode);

    nEdges--;
    // TODO double check this
    nArcs -= 2;
  }

  // remove the node from degree vector
  degreeVector[degree[node]].erase(node);

  // update mapToPre
  const NodeSet& nodes = preOrigNodes[node];
  for (NodeSetIt nodeIt = nodes.begin(); nodeIt != nodes.end(); nodeIt++)
  {
    // unmap
    mapToPre[*nodeIt] = lemon::INVALID;
  }

  // remove the node from the graph
  g.erase(node);
  nNodes--;
}

template<typename GR, typename WGHT>
inline typename MwcsPreprocessRuleBase<GR, WGHT>::Node
MwcsPreprocessRuleBase<GR, WGHT>::merge(Graph& g,
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
                                        Node node1,
                                        Node node2,
                                        double& LB)
{
  // TODO: prevent multiple edges from occurring
  if (node1 == node2)
    return node1;

  // find the node with minimum degree, this is the node to delete
  Node minNode, maxNode;
  if (degree[node1] < degree[node2])
  {
    minNode = node1;
    maxNode = node2;
  }
  else
  {
    minNode = node2;
    maxNode = node1;
  }

  // erase the degree of maxNode and minNode
  degreeVector[degree[minNode]].erase(minNode);
  degreeVector[degree[maxNode]].erase(maxNode);

  // now rewire the edges incident to minNode to maxNode
  degree[maxNode]--;
  for (IncEdgeIt e(g, minNode); e != lemon::INVALID;)
  {
    Node node = g.oppositeNode(minNode, e);
    if (node != maxNode)// && arcLookUp(maxNode, node) == lemon::INVALID)
    {

      if (arcLookUp(maxNode, node) == lemon::INVALID)
      {
        // introduce new edge between maxNode and node
        g.addEdge(maxNode, node);

        // update degree of maxNode
        degree[maxNode]++;
      }
      else
      {
        // adjust degree of node
        int d = degree[node];
        degreeVector[d].erase(node);
        degreeVector[d-1].insert(node);
        degree[node]--;

        nEdges--;
        nArcs -= 2;
      }

      // remove edge between minNode and node
      Edge toDelete = e;
      ++e;
      g.erase(toDelete);
    }
    else
    {
      ++e;
    }
  }

  // update degree of maxNode
  int d = degree[maxNode];
  if (degreeVector.size() <= static_cast<size_t>(d))
  {
    // add node sets to degreeVector
    int len = d - degreeVector.size() + 1;
    for (int i = 0; i < len; i++)
      degreeVector.push_back(NodeSet());
  }

  degreeVector[d].insert(maxNode);

  // update score of maxNode
  score[maxNode] += score[minNode];

  // update set of original nodes corresponding to maxNode
  const NodeSet& minOrgNodeSet = preOrigNodes[minNode];
  for (NodeSetIt nodeIt = minOrgNodeSet.begin(); nodeIt != minOrgNodeSet.end(); nodeIt++)
  {
    preOrigNodes[maxNode].insert(*nodeIt);
    mapToPre[*nodeIt] = maxNode;
  }

  // merge the labels
  label[maxNode] += "\t" + label[minNode];

  // erase minNode
  g.erase(minNode);
  nNodes--;
  nEdges--;
  nArcs -= 2;
  
  // update LB if necessary
  if (LB < score[maxNode])
  {
    LB = score[maxNode];
  }

  return maxNode;
}

} // namespace mwcs
} // namespace nina

#endif // MWCSPREPROCESSRULEBASE_H
