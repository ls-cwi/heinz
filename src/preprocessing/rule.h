/*
 * rule.h
 *
 *  Created on: 12-jan-2013
 *      Author: M. El-Kebir
 */

#ifndef RULE_H
#define RULE_H

#include <lemon/core.h>
#include <string>
#include <vector>
#include <set>

namespace nina {
namespace mwcs {

template<typename GR,
         typename WGHT = typename GR::template NodeMap<double> >
class Rule
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

public:
  Rule()
  {
  }
  
  virtual ~Rule()
  {
  }
  
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
                    double& LB) = 0;
  
  virtual std::string name() const = 0;
  
protected:
  void remove(Graph& g,
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
              Node node)
  {
    // decrease the degrees of adjacent nodes and update neighbors
    for (IncEdgeIt e(g, node); e != lemon::INVALID; ++e)
    {
      Node adjNode = g.oppositeNode(node, e);
      int d = degree[adjNode]--;
      degreeVector[d].erase(adjNode);
      degreeVector[d-1].insert(adjNode);
      neighbors[adjNode].erase(node);
      
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
      mapToPre[*nodeIt].erase(node);
    }
    
    // remove the node from the graph
    if (degree[node] == 0)
    {
      --nComponents;
    }
    
    g.erase(node);
    --nNodes;
  }
  
  Node extract(Graph& g,
               LabelNodeMap& label,
               WeightNodeMap& score,
               IntNodeMap& comp,
               NodeSetMap& mapToPre,
               NodeSetMap& preOrigNodes,
               int& nNodes,
               int& nArcs,
               int& nEdges,
               int& nComponents,
               DegreeNodeMap& degree,
               DegreeNodeSetVector& degreeVector,
               Node node)
  {
    const NodeSet& orgNodes = preOrigNodes[node];
    Node newNode = g.addNode();
    
    label[newNode] = label[node];
    score[newNode] = score[node];
    preOrigNodes[newNode] = orgNodes;
    
    comp[newNode] = nComponents++;
    degree[newNode] = 0;
    
    ++nNodes;
    
    for (NodeSetIt orgNodeIt = orgNodes.begin(); orgNodeIt != orgNodes.end(); ++orgNodeIt)
    {
      mapToPre[*orgNodeIt].insert(newNode);
    }

    degreeVector[0].insert(newNode);
    
    return newNode;
  }
  
  Node merge(Graph& g,
             const ArcLookUpType& arcLookUp,
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
      neighbors[node].erase(minNode);
      if (node != maxNode)// && arcLookUp(maxNode, node) == lemon::INVALID)
      {
        if (arcLookUp(maxNode, node) == lemon::INVALID)
        {
          neighbors[node].insert(maxNode);
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
      mapToPre[*nodeIt].erase(*nodeIt);
      mapToPre[*nodeIt].insert(maxNode);
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
};

} // namespace mwcs
} // namespace nina

#endif // RULE_H
