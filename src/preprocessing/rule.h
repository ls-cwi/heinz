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

public:
  Rule()
  {
  }
  
  virtual ~Rule()
  {
  }
  
  virtual int apply(Graph& g,
                    const NodeSet& rootNodes,
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
    
//    assert(isValid(g, mapToPre, preOrigNodes, neighbors, nNodes, nArcs, nEdges, degree, degreeVector));
  }
  
  Node extract(Graph& g,
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
               Node node)
  {
    Node newNode = g.addNode();
    
    label[newNode] = label[node];
    score[newNode] = score[node];
    
    const NodeSet& orgNodes = preOrigNodes[node];
    preOrigNodes[newNode].insert(orgNodes.begin(), orgNodes.end());
    
    comp[newNode] = nComponents++;
    degree[newNode] = 0;
    
    ++nNodes;
    
    for (NodeSetIt orgNodeIt = orgNodes.begin(); orgNodeIt != orgNodes.end(); ++orgNodeIt)
    {
      mapToPre[*orgNodeIt].insert(newNode);
    }

    degreeVector[0].insert(newNode);
    
//    assert(isValid(g, mapToPre, preOrigNodes, neighbors, nNodes, nArcs, nEdges, degree, degreeVector));
    return newNode;
  }
  
  Node merge(Graph& g,
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
    if (node1 == node2)
      return node1;
    
    // node1 is deleted, node2 is kept
    Node minNode = node1, maxNode = node2;
    
    // erase maxNode and minNode from degreeVector, we'll reinsert them later
    degreeVector[degree[minNode]].erase(minNode);
    degreeVector[degree[maxNode]].erase(maxNode);
    
    // now rewire the edges incident to minNode to maxNode
    NodeSet& maxNodeNeighbors = neighbors[maxNode];
    NodeSet& minNodeNeighbors = neighbors[minNode];
    
    minNodeNeighbors.erase(maxNode);
    maxNodeNeighbors.erase(minNode);
    
    NodeSet intersection;
    std::set_intersection(minNodeNeighbors.begin(), minNodeNeighbors.end(),
                          maxNodeNeighbors.begin(), maxNodeNeighbors.end(),
                          std::inserter(intersection, intersection.begin()));
    
    // remove edges incident to minNode and intersection => prevent multiple edges
    for (IncEdgeIt e(g, minNode); e != lemon::INVALID;)
    {
      Node node = g.oppositeNode(minNode, e);
      
      if (node == maxNode)
      {
        ++e;
        continue;
      }
      
      assert(neighbors[node].find(minNode) != neighbors[node].end());
      neighbors[node].erase(minNode);
      neighbors[node].insert(maxNode);
      
      if (intersection.find(node) != intersection.end())
      {
        // adjust degree of node
        int d = degree[node];
        degreeVector[d].erase(node);
        degreeVector[d-1].insert(node);
        --degree[node];
        
        // remove edge
        Edge toDelete = e;
        ++e;
        g.erase(toDelete);
        
        // update edge and arc count
        --nEdges;
        nArcs -= 2;
      }
      else
      {
        ++e;
      }
    }
    
    // now update the neighbor sets
    maxNodeNeighbors.insert(minNodeNeighbors.begin(), minNodeNeighbors.end());
    
    // update degree of maxNode
    int d = degree[maxNode] = static_cast<int>(maxNodeNeighbors.size());
    
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
      mapToPre[*nodeIt].erase(minNode);
      mapToPre[*nodeIt].insert(maxNode);
    }
    
    // merge the labels
    label[maxNode] += "\t" + label[minNode];
    
    // erase minNode
    g.contract(maxNode, minNode, true);
    nNodes--;
    nEdges--;
    nArcs -= 2;
    
    // update LB if necessary
    if (LB < score[maxNode])
    {
      LB = score[maxNode];
    }
    
#ifdef DEBUG
    int d2 = 0;
    for (IncEdgeIt e(g, maxNode); e != lemon::INVALID; ++e)
    {
      ++d2;
    }
    assert(neighbors[maxNode].size() == static_cast<size_t>(d2));
    assert(degree[maxNode] == static_cast<int>(neighbors[maxNode].size()));
    assert(degree[maxNode] >= 0);
    assert(lemon::simpleGraph(g));
//    assert(isValid(g, mapToPre, preOrigNodes, neighbors, nNodes, nArcs, nEdges, degree, degreeVector));
#endif
    
    return maxNode;
  }
  
  bool isValid(Graph& g,
               NodeSetMap& mapToPre,
               NodeSetMap& preOrigNodes,
               NodeSetMap& neighbors,
               int& nNodes,
               int& nArcs,
               int& nEdges,
               DegreeNodeMap& degree,
               DegreeNodeSetVector& degreeVector)
  {
    if (!lemon::simpleGraph(g))
      return false;
    
    IntNodeMap newDeg(g, 0);
    NodeSetMap newNeighbors(g);
    for (NodeIt v(g); v != lemon::INVALID; ++v)
    {
      for (IncEdgeIt e(g, v); e != lemon::INVALID; ++e)
      {
        ++newDeg[v];
        newNeighbors[v].insert(g.oppositeNode(v, e));
      }
    }
    
    for (NodeIt v(g); v != lemon::INVALID; ++v)
    {
      if (degree[v] != newDeg[v])
        return false;
      
      if (neighbors[v].size() != static_cast<size_t>(newDeg[v]))
        return false;
      
      if (degreeVector[degree[v]].find(v) == degreeVector[degree[v]].end())
        return false;
      
      if (newNeighbors[v] != neighbors[v])
        return false;
    }
    
    return true;
  }
};

} // namespace mwcs
} // namespace nina

#endif // RULE_H
