/*
 * mwcspreprocessedgraph.h
 *
 *  Created on: 11-jan-2013
 *      Author: M. El-Kebir
 */

#ifndef MWCSPREPROCESSEDGRAPH_H
#define MWCSPREPROCESSEDGRAPH_H

#include "mwcsgraphparser.h"
#include "preprocessing/rule.h"
#include <set>
#include <vector>
#include <algorithm>
#include <lemon/core.h>

#include "preprocessing/negdeg01.h"
#include "preprocessing/posedge.h"
#include "preprocessing/negedge.h"
#include "preprocessing/negcircuit.h"
#include "preprocessing/negdiamond.h"
#include "preprocessing/negmirroredhubs.h"
#include "preprocessing/negdominatedhubs.h"
#include "preprocessing/posdeg01.h"
#include "preprocessing/posdiamond.h"
#include "preprocessing/shortestpath.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class MwcsPreprocessedGraph : public MwcsGraphParser<GR, NWGHT, NLBL, EWGHT>
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;

  typedef MwcsGraphParser<GR, NWGHT, NLBL, EWGHT> Parent;
  typedef Rule<GR, NWGHT> RuleType;
  typedef typename Parent::ParserType ParserType;
  typedef typename Parent::InvLabelNodeMap InvLabelNodeMap;
  typedef typename Parent::InvLabelNodeMapIt InvLabelNodeMapIt;
  typedef typename RuleType::DegreeNodeMap DegreeNodeMap;
  typedef typename RuleType::DegreeNodeSetVector DegreeNodeSetVector;
  typedef typename RuleType::NodeMap NodeMap;
  typedef typename RuleType::NodeSet NodeSet;
  typedef typename RuleType::NodeSetIt NodeSetIt;
  typedef typename RuleType::NodeSetMap NodeSetMap;
  
  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  
  typedef std::set<Edge> EdgeSet;
  typedef typename EdgeSet::const_iterator EdgeSetIt;

  using Parent::getGraph;
  using Parent::getScores;
  using Parent::getOrgArcCount;
  using Parent::getOrgComponent;
  using Parent::getOrgComponentCount;
  using Parent::getOrgComponentMap;
  using Parent::getOrgEdgeCount;
  using Parent::getOrgGraph;
  using Parent::getOrgLabel;
  using Parent::getOrgLabels;
  using Parent::getOrgNodeByLabel;
  using Parent::getOrgNodeCount;
  using Parent::getOrgScore;
  using Parent::getOrgScores;
  using Parent::_parserInit;

private:
  typedef std::vector<RuleType*> RuleVector;
  typedef typename RuleVector::const_iterator RuleVectorIt;
  typedef typename RuleVector::iterator RuleVectorNonConstIt;
  typedef std::vector<RuleVector> RuleMatrix;

public:
  MwcsPreprocessedGraph();
  virtual ~MwcsPreprocessedGraph();
  virtual bool init(ParserType* pParser, bool pval);
  void preprocess(const NodeSet& rootNodes);
  void updateComponentMap()
  {
    _pGraph->_nComponents = lemon::connectedComponents(*_pGraph->_pG, *_pGraph->_pComp);
  }

protected:
  typedef NegDeg01<Graph> NegDeg01Type;
  typedef PosEdge<Graph> PosEdgeType;
  typedef NegEdge<Graph> NegEdgeType;
  typedef NegCircuit<Graph> NegCircuitType;
  typedef NegDiamond<Graph> NegDiamondType;
  typedef NegMirroredHubs<Graph> NegMirroredHubsType;
  typedef NegDominatedHubs<Graph> NegDominatedHubsType;
  typedef PosDeg01<Graph> PosDeg01Type;
  typedef PosDiamond<Graph> PosDiamondType;
  typedef ShortestPath<Graph> ShortestPathType;

private:
  typedef struct GraphStruct
  {
    Graph* _pG;
    LabelNodeMap* _pLabel;
    WeightNodeMap* _pScore;
    IntNodeMap* _pComp;
    NodeSetMap* _pPreOrigNodes;
    NodeSetMap* _pMapToPre;
    int _nNodes;
    int _nEdges;
    int _nArcs;
    int _nComponents;

    // constructor
    GraphStruct(const Graph& orgG)
      : _pG(new Graph())
      , _pLabel(new LabelNodeMap(*_pG))
      , _pScore(new WeightNodeMap(*_pG))
      , _pComp(new IntNodeMap(*_pG))
      , _pPreOrigNodes(new NodeSetMap(*_pG))
      , _pMapToPre(new NodeSetMap(orgG))
      , _nNodes(0)
      , _nEdges(0)
      , _nArcs(0)
      , _nComponents(0)
    {
    }

    // destructor
    ~GraphStruct()
    {
      delete _pPreOrigNodes;
      delete _pMapToPre;
      delete _pComp;
      delete _pScore;
      delete _pLabel;
      delete _pG;
    }
  } GraphStruct;

private:
  GraphStruct* _pGraph;
  GraphStruct* _pBackupGraph;
  RuleMatrix _rules;

protected:
  virtual void initParserMembers(Graph*& pG,
                                 LabelNodeMap*& pLabel,
                                 WeightNodeMap*& pScore,
                                 WeightNodeMap*& pPVal)
  {
    Parent::initParserMembers(pG, pLabel, pScore, pPVal);
    _pGraph = new GraphStruct(*pG);
  }

public:
  virtual const Graph& getGraph() const
  {
    return *_pGraph->_pG;
  }

  virtual Graph& getGraph()
  {
    return *_pGraph->_pG;
  }

  virtual NodeSet getOrgNodes(Node node) const
  {
    return (*_pGraph->_pPreOrigNodes)[node];
  }

  virtual NodeSet getOrgNodes(const NodeSet& nodes) const
  {
    NodeSet result;
    const NodeSetMap& preOrigNodes = *_pGraph->_pPreOrigNodes;

    for (NodeSetIt nodeIt = nodes.begin(); nodeIt != nodes.end(); nodeIt++)
    {
      result.insert(preOrigNodes[*nodeIt].begin(), preOrigNodes[*nodeIt].end());
    }

    return result;
  }

  virtual const LabelNodeMap& getLabels() const
  {
    return *_pGraph->_pLabel;
  }

  virtual LabelNodeMap& getLabels()
  {
    return *_pGraph->_pLabel;
  }

  virtual const WeightNodeMap& getScores() const
  {
    return *_pGraph->_pScore;
  }

  virtual WeightNodeMap& getScores()
  {
    return *_pGraph->_pScore;
  }

  virtual int getNodeCount() const
  {
    return _pGraph->_nNodes;
  }

  virtual int getEdgeCount() const
  {
    return _pGraph->_nEdges;
  }

  virtual int getArcCount() const
  {
    return _pGraph->_nArcs;
  }

  virtual int getComponentCount() const
  {
    return _pGraph->_nComponents;
  }

  virtual int getComponent(Node n) const
  {
    assert(n != lemon::INVALID);
    return (*_pGraph->_pComp)[n];
  }

  virtual const IntNodeMap& getComponentMap() const
  {
    return *_pGraph->_pComp;
  }

  virtual double getScore(Node n) const
  {
    assert(n != lemon::INVALID);
    return (*_pGraph->_pScore)[n];
  }

  virtual NodeSet getNodeByLabel(const std::string& label) const
  {
    Node orgNode = getOrgNodeByLabel(label);
    if (orgNode != lemon::INVALID)
    {
      return (*_pGraph->_pMapToPre)[orgNode];
    }
    
    return NodeSet();
  }

  void addPreprocessRule(int phase, RuleType* pRule)
  {
    while (static_cast<int>(_rules.size()) < phase)
      _rules.push_back(RuleVector());
    
    _rules[phase - 1].push_back(pRule);
  }

  virtual std::string getLabel(Node n) const
  {
    assert(n != lemon::INVALID);
    return (*_pGraph->_pLabel)[n];
  }

  virtual bool init(Graph* pG,
                    LabelNodeMap* pLabel,
                    WeightNodeMap* pScore,
                    WeightNodeMap* pPVal);

  virtual void computeScores(double lambda, double a, double FDR);

  virtual void computeScores(double tau);

  virtual NodeSet getPreNodes(Node orgNode) const
  {
    assert(_pGraph->_pMapToPre);
    return (*_pGraph->_pMapToPre)[orgNode];
  }
  
  virtual NodeSet getPreNodes(const NodeSet orgNodes) const
  {
    assert(_pGraph->_pMapToPre);
    NodeSet res;
    for (NodeSetIt nodeIt = orgNodes.begin(); nodeIt != orgNodes.end(); ++nodeIt)
    {
      const NodeSet& preNodes = (*_pGraph->_pMapToPre)[*nodeIt];
      res.insert(preNodes.begin(), preNodes.end());
    }
    return res;
  }
  
  
  void clear()
  {
    if (!_pGraph)
    {
      _pGraph = new GraphStruct(getOrgGraph());
    }
    
    _pGraph->_pG->clear();
    _pGraph->_nNodes = getOrgNodeCount();
    _pGraph->_nEdges = getOrgEdgeCount();
    _pGraph->_nArcs = getOrgArcCount();
    _pGraph->_nComponents = getOrgComponentCount();
    
    NodeMap nodeRef(getOrgGraph());

    lemon::graphCopy(getOrgGraph(), *_pGraph->_pG)
        .nodeMap(getOrgScores(), *_pGraph->_pScore)
        .nodeMap(getOrgLabels(), *_pGraph->_pLabel)
        .nodeMap(getOrgComponentMap(), *_pGraph->_pComp)
        .nodeRef(nodeRef)
        .run();

    for (NodeIt n(getOrgGraph()); n != lemon::INVALID; ++n)
    {
      Node preNode = nodeRef[n];
      (*_pGraph->_pPreOrigNodes)[preNode].clear();
      (*_pGraph->_pPreOrigNodes)[preNode].insert(n);
      (*_pGraph->_pMapToPre)[n].clear();
      (*_pGraph->_pMapToPre)[n].insert(preNode);
    }
  }
  
  void remove(Node node)
  {
    Graph& g = *_pGraph->_pG;
    
    // update edge and arc counts
    bool isolated = true;
    for (IncEdgeIt e(g, node); e != lemon::INVALID; ++e)
    {
      --_pGraph->_nEdges;
      _pGraph->_nArcs -= 2;
      isolated = false;
    }
    
    // update mapToPre
    const NodeSet& orgNodes = (*_pGraph->_pPreOrigNodes)[node];
    for (NodeSetIt orgNodeIt = orgNodes.begin(); orgNodeIt != orgNodes.end(); ++orgNodeIt)
    {
      // unmap
      Node orgNode = *orgNodeIt;
      (*_pGraph->_pMapToPre)[orgNode].erase(node);
    }
    
    g.erase(node);
    --_pGraph->_nNodes;
  }
  
  Edge addEdge(Node u, Node v)
  {
    assert(u != lemon::INVALID);
    assert(v != lemon::INVALID);
    
    Graph& g = *_pGraph->_pG;
    
    ++_pGraph->_nEdges;
    _pGraph->_nArcs += 2;
    
    return g.addEdge(u, v);
  }
  
  void remove(Edge e)
  {
    assert(e != lemon::INVALID);
    
    Graph& g = *_pGraph->_pG;
    
    _pGraph->_nEdges--;
    _pGraph->_nArcs -= 2;
    
    g.erase(e);
  }
  
  void remove(const NodeSet& nodes)
  {
    for (NodeSetIt nodeIt = nodes.begin(); nodeIt != nodes.end(); ++nodeIt)
    {
      remove(*nodeIt);
    }
  }
  
  Node merge(Node u, Node v)
  {
    Graph& g = *_pGraph->_pG;
    
    NodeSet neighborsU;
    for (IncEdgeIt e(g, u); e != lemon::INVALID; ++e)
    {
      Node node = g.oppositeNode(u, e);
      neighborsU.insert(node);
    }
    
    NodeSet neighborsV;
    for (IncEdgeIt e(g, v); e != lemon::INVALID; ++e)
    {
      Node node = g.oppositeNode(v, e);
      neighborsV.insert(node);
    }
    
    NodeSet intersection;
    std::set_intersection(neighborsU.begin(), neighborsU.end(),
                          neighborsV.begin(), neighborsV.end(),
                          std::inserter(intersection, intersection.begin()));
    
    // remove edges incident to minNode and intersection => prevent multiple edges
    for (IncEdgeIt e(g, u); e != lemon::INVALID;)
    {
      Node node = g.oppositeNode(u, e);      
      if (node == v || intersection.find(node) != intersection.end())
      {
        // remove edge
        Edge toDelete = e;
        ++e;
        g.erase(toDelete);
        
        // update edge and arc count
        --_pGraph->_nEdges;
        _pGraph->_nArcs -= 2;
      }
      else
      {
        ++e;
      }
    }
    
    // update score of v
    (*_pGraph->_pScore)[v] += (*_pGraph->_pScore)[u];
    
    // update set of original nodes corresponding to v
    const NodeSet& uOrgNodeSet = (*_pGraph->_pPreOrigNodes)[u];
    for (NodeSetIt nodeIt = uOrgNodeSet.begin(); nodeIt != uOrgNodeSet.end(); ++nodeIt)
    {
      (*_pGraph->_pPreOrigNodes)[v].insert(*nodeIt);
      (*_pGraph->_pMapToPre)[*nodeIt].erase(u);
      (*_pGraph->_pMapToPre)[*nodeIt].insert(v);
    }
    
    // merge the labels
    (*_pGraph->_pLabel)[v] += "\t" + (*_pGraph->_pLabel)[u];
    
    // erase minNode
    g.contract(v, u, true);
    --_pGraph->_nNodes;
    
    assert(lemon::simpleGraph(g));
    
    return v;
  }
  
  Node merge(NodeSet nodes)
  {
    if (nodes.empty())
    {
      return lemon::INVALID;
    }
    
    NodeSetIt nodeIt = nodes.begin();
    Node res = *nodeIt;
    for (++nodeIt; nodeIt != nodes.end(); ++nodeIt)
    {
      res = merge(res, *nodeIt);
    }
    
    return res;
  }
  
  Node merge(Node nodeToKeep, NodeSet nodes)
  {
    assert(nodes.find(nodeToKeep) != nodes.end());
    
    for (NodeSetIt nodeIt = nodes.begin(); nodeIt != nodes.end(); ++nodeIt)
    {
      if (*nodeIt != nodeToKeep)
      {
        Node res = merge(*nodeIt, nodeToKeep);
        assert(res == nodeToKeep);
      }
    }
    
    return nodeToKeep;
  }
  
  Node extract(NodeSet nodes)
  {
    Graph& g = *_pGraph->_pG;
    Node res = g.addNode();
    ++_pGraph->_nNodes;
    
    (*_pGraph->_pScore)[res] = 0;
    
    bool first = true;
    for (NodeSetIt nodeIt = nodes.begin(); nodeIt != nodes.end(); ++nodeIt)
    {
      Node v = *nodeIt;
      const NodeSet& orgNodes = (*_pGraph->_pPreOrigNodes)[v];
      
      (*_pGraph->_pScore)[res] += (*_pGraph->_pScore)[v];
      
      if (!first)
      {
        (*_pGraph->_pLabel)[res] += "\t";
      }
      else
      {
        first = false;
      }
      (*_pGraph->_pLabel)[res] += (*_pGraph->_pLabel)[v];
      (*_pGraph->_pPreOrigNodes)[res].insert(orgNodes.begin(), orgNodes.end());
      
      for (NodeSetIt orgNodeIt = orgNodes.begin(); orgNodeIt != orgNodes.end(); ++orgNodeIt)
      {
        (*_pGraph->_pMapToPre)[*orgNodeIt].insert(res);
      }
    }
    
    return res;
  }

protected:
  void constructDegreeMap(DegreeNodeMap& degree,
                          DegreeNodeSetVector& degreeVector) const;
  void constructNeighborMap(NodeSetMap& neighbors) const;
};

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline MwcsPreprocessedGraph<GR, NWGHT, NLBL, EWGHT>::MwcsPreprocessedGraph()
  : Parent()
  , _pGraph(NULL)
  , _pBackupGraph(NULL)
  , _rules()
{
  addPreprocessRule(1, new NegDeg01Type());
  addPreprocessRule(1, new PosEdgeType());
  addPreprocessRule(1, new NegEdgeType());
//  addPreprocessRule(1, new NegCircuitType());
//  addPreprocessRule(1, new NegDiamondType());
  addPreprocessRule(1, new PosDeg01Type());
  
  addPreprocessRule(2, new PosDiamondType());
  addPreprocessRule(2, new NegMirroredHubsType());
  //addPreprocessRule(2, new NegDominatedHubsType());
  
  addPreprocessRule(3, new ShortestPathType());
  
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline MwcsPreprocessedGraph<GR, NWGHT, NLBL, EWGHT>::~MwcsPreprocessedGraph()
{
  for (size_t i = 0; i < _rules.size(); ++i)
  {
    for (RuleVectorNonConstIt it = _rules[i].begin(); it != _rules[i].end(); it++)
    {
      delete *it;
    }
  }
  
  delete _pGraph;
  delete _pBackupGraph;
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsPreprocessedGraph<GR, NWGHT, NLBL, EWGHT>::preprocess(const NodeSet& rootNodes)
{
  DegreeNodeMap degree(*_pGraph->_pG);
  DegreeNodeSetVector degreeVector;
  NodeSetMap neighbors(*_pGraph->_pG);
  
  constructDegreeMap(degree, degreeVector);
  constructNeighborMap(neighbors);

  // determine max score
  double LB = std::max((*_pGraph->_pScore)[lemon::mapMax(*_pGraph->_pG, *_pGraph->_pScore)], 0.);
  
  // now let's preprocess the graph
  // in phases: first do phase 0 until no more change
  // then move on to phase 1 upon change fallback to phase 0
  //
  int uberTotRemovedNodes;
  do
  {
    uberTotRemovedNodes = 0;
    for (size_t phase = 0; phase < _rules.size(); ++phase)
    {
      int totRemovedNodes;
      do
      {
        totRemovedNodes = 0;
        for (RuleVectorIt ruleIt = _rules[phase].begin(); ruleIt != _rules[phase].end(); ruleIt++)
        {
          int removedNodes = (*ruleIt)->apply(*_pGraph->_pG, rootNodes,
                                              *_pGraph->_pLabel,
                                              *_pGraph->_pScore, *_pGraph->_pMapToPre,
                                              *_pGraph->_pPreOrigNodes, neighbors,
                                              _pGraph->_nNodes, _pGraph->_nArcs, _pGraph->_nEdges,
                                              degree, degreeVector, LB);
          
          assert(lemon::countNodes(*_pGraph->_pG) == _pGraph->_nNodes);
          assert(lemon::countEdges(*_pGraph->_pG) == _pGraph->_nEdges);
          
          totRemovedNodes += removedNodes;

          if (g_verbosity >= VERBOSE_DEBUG && removedNodes > 0)
          {
            std::cout << "// Phase " << phase + 1
                      << ": applied rule '" << (*ruleIt)->name()
                      << "' and removed " << removedNodes
                      << " node(s)" << std::endl;
          }
        }
        
        if (totRemovedNodes > 0)
        {
          phase = 0;
          uberTotRemovedNodes += totRemovedNodes;
        }
      } while (totRemovedNodes > 0);
    }
  } while (uberTotRemovedNodes > 0);

  // determine the connected components
  updateComponentMap();

  if (g_verbosity >= VERBOSE_ESSENTIAL)
  {
    std::cout << "// Preprocessing successfully applied"
              << ": " << _pGraph->_nNodes << " nodes, "
              << _pGraph->_nEdges << " edges and "
              << _pGraph->_nComponents << " component(s) remaining" << std::endl;
  }
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline bool MwcsPreprocessedGraph<GR, NWGHT, NLBL, EWGHT>::init(ParserType* pParser, bool pval)
{
  if (!Parent::init(pParser, pval))
    return false;

  // start by making a copy of the graph
  clear();

  return true;
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline bool MwcsPreprocessedGraph<GR, NWGHT, NLBL, EWGHT>::init(Graph* pG,
                                                                LabelNodeMap* pLabel,
                                                                WeightNodeMap* pScore,
                                                                WeightNodeMap* pPVal)
{
  if (!Parent::init(pG, pLabel, pScore, pPVal))
    return false;

  // start by making a copy of the graph
  clear();

  return true;
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsPreprocessedGraph<GR, NWGHT, NLBL, EWGHT>::constructNeighborMap(NodeSetMap& neighbors) const
{
  const Graph& g = *_pGraph->_pG;
  for (NodeIt n(g); n != lemon::INVALID; ++n)
  {
    NodeSet& neighborSet = neighbors[n];
    neighborSet.clear();
    for (IncEdgeIt e(g, n); e != lemon::INVALID; ++e)
    {
      neighborSet.insert(g.oppositeNode(n, e));
    }
  }
}
  
template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsPreprocessedGraph<GR, NWGHT, NLBL, EWGHT>::constructDegreeMap(
    DegreeNodeMap& degree,
    DegreeNodeSetVector& degreeVector) const
{
  for (NodeIt n(*_pGraph->_pG); n != lemon::INVALID; ++n)
  {
    int d = 0;
    for (IncEdgeIt e(*_pGraph->_pG, n); e != lemon::INVALID; ++e, d++) ;

    degree[n] = d;
    if (degreeVector.size() <= static_cast<size_t>(d))
    {
      // add node sets to degreeVector
      int len = d - degreeVector.size() + 1;
      for (int i = 0; i < len; i++)
        degreeVector.push_back(NodeSet());
    }

    degreeVector[d].insert(n);
  }
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsPreprocessedGraph<GR, NWGHT, NLBL, EWGHT>::computeScores(double lambda,
                                                                         double a,
                                                                         double FDR)
{
  Parent::computeScores(lambda, a, FDR);
  clear();
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsPreprocessedGraph<GR, NWGHT, NLBL, EWGHT>::computeScores(double tau)
{
  Parent::computeScores(tau);
  clear();
}

} // namespace mwcs
} // namespace nina

#endif // MWCSPREPROCESSEDGRAPH_H
