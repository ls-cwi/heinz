/*
 * mwcspreprocessedgraph.h
 *
 *  Created on: 11-jan-2013
 *      Author: M. El-Kebir
 */

#ifndef MWCSPREPROCESSEDGRAPH_H
#define MWCSPREPROCESSEDGRAPH_H

#include "mwcsgraphparser.h"
#include "preprocessing/unrootedrule.h"
#include "preprocessing/rootedrule.h"
#include <set>
#include <vector>
#include <algorithm>
#include <lemon/core.h>

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
  typedef UnrootedRule<GR, NWGHT> UnrootedRuleType;
  typedef RootedRule<GR, NWGHT> RootedRuleType;
  typedef typename Parent::ParserType ParserType;
  typedef typename Parent::InvLabelNodeMap InvLabelNodeMap;
  typedef typename Parent::InvLabelNodeMapIt InvLabelNodeMapIt;
  typedef typename UnrootedRuleType::DegreeNodeMap DegreeNodeMap;
  typedef typename UnrootedRuleType::DegreeNodeSetVector DegreeNodeSetVector;
  typedef typename UnrootedRuleType::NodeMap NodeMap;
  typedef typename UnrootedRuleType::NodeSet NodeSet;
  typedef typename UnrootedRuleType::NodeSetIt NodeSetIt;
  typedef typename UnrootedRuleType::NodeSetMap NodeSetMap;

  using Parent::getGraph;
  using Parent::getScores;
  using Parent::getOrgArcCount;
  using Parent::getOrgArcLookUp;
  using Parent::getOrgComponent;
  using Parent::getOrgComponentCount;
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
  TEMPLATE_GRAPH_TYPEDEFS(Graph);

  typedef typename Graph::Snapshot Snapshot;
  typedef std::vector<UnrootedRuleType*> UnrootedRuleVector;
  typedef typename UnrootedRuleVector::const_iterator UnrootedRuleVectorIt;
  typedef typename UnrootedRuleVector::iterator UnrootedRuleVectorNonConstIt;
  typedef std::vector<UnrootedRuleVector> UnrootedRuleMatrix;

  typedef std::vector<RootedRuleType*> RootedRuleVector;
  typedef typename RootedRuleVector::const_iterator RootedRuleVectorIt;
  typedef typename RootedRuleVector::iterator RootedRuleVectorNonConstIt;
  typedef std::vector<RootedRuleVector> RootedRuleMatrix;
public:
  MwcsPreprocessedGraph();
  virtual ~MwcsPreprocessedGraph();
  virtual bool init(ParserType* pParser, bool pval);
  virtual Node init(Node root);
  virtual bool deinit();

protected:
  typedef typename Parent::ArcLookUpType ArcLookUpType;

private:
  typedef struct GraphStruct
  {
    Graph* _pG;
    LabelNodeMap* _pLabel;
    WeightNodeMap* _pScore;
    IntNodeMap* _pComp;
    NodeSetMap* _pPreOrigNodes;
    int _nNodes;
    int _nEdges;
    int _nArcs;
    int _nComponents;
    ArcLookUpType* _pArcLookUp;

    // constructor
    GraphStruct()
      : _pG(new Graph())
      , _pLabel(new LabelNodeMap(*_pG))
      , _pScore(new WeightNodeMap(*_pG))
      , _pComp(new IntNodeMap(*_pG))
      , _pPreOrigNodes(new NodeSetMap(*_pG))
      , _nNodes(0)
      , _nEdges(0)
      , _nArcs(0)
      , _nComponents(0)
      , _pArcLookUp(new ArcLookUpType(*_pG))
    {
    }

    // destructor
    ~GraphStruct()
    {
      delete _pArcLookUp;
      delete _pPreOrigNodes;
      delete _pComp;
      delete _pScore;
      delete _pLabel;
      delete _pG;
    }
  } GraphStruct;

private:
  GraphStruct* _pGraph;
  GraphStruct* _pBackupGraph;
  NodeMap* _pMapToPre;
  UnrootedRuleMatrix _rules;
  RootedRuleMatrix _rootRules;

protected:
  virtual void initParserMembers(Graph*& pG,
                                 LabelNodeMap*& pLabel,
                                 WeightNodeMap*& pScore,
                                 WeightNodeMap*& pPVal)
  {
    Parent::initParserMembers(pG, pLabel, pScore, pPVal);
    _pMapToPre = new NodeMap(*pG);
  }

  void clear()
  {
    delete _pMapToPre;
    _pMapToPre = new NodeMap(getOrgGraph());

    _pGraph->_nNodes = getOrgNodeCount();
    _pGraph->_nEdges = getOrgEdgeCount();
    _pGraph->_nArcs = getOrgArcCount();

    lemon::graphCopy(getOrgGraph(), *_pGraph->_pG)
        .nodeMap(getOrgScores(), *_pGraph->_pScore)
        .nodeMap(getOrgLabels(), *_pGraph->_pLabel)
        .nodeRef(*_pMapToPre)
        .run();

    for (NodeIt n(getOrgGraph()); n != lemon::INVALID; ++n)
    {
      Node preNode = (*_pMapToPre)[n];
      (*_pGraph->_pPreOrigNodes)[preNode].clear();
      (*_pGraph->_pPreOrigNodes)[preNode].insert(n);
    }
  }

  void preprocess();

  virtual const ArcLookUpType& getArcLookUp() const { return *_pGraph->_pArcLookUp; }

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

  virtual Node getNodeByLabel(const std::string& label) const
  {
    Node node = Parent::getNodeByLabel(label);

    if (node == lemon::INVALID)
      return node;
    else
      return (*_pMapToPre)[node];
  }

  void addPreprocessRootRule(int phase, RootedRuleType* pRule)
  {
    while (static_cast<int>(_rootRules.size()) < phase)
      _rootRules.push_back(RootedRuleVector());
    
    _rootRules[phase - 1].push_back(pRule);
  }

  void addPreprocessRule(int phase, UnrootedRuleType* pRule)
  {
    while (static_cast<int>(_rules.size()) < phase)
      _rules.push_back(UnrootedRuleVector());
    
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

  virtual Node getPreNode(Node orgNode) const
  {
    assert(_pMapToPre);
    return (*_pMapToPre)[orgNode];
  }

protected:
  void constructDegreeMap(DegreeNodeMap& degree,
                          DegreeNodeSetVector& degreeVector) const;
  void constructNeighborMap(NodeSetMap& neighbors) const;
  void backup();
  void restore();
};

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline MwcsPreprocessedGraph<GR, NWGHT, NLBL, EWGHT>::MwcsPreprocessedGraph()
  : Parent()
  , _pGraph(new GraphStruct())
  , _pBackupGraph(NULL)
  , _pMapToPre(NULL)
  , _rules()
  , _rootRules()
{
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline MwcsPreprocessedGraph<GR, NWGHT, NLBL, EWGHT>::~MwcsPreprocessedGraph()
{
  for (size_t i = 0; i < _rules.size(); ++i)
  {
    for (UnrootedRuleVectorNonConstIt it = _rules[i].begin(); it != _rules[i].end(); it++)
    {
      delete *it;
    }
  }

  for (size_t i = 0; i < _rootRules.size(); ++i)
  {
    for (RootedRuleVectorNonConstIt it = _rootRules[i].begin(); it != _rootRules[i].end(); it++)
    {
      delete *it;
    }
  }
  
  delete _pMapToPre;
  delete _pGraph;
  delete _pBackupGraph;
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsPreprocessedGraph<GR, NWGHT, NLBL, EWGHT>::backup()
{
  delete _pBackupGraph;
  _pBackupGraph = new GraphStruct();

  lemon::graphCopy(*_pGraph->_pG, *_pBackupGraph->_pG)
      .nodeMap(*_pGraph->_pLabel, *_pBackupGraph->_pLabel)
      .nodeMap(*_pGraph->_pScore, *_pBackupGraph->_pScore)
      .nodeMap(*_pGraph->_pComp, *_pBackupGraph->_pComp)
      .nodeMap(*_pGraph->_pPreOrigNodes, *_pBackupGraph->_pPreOrigNodes)
      .run();

  _pBackupGraph->_nNodes = _pGraph->_nNodes;
  _pBackupGraph->_nEdges = _pGraph->_nEdges;
  _pBackupGraph->_nArcs = _pGraph->_nArcs;
  _pBackupGraph->_nComponents = _pGraph->_nComponents;

  for (NodeIt preNode(*_pGraph->_pG); preNode != lemon::INVALID; ++preNode)
  {
    const NodeSet& orgNodes = (*_pGraph->_pPreOrigNodes)[preNode];

    for (NodeSetIt orgNodeIt = orgNodes.begin();
         orgNodeIt != orgNodes.end(); orgNodeIt++)
    {
      (*_pMapToPre)[*orgNodeIt] = preNode;
    }
  }
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsPreprocessedGraph<GR, NWGHT, NLBL, EWGHT>::restore()
{
  delete _pGraph;
  _pGraph = _pBackupGraph;
  _pBackupGraph = NULL;

  // recompute _pMapToPre
  assert(_pMapToPre);

  for (NodeIt preNode(*_pGraph->_pG); preNode != lemon::INVALID; ++preNode)
  {
    const NodeSet& orgNodes = (*_pGraph->_pPreOrigNodes)[preNode];

    for (NodeSetIt orgNodeIt = orgNodes.begin();
         orgNodeIt != orgNodes.end(); orgNodeIt++)
    {
      (*_pMapToPre)[*orgNodeIt] = preNode;
    }
  }
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsPreprocessedGraph<GR, NWGHT, NLBL, EWGHT>::preprocess()
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
        for (UnrootedRuleVectorIt ruleIt = _rules[phase].begin(); ruleIt != _rules[phase].end(); ruleIt++)
        {
          int removedNodes = (*ruleIt)->apply(*_pGraph->_pG, getArcLookUp(), *_pGraph->_pLabel,
                                              *_pGraph->_pScore, (*_pMapToPre),
                                              *_pGraph->_pPreOrigNodes, neighbors,
                                              _pGraph->_nNodes, _pGraph->_nArcs,
                                              _pGraph->_nEdges, degree, degreeVector, LB);
          totRemovedNodes += removedNodes;

          if (g_verbosity >= VERBOSE_DEBUG && removedNodes > 0)
          {
            std::cout << "// Phase " << phase
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
  _pGraph->_nComponents = lemon::connectedComponents(*_pGraph->_pG, *_pGraph->_pComp);

  if (g_verbosity >= VERBOSE_ESSENTIAL)
  {
    std::cout << "// Preprocessing successfully applied (stage 1)"
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

  if (!pval)
    preprocess();

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

  preprocess();

  return true;
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline bool MwcsPreprocessedGraph<GR, NWGHT, NLBL, EWGHT>::deinit()
{
  restore();
  return true;
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline typename MwcsPreprocessedGraph<GR, NWGHT, NLBL, EWGHT>::Node
MwcsPreprocessedGraph<GR, NWGHT, NLBL, EWGHT>::init(Node root)
{
  backup();

  DegreeNodeMap degree(*_pGraph->_pG);
  DegreeNodeSetVector degreeVector;
  NodeSetMap neighbors(*_pGraph->_pG);

  constructDegreeMap(degree, degreeVector);
  constructNeighborMap(neighbors);

  // now let's preprocess the graph
  int totRemovedNodes;
  do
  {
    totRemovedNodes = 0;
    for (size_t phase = 0; phase < _rules.size(); ++phase)
    {
      for (UnrootedRuleVectorIt ruleIt = _rules[phase].begin();
           ruleIt != _rules[phase].end(); ruleIt++)
      {
        int removedNodes = (*ruleIt)->apply(*_pGraph->_pG, getArcLookUp(),
                                            *_pGraph->_pLabel, *_pGraph->_pScore, *_pMapToPre,
                                            *_pGraph->_pPreOrigNodes, neighbors,
                                            _pGraph->_nNodes, _pGraph->_nArcs,
                                            _pGraph->_nEdges, degree, degreeVector, (*_pGraph->_pScore)[root]);
        totRemovedNodes += removedNodes;

        if (g_verbosity >= VERBOSE_ESSENTIAL && removedNodes > 0)
        {
          std::cout << "// Applied rule '" << (*ruleIt)->name()
                    << "' and removed " << removedNodes
                    << " node(s)" << std::endl;
        }
      }
    }

    for (size_t phase = 0; phase < _rootRules.size(); ++phase)
    {
      for (RootedRuleVectorIt ruleIt = _rootRules[phase].begin();
           ruleIt != _rootRules[phase].end(); ruleIt++)
      {
        int removedNodes = (*ruleIt)->apply(*_pGraph->_pG, root, getArcLookUp(),
                                            *_pGraph->_pLabel, *_pGraph->_pScore, *_pMapToPre,
                                            *_pGraph->_pPreOrigNodes, neighbors,
                                            _pGraph->_nNodes, _pGraph->_nArcs,
                                            _pGraph->_nEdges, degree, degreeVector);
        totRemovedNodes += removedNodes;

        if (g_verbosity >= VERBOSE_ESSENTIAL && removedNodes > 0)
        {
          std::cout << "// Applied rule '" << (*ruleIt)->name()
                    << "' and removed " << removedNodes
                    << " node(s)" << std::endl;
        }
      }
    }
  } while (totRemovedNodes > 0);

  // determine the connected components
  _pGraph->_nComponents = lemon::connectedComponents(*_pGraph->_pG, *_pGraph->_pComp);

  if (g_verbosity >= VERBOSE_ESSENTIAL)
  {
    std::cout << "// Succesfully preprocessed (stage 2)"
              << ": contains " << _pGraph->_nNodes << " nodes, "
              << _pGraph->_nEdges << " edges and "
              << _pGraph->_nComponents << " component(s)" << std::endl;
  }

  return (*_pMapToPre)[root];
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
  preprocess();
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsPreprocessedGraph<GR, NWGHT, NLBL, EWGHT>::computeScores(double tau)
{
  Parent::computeScores(tau);
  clear();
  preprocess();
}

} // namespace mwcs
} // namespace nina

#endif // MWCSPREPROCESSEDGRAPH_H
