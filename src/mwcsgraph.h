/*
 * mwcsgraph.h
 *
 *  Created on: 6-aug-2012
 *     Authors: C.I. Bucur and M. El-Kebir
 */

#ifndef MWCSGRAPH_H
#define MWCSGRAPH_H

#include <map>
#include <assert.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <lemon/core.h>
#include <lemon/lgf_writer.h>
#include <lemon/connectivity.h>
#include "verbose.h"
#include "parser/parser.h"
#include "solver/spqrtree.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class MwcsGraph
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;

  TEMPLATE_GRAPH_TYPEDEFS(Graph);

public:
  /// Parser type
  typedef Parser<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> ParserType;
  typedef typename ParserType::InvIdNodeMap InvLabelNodeMap;
  typedef typename InvLabelNodeMap::const_iterator InvLabelNodeMapIt;

  typedef typename std::set<Node> NodeSet;
  typedef typename NodeSet::const_iterator NodeSetIt;
  typedef typename std::vector<NodeSet> NodeSetVector;
  typedef typename NodeSetVector::const_iterator NodeSetVectorIt;
  typedef lemon::DynArcLookUp<Graph> ArcLookUpType;
  
public:
  MwcsGraph();
  virtual ~MwcsGraph()
  {
    delete _pArcLookUp;
    delete _pComp;

    if (_parserInit)
    {
      delete _pLabel;
      delete _pPVal;
      delete _pScore;
      delete _pG;
    }
  }
  virtual bool init(ParserType* pParser, bool pval);
  virtual bool init(Graph* pG,
                    LabelNodeMap* pLabel,
                    WeightNodeMap* pScore,
                    WeightNodeMap* pPval);
  virtual Node init(Node root) { return root; }
  virtual bool deinit() { return false; }
  void resetCounts();



private:
  Graph* _pG;
  LabelNodeMap* _pLabel;
  WeightNodeMap* _pPVal;
  WeightNodeMap* _pScore;
  InvLabelNodeMap _invLabel;
  IntNodeMap* _pComp;
  int _nNodes;
  int _nEdges;
  int _nArcs;
  int _nComponents;
  ArcLookUpType* _pArcLookUp;

protected:
  bool _parserInit;

protected:
  virtual void initParserMembers(Graph*& pG,
                                 LabelNodeMap*& pLabel,
                                 WeightNodeMap*& pScore,
                                 WeightNodeMap*& pPVal)
  {
  }

public:
  const ArcLookUpType& getOrgArcLookUp() const { return *_pArcLookUp; }
  
  virtual const ArcLookUpType& getArcLookUp() const { return *_pArcLookUp; }

  virtual const Graph& getGraph() const
  {
    return getOrgGraph();
  }

  virtual Graph& getGraph()
  {
    return getOrgGraph();
  }

  virtual const LabelNodeMap& getLabels() const
  {
    return getOrgLabels();
  }

  virtual LabelNodeMap& getLabels()
  {
    return getOrgLabels();
  }

  virtual NodeSet getOrgNodes(Node node) const
  {
    NodeSet result;
    result.insert(node);
    return result;
  }

  virtual NodeSet getOrgNodes(const NodeSet& nodes) const
  {
    NodeSet result;
    for (NodeSetIt nodeIt = nodes.begin(); nodeIt != nodes.end(); nodeIt++)
    {
      NodeSet tmp = getOrgNodes(*nodeIt);
      result.insert(tmp.begin(), tmp.end());
    }
    return result;
  }

  virtual const WeightNodeMap& getScores() const
  {
    return getOrgScores();
  }

  virtual WeightNodeMap& getScores()
  {
    return getOrgScores();
  }

  virtual int getNodeCount() const
  {
    return getOrgNodeCount();
  }

  virtual int getEdgeCount() const
  {
    return getOrgEdgeCount();
  }

  virtual int getArcCount() const
  {
    return getOrgArcCount();
  }

  virtual const IntNodeMap& getComponentMap() const
  {
    return *_pComp;
  }

  virtual int getComponentCount() const
  {
    return getOrgComponentCount();
  }

  virtual int getComponent(Node n) const
  {
    return getOrgComponent(n);
  }

  virtual std::string getLabel(Node n) const
  {
    return getOrgLabel(n);
  }

  virtual double getScore(Node n) const
  {
    return getOrgScore(n);
  }

  virtual Node getNodeByLabel(const std::string& label) const
  {
    return getOrgNodeByLabel(label);
  }

  void writeLGF(std::ostream& out) const
  {
    graphWriter(_pG, out)
      .nodeMap("pval", *_pPVal)
      .nodeMap("score", *_pScore)
      .nodeMap("id", _pLabel)
      .run();
  }

  const Graph& getOrgGraph() const
  {
    return *_pG;
  }

  Graph& getOrgGraph()
  {
    return *_pG;
  }

  const LabelNodeMap& getOrgLabels() const
  {
    return *_pLabel;
  }

  LabelNodeMap& getOrgLabels()
  {
    return *_pLabel;
  }

  const WeightNodeMap& getOrgScores() const
  {
    return *_pScore;
  }

  WeightNodeMap& getOrgScores()
  {
    return *_pScore;
  }

  WeightNodeMap* getOrgPValues()
  {
    return _pPVal;
  }

  const WeightNodeMap* getOrgPValues() const
  {
    return _pPVal;
  }

  int getOrgNodeCount() const
  {
    return _nNodes;
  }

  int getOrgEdgeCount() const
  {
    return _nEdges;
  }

  int getOrgArcCount() const
  {
    return _nArcs;
  }

  int getOrgComponentCount() const
  {
    return _nComponents;
  }

  int getOrgComponent(Node n) const
  {
    assert(n != lemon::INVALID);
    return (*_pComp)[n];
  }

  std::string getOrgLabel(Node n) const
  {
    assert(n != lemon::INVALID);
    return (*_pLabel)[n];
  }

  double getOrgScore(Node n) const
  {
    assert(n != lemon::INVALID);
    return (*_pScore)[n];
  }

  double getOrgPValue(Node n) const
  {
    assert(n != lemon::INVALID);
    assert(_pPVal);

    return (*_pPVal)[n];
  }

  virtual Node getPreNode(Node orgNode) const
  {
    return orgNode;
  }

  Node getOrgNodeByLabel(const std::string& label) const
  {
    InvLabelNodeMapIt it = _invLabel.find(label);

    if (it != _invLabel.end())
      return it->second;
    else
      return lemon::INVALID;
  }
  
  virtual void printNodeList(std::ostream& out,
                             bool orig = false) const;
  virtual void printEdgeList(std::ostream& out,
                             bool orig = false) const;

  virtual void print(std::ostream& out,
                     bool orig = false) const;
  virtual void printModule(const NodeSet& module,
                           std::ostream& out,
                           bool orig = false) const;
  virtual void printModule(const BoolNodeMap& module,
                           std::ostream& out,
                           bool orig = false) const;
  virtual void printModules(const NodeSetVector& modules,
                            std::ostream& out,
                            bool orig = false) const;
  virtual void printHeinz(const NodeSet& module,
                          std::ostream& out) const;
  virtual void printHeinzOrg(const NodeSet& module,
                          std::ostream& out) const;
  virtual void computeScores(double lambda, double a, double FDR) {}
  virtual void computeScores(double tau) {}
};

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline MwcsGraph<GR, NWGHT, NLBL, EWGHT>::MwcsGraph()
  : _pG(NULL)
  , _pLabel(NULL)
  , _pPVal(NULL)
  , _pScore(NULL)
  , _invLabel()
  , _pComp(NULL)
  , _nNodes(0)
  , _nEdges(0)
  , _nArcs(0)
  , _nComponents(0)
  , _pArcLookUp(NULL)
  , _parserInit(false)
{
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline bool MwcsGraph<GR, NWGHT, NLBL, EWGHT>::init(ParserType* pParser, bool pval)
{
  assert(pParser);

  if (_parserInit)
  {
    delete _pLabel;
    delete _pPVal;
    delete _pScore;
    delete _pG;
  }

  delete _pComp;
  _parserInit = true;
  initParserMembers(_pG, _pLabel, _pPVal, _pScore);

  pParser->setGraph(_pG);
  pParser->setIdNodeMap(_pLabel);

  if (pval)
    pParser->setWeightNodeMap(_pPVal);
  else
    pParser->setWeightNodeMap(_pScore);

  pParser->setInvIdNodeMap(&_invLabel);

  if (pParser->parse())
  {
    _nNodes = pParser->getNodeCount();
    _nEdges = pParser->getEdgeCount();

    // TODO: this is not the way to do it
    _nArcs = lemon::countArcs(*_pG);

    // determine the components
    _pComp = new IntNodeMap(*_pG, -1);
    _nComponents = lemon::connectedComponents(*_pG, *_pComp);

    _pArcLookUp = new ArcLookUpType(*_pG);

    if (g_verbosity >= VERBOSE_ESSENTIAL)
    {
      std::cout << "// Successfully parsed '"
                << pParser->getFilename()
                << "': contains " << _nNodes << " nodes, "
                << _nEdges << " edges and "
                << _nComponents << " component(s)" << std::endl;
    }

    return true;
  }
  else
  {
    return false;
  }
}


template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsGraph<GR, NWGHT, NLBL, EWGHT>::printModule(const BoolNodeMap& module,
                                                           std::ostream& out,
                                                           bool orig) const
{
  const Graph& g = orig ? getOrgGraph() : getGraph();

  NodeSet nodes;
  for (NodeIt n(g); n != lemon::INVALID; ++n)
  {
    if (module[n])
    {
      nodes.insert(n);
    }
  }

  printModule(nodes, out);
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsGraph<GR, NWGHT, NLBL, EWGHT>::printNodeList(std::ostream& out,
                                                             bool orig) const
{
  const Graph& g = orig ? getOrgGraph() : getGraph();
  const WeightNodeMap& weight = orig ? getOrgScores() : getScores();
  const LabelNodeMap& label = orig ? getOrgLabels() : getLabels();
  const WeightNodeMap* pPVal = orig ? getOrgPValues() : NULL;
  
  for (NodeIt n(g); n != lemon::INVALID; ++n)
  {
    out << g.id(n) << " " << weight[n];
    
    out << " " << label[n];
    
    if (pPVal)
    {
      out << " " << (*_pPVal)[n];
    }
    
    out << std::endl;
  }
}
  
template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsGraph<GR, NWGHT, NLBL, EWGHT>::printEdgeList(std::ostream& out,
                                                             bool orig) const
{
  const Graph& g = orig ? getOrgGraph() : getGraph();
//  const WeightNodeMap& weight = orig ? getOrgScores() : getScores();
//  const LabelNodeMap& label = orig ? getOrgLabels() : getLabels();
  
  SpqrTree<Graph> spqr(g);
  spqr.run();
  
  for (EdgeIt e(g); e != lemon::INVALID; ++e)
  {
    out << g.id(g.u(e)) << " (pp) " << g.id(g.v(e)) << "\t" << spqr.getSpqrTree().id(spqr.toSpqrNode(e)) << "\t";
    out << spqr.toChar(spqr.getSpqrNodeType(e));
   
    out << std::endl;
  }
}
  
template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsGraph<GR, NWGHT, NLBL, EWGHT>::print(std::ostream& out,
                                                     bool orig) const
{
  const Graph& g = orig ? getOrgGraph() : getGraph();
  const WeightNodeMap& weight = orig ? getOrgScores() : getScores();
  const LabelNodeMap& label = orig ? getOrgLabels() : getLabels();
  const WeightNodeMap* pPVal = orig ? getOrgPValues() : NULL;

  // header
  out << "graph G {" << std::endl;
  out << "\toverlap=scale" << std::endl;
  out << "\tlayout=neato" << std::endl;

  // nodes
  for (NodeIt n(g); n != lemon::INVALID; ++n)
  {
    out << "\t" << g.id(n) << " [label=\""
       << label[n] << "\\n"
       << weight[n] << "\\n";

    if (pPVal)
    {
      out << (*_pPVal)[n] << "\\n";
    }

    out << g.id(n)
       << "\""
       << (weight[n] < 0 ? ",shape=box" : "")
       << "]" << std::endl;
  }

  // edges
  for (EdgeIt e(g); e != lemon::INVALID; ++e)
  {
    out << "\t" << g.id(g.u(e)) << " -- " << g.id(g.v(e)) << std::endl;
  }

  // footer
  out << "}" << std::endl;
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsGraph<GR, NWGHT, NLBL, EWGHT>::printModule(const NodeSet& module,
                                                           std::ostream& out,
                                                           bool orig) const
{
  const Graph& g = orig ? getOrgGraph() : getGraph();
  const WeightNodeMap& weight = orig ? getOrgScores() : getScores();
  const LabelNodeMap& label = orig ? getOrgLabels() : getLabels();
  const ArcLookUpType& arcLookUp = orig ? getOrgArcLookUp() : getArcLookUp();
  const WeightNodeMap* pPVal = orig ? getOrgPValues() : NULL;

  std::vector<Edge> edges;

  // print header
  out << "graph G {" << std::endl;
  out << "\toverlap=scale" << std::endl;
  out << "\tlayout=neato" << std::endl;

  double totalWeight = 0;
  for (NodeSetIt nodeIt1 = module.begin(); nodeIt1 != module.end(); nodeIt1++)
  {
    totalWeight += weight[*nodeIt1];
    out << "\t" << g.id(*nodeIt1) << " [label=\""
        << label[*nodeIt1] << "\\n"
        << weight[*nodeIt1] << "\\n";

    if (pPVal)
    {
      out << (*_pPVal)[*nodeIt1] << "\\n";
    }

    out << g.id(*nodeIt1)
        << "\""
        << (weight[*nodeIt1] < 0 ? ",shape=box" : "")
        << "]" << std::endl;

    // determine incident edges
    for (NodeSetIt nodeIt2 = nodeIt1; nodeIt2 != module.end(); nodeIt2++)
    {
      Edge e = arcLookUp(*nodeIt1, *nodeIt2);
      if (e != lemon::INVALID)
      {
        edges.push_back(e);
      }
    }
  }

  out << "\tlabel=\"Total weight: " << totalWeight << '"' << std::endl;

  // print edges
  for (typename std::vector<Edge>::const_iterator edgeIt = edges.begin();
       edgeIt != edges.end(); edgeIt++)
  {
    Node u = g.u(*edgeIt);
    Node v = g.v(*edgeIt);

    out << "\t" << g.id(u) << " -- " << g.id(v) << std::endl;
  }

  // print footer
  out << "}" << std::endl;
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsGraph<GR, NWGHT, NLBL, EWGHT>::printHeinzOrg(const NodeSet& module,
                                                             std::ostream& out) const
{
  const Graph& g = getOrgGraph();

  out << "#label\tscore" << std::endl;
  double score = 0;
  for (NodeIt v(g); v != lemon::INVALID; ++v)
  {
    bool inSolution = (module.find(v) != module.end());
    if (inSolution)
    {
      score += getOrgScore(v);
      out << getOrgLabel(v) << "\t" << getOrgScore(v) << std::endl;
    }
    else
    {
      out << getOrgLabel(v) << "\tNaN" << std::endl;
    }
  }
  out << "#total score\t" << score << std::endl;
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsGraph<GR, NWGHT, NLBL, EWGHT>::printHeinz(const NodeSet& module,
                                                          std::ostream& out) const
{
  const Graph& g = getOrgGraph();

  out << "#label\tscore" << std::endl;
  double score = 0;
  for (NodeIt v(g); v != lemon::INVALID; ++v)
  {
    Node pre_v = this->getPreNode(v);
    bool inSolution = (module.find(pre_v) != module.end());
    if (inSolution)
    {
      score += getOrgScore(v);
      out << getOrgLabel(v) << "\t" << getOrgScore(v) << std::endl;
    }
    else
    {
      out << getOrgLabel(v) << "\tNaN" << std::endl;
    }
  }
  out << "#total score\t" << score << std::endl;
}


template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsGraph<GR, NWGHT, NLBL, EWGHT>::printModules(const NodeSetVector& modules,
                                                            std::ostream& out,
                                                            bool orig) const
{
  const Graph& g = orig ? getOrgGraph() : getGraph();
  const WeightNodeMap& weight = orig ? getOrgScores() : getScores();
  const LabelNodeMap& label = orig ? getOrgLabels() : getLabels();
  const ArcLookUpType& arcLookUp = orig ? getOrgArcLookUp() : getArcLookUp();

  NodeSet nodes;

  // print header
  out << "graph G {" << std::endl;
  out << "\toverlap=scale" << std::endl;
  out << "\tlayout=neato" << std::endl;

  // print nodes
  int i = 0;
  for (NodeSetVectorIt moduleIt = modules.begin(); moduleIt != modules.end(); moduleIt++, i++)
  {
    if (moduleIt->size() == 0)
      continue;

    out << "\tsubgraph cluster_" << i << " {" << std::endl;

    double totalWeight = 0;
    for (NodeSetIt nodeIt1 = moduleIt->begin(); nodeIt1 != moduleIt->end(); nodeIt1++)
    {
      totalWeight += weight[*nodeIt1];
      out << "\t\t" << g.id(*nodeIt1) << " [label=\""
          << label[*nodeIt1] << "\\n"
          << weight[*nodeIt1] << "\\n"
          << g.id(*nodeIt1)
          << "\"]" << std::endl;

      nodes.insert(*nodeIt1);
    }

    out << "\t\tlabel=\"Total weight: " << totalWeight << "\"" << std::endl;
    out << "\t}" << std::endl;
  }

  // print edges
  for (NodeSetIt nodeIt1 = nodes.begin(); nodeIt1 != nodes.end(); nodeIt1++)
  {
    for (NodeSetIt nodeIt2 = nodeIt1; nodeIt2 != nodes.end(); nodeIt2++)
    {
      if (arcLookUp(*nodeIt1, *nodeIt2) != lemon::INVALID)
      {
        out << "\t" << g.id(*nodeIt1) << " -- " << g.id(*nodeIt2) << std::endl;
      }
    }
  }

  // print footer
  out << "}" << std::endl;
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsGraph<GR, NWGHT, NLBL, EWGHT>::resetCounts()
{
  _nNodes = lemon::countNodes(*_pG);
  _nEdges = lemon::countEdges(*_pG);
  _nArcs = lemon::countArcs(*_pG);
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline bool MwcsGraph<GR, NWGHT, NLBL, EWGHT>::init(Graph* pG,
                                                    LabelNodeMap* pLabel,
                                                    WeightNodeMap* pScore,
                                                    WeightNodeMap* pPVal)
{
  if (_parserInit)
  {
    delete _pLabel;
    delete _pPVal;
    delete _pScore;
    delete _pG;
  }

  delete _pComp;
  _pG = pG;
  _pLabel = pLabel;
  _pPVal = pPVal;
  _pScore = pScore;

  _nNodes = lemon::countNodes(*_pG);
  _nEdges = lemon::countEdges(*_pG);
  _nArcs = lemon::countArcs(*_pG);

  // determine the components
  _pComp = new IntNodeMap(*_pG, -1);
  _nComponents = lemon::connectedComponents(*_pG, *_pComp);

  _parserInit = false;

  return true;
}

} // namespace mwcs
} // namespace nina

#endif // MWCSGRAPH_H
