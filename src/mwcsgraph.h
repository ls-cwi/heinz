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
#include "utils.h"
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
  
public:
  MwcsGraph();
  virtual ~MwcsGraph()
  {
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

  virtual NodeSet getNodeByLabel(const std::string& label) const
  {
    NodeSet res;
    
    Node node = getOrgNodeByLabel(label);
    if (node != lemon::INVALID)
    {
      res.insert(node);
    }
    
    return res;
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
  
  const IntNodeMap& getOrgComponentMap() const
  {
    return *_pComp;
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

  virtual NodeSet getPreNodes(Node orgNode) const
  {
    return NodeSet();
  }

  Node getOrgNodeByLabel(const std::string& label) const
  {
    InvLabelNodeMapIt it = _invLabel.find(label);

    if (it != _invLabel.end())
      return it->second;
    else
      return lemon::INVALID;
  }
  
  bool allNodesNegative() const
  {
    const Graph& g = getGraph();
    for (NodeIt v(g); v != lemon::INVALID; ++v)
    {
      if (getScore(v) > 0)
      {
        return false;
      }
    }
    return true;
  }
  
  double getTotalNodeProfitPCST() const;
  
  virtual void printNodeList(std::ostream& out,
                             bool orig = false) const;
  virtual void printEdgeList(std::ostream& out,
                             bool orig = false) const;
  virtual void print(std::ostream& out,
                     bool orig = false) const;
  virtual void printSTP(const std::string& name,
                        const std::string& creator,
                        std::ostream& out,
                        bool orig = false) const;
  
  virtual void printModule(const NodeSet& module,
                           std::ostream& out,
                           bool orig = false) const;
  virtual void printModule(const BoolNodeMap& module,
                           std::ostream& out,
                           bool orig = false) const;
  virtual void printHeinz(const NodeSet& module,
                          std::ostream& out) const;
  virtual void printHeinzOrg(const NodeSet& module,
                          std::ostream& out) const;
  virtual void printMwcsDimacs(const NodeSet& module,
                               std::ostream& out) const;
  virtual void printPcstDimacs(const NodeSet& module,
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
  }

  out << "\tlabel=\"Total weight: " << totalWeight << '"' << std::endl;

  // print edges
  for (EdgeIt e(g); e != lemon::INVALID; ++e)
  {
    Node u = g.u(e);
    Node v = g.v(e);
    
    if (module.find(u) != module.end() && module.find(v) != module.end())
    {
      out << "\t" << g.id(u) << " -- " << g.id(v) << std::endl;
    }
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
    NodeSet preNodes = this->getPreNodes(v);
    bool inSolution = false;
    for (NodeSetIt preNodeIt = preNodes.begin(); preNodeIt != preNodes.end(); ++preNodeIt)
    {
      inSolution |= module.find(*preNodeIt) != module.end();
    }
    
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
inline void MwcsGraph<GR, NWGHT, NLBL, EWGHT>::printMwcsDimacs(const NodeSet& module,
                                                               std::ostream& out) const
{
  // determine vertices
  int n = 0;
  for (NodeSetIt nodeIt = module.begin(); nodeIt != module.end(); ++nodeIt)
  {
    n += getOrgNodes(*nodeIt).size();
  }
  out << "Vertices " << n << std::endl;
  
  std::set<Node> nodes;
  for (NodeSetIt nodeIt = module.begin(); nodeIt != module.end(); ++nodeIt)
  {
    NodeSet orgNodes = getOrgNodes(*nodeIt);
    for (NodeSetIt nodeIt2 = orgNodes.begin(); nodeIt2 != orgNodes.end(); ++nodeIt2)
    {
      nodes.insert(*nodeIt2);
      out << "V " << getOrgLabel(*nodeIt2) << std::endl;
    }
  }
  
  // determine edges
  std::set<Edge> edges;
  const Graph& g = getOrgGraph();
  for (EdgeIt e(g); e != lemon::INVALID; ++e)
  {
    Node u = g.u(e);
    Node v = g.v(e);
    
    if (nodes.find(u) != nodes.end() && nodes.find(v) != nodes.end())
    {
      edges.insert(e);
    }
  }
  
  out << "Edges " << edges.size() << std::endl;
  for (typename std::set<Edge>::const_iterator edgeIt = edges.begin();
       edgeIt != edges.end(); ++edgeIt)
  {
    out << "E " << getOrgLabel(g.u(*edgeIt)) << " " << getOrgLabel(g.v(*edgeIt)) << std::endl;
  }
}
  
template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsGraph<GR, NWGHT, NLBL, EWGHT>::printPcstDimacs(const NodeSet& module,
                                                               std::ostream& out) const
{
  typedef std::pair<int, int> IntPair;
  typedef std::vector<IntPair> IntPairVector;
  typedef IntPairVector::const_iterator IntPairVectorIt;
  typedef std::vector<int> IntVector;
  typedef IntVector::const_iterator IntVectorIt;
  
  IntVector vertices;
  IntPairVector edges;
  
  // determine vertices and edges
  for (NodeSetIt nodeIt = module.begin(); nodeIt != module.end(); ++nodeIt)
  {
    NodeSet orgNodes = getOrgNodes(*nodeIt);
    for (NodeSetIt nodeIt2 = orgNodes.begin(); nodeIt2 != orgNodes.end(); ++nodeIt2)
    {
      const std::string& label = getOrgLabel(*nodeIt2);
      int u = -1, v = -1;
      if (sscanf(label.c_str(), "%d--%d", &u, &v) == 2)
      {
        edges.push_back(std::make_pair(u, v));
      }
      else if (sscanf(label.c_str(), "%d", &u) == 1)
      {
        vertices.push_back(u);
      }
      else
      {
        assert(false);
      }
    }
  }
  
  out << "Vertices " << vertices.size() << std::endl;
  for (IntVectorIt nodeIt = vertices.begin(); nodeIt != vertices.end(); ++nodeIt)
  {
    out << "V " << *nodeIt << std::endl;
  }
  
  out << "Edges " << edges.size() << std::endl;
  for (IntPairVectorIt edgeIt = edges.begin(); edgeIt != edges.end(); ++edgeIt)
  {
    out << "E " << edgeIt->first << " " << edgeIt->second << std::endl;
  }
}
  
template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsGraph<GR, NWGHT, NLBL, EWGHT>::printSTP(const std::string& name,
                                                        const std::string& creator,
                                                        std::ostream& out,
                                                        bool orig) const
{
  const Graph& g = orig ? getOrgGraph() : getGraph();
  const WeightNodeMap& weight = orig ? getOrgScores() : getScores();

  out << "33D32945 STP File, STP Format Version 1.0" << std::endl;
  out << std::endl;
  
  out << "SECTION Comment" << std::endl;
  out << "Name \"" << name << '"' << std::endl;
  out << "Creator \"" << creator << '"' << std::endl;
  out << "Problem \"MWCS\"" << std::endl;
  out << "END" << std::endl;
  out << std::endl;
  
  out << "SECTION Graph" << std::endl;
  out << "Nodes " << getNodeCount() << std::endl;
  out << "Edges " << getEdgeCount() << std::endl;
  
  IntNodeMap id(g);
  int idx = 1;
  for (NodeIt v(g); v != lemon::INVALID; ++v, ++idx)
  {
    id[v] = idx;
  }
  
  for (EdgeIt e(g); e != lemon::INVALID; ++e)
  {
    out << "E " << id[g.u(e)] << " " << id[g.v(e)] << std::endl;
  }
  out << "END" << std::endl;
  out << std::endl;
  
  out << "SECTION Terminals" << std::endl;
  out << "Terminals " << getNodeCount() << std::endl;
  for (NodeIt v(g); v != lemon::INVALID; ++v)
  {
    out << "T " << id[v] << " " << weight[v] << std::endl;
  }
  out << "END" << std::endl;
  out << std::endl;
  out << "EOF" << std::endl;
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

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline double MwcsGraph<GR, NWGHT, NLBL, EWGHT>::getTotalNodeProfitPCST() const
{
  double res = 0;
  
  const Graph& orgG = getOrgGraph();
  for (NodeIt v(orgG); v != lemon::INVALID; ++v)
  {
    const std::string& label = getOrgLabel(v);
    
    int idU = -1;
    char c = '\0';
    if (sscanf(label.c_str(), "%d%c", &idU, &c) == 1)
    {
      res += getOrgScore(v);
    }
  }
  
  return res;
}
  
} // namespace mwcs
} // namespace nina

#endif // MWCSGRAPH_H
