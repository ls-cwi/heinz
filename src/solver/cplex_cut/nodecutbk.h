/*
 * nodecutbk.h
 *
 *  Created on: 17-feb-2014
 *      Author: M. El-Kebir
 */

#ifndef NODECUTBK_H
#define NODECUTBK_H

#include <ilcplex/ilocplex.h>
#include <ilcplex/ilocplexi.h>
#include <ilconcert/ilothread.h>
#include <lemon/tolerance.h>
#include <lemon/smart_graph.h>
#include <lemon/adaptors.h>
#include <vector>
#include <set>
#include <queue>
#include <list>
#include "bk_alg.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class NodeCutBk
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;

  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  typedef lemon::SmartDigraph Digraph;
  typedef typename Digraph::Arc DiArc;
  typedef typename Digraph::Node DiNode;
  typedef typename Digraph::NodeIt DiNodeIt;
  typedef typename Digraph::ArcIt DiArcIt;
  typedef typename Digraph::InArcIt DiInArcIt;
  typedef typename Digraph::OutArcIt DiOutArcIt;

  typedef typename Digraph::template NodeMap<Node> DiNodeNodeMap;
  typedef typename Graph::template NodeMap<DiNode> NodeDiNodeMap;
  typedef typename Graph::template NodeMap<DiArc> NodeDiArcMap;
  typedef typename Digraph::ArcMap<double> CapacityMap;

  typedef std::vector<Node> NodeVector;
  typedef typename NodeVector::const_iterator NodeVectorIt;
  typedef std::vector<NodeVector> NodeMatrix;
  typedef std::vector<double> DoubleVector;
  typedef std::set<Node> NodeSet;
  typedef typename NodeSet::const_iterator NodeSetIt;
  typedef std::set<DiArc> DiArcSet;
  typedef typename DiArcSet::const_iterator DiArcSetIt;

  typedef std::queue<Node> NodeQueue;
  typedef std::queue<DiNode> DiNodeQueue;
  typedef std::set<DiNode> DiNodeSet;
  typedef std::list<DiNode> DiNodeList;
  typedef DiNodeList::const_iterator DiNodeListIt;

  typedef lemon::FilterNodes<const Graph, const BoolNodeMap> SubGraph;
  typedef typename SubGraph::NodeIt SubNodeIt;
  typedef nina::BkFlowAlg<Digraph> BkAlg;
  typedef typename Digraph::NodeMap<bool> DiBoolNodeMap;

protected:
  IloBoolVarArray _x;
  const Graph& _g;
  const WeightNodeMap& _weight;
  const Node _root;
  const IntNodeMap& _nodeMap;
  const int _n;
  const int _m;
  const int _maxNumberOfCuts;
  const IntNodeMap& _comp;
  const lemon::Tolerance<double> _tol;
  Digraph _h;
  CapacityMap _cap;
  NodeDiNodeMap* _pG2h1;
  NodeDiNodeMap* _pG2h2;
  DiNodeNodeMap _h2g;
  DiNode _diRoot;
  BkAlg* _pBK;
  DiBoolNodeMap _marked;
  IloFastMutex* _pMutex;

  // used for determining non-zero components
  BoolNodeMap* _pNodeFilterMap;
  const SubGraph* _pSubG;
  IntNodeMap* _pComp;

public:
  NodeCutBk(IloBoolVarArray x,
            const Graph& g,
            const WeightNodeMap& weight,
            Node root,
            const IntNodeMap& nodeMap,
            int n,
            int m,
            int maxNumberOfCuts,
            const IntNodeMap& comp,
            IloFastMutex* pMutex)
    : _x(x)
    , _g(g)
    , _weight(weight)
    , _root(root)
    , _nodeMap(nodeMap)
    , _n(n)
    , _m(m)
    , _maxNumberOfCuts(maxNumberOfCuts)
    , _comp(comp)
    , _tol()
    , _h()
    , _cap(_h)
    , _pG2h1(NULL)
    , _pG2h2(NULL)
    , _h2g(_h)
    , _diRoot(lemon::INVALID)
    , _pBK(NULL)
    , _marked(_h, false)
    , _pMutex(pMutex)
    , _pNodeFilterMap(NULL)
    , _pSubG(NULL)
    , _pComp(NULL)
  {
    lock();
    _pG2h1 = new NodeDiNodeMap(_g);
    _pG2h2 = new NodeDiNodeMap(_g);
    _pNodeFilterMap = new BoolNodeMap(_g);
    _pSubG = new SubGraph(_g, *_pNodeFilterMap);
    _pComp = new IntNodeMap(_g);
    unlock();
  }

  NodeCutBk(const NodeCutBk& other)
    : _x(other._x)
    , _g(other._g)
    , _weight(other._weight)
    , _root(other._root)
    , _nodeMap(other._nodeMap)
    , _n(other._n)
    , _m(other._m)
    , _maxNumberOfCuts(other._maxNumberOfCuts)
    , _comp(other._comp)
    , _tol(other._tol)
    , _h()
    , _cap(_h)
    , _pG2h1(NULL)
    , _pG2h2(NULL)
    , _h2g(_h)
    , _diRoot(lemon::INVALID)
    , _pBK(NULL)
    , _marked(_h, false)
    , _pMutex(other._pMutex)
    , _pNodeFilterMap(NULL)
    , _pSubG(NULL)
    , _pComp(NULL)
  {
  }

  virtual ~NodeCutBk()
  {
    delete _pBK;

    lock();
    delete _pG2h1;
    delete _pG2h2;
    delete _pSubG;
    delete _pNodeFilterMap;
    delete _pComp;
    unlock();
  }

protected:
  void lock()
  {
    if (_pMutex)
      _pMutex->lock();
  }

  void unlock()
  {
    if (_pMutex)
      _pMutex->unlock();
  }

  virtual void init() = 0;

  virtual bool isUser() const = 0;

  virtual void addConstraint(IloConstraint con) = 0;

  void printNonZeroVars(IloCplex::ControlCallbackI& cbk,
                        IloBoolVarArray variables,
                        IloNumArray values) const
  {
    std::cerr << cbk.getNnodes() << ":";
    for (NodeIt i(_g); i != lemon::INVALID; ++i)
    {
      double i_value = values[_nodeMap[i]];
      if (!_tol.nonZero(i_value)) continue;
      std::cerr << " " << variables[_nodeMap[i]].getName()
                << " (" << _g.id(i) << ", " << _weight[i] << ", " << i_value << ") " ;

      if (cbk.getDirection(variables[_nodeMap[i]]) == CPX_BRANCH_UP)
        std::cerr << "*";
    }
    std::cerr << std::endl;
  }

  void printNodeSet(const NodeSet& nodes,
                    IloBoolVarArray variables,
                    IloNumArray values) const
  {
    bool first = true;
    for (NodeSetIt it = nodes.begin(); it != nodes.end(); it++)
    {
      if (!first)
        std::cout << " ";
      else
        first = false;

      std::cout << values[_nodeMap[*it]];
    }
    std::cout << std::endl;
  }
};

} // namespace mwcs
} // namespace nina


#endif // NODECUTROOTEDBK_H
