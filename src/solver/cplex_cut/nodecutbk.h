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

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class NodeCutRootedBkCallback
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

  NodeCutBkCallback(IloBoolVarArray x,
                    const Graph& g,
                    const WeightNodeMap& weight,
                    Node root,
                    const IntNodeMap& nodeMap,
                    int n,
                    int m,
                    int maxNumberOfCuts,
                    const IntNodeMap& comp,
                    IloFastMutex* pMutex)
    : IloCplex::LazyConstraintCallbackI(env)
    , _x(x)
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

    init();
    _pBK = new BkAlg(_h, _cap);
  }

protected:

};

} // namespace mwcs
} // namespace nina


#endif // NODECUTROOTEDBK_H
