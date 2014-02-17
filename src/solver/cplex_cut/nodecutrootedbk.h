/*
 * nodecutrootedbk.h
 *
 *  Created on: 12-may-2013
 *      Author: M. El-Kebir
 */


#ifndef NODECUTROOTEDBK_H
#define NODECUTROOTEDBK_H

#include <ilcplex/ilocplex.h>
#include <ilcplex/ilocplexi.h>
#include <ilconcert/ilothread.h>
#include <lemon/tolerance.h>
#include <lemon/time_measure.h>
#include <lemon/smart_graph.h>
#include <lemon/connectivity.h>
#include <vector>
#include <set>
#include <queue>
#include <stack>
#include <list>
#include "bk_alg.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class NodeCutRootedBkCallback : public IloCplex::LazyConstraintCallbackI
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

  NodeCutRootedBkCallback(IloEnv env,
                          IloBoolVarArray x,
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

  NodeCutRootedBkCallback(const NodeCutRootedBkCallback& other)
    : IloCplex::LazyConstraintCallbackI(other)
    , _x(other._x)
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
    typename Digraph::template NodeMap<DiNode> nodeMap(other._h);
    typename Digraph::template ArcMap<DiArc> arcMap(other._h);

    lemon::digraphCopy(other._h, _h)
        .nodeRef(nodeMap)
        .arcRef(arcMap)
        .nodeMap(other._h2g, _h2g)
        .arcMap(other._cap, _cap)
        .run();

    lock();
    _pG2h1 = new NodeDiNodeMap(_g);
    _pG2h2 = new NodeDiNodeMap(_g);
    _pNodeFilterMap = new BoolNodeMap(_g);
    _pSubG = new SubGraph(_g, *_pNodeFilterMap);
    _pComp = new IntNodeMap(_g);
    unlock();

    for (NodeIt v(_g); v != lemon::INVALID; ++v)
    {
      DiNode v1 = (*other._pG2h1)[v];
      DiNode v2 = (*other._pG2h2)[v];

      _pG2h1->set(v, nodeMap[v1]);
      if (v2 != lemon::INVALID)
        _pG2h2->set(v, nodeMap[v2]);
    }

    _diRoot = nodeMap[other._diRoot];
    _pBK = new BkAlg(_h, _cap);
  }

  virtual ~NodeCutRootedBkCallback()
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
  virtual void main();
  virtual IloCplex::CallbackI* duplicateCallback() const
  {
    return (new (getEnv()) NodeCutRootedBkCallback(*this));
  }

private:
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

  typedef nina::BkFlowAlg<Digraph> BkAlg;
  typedef typename Digraph::NodeMap<bool> DiBoolNodeMap;
  Digraph _h;
  CapacityMap _cap;
  NodeDiNodeMap* _pG2h1;
  NodeDiNodeMap* _pG2h2;
  DiNodeNodeMap _h2g;
  DiNode _diRoot;
  BkAlg* _pBK;
  DiBoolNodeMap _marked;
  IloFastMutex* _pMutex;

  typedef lemon::FilterNodes<const Graph, const BoolNodeMap> SubGraph;
  typedef typename SubGraph::NodeIt SubNodeIt;
  BoolNodeMap* _pNodeFilterMap;
  const SubGraph* _pSubG;
  IntNodeMap* _pComp;

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

  void init()
  {
    // we initialize _h:
    // - for every node i, there will be two nodes i1 and i2
    //   connected by an arc from i1 to i2
    // - for every edge (i,j) there are two arcs
    //   in h: (i2,j1) and (j2,i1) with capacties 1
    _h.clear();
    for (NodeIt i(_g); i != lemon::INVALID; ++i)
    {
      if (i != _root)
      {
        DiNode i1 = _h.addNode();
        DiNode i2 = _h.addNode();
        _pG2h1->set(i, i1);
        _pG2h2->set(i, i2);
        _h2g[i1] = i;
        _h2g[i2] = i;

        DiArc i1i2 = _h.addArc(i1, i2);
        _cap[i1i2] = 0;
      }
      else
      {
        _diRoot = _h.addNode();
        _pG2h1->set(i, _diRoot);
        _pG2h2->set(i, _diRoot);
        _h2g[_diRoot] = _root;
      }
    }

    for (EdgeIt e(_g); e != lemon::INVALID; ++e)
    {
      Node i = _g.u(e);
      Node j = _g.v(e);
      DiNode i1 = (*_pG2h1)[i];
      DiNode i2 = (*_pG2h2)[i];
      DiNode j1 = (*_pG2h1)[j];
      DiNode j2 = (*_pG2h2)[j];

      if (i != _root && j != _root)
      {
        DiArc i2j1 = _h.addArc(i2, j1);
        DiArc j2i1 = _h.addArc(j2, i1);
        _cap[i2j1] = _cap[j2i1] = 1;
      }
      else if (i == _root)
      {
        DiArc ij1 = _h.addArc(_diRoot, j1);
        _cap[ij1] = 1;
      }
      else if (j == _root)
      {
        DiArc ji1 = _h.addArc(_diRoot, i1);
        _cap[ji1] = 1;
      }
    }
  }

  void computeCapacities(CapacityMap& capacity,
                         IloNumArray x_values)
  {
    // cap((i,j)) = x_i
    for (NodeIt v(_g); v != lemon::INVALID; ++v)
    {
      double val = x_values[_nodeMap[v]];
      if (!_tol.nonZero(val)) val = 1e-9;

      DiNode v1 = (*_pG2h1)[v];
      DiOutArcIt a(_h, v1);

      if (a != lemon::INVALID)
        capacity[a] = val;
      else
        assert(v == _root);
    }
  }

  void addViolatedConstraint(Node target, const NodeSet& dS)
  {
    if (dS.empty())
    {
      //std::cout << getNnodes() << ": " << _x[_nodeMap[target]].getName() << " <= 0" << std::endl;
      add(_x[_nodeMap[target]] <= 0);
    }
    else
    {
      IloExpr expr(getEnv());

      //bool first = true;
      //std::cout << getNnodes() << ": " << _x[_nodeMap[target]].getName() << " <=";
      for (NodeSetIt nodeIt = dS.begin(); nodeIt != dS.end(); nodeIt++)
      {
        expr += _x[_nodeMap[*nodeIt]];

        //std::cout << (first ? " " : " + ") << _x[_nodeMap[*nodeIt]].getName();
        //first = false;
      }
      //std::cout << std::endl;

      IloConstraint constraint = _x[_nodeMap[target]] <= expr;
      add(constraint);
      constraint.end();

      expr.end();
    }
  }

  void determineBwdCutSet(const Digraph& h,
                          const BkAlg& bk,
                          NodeSet& dS)
  {
    DiNodeList diS;
    DiNode target = bk.getTarget();

    lemon::mapFill(_h, _marked, false);

    DiNodeQueue queue;
    queue.push(target);
    _marked[target] = true;

    while (!queue.empty())
    {
      DiNode v = queue.front();
      queue.pop();
      diS.push_back(v);

      for (DiInArcIt a(h, v); a != lemon::INVALID; ++a)
      {
        DiNode u = _h.source(a);

        if (!_marked[u] && _tol.nonZero(bk.resCap(a)))
        {
          queue.push(u);
          _marked[u] = true;
        }
      }
    }

    for (DiNodeListIt nodeIt = diS.begin(); nodeIt != diS.end(); nodeIt++)
    {
      DiNode v = *nodeIt;
      assert(_marked[v]);
      if (v == target) continue;

      for (DiInArcIt a(h, v); a != lemon::INVALID; ++a)
      {
        DiNode u = h.source(a);
        if (!_marked[u])
        {
          //std::cout << _h.id(u) << " -> "
          //          << _h.id(v) << " "
          //          << bk.flow(a) << "/" << bk.cap(a) << std::endl;
          dS.insert(_h2g[v]);
        }
      }
    }
  }

  void determineFwdCutSet(const Digraph& h,
                          const BkAlg& bk,
                          NodeSet& dS)
  {
    DiNode target = bk.getTarget();
    DiNodeList diS;

    lemon::mapFill(_h, _marked, false);

    DiNodeQueue queue;
    queue.push(_diRoot);
    _marked[_diRoot] = true;

    while (!queue.empty())
    {
      DiNode v = queue.front();
      queue.pop();
      diS.push_back(v);

      for (DiOutArcIt a(h, v); a != lemon::INVALID; ++a)
      {
        DiNode w = _h.target(a);

        if (!_marked[w] && _tol.nonZero(bk.resCap(a)))
        {
          queue.push(w);
          _marked[w] = true;
        }
      }
    }

    for (DiNodeListIt nodeIt = diS.begin(); nodeIt != diS.end(); nodeIt++)
    {
      DiNode v = *nodeIt;
      assert(_marked[v]);

      for (DiOutArcIt a(h, v); a != lemon::INVALID; ++a)
      {
        DiNode w = h.target(a);
        if (!_marked[w] && w != target)
        {
          //std::cout << _h.id(v) << " -> "
          //          << _h.id(w) << " "
          //          << bk.flow(a) << "/" << bk.cap(a) << std::endl;
          dS.insert(_h2g[w]);
        }
      }
    }
  }
};

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void NodeCutRootedBkCallback<GR, NWGHT, NLBL, EWGHT>::main()
{
  //lemon::Timer t;
  IloNumArray x_values(getEnv(), _n);
  getValues(x_values, _x);

  computeCapacities(_cap, x_values);

  int nCuts = 0;
  int nBackCuts = 0;
  int nNestedCuts = 0;

  lemon::mapFill(_g, *_pNodeFilterMap, false);

  for (NodeIt i(_g); i != lemon::INVALID; ++i)
  {
    double x_i_value = x_values[_nodeMap[i]];
    if (_tol.nonZero(x_i_value))
      _pNodeFilterMap->set(i, true);
  }

  int nComp = lemon::connectedComponents(*_pSubG, *_pComp);

  typedef std::pair<double, Node> NodeWeightPair;
  typedef std::vector<NodeWeightPair> NodeWeightPairVector;
  typedef std::vector<NodeWeightPairVector> NodeWeightPairMatrix;

  NodeWeightPairMatrix compMatrix(nComp, NodeWeightPairVector());

  for (SubNodeIt i(*_pSubG); i != lemon::INVALID; ++i)
  {
    if (i == _root)
      continue;

    double x_i_value = x_values[_nodeMap[i]];
    int compIdx = (*_pComp)[i];

    compMatrix[compIdx].push_back(std::make_pair(x_i_value, i));
  }

  // sort compMatrix
  for (int compIdx = 0; compIdx < nComp; compIdx++)
  {
    std::sort(compMatrix[compIdx].begin(), compMatrix[compIdx].end());
  }

  for (int compIdx = 0; compIdx < nComp; compIdx++)
  {
    bool foundCut = false;
    const NodeWeightPairVector& compVector = compMatrix[compIdx];
    for (typename NodeWeightPairVector::const_iterator it = compVector.begin();
         !foundCut && it != compVector.end(); ++it)
    {
      const double x_i_value = it->first;
      const Node i = it->second;

      _pBK->setSource(_diRoot);
      _pBK->setTarget((*_pG2h2)[i]);
      _pBK->setCap(_cap);

      bool first = true;
      bool nestedCut = false;
      while (true)
      {
        if (first)
        {
          _pBK->run();
          first = false;
        }
        else
        {
          _pBK->run(true);
        }

        // let's see if there's a violated constraint
        double minCutValue = _pBK->maxFlow();
        if (_tol.less(minCutValue, x_i_value))
        {
          foundCut = true;

          // determine N (forward)
          NodeSet fwdDS;
          determineFwdCutSet(_h, *_pBK, fwdDS);

          NodeSet bwdDS;
          determineBwdCutSet(_h, *_pBK, bwdDS);

          // add violated constraints
          for (typename NodeWeightPairVector::const_iterator it2 = compVector.begin(); it2 != compVector.end(); ++it2)
          {
            //const double x_j_value = it2->first;
            const Node j = it2->second;
            assert(_tol.less(minCutValue, it2->first));
            addViolatedConstraint(j, fwdDS);

            nCuts++;
            if (nestedCut)
            {
              nNestedCuts++;
            }
          }

          if (fwdDS.size() != bwdDS.size() || fwdDS != bwdDS)
          {
            for (typename NodeWeightPairVector::const_iterator it2 = compVector.begin(); it2 != compVector.end(); ++it2)
            {
              //const double x_j_value = it2->first;
              const Node j = it2->second;
              assert(_tol.less(minCutValue, it2->first));
              addViolatedConstraint(j, bwdDS);
              nBackCuts++;
              nCuts++;
            }
          }

          // generate nested-cuts
          for (NodeSetIt nodeIt = fwdDS.begin(); nodeIt != fwdDS.end(); nodeIt++)
          {
            nestedCut = true;
            // update the capactity to generate nested-cuts
            _pBK->incCap(DiOutArcIt(_h, (*_pG2h1)[*nodeIt]), 1);
          }

          if (fwdDS.empty()) break;
        }
        else
        {
          break;
        }
      }
    }
  }

  x_values.end();
  //std::cerr << "Generated " << nCuts
  //          << " cuts of which " << nBackCuts << " are back-cuts and "
  //          << nNestedCuts << " are nested cuts" << std::endl;
  //std::cerr << "Time: " << t.realTime() << "s" << std::endl;

  //std::cerr << "[";
  //for (int idx = 0; idx < nComp; idx++)
  //{
  //  std::cerr << " " << compMatrix[idx].size();
  //}
  //std::cerr << " ]" << std::endl;
}

} // namespace mwcs
} // namespace nina

#endif // NODECUTROOTEDBK_H
