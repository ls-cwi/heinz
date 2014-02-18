/*
 * nodecutunrootedbk.h
 *
 *  Created on: 12-may-2013
 *      Author: M. El-Kebir
 */

#ifndef NODECUTUNROOTEDBK_H
#define NODECUTUNROOTEDBK_H

#include <ilcplex/ilocplex.h>
#include <ilcplex/ilocplexi.h>
#include <ilconcert/ilothread.h>
#include <lemon/time_measure.h>
#include <lemon/connectivity.h>
#include "nodecutbk.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class NodeCutUnrootedBk : public NodeCutBk<GR, NWGHT, NLBL, EWGHT>
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;
  typedef NodeCutBk<GR, NWGHT, NLBL, EWGHT> Parent;

protected:
  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  typedef typename Parent::Digraph Digraph;
  typedef typename Parent::DiArc DiArc;
  typedef typename Parent::DiNode DiNode;
  typedef typename Parent::DiNodeIt DiNodeIt;
  typedef typename Parent::DiArcIt DiArcIt;
  typedef typename Parent::DiInArcIt DiInArcIt;
  typedef typename Parent::DiOutArcIt DiOutArcIt;
  typedef typename Parent::DiNodeNodeMap DiNodeNodeMap;
  typedef typename Parent::NodeDiNodeMap NodeDiNodeMap;
  typedef typename Parent::NodeDiArcMap NodeDiArcMap;
  typedef typename Parent::CapacityMap CapacityMap;
  typedef typename Parent::NodeVector NodeVector;
  typedef typename Parent::NodeVectorIt NodeVectorIt;
  typedef typename Parent::NodeMatrix NodeMatrix;
  typedef typename Parent::DoubleVector DoubleVector;
  typedef typename Parent::NodeSet NodeSet;
  typedef typename Parent::NodeSetIt NodeSetIt;
  typedef typename Parent::DiArcSet DiArcSet;
  typedef typename Parent::DiArcSetIt DiArcSetIt;
  typedef typename Parent::NodeQueue NodeQueue;
  typedef typename Parent::DiNodeQueue DiNodeQueue;
  typedef typename Parent::DiNodeSet DiNodeSet;
  typedef typename Parent::DiNodeList DiNodeList;
  typedef typename Parent::DiNodeListIt DiNodeListIt;
  typedef typename Parent::SubGraph SubGraph;
  typedef typename Parent::SubNodeIt SubNodeIt;
  typedef typename Parent::BkAlg BkAlg;
  typedef typename Parent::DiBoolNodeMap DiBoolNodeMap;

  using Parent::_x;
  using Parent::_g;
  using Parent::_weight;
  using Parent::_root;
  using Parent::_nodeMap;
  using Parent::_n;
  using Parent::_m;
  using Parent::_maxNumberOfCuts;
  using Parent::_comp;
  using Parent::_tol;
  using Parent::_h;
  using Parent::_cap;
  using Parent::_pG2h1;
  using Parent::_pG2h2;
  using Parent::_h2g;
  using Parent::_diRoot;
  using Parent::_pBK;
  using Parent::_marked;
  using Parent::_pMutex;
  using Parent::_pNodeFilterMap;
  using Parent::_pSubG;
  using Parent::_pComp;
  using Parent::lock;
  using Parent::unlock;
  using Parent::printNonZeroVars;
  using Parent::addConstraint;

  IloBoolVarArray _y;
  NodeDiArcMap* _pG2hRootArc;

public:
  NodeCutUnrootedBk(IloBoolVarArray x,
                    IloBoolVarArray y,
                    const Graph& g,
                    const WeightNodeMap& weight,
                    const IntNodeMap& nodeMap,
                    int n,
                    int m,
                    int maxNumberOfCuts,
                    const IntNodeMap& comp,
                    IloFastMutex* pMutex)
    : Parent(x, g, weight, lemon::INVALID, nodeMap, n, m, maxNumberOfCuts, comp, pMutex)
    , _y(y)
    , _pG2hRootArc(NULL)
  {
    lock();
    _pG2hRootArc = new NodeDiArcMap(_g);
    unlock();

    init();
    _pBK = new BkAlg(_h, _cap);
  }

  NodeCutUnrootedBk(const NodeCutUnrootedBk& other)
    : Parent(other)
    , _y(other._y)
    , _pG2hRootArc(NULL)
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
    _pG2hRootArc = new NodeDiArcMap(_g);
    _pNodeFilterMap = new BoolNodeMap(_g);
    _pSubG = new SubGraph(_g, *_pNodeFilterMap);
    _pComp = new IntNodeMap(_g);
    unlock();

    for (NodeIt v(_g); v != lemon::INVALID; ++v)
    {
      DiNode v1 = (*other._pG2h1)[v];
      DiNode v2 = (*other._pG2h2)[v];
      DiArc root_arc_v = (*other._pG2hRootArc)[v];

      _pG2h1->set(v, nodeMap[v1]);
      _pG2h2->set(v, nodeMap[v2]);
      _pG2hRootArc->set(v, arcMap[root_arc_v]);
    }

    _diRoot = nodeMap[other._diRoot];
    _pBK = new BkAlg(_h, _cap);
  }

  virtual ~NodeCutUnrootedBk()
  {
    lock();
    delete _pG2hRootArc;
    unlock();
  }

protected:
  template<class CBK>
  void separate(CBK& cbk)
  {
    //lemon::Timer t;
    IloNumArray x_values(cbk.getEnv(), _n);
    cbk.getValues(x_values, _x);

    IloNumArray y_values(cbk.getEnv(), _n);
    cbk.getValues(y_values, _y);

    computeCapacities(_cap, x_values, y_values);

    int nCuts = 0;
    int nBackCuts = 0;
    int nNestedCuts = 0;

    //printNonZeroX(x_values);
    // COMMENTED OUT: 22-10-2013
    printNonZeroVars(cbk, _y, y_values);
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
            NodeSet fwdS, fwdDS;
            determineFwdCutSet(_h, *_pBK, fwdDS, fwdS);

            NodeSet bwdS, bwdDS;
            determineBwdCutSet(_h, *_pBK, bwdDS, bwdS);

            // add violated constraints
            for (typename NodeWeightPairVector::const_iterator it2 = it; it2 != compVector.end(); ++it2)
            {
              //const double x_j_value = it2->first;
              const Node j = it2->second;
              assert(_tol.less(minCutValue, it2->first));
              addViolatedConstraint(cbk, j, fwdDS, fwdS);

              nCuts++;
              if (nestedCut)
              {
                nNestedCuts++;
              }
            }

            if (fwdDS.size() != bwdDS.size() || fwdS.size() != bwdS.size() ||
                fwdDS != bwdDS || bwdS != fwdS)
            {
              for (typename NodeWeightPairVector::const_iterator it2 = it; it2 != compVector.end(); ++it2)
              {
                //const double x_j_value = it2->first;
                const Node j = it2->second;
                assert(_tol.less(minCutValue, it2->first));
                addViolatedConstraint(cbk, j, bwdDS, bwdS);
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
    y_values.end();

    //std::cerr << "Generated " << nCuts
    //          << " cuts of which " << nBackCuts << " are back-cuts and "
    //          << nNestedCuts << " are nested cuts" << std::endl;


    // COMMENTED OUT: 22-10-2013
    std::cerr << "[";
    for (int idx = 0; idx < nComp; idx++)
    {
      std::cerr << " " << compMatrix[idx].size();
    }
    std::cerr << " ]" << std::endl;

    //std::cerr << "Time: " << t.realTime() << "s" << std::endl;

  }

  void init()
  {
    // we initialize _h:
    // - for every node i, there will be two nodes i1 and i2
    //   connected by an arc from i1 to i2
    // - for every edge (i,j) there are two arcs
    //   in h: (i2,j1) and (j2,i1) with capacties 1
    _h.clear();
    _diRoot = _h.addNode();
    for (NodeIt i(_g); i != lemon::INVALID; ++i)
    {
      DiNode i1 = _h.addNode();
      DiNode i2 = _h.addNode();
      _pG2h1->set(i, i1);
      _pG2h2->set(i, i2);
      _h2g[i1] = i;
      _h2g[i2] = i;

      DiArc i1i2 = _h.addArc(i1, i2);
      _cap[i1i2] = 0;

      DiArc ri1 = _h.addArc(_diRoot, i1);
      _pG2hRootArc->set(i, ri1);
      _cap[ri1] = 1;
    }

    for (EdgeIt e(_g); e != lemon::INVALID; ++e)
    {
      Node i = _g.u(e);
      Node j = _g.v(e);
      DiNode i1 = (*_pG2h1)[i];
      DiNode i2 = (*_pG2h2)[i];
      DiNode j1 = (*_pG2h1)[j];
      DiNode j2 = (*_pG2h2)[j];

      DiArc i2j1 = _h.addArc(i2, j1);
      DiArc j2i1 = _h.addArc(j2, i1);
      _cap[i2j1] = _cap[j2i1] = 1;
    }
  }

  void computeCapacities(CapacityMap& capacity,
                         IloNumArray x_values,
                         IloNumArray y_values)
  {
    // cap((i,j)) = x_i
    for (NodeIt v(_g); v != lemon::INVALID; ++v)
    {
      double val = x_values[_nodeMap[v]];
      if (!_tol.nonZero(val)) val = 1e-9;
      DiNode v1 = (*_pG2h1)[v];
      capacity[DiOutArcIt(_h, v1)] = val;

      // cap((r,i)) = y_i
      val = y_values[_nodeMap[v]];
      if (!_tol.nonZero(val)) val = 1e-9;
      capacity[(*_pG2hRootArc)[v]] = val;
    }
  }

  template<class CBK>
  void addViolatedConstraint(CBK& cbk, Node target, const NodeSet& dS, const NodeSet& S)
  {
    if (dS.empty() && S.empty())
    {
      //std::cout << getNnodes() << ": " << _x[_nodeMap[target]].getName() << " <= 0" << std::endl;
      addConstraint(_x[_nodeMap[target]] <= 0);
    }
    else
    {
      IloExpr expr(cbk.getEnv());

      //bool first = true;
      //std::cout << getNnodes() << ": " << _x[_nodeMap[target]].getName() << " <=";
      for (NodeSetIt nodeIt = dS.begin(); nodeIt != dS.end(); nodeIt++)
      {
        expr += _x[_nodeMap[*nodeIt]];

        //std::cout << (first ? " " : " + ") << _x[_nodeMap[*nodeIt]].getName();
        //first = false;
      }

      for (NodeSetIt nodeIt = S.begin(); nodeIt != S.end(); nodeIt++)
      {
        expr += _y[_nodeMap[*nodeIt]];

        //std::cout << (first ? " " : " + ") << _y[_nodeMap[*nodeIt]].getName();
        //first = false;
      }

      //std::cout << std::endl;

      IloConstraint constraint = _x[_nodeMap[target]] <= expr;
      addConstraint(constraint);
      constraint.end();

      expr.end();
    }
  }

  void determineBwdCutSet(const Digraph& h,
                          const BkAlg& bk,
                          NodeSet& dS,
                          NodeSet& S)
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
        if (u == _diRoot)
        {
          assert(!_marked[u]);
          //std::cout << _h.id(u) << " -> "
          //          << _h.id(v) << " "
          //          << bk.flow(a) << "/" << bk.cap(a) << std::endl;
          S.insert(_h2g[v]);
        }
        else if (!_marked[u])
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
                          NodeSet& dS,
                          NodeSet& S)
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

      if (v == _diRoot)
      {
        for (DiOutArcIt a(h, v); a != lemon::INVALID; ++a)
        {
          DiNode w = h.target(a);
          if (!_marked[w])
          {
            //std::cout << _h.id(v) << " -> "
            //          << _h.id(w) << " "
            //          << bk.flow(a) << "/" << bk.cap(a) << std::endl;
            S.insert(_h2g[w]);
          }
        }
      }
      else
      {
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
  }
};

} // namespace mwcs
} // namespace nina

#endif // NODECUTUNROOTEDBK_H
