/*
 * edgecut.h
 *
 *  Created on: 10-sep-2014
 *      Author: M. El-Kebir
 */

#ifndef EDGECUT_H
#define EDGECUT_H

#include <ilcplex/ilocplex.h>
#include <ilcplex/ilocplexi.h>
#include <ilconcert/ilothread.h>
#include <lemon/adaptors.h>
#include <lemon/tolerance.h>
#include <set>
#include <queue>
#include <list>

namespace nina {
namespace mwcs {

template<typename DGR,
         typename NWGHT = typename DGR::template NodeMap<double> >
class EdgeCut
{
public:
  typedef DGR Graph;
  typedef NWGHT WeightNodeMap;

protected:
  TEMPLATE_DIGRAPH_TYPEDEFS(Graph);
  typedef std::set<Node> NodeSet;
  typedef typename NodeSet::const_iterator NodeSetIt;
  typedef std::set<Arc> ArcSet;
  typedef typename ArcSet::const_iterator ArcSetIt;
  typedef std::vector<NodeSet> NodeSetVector;
  typedef typename NodeSetVector::const_iterator NodeSetVectorIt;
  typedef lemon::FilterNodes<const Graph, const BoolNodeMap> SubGraph;
  typedef typename SubGraph::NodeIt SubNodeIt;
  typedef typename SubGraph::ArcIt SubArcIt;
  typedef std::queue<Node> NodeQueue;
  typedef std::list<Node> NodeList;
  typedef typename NodeList::const_iterator NodeListIt;

protected:
  IloBoolVarArray _x;
  IloBoolVarArray _z;
  const Graph& _d;
  const WeightNodeMap& _weight;
  const IntNodeMap& _nodeMap;
  const IntArcMap& _arcMap;
  const int _n;
  const int _m;
  const int _maxNumberOfCuts;
  const lemon::Tolerance<double> _tol;
  BoolNodeMap* _pNodeBoolMap;
  const SubGraph* _pSubG;
  IntNodeMap* _pComp;
  IloFastMutex* _pMutex;

  // 1e-5 is the epsilon that CPLEX uses (for deciding integrality),
  // i.e. if |x| < 1e-5 it's considered to be 0 by CPLEX.
  // if 1 - |x| < 1e-5 it's considered to be 1 by CPLEX.
  // we use the same epsilon to separate violated cuts.
  static const double _epsilon = 1e-5;

public:
  EdgeCut(IloBoolVarArray x,
          IloBoolVarArray z,
          const Graph& d,
          const WeightNodeMap& weight,
          const IntNodeMap& nodeMap,
          const IntArcMap& arcMap,
          int n,
          int m,
          int maxNumberOfCuts,
          IloFastMutex* pMutex)
    : _x(x)
    , _z(z)
    , _d(d)
    , _weight(weight)
    , _nodeMap(nodeMap)
    , _arcMap(arcMap)
    , _n(n)
    , _m(m)
    , _maxNumberOfCuts(maxNumberOfCuts)
    , _tol(_epsilon)
    , _pNodeBoolMap(NULL)
    , _pSubG(NULL)
    , _pComp(NULL)
    , _pMutex(pMutex)
  {
    lock();
    _pNodeBoolMap = new BoolNodeMap(_d);
    _pSubG = new SubGraph(_d, *_pNodeBoolMap);
    _pComp = new IntNodeMap(_d);
    unlock();
  }

  EdgeCut(const EdgeCut& other)
    : _x(other._x)
    , _z(other._z)
    , _d(other._d)
    , _weight(other._weight)
    , _nodeMap(other._nodeMap)
    , _arcMap(other._arcMap)
    , _n(other._n)
    , _m(other._m)
    , _maxNumberOfCuts(other._maxNumberOfCuts)
    , _tol(other._tol)
    , _pNodeBoolMap(NULL)
    , _pSubG(NULL)
    , _pComp(NULL)
    , _pMutex(other._pMutex)
  {
    lock();
    _pNodeBoolMap = new BoolNodeMap(_d);
    _pSubG = new SubGraph(_d, *_pNodeBoolMap);
    _pComp = new IntNodeMap(_d);
    unlock();
  }

  virtual ~EdgeCut()
  {
    lock();
    delete _pNodeBoolMap;
    delete _pSubG;
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
  
  NodeSetVector determineStronglyConnectedComponents(const IloNumArray& x_values)
  {
    // update _pSubG
    for (NodeIt v(_d); v != lemon::INVALID; ++v)
    {
      double val = x_values[_nodeMap[v]];
      _pNodeBoolMap->set(v, _tol.nonZero(val));
    }
    
    int nComp = lemon::stronglyConnectedComponents(*_pSubG, *_pComp);
    
    NodeSetVector nonZeroComponents(nComp, NodeSet());
    for (SubNodeIt i(*_pSubG); i != lemon::INVALID; ++i)
    {
      int compIdx = (*_pComp)[i];
      nonZeroComponents[compIdx].insert(i);
    }
    
    return nonZeroComponents;
  }
  
  template<typename CBK>
  void separateStronglyConnectedComponent(const NodeSet& S,
                                          const NodeSet rootNodes,
                                          const IloNumArray& x_values,
                                          const IloNumArray& z_values,
                                          CBK& cbk,
                                          int& nCuts)
  {
#ifdef DEBUG
    NodeSet intersection;
    std::set_intersection(rootNodes.begin(), rootNodes.end(),
                          S.begin(), S.end(),
                          std::inserter(intersection, intersection.begin()));
    assert(intersection.size() == 0);
#endif

    IloExpr rhs(cbk.getEnv());
    
    // determine dS
    ArcSet dS;
    for (NodeSetIt it = S.begin(); it != S.end(); ++it)
    {
      const Node i = *it;
      for (InArcIt a(_d, i); a != lemon::INVALID; ++a)
      {
        const Node j = _d.source(a);
        if (S.find(j) == S.end())
        {
          bool val_a = _tol.nonZero(z_values[_arcMap[a]]);
          if (!val_a)
          {
            dS.insert(a);
          }
          else
          {
            for (InArcIt a2(_d, j); a2 != lemon::INVALID; ++a2)
            {
              bool val_a2 = _tol.nonZero(z_values[_arcMap[a2]]);
              assert(!val_a2);
              dS.insert(a2);
            }
          }
        }
      }
    }
    
    constructRHS(rhs, dS);
//    std::cout << std::endl << "dS" << std::endl;
//    printArcSet(dS, _z, z_values);
//    std::cout << "S" << std::endl;
//    printNodeSet(S, _x, x_values);
//    std::cout << std::endl << std::endl;
    for (NodeSetIt it = S.begin(); it != S.end(); ++it)
    {
      assert(isValid(*it, dS));
      cbk.add(_x[_nodeMap[*it]] <= rhs, IloCplex::UseCutPurge).end();
      ++nCuts;
    }

    rhs.end();
  }
  
  template<typename CBK>
  void separateRootedStronglyConnectedComponent(const NodeSet& S,
                                                const Node root,
                                                const IloNumArray& x_values,
                                                CBK& cbk,
                                                int& nCuts)
  {
    assert(S.find(root) == S.end());

    IloExpr rhs(cbk.getEnv());
    
    // determine dS
    NodeSet dS;
    for (NodeSetIt it = S.begin(); it != S.end(); ++it)
    {
      const Node i = *it;
      for (OutArcIt a(_d, i); a != lemon::INVALID; ++a)
      {
        const Node j = _d.target(a);
        if (S.find(j) == S.end())
        {
          dS.insert(j);
        }
      }
    }
    
    constructRHS(rhs, dS);
    for (NodeSetIt it = S.begin(); it != S.end(); ++it)
    {
      assert(isValid(*it, dS, S));
      cbk.add(_x[_nodeMap[*it]] <= rhs, IloCplex::UseCutPurge).end();
      ++nCuts;
    }
    
    rhs.end();
  }

  void printNonZeroVars(IloCplex::ControlCallbackI& cbk,
                        IloBoolVarArray variables,
                        IloNumArray values) const
  {
    std::cerr << cbk.getNnodes() << ":";
    for (NodeIt i(_d); i != lemon::INVALID; ++i)
    {
      double i_value = values[_nodeMap[i]];
      if (!_tol.nonZero(i_value)) continue;
      std::cerr << " " << variables[_nodeMap[i]].getName()
                << " (" << _d.id(i) << ", " << _weight[i] << ", " << i_value << ") " ;

      if (cbk.getDirection(variables[_nodeMap[i]]) == CPX_BRANCH_UP)
        std::cerr << "*";
      std::cerr << std::endl;
    }
    std::cerr << std::endl;
  }
  
  void printArcSet(const ArcSet& arcs,
                   IloBoolVarArray variables,
                   IloNumArray values) const
  {
    std::cout.precision(std::numeric_limits<double>::digits10);
    bool first = true;
    for (ArcSetIt it = arcs.begin(); it != arcs.end(); it++)
    {
      if (!first)
        std::cout << " ";
      else
        first = false;
      
      std::cout << _d.id(*it) << "(" << variables[_arcMap[*it]].getName()
                << " = " << std::fixed << values[_arcMap[*it]]
                << (_tol.nonZero(values[_arcMap[*it]]) ? "*" : "") << ")";
    }
    std::cout << std::endl;
  }
  
  void printArcSet(const ArcSet& arcs,
                    IloBoolVarArray variables) const
  {
    bool first = true;
    for (ArcSetIt it = arcs.begin(); it != arcs.end(); it++)
    {
      if (!first)
        std::cout << " ";
      else
        first = false;
      
      std::cout << _d.id(*it) << "(" << variables[_arcMap[*it]].getName() << ")";
    }
    std::cout << std::endl;
  }

  void printNodeSet(const NodeSet& nodes,
                    IloBoolVarArray variables,
                    IloNumArray values) const
  {
    std::cout.precision(std::numeric_limits<double>::digits10);
    bool first = true;
    for (NodeSetIt it = nodes.begin(); it != nodes.end(); it++)
    {
      if (!first)
        std::cout << " ";
      else
        first = false;

      std::cout << _d.id(*it) << "(" << variables[_nodeMap[*it]].getName()
                << " = " << std::fixed << values[_nodeMap[*it]]
                << (_tol.nonZero(values[_nodeMap[*it]]) ? "*" : "") << ")";
    }
    std::cout << std::endl;
  }
  
  void printNodeSet(const NodeSet& nodes,
                    IloBoolVarArray variables) const
  {
    bool first = true;
    for (NodeSetIt it = nodes.begin(); it != nodes.end(); it++)
    {
      if (!first)
        std::cout << " ";
      else
        first = false;
      
      std::cout << _d.id(*it) << "(" << variables[_nodeMap[*it]].getName() << ")";
    }
    std::cout << std::endl;
  }
  
  template<typename CBK>
  void addViolatedConstraint(CBK& cbk,
                             Node target,
                             const NodeSet& dS)
  {
    assert(isValid(target, dS));
    if (dS.empty())
    {
      //std::cout << cbk.getNnodes() << ": " << _x[_nodeMap[target]].getName() << " <= 0" << std::endl;
      cbk.add(_x[_nodeMap[target]] <= 0);
      // there should only be one component!
      assert(false);
    }
    else
    {
      IloExpr expr(cbk.getEnv());
      
      //bool first = true;
      //std::cout << cbk.getNnodes() << ": " << _x[_nodeMap[target]].getName() << " <=";
      for (NodeSetIt nodeIt = dS.begin(); nodeIt != dS.end(); nodeIt++)
      {
        expr += _x[_nodeMap[*nodeIt]];
        
        //std::cout << (first ? " " : " + ") << _x[_nodeMap[*nodeIt]].getName();
        //first = false;
      }
      
      IloConstraint constraint = _x[_nodeMap[target]] <= expr;
      cbk.add(constraint, IloCplex::UseCutPurge);
      constraint.end();
      
      expr.end();
    }
  }
  
  bool isValid(Node target, const ArcSet& dS) const
  {
    return true;
  }
  
  template<typename CBK>
  void addViolatedConstraint(CBK& cbk,
                             Node target,
                             const NodeSet& dS,
                             const NodeSet& S)
  {
    assert(isValid(target, dS, S));
    if (dS.empty() && S.empty())
    {
      //std::cout << cbk.getNnodes() << ": " << _x[_nodeMap[target]].getName() << " <= 0" << std::endl;
      cbk.add(_x[_nodeMap[target]] <= 0);
      // target should be in S!
      assert(false);
    }
    else
    {
      IloExpr expr(cbk.getEnv());
      
      //bool first = true;
      //std::cout << cbk.getNnodes() << ": " << _x[_nodeMap[target]].getName() << " <=";
      for (NodeSetIt nodeIt = dS.begin(); nodeIt != dS.end(); nodeIt++)
      {
        expr += _x[_nodeMap[*nodeIt]];
        
        //std::cout << (first ? " " : " + ") << _x[_nodeMap[*nodeIt]].getName();
        //first = false;
      }
      
//      for (NodeSetIt nodeIt = S.begin(); nodeIt != S.end(); nodeIt++)
//      {
//        expr += _y[_nodeMap[*nodeIt]];
//        
//        //std::cout << (first ? " " : " + ") << _y[_nodeMap[*nodeIt]].getName();
//        //first = false;
//      }
      
      //std::cout << std::endl;
      
      IloConstraint constraint = _x[_nodeMap[target]] <= expr;
      cbk.add(constraint, IloCplex::UseCutPurge);
      constraint.end();
      
      expr.end();
    }
  }

  void constructRHS(IloExpr& rhs,
                    const ArcSet& dS)
  {
    rhs.clear();
    for (ArcSetIt arcIt = dS.begin(); arcIt != dS.end(); ++arcIt)
    {
      rhs += _z[_arcMap[*arcIt]];
    }
  }
  
//  bool isValid(Node target,
//               const ArcSet& dS,
//               const NodeSet& S) const
//  {
//    ArcSet ddS;
//    for (NodeSetIt nodeIt = S.begin(); nodeIt != S.end(); ++nodeIt)
//    {
//      Node u = *nodeIt;
//      for (InArcIt a(_d, u); a != lemon::INVALID; ++a)
//      {
//        Node v = _d.source(a);
//        if (S.find(v) == S.end())
//        {
//          ddS.insert(a);
//        }
//      }
//    }
//    
//    // let's do a bfs from target
//    NodeSet SS;
//    BoolNodeMap visited(_d, false);
//
//    NodeQueue Q;
//    Q.push(target);
//    while (!Q.empty())
//    {
//      Node v = Q.front();
//      Q.pop();
//      visited[v] = true;
//      SS.insert(v);
//      
////      for (IncEdgeIt e(_d, v); e != lemon::INVALID; ++e)
////      {
////        Node u = _d.oppositeNode(v, e);
////        if (!visited[u] && dS.find(u) == dS.end())
////        {
////          Q.push(u);
////        }
////      }
//    }
//    
//    if (S != SS)
//    {
//      std::cerr << std::endl << "SS" << std::endl;
//      printNodeSet(SS, _y);
//      return false;
//    }
//    
//    if (dS != ddS)
//    {
//      std::cerr << std::endl << "ddS" << std::endl;
//      printArcSet(ddS, _z);
//      return false;
//    }
//    else
//    {
//      // target must be in S
//      if (S.find(target) == S.end())
//      {
//        return false;
//      }
//      
//      // but must not be in dS
//      for (ArcSetIt arcIt = dS.begin(); arcIt != dS.end(); ++arcIt)
//      {
//        if (target == _d.source(*arcIt))
//        {
//          return false;
//        }
//      }
//      
//      return true;
//    }
//  }
};

} // namespace mwcs
} // namespace nina

#endif // EDGECUT_H
