/*
 * nodecut.h
 *
 *  Created on: 17-feb-2014
 *      Author: M. El-Kebir
 */

#ifndef NODECUT_H
#define NODECUT_H

#include <ilcplex/ilocplex.h>
#include <ilcplex/ilocplexi.h>
#include <ilconcert/ilothread.h>
#include <lemon/tolerance.h>
#include <set>

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class NodeCut
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;

  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  typedef std::vector<Node> NodeVector;
  typedef typename NodeVector::const_iterator NodeVectorIt;
  typedef std::vector<NodeVector> NodeMatrix;
  typedef std::set<Node> NodeSet;
  typedef typename NodeSet::const_iterator NodeSetIt;

protected:
  IloBoolVarArray _x;
  IloBoolVarArray _y;
  const Graph& _g;
  const WeightNodeMap& _weight;
  const Node _root;
  const IntNodeMap& _nodeMap;
  const int _n;
  const int _m;
  const int _maxNumberOfCuts;
  const lemon::Tolerance<double> _tol;
  BoolNodeMap* _pNodeBoolMap;
  IloFastMutex* _pMutex;

  // 1e-5 is the epsilon that CPLEX uses (for deciding integrality),
  // i.e. if |x| < 1e-5 it's considered to be 0 by CPLEX.
  // if 1 - |x| < 1e-5 it's considered to be 1 by CPLEX.
  // we use the same epsilon to separate violated cuts.
  static const double _epsilon = 1e-5;

public:
  NodeCut(IloBoolVarArray x,
          IloBoolVarArray y,
          const Graph& g,
          const WeightNodeMap& weight,
          Node root,
          const IntNodeMap& nodeMap,
          int n,
          int m,
          int maxNumberOfCuts,
          IloFastMutex* pMutex)
    : _x(x)
    , _y(y)
    , _g(g)
    , _weight(weight)
    , _root(root)
    , _nodeMap(nodeMap)
    , _n(n)
    , _m(m)
    , _maxNumberOfCuts(maxNumberOfCuts)
    , _tol(_epsilon)
    , _pNodeBoolMap(NULL)
    , _pMutex(pMutex)
  {
    lock();
    _pNodeBoolMap = new BoolNodeMap(_g);
    unlock();
  }

  NodeCut(const NodeCut& other)
    : _x(other._x)
    , _y(other._y)
    , _g(other._g)
    , _weight(other._weight)
    , _root(other._root)
    , _nodeMap(other._nodeMap)
    , _n(other._n)
    , _m(other._m)
    , _maxNumberOfCuts(other._maxNumberOfCuts)
    , _tol(other._tol)
    , _pNodeBoolMap(NULL)
    , _pMutex(other._pMutex)
  {
    lock();
    _pNodeBoolMap = new BoolNodeMap(_g);
    unlock();
  }

  virtual ~NodeCut()
  {
    lock();
    delete _pNodeBoolMap;
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
    std::cout.precision(std::numeric_limits<double>::digits10);
    bool first = true;
    for (NodeSetIt it = nodes.begin(); it != nodes.end(); it++)
    {
      if (!first)
        std::cout << " ";
      else
        first = false;

      std::cout << _g.id(*it) << "(" << variables[_nodeMap[*it]].getName()
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
      
      std::cout << _g.id(*it) << "(" << variables[_nodeMap[*it]].getName() << ")";
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
      cbk.add(constraint);
      constraint.end();
      
      expr.end();
    }
  }
  
  bool isValid(Node target, const NodeSet& dS) const
  {
    return dS.find(target) == dS.end();
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
      
      for (NodeSetIt nodeIt = S.begin(); nodeIt != S.end(); nodeIt++)
      {
        expr += _y[_nodeMap[*nodeIt]];
        
        //std::cout << (first ? " " : " + ") << _y[_nodeMap[*nodeIt]].getName();
        //first = false;
      }
      
      //std::cout << std::endl;
      
      IloConstraint constraint = _x[_nodeMap[target]] <= expr;
      cbk.add(constraint);
      constraint.end();
      
      expr.end();
    }
  }

  bool isValid(Node target,
               const NodeSet& dS,
               const NodeSet& S) const
  {
    NodeSet ddS;
    for (NodeSetIt nodeIt = S.begin(); nodeIt != S.end(); ++nodeIt)
    {
      Node u = *nodeIt;
      for (OutArcIt a(_g, u); a != lemon::INVALID; ++a)
      {
        Node v = _g.target(a);
        if (S.find(v) == S.end())
        {
          ddS.insert(v);
        }
      }
    }
    
    if (dS != ddS)
    {
      return false;
    }
    else
    {
      // target must be in S
      if (S.find(target) == S.end())
      {
        return false;
      }
      
      // but must not be in dS
      if (dS.find(target) != dS.end())
      {
        return false;
      }
      
      return true;
    }
  }
};

} // namespace mwcs
} // namespace nina


#endif // NODECUT_H
