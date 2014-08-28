/*
 * branch.h
 *
 *  Created on: 7-mar-2014
 *      Author: M. El-Kebir
 */

#ifndef BRANCH_H
#define BRANCH_H

#include <ilcplex/ilocplex.h>
#include <ilcplex/ilocplexi.h>
#include <lemon/core.h>
#include <set>

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class Branch : public IloCplex::BranchCallbackI
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;
  typedef NodeCut<GR, NWGHT, NLBL, EWGHT> Parent;
  
protected:
  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  
protected:
  IloBoolVarArray _x;
  const Graph& _g;
  const WeightNodeMap& _weight;
  const Node _root;
  const IntNodeMap& _nodeMap;
  const int _n;
  const IntNodeMap& _deg;
  const IntNodeMap& _posDeg;
  const lemon::Tolerance<double> _tol;

  // 1e-5 is the epsilon that CPLEX uses (for deciding integrality),
  // i.e. if |x| < 1e-5 it's considered to be 0 by CPLEX.
  // if 1 - |x| < 1e-5 it's considered to be 1 by CPLEX.
  // we use the same epsilon to separate violated cuts.
  static const double _epsilon = 1e-5;
  
public:
  Branch(IloEnv env,
         IloBoolVarArray x,
         const Graph& g,
         const WeightNodeMap& weight,
         const IntNodeMap& nodeMap,
         int n,
         const IntNodeMap& deg,
         const IntNodeMap& posDeg)
    : IloCplex::BranchCallbackI(env)
    , _x(x)
    , _g(g)
    , _weight(weight)
    , _nodeMap(nodeMap)
    , _n(n)
    , _deg(deg)
    , _posDeg(posDeg)
    , _tol(_epsilon)
  {
  }
  
  Branch(const Branch& other)
    : IloCplex::BranchCallbackI(other)
    , _x(other._x)
    , _g(other._g)
    , _weight(other._weight)
    , _nodeMap(other._nodeMap)
    , _n(other._n)
    , _deg(other._deg)
    , _posDeg(other._deg)
    , _tol(other._tol)
  {
  }
  
  virtual ~Branch()
  {
  }
  
protected:
  void main()
  {
    if (getBranchType() != BranchOnVariable)
      return;
    
    // do some stuff here
//    std::cout << "Hello, this is node " //  << this->getNodeId()
//              << " and I'm about to create " << this->getNbranches()
//              << " branches" << std::endl;

    IloNumArray x_values(getEnv(), _n);
    getValues(x_values, _x);
    
    // find a fractional node with heighest pos deg
    Node branchOnMe = lemon::INVALID;
    
    std::set<Node> fracNodes;
    for (NodeIt i(_g); i != lemon::INVALID; ++i) {
      double x_i_value = x_values[_nodeMap[i]];
      if (_tol.nonZero(x_i_value) && _tol.different(x_i_value, 1))
      {
        fracNodes.insert(i);
//        std::cout << _x[_nodeMap[i]].getName() << " = " << x_i_value << ", deg: " << _deg[i] << ", posDeg: " << _posDeg[i] <<std::endl;
        
        if (branchOnMe == lemon::INVALID || _posDeg[i] > _posDeg[branchOnMe])
        {
          branchOnMe = i;
        }
      }
    }
    
    if (branchOnMe != lemon::INVALID)
    {
//      branch(branchOnMe);
      std::cout << "Suggestion: " << _x[_nodeMap[branchOnMe]].getName() << ", posDeg: " << _posDeg[branchOnMe] << std::endl;
    }
    else
    {
      std::cout << "Nothing found" << std::endl;
    }
    
//    IloNumVarArray vars(this->getEnv());
//    IloNumArray bounds(this->getEnv());
//    IloCplex::BranchDirectionArray dirs(this->getEnv());
//    
//    for (int i = 0; i < this->getNbranches(); ++i)
//    {
//      vars.clear();
//      bounds.clear();
//      std::cout << "Branch " << i << std::endl;
//      
//      this->getBranch(vars, bounds, dirs, i);
//      for (int j = 0; j < vars.getSize(); ++j)
//      {
//        std::cout << vars[j].getName() << " = " << bounds[j] << ", frac value: " << getValue(vars[j]) << std::endl;
//      }
//    }
//    
//    vars.end();
//    bounds.end();
//    dirs.end();
//    x_values.end();
  }
  
  IloCplex::CallbackI* duplicateCallback() const
  {
    return (new (getEnv()) Branch(*this));
  }
  
  void branch(Node i)
  {
    double obj_value = getObjValue();
    // up
    IloConstraintArray cons(getEnv());
    cons.add(_x[_nodeMap[i]] == 1);
    for (IncEdgeIt e(_g, i); e != lemon::INVALID; ++e)
    {
      Node j = _g.oppositeNode(i, e);
      if (_weight[j] >= 0)
      {
        cons.add(_x[_nodeMap[j]] == 1);
      }
    }
    makeBranch(cons, obj_value);
    
    // down
    makeBranch(_x[_nodeMap[i]] == 0, obj_value);
  }
};
  
} // namespace mwcs
} // namespace nina

#endif // BRANCH_H
