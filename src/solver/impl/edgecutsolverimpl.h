/*
 * edgecutsolverimpl.h
 *
 *  Created on: 09-sep-2014
 *     Authors: M. El-Kebir
 */

#ifndef EDGECUTSOLVERIMPL_H
#define EDGECUTSOLVERIMPL_H

#include "cplexsolverimpl.h"

#include <lemon/smart_graph.h>

// ILOG stuff
#include <ilconcert/iloalg.h>
#include <ilcplex/ilocplex.h>

#include <string>

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class EdgeCutSolverImpl : public CplexSolverImpl<GR, NWGHT, NLBL, EWGHT>
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;
  
  typedef CplexSolverImpl<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> Parent;
  
  typedef typename Parent::MwcsGraphType MwcsGraphType;
  typedef typename Parent::MwcsAnalyzeType MwcsAnalyzeType;
  
  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  
  typedef typename Parent::NodeSet NodeSet;
  typedef typename Parent::NodeSetIt NodeSetIt;
  typedef typename Parent::NodeVector NodeVector;
  typedef typename Parent::NodeVectorIt NodeVectorIt;
  typedef typename Parent::InvNodeIntMap InvNodeIntMap;
  typedef typename Parent::InvArcIntMap InvArcIntMap;
  typedef typename Parent::Options Options;
  
  using Parent::_options;
  using Parent::_pAnalysis;
  using Parent::_n;
  using Parent::_m;
  using Parent::_env;
  using Parent::_model;
  using Parent::_cplex;
  using Parent::_x;
  using Parent::initVariables;
  using Parent::initConstraints;
  using Parent::clean;
  using Parent::solveCplex;
  using Parent::exportModel;
  
protected:
  EdgeCutSolverImpl(const Options& options)
    : Parent(options)
    , _z()
    , _d()
    , _pToD(NULL)
    , _toG(_d)
    , _weight(_d)
    , _arcCost(0)
    , _label(_d)
    , _diNode(_d)
    , _diArc(_d)
    , _invDiNode()
    , _invDiArc()
  {
  }
  
  virtual ~EdgeCutSolverImpl()
  {
    delete _pToD;
  }
  
  virtual void initGraph(const MwcsGraphType& mwcsGraph) = 0;
  virtual void initVariables(const MwcsGraphType& mwcsGraph);
  virtual void initConstraints(const MwcsGraphType& mwcsGraph);
  
  virtual void extractSolution(double& score,
                               BoolNodeMap& solutionMap,
                               NodeSet& solutionSet)
  {
    lemon::Tolerance<double> tol(1e-5);
    // solution
    solutionSet.clear();
    for (int i = 0; i < _n ; i++)
    {
      Node node = _toG[_invDiNode[i]];
      if (node != lemon::INVALID)
      {
        bool inSol = tol.nonZero(_cplex.getValue(_x[i]));
        solutionMap[node] = inSol;
        if (inSol)
        {
          solutionSet.insert(node);
        }
      }
    }
    
    score = _cplex.getObjValue();
  }
  
public:
  virtual void printVariables(const MwcsGraphType&, std::ostream& out)
  {
    for (int id_v = 0; id_v < _n; ++id_v)
    {
      out << "Node '" << _label[_invDiNode[id_v]]
          << "', x_" << id_v << " = " << _cplex.getValue(_x[id_v])
          << ", w = " << _weight[_invDiNode[id_v]]
          << std::endl;
    }

    for (int id_a = 0; id_a < _m; ++id_a)
    {
      DiArc a = _invDiArc[id_a];
      DiNode u = _d.source(a);
      DiNode v = _d.target(a);
      
      out << "Arc '(" << _diNode[u] << ", " << _diNode[v] << ")"
          << "', z_" << id_a << " = " << _cplex.getValue(_z[id_a])
          << std::endl;
    }
  }
  
protected:
  typedef lemon::SmartDigraph DiGraph;
  typedef typename DiGraph::Node DiNode;
  typedef typename DiGraph::NodeIt DiNodeIt;
  typedef typename DiGraph::Arc DiArc;
  typedef typename DiGraph::ArcIt DiArcIt;
  typedef typename DiGraph::InArcIt DiInArcIt;
  typedef typename DiGraph::OutArcIt DiOutArcIt;
  
  typedef std::set<DiNode> DiNodeSet;
  typedef std::vector<DiNode> DiNodeVector;
  typedef typename DiNodeVector::const_iterator DiNodeVectorIt;
  typedef std::vector<DiArc> DiArcVector;
  typedef typename DiArcVector::const_iterator DiArcVectorIt;
  typedef std::vector<DiNode> InvDiNodeIntMap;
  typedef std::vector<DiArc> InvDiArcIntMap;
  typedef DiGraph::template NodeMap<int> IntDiNodeMap;
  typedef DiGraph::template ArcMap<int> IntDiArcMap;
  typedef DiGraph::template NodeMap<int> DoubleDiNodeMap;
  typedef DiGraph::template NodeMap<std::string> LabelDiNodeMap;
  typedef DiGraph::template NodeMap<Node> NodeDiNodeMap;
  typedef typename Graph::template NodeMap<DiNode> DiNodeNodeMap;
  
  IloBoolVarArray _z;
  DiGraph _d;
  DiNodeNodeMap* _pToD;
  NodeDiNodeMap _toG;
  DoubleDiNodeMap _weight;
  double _arcCost;
  LabelDiNodeMap _label;
  IntDiNodeMap _diNode;
  IntDiArcMap _diArc;
  InvDiNodeIntMap _invDiNode;
  InvDiArcIntMap _invDiArc;
  
  struct DiNodesDegComp
  {
  private:
    const IntDiNodeMap& _deg;
    
  public:
    DiNodesDegComp(const IntDiNodeMap& deg)
      : _deg(deg)
    {
    }
    
    bool operator ()(DiNode u, DiNode v)
    {
      return _deg[u] > _deg[v];
    }
  };
  
  struct DiNodesWeightComp
  {
  private:
    const DoubleDiNodeMap& _weight;
    
  public:
    DiNodesWeightComp(const DoubleDiNodeMap& weight)
      : _weight(weight)
    {
    }
    
    bool operator ()(DiNode u, DiNode v)
    {
      return _weight[u] > _weight[v];
    }
  };
};

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void EdgeCutSolverImpl<GR, NWGHT, NLBL, EWGHT>::initVariables(const MwcsGraphType&)
{
  _x = IloBoolVarArray(_env, _n);
  
  // let's sort the nodes on degree
  IntDiNodeMap deg(_d);
  _invDiNode.clear();
  _invDiNode.reserve(_n);
  for (DiNodeIt v(_d); v != lemon::INVALID; ++v)
  {
    int d = 0;
    for (DiOutArcIt a(_d, v); a != lemon::INVALID; ++a) ++d;
    deg[v] = d;
    _invDiNode.push_back(v);
  }
  
  DiNodesDegComp comp(deg);
  std::sort(_invDiNode.begin(), _invDiNode.end(), comp);
  
  char buf[1024];
  
  int i = 0;
  for (DiNodeVectorIt it = _invDiNode.begin(); it != _invDiNode.end(); ++it, ++i)
  {
    DiNode v = *it;
    
    // x_i = 0 if node i is not in the subgraph
    // x_i = 1 if node i is the subgraph
    snprintf(buf, 1024, "x_%s", _label[v].c_str());
    //    snprintf(buf, 1024, "x_%d", g.id(v));
    _x[i].setName(buf);
    
    _diNode[v] = i;
  }
  
  _z = IloBoolVarArray(_env, _m);
  
  _invDiArc.clear();
  _invDiArc.reserve(_m);
  for (DiArcIt a(_d); a != lemon::INVALID; ++a)
  {
    _invDiArc.push_back(a);
  }
  
  i = 0;
  for (DiArcVectorIt it = _invDiArc.begin(); it != _invDiArc.end(); ++it, ++i)
  {
    DiArc a = *it;
    DiNode u = _d.source(a);
    DiNode v = _d.target(a);
    
    snprintf(buf, 1024, "z_%s_%s", _label[u].c_str(), _label[v].c_str());
    _z[i].setName(buf);
    
    _diArc[a] = i;
  }
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void EdgeCutSolverImpl<GR, NWGHT, NLBL, EWGHT>::initConstraints(const MwcsGraphType& mwcsGraph)
{
  assert(_pToD);
  
  const Graph& g = mwcsGraph.getGraph();
  
  // objective function
  IloExpr expr(_env);
  for (int i = 0; i < _n ; i++)
  {
    expr += _x[i] * _weight[_invDiNode[i]];
  }
  for (int i = 0; i < _m ; i++)
  {
    expr += _z[i] * _arcCost;
  }
  expr += _arcCost;
  
  _model.add(IloObjective(_env, expr, IloObjective::Maximize));

  // add equality constraints
  if (_pAnalysis)
  {
    int nAnalyzeConstraints = 0;
    for (NodeIt i(g); i != lemon::INVALID; ++i)
    {
      if (mwcsGraph.getScore(i) > 0)
      {
        for (NodeIt j(g); j != lemon::INVALID; ++j)
        {
          if (mwcsGraph.getScore(j) > 0 && _pAnalysis->ok(i, j))
          {
            _model.add(_x[_diNode[(*_pToD)[i]]] <= _x[_diNode[(*_pToD)[j]]]);
            nAnalyzeConstraints++;
          }
        }
      }
    }
    
    if (g_verbosity >= VERBOSE_ESSENTIAL)
    {
      std::cout << "// Added " << nAnalyzeConstraints
                << " analyze constraints"
                << std::endl;
    }
  }
}

} // namespace mwcs
} // namespace nina

#endif // EDGECUTSOLVERIMPL_H