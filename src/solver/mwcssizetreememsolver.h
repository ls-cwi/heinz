/*
 * mwcssizetreememsolver.h
 *
 *  Created on: 3-jul-2013
 *      Author: M. El-Kebir
 */

#ifndef MWCSSIZETREEMEMSOLVER_H
#define MWCSSIZETREEMEMSOLVER_H

#include "mwcssolver.h"
#include <set>
#include <map>
#include <tr1/unordered_map>
#include <limits>
#include <vector>
#include <assert.h>
#include <lemon/bfs.h>
#include <lemon/connectivity.h>

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class MwcsSizeTreeMemSolver : public MwcsSolver<GR, NWGHT, NLBL, EWGHT>
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;

  typedef MwcsSolver<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> Parent;
  typedef typename Parent::MwcsGraphType MwcsGraphType;

  TEMPLATE_GRAPH_TYPEDEFS(Graph);

  using Parent::_mwcsGraph;
  using Parent::_root;
  using Parent::_score;
  using Parent::_solutionMap;
  using Parent::_solutionSet;
  using Parent::printSolution;
  using Parent::getSolutionWeight;
  using Parent::getSolutionNodeMap;
  using Parent::getSolutionModule;
  using Parent::isNodeInSolution;
  using Parent::init;

private:
  typedef std::set<Node> NodeSet;
  typedef typename NodeSet::const_iterator NodeSetIt;

  typedef std::vector<Node> NodeVector;
  typedef typename NodeVector::const_iterator NodeVectorIt;

  typedef typename Graph::template NodeMap<NodeVector> NodeVectorMap;
  typedef typename lemon::Bfs<Graph>::PredMap PredMap;
  typedef typename lemon::Bfs<Graph>::DistMap DistMap;

  typedef std::tr1::unordered_map<int, double> DoubleDict;

  typedef std::vector<int> IntVector;

  typedef IntVector::const_iterator IntVectorIt;
  typedef IntVector::const_reverse_iterator IntVectorRevIt;

  typedef std::tr1::unordered_map<int, IntVector> IntVectorDict;

  typedef std::vector<NodeVector> NodeMatrix;
  typedef typename NodeMatrix::const_iterator NodeMatrixIt;
  typedef typename NodeMatrix::const_reverse_iterator NodeMatrixRevIt;

  struct DpEntry {
    NodeVector _children;
    IntVectorDict _perm;
    DoubleDict _weight;

    DpEntry()
      : _children()
      , _perm()
      , _weight()
    {
    }

    void clear()
    {
      _weight.clear();// = -std::numeric_limits<double>::max();
      _perm.clear();
      _children.clear();
    }
  };

  typedef std::tr1::unordered_map<int, DpEntry> DpEntryDict;

public:
  MwcsSizeTreeMemSolver(const MwcsGraphType& mwcsGraph, int k);
  virtual void init(Node root);
  virtual bool solve();
  double solve(const Node v, const int k);
  virtual void setLowerBound(double LB) {}

private:
  const int _k;
  DpEntryDict _dpMap;

  bool visited(const Node v) const
  {
    const Graph& g = _mwcsGraph.getGraph();
    return _dpMap.find(g.id(v)) != _dpMap.end();
  }

  bool visited(const Node v, const int k) const
  {
    if (!visited(v))
      return false;

    const Graph& g = _mwcsGraph.getGraph();
    const DpEntry& entry = _dpMap.find(g.id(v))->second;

    return entry._weight.find(k) != entry._weight.end();
  }

  IntVector getInitSplit(const DpEntry& entry, const Node v) const
  {
    const int n = static_cast<int>(entry._children.size() - 1);
    return IntVector(n, 0);
  }

  bool nextSplit(const Node v, const int k, IntVector& split) const
  {
    const int n = static_cast<int>(split.size());
    assert(visited(v));
    assert(n == static_cast<int>(_dpMap.find(_mwcsGraph.getGraph().id(v))->second._children.size()) - 1);

    for (int i = n - 1; i >= 0; i--)
    {
      if (split[i] < k)
      {
        split[i]++;

        // set everything to the right to split[i]
        for (int j = i + 1; j < n; j++) split[j] = split[i];

        return true;
      }
    }
    return false;
  }

  IntVector convertSplit(const Node v, const int k, const IntVector& split) const
  {
    const int n = static_cast<int>(split.size());
    IntVector result(n+1, 0);

    if (n >= 1)
    {
      result[0] = split[0];
      for (int i = 1; i < n; i++)
      {
        result[i] = split[i] - split[i-1];
      }
      result[n] = k - split[n-1];
    }
    else
    {
      result[0] = k;
    }

    return result;
  }

  double getPermWeight(const DpEntry& entry, const Node v, const IntVector& perm)
  {
    double res = 0;

    const int n = static_cast<int>(perm.size());
    for (int i = 0; i < n; i++)
    {
      res += _dpMap[entry._children[i]]._weight[perm[i]-1];
    }

    return res;
  }

  void reconstruct(const Node v, const int k)
  {
    const Graph& g = _mwcsGraph.getGraph();
    DpEntry& entry = _dpMap[g.id(v)];
    const IntVector& perm = entry._perm[k];

    _solutionSet.insert(v);
    _solutionMap[v] = true;

    const int n = static_cast<int>(entry._children.size());
    for (int i = 0; i < n; i++)
    {
      Node u = entry._children[i];
      if (perm[i] > 0)
      {
        reconstruct(u, perm[i]);
      }
    }
  }
};

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline MwcsSizeTreeMemSolver<GR, NWGHT, NLBL, EWGHT>::MwcsSizeTreeMemSolver(const MwcsGraphType& mwcsGraph,
                                                                            int k)
  : Parent(mwcsGraph)
  , _k(k)
  , _dpMap(k)
{
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsSizeTreeMemSolver<GR, NWGHT, NLBL, EWGHT>::init(Node root)
{
  const Graph& g = _mwcsGraph.getGraph();

  _score = 0;
  _solutionSet.clear();

  lemon::mapFill(g, _solutionMap, false);
  _dpMap.clear();

  _root = root;

  assert(lemon::acyclic(g));
  assert(root != lemon::INVALID);
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline double MwcsSizeTreeMemSolver<GR, NWGHT, NLBL, EWGHT>::solve(const Node v, const int k)
{
  const Graph& g = _mwcsGraph.getGraph();
  const WeightNodeMap& weight = _mwcsGraph.getScores();

  assert(k >= 1);

  // check if already computed
  if (visited(v, k))
  {
    return _dpMap[g.id(v)]._weight[k];
  }

  // init children if necessary
  if (!visited(v))
  {
    DpEntry& entry = _dpMap[g.id(v)];

    for (IncEdgeIt e(g, v); e != lemon::INVALID; ++e)
    {
      Node u = g.oppositeNode(v, e);
      if (!visited(u))
      {
        entry._children.push_back(u);
      }
    }
  }

  DpEntry& entry = _dpMap[g.id(v)];

  // set k = 1
  if (k == 1)
  {
    entry._weight[k] = weight[v];
    entry._perm[k] = IntVector(entry._children.size(), 0);

    return entry._weight[k];
  }
  else if (k - 1 <= static_cast<int>(entry._children.size()))
  {
    // for all permutations summing up to k-1
    IntVector bestPerm;
    double bestPermWeight = -std::numeric_limits<double>::max();

    IntVector split = getInitSplit(entry, v);
    do {
      IntVector perm = convertSplit(v, k-1, split);
      double permWeight = 0;
      for (size_t i = 0; i < perm.size(); i++)
      {
        if (perm[i] > 0)
        {
          double w = solve(entry._children[i], perm[i]);
          if (w == -std::numeric_limits<double>::max())
          {
            permWeight = w;
            break;
          }
          else
          {
            permWeight += w;
          }
        }
      }

      if (permWeight > bestPermWeight)
      {
        bestPerm = perm;
        bestPermWeight = permWeight;
      }
    } while (nextSplit(v, k-1, split));

    entry._perm[k] = bestPerm;
    entry._weight[k] = weight[v] + bestPermWeight;

    return entry._weight[k];
  }
  else
  {
    return -std::numeric_limits<double>::max();
  }
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline bool MwcsSizeTreeMemSolver<GR, NWGHT, NLBL, EWGHT>::solve()
{
  _score = solve(_root, _k);
  if (_score == -std::numeric_limits<double>::max())
  {
    return false;
  }
  else
  {
    reconstruct(_root, _k);
    return true;
  }
}

} // namespace mwcs
} // namespace nina

#endif // MWCSSIZETREEMEMSOLVER_H
