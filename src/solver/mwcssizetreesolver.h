/*
 * mwcssizetreesolver.h
 *
 *  Created on: 2-jul-2013
 *      Author: M. El-Kebir
 */

#ifndef MWCSSIZETREESOLVER_H
#define MWCSSIZETREESOLVER_H

#include "mwcssolverunrooted.h"
#include <set>
#include <map>
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
class MwcsSizeTreeSolver : public MwcsSolverUnrooted<GR, NWGHT, NLBL, EWGHT>
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;

  typedef MwcsSolverUnrooted<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> Parent;
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

  typedef std::vector<double> DoubleVector;
  typedef std::vector<int> IntVector;
  typedef std::vector<IntVector> IntMatrix;
  typedef IntVector::const_iterator IntVectorIt;
  typedef IntVector::const_reverse_iterator IntVectorRevIt;

  typedef std::vector<NodeVector> NodeMatrix;
  typedef typename NodeMatrix::const_iterator NodeMatrixIt;
  typedef typename NodeMatrix::const_reverse_iterator NodeMatrixRevIt;

  struct DpEntry {
    int _subTreeSize;
    NodeVector _children;
    IntMatrix _solutionPerm;
    DoubleVector _weight;

    DpEntry()
      : _subTreeSize(1)
      , _children()
      , _solutionPerm()
      , _weight()
    {
    }

    void clear()
    {
      _subTreeSize = 1;
      _weight.clear();// = -std::numeric_limits<double>::max();
      _solutionPerm.clear();
      _children.clear();
    }
  };

  typedef typename Graph::template NodeMap<DpEntry> DpEntryMap;

public:
  MwcsSizeTreeSolver(const MwcsGraphType& mwcsGraph, int k);
  virtual void init(Node root);
  int initSubTreeSize(Node node);
  virtual bool solve();
  virtual void setLowerBound(double LB) {}

private:
  const int _k;
  DpEntryMap _dpMap;
  NodeMatrix _nodesPerLevel;

  IntVector getInitSplit(const Node v) const
  {
    const int n = static_cast<int>(_dpMap[v]._children.size() - 1);
    return IntVector(n, 0);
  }

  bool nextSplit(const Node v, const int k, IntVector& split) const
  {
    const int n = static_cast<int>(split.size());
    assert(n == static_cast<int>(_dpMap[v]._children.size()) - 1);

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

  bool isFeasiblePerm(const Node v, const int k, const IntVector& perm) const
  {
    const int n = static_cast<int>(perm.size());
    const DpEntry& entry = _dpMap[v];

    assert(perm.size() == entry._children.size());

#ifndef NDEBUG
    // check sum
    int sum = 0;
    //std::cout << "Perm " << k << ": ";
    for (int i = 0; i < n; i++)
    {
      sum += perm[i];
      //std::cout << perm[i] << " ";
    }
    //std::cout << std::endl;
    assert(sum == k);
#endif

    // check whether sizes make sense
    for (int i = 0; i < n; i++)
    {
      if (perm[i] > _dpMap[entry._children[i]]._subTreeSize) return false;
    }

    return true;
  }

  double getPermWeight(const Node v, const IntVector& perm)
  {
    const DpEntry& entry = _dpMap[v];
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
    const DpEntry& entry = _dpMap[v];
    const IntVector& perm = entry._solutionPerm[k-1];

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
inline MwcsSizeTreeSolver<GR, NWGHT, NLBL, EWGHT>::MwcsSizeTreeSolver(const MwcsGraphType& mwcsGraph,
                                                                      int k)
  : Parent(mwcsGraph)
  , _k(k)
  , _dpMap(mwcsGraph.getGraph())
  , _nodesPerLevel()
{
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline int MwcsSizeTreeSolver<GR, NWGHT, NLBL, EWGHT>::initSubTreeSize(Node node)
{
  DpEntry& entry = _dpMap[node];
  entry._subTreeSize = 1;

  if (entry._children.size() > 0)
  {
    // recurse on children
    for (NodeVectorIt nodeIt = entry._children.begin(); nodeIt != entry._children.end(); nodeIt++)
    {
      entry._subTreeSize += initSubTreeSize(*nodeIt);
    }
  }

  entry._weight = DoubleVector(std::min(_k, entry._subTreeSize), 0);
  entry._solutionPerm = IntMatrix(std::min(_k, entry._subTreeSize));

  return entry._subTreeSize;
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline void MwcsSizeTreeSolver<GR, NWGHT, NLBL, EWGHT>::init(Node root)
{
  const Graph& g = _mwcsGraph.getGraph();

  PredMap pred(g);
  DistMap level(g);

  _nodesPerLevel.clear();
  _score = 0;
  _solutionSet.clear();
  for (NodeIt v(g); v != lemon::INVALID; ++v)
  {
    _solutionMap[v] = false;
    _dpMap[v].clear();
  }

  _root = root;

  assert(lemon::acyclic(g));
  assert(root != lemon::INVALID);

  // start by doing a bfs to determine levels
  lemon::Bfs<Graph> bfs(g);
  bfs.init();
  bfs.predMap(pred);
  bfs.distMap(level);

  bfs.addSource(root);

  while (!bfs.emptyQueue())
  {
    Node node = bfs.processNextNode();
    int bfsLevel = level[node];

    if (bfsLevel == static_cast<int>(_nodesPerLevel.size()))
      _nodesPerLevel.push_back(NodeVector());

    _nodesPerLevel[bfsLevel].push_back(node);

    Node parent = pred[node] != lemon::INVALID ?
          g.oppositeNode(node, pred[node]) : lemon::INVALID;

    if (parent != lemon::INVALID)
    {
      //std::cout << "Child: " << _mwcsGraph.getLabel(node)
      //          << " parent: " << _mwcsGraph.getLabel(parent) << std::endl;
      _dpMap[parent]._children.push_back(node);
    }
  }

  initSubTreeSize(root);
}

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
inline bool MwcsSizeTreeSolver<GR, NWGHT, NLBL, EWGHT>::solve()
{
  if (_k > _dpMap[_root]._subTreeSize)
    return false;

  const WeightNodeMap& weight = _mwcsGraph.getScores();

  // work bottom-up
  for (NodeMatrixRevIt nodeSetIt = _nodesPerLevel.rbegin();
       nodeSetIt != _nodesPerLevel.rend(); nodeSetIt++)
  {
    for (NodeVectorIt nodeIt = nodeSetIt->begin(); nodeIt != nodeSetIt->end(); nodeIt++)
    {
      Node node = *nodeIt;
      DpEntry& entry = _dpMap[node];

      entry._weight[0] = weight[node];
      entry._solutionPerm[0] = IntVector(entry._children.size(), 0);

      // for all k
      for (int i = 2; i <= std::min(entry._subTreeSize, _k); i++)
      {
        IntVector bestPerm;
        double bestPermWeight = -std::numeric_limits<double>::max();

        // for all permutations summing up to i-1
        IntVector split = getInitSplit(node);
        do {
          // print split
          //std::cout << "Split: ";
          //for (size_t l = 0; l < split.size(); l++)
          //{
          //  std::cout << split[l] << " ";
          //}
          //std::cout << std::endl;

          IntVector perm = convertSplit(node, i-1, split);
          if (isFeasiblePerm(node, i-1, perm))
          {
            double permWeight = getPermWeight(node, perm);
            if (permWeight > bestPermWeight)
            {
              bestPerm = perm;
              bestPermWeight = permWeight;
            }
          }
        } while (nextSplit(node, i-1, split));

        entry._solutionPerm[i-1] = bestPerm;
        entry._weight[i-1] = weight[node] + bestPermWeight;
      }
    }
  }

  // construct the solution
  reconstruct(_root, _k);

  return true;
}

} // namespace mwcs
} // namespace nina

#endif // MWCSSIZETREESOLVER_H
