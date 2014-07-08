/*
 * analysis.h
 *
 *  Created on: 25-jan-2013
 *      Author: M. El-Kebir
 */

#ifndef MWCSANALYZE_H
#define MWCSANALYZE_H

#include <lemon/core.h>
#include <lemon/dijkstra.h>
#include <lemon/suurballe.h>
#include <limits>
#include <set>
#include <vector>
#include <list>
#include "mwcsgraph.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename WGHT = typename GR::template NodeMap<double> >
class MwcsAnalyze
{
public:
  typedef GR Graph;
  typedef WGHT WeightNodeMap;
  typedef MwcsGraph<GR, WGHT> MwcsGraphType;
  typedef typename GR::template ArcMap<double> WeightArcMap;
  typedef typename GR::template NodeMap<WeightNodeMap*> WeightNodeMatrix;

  TEMPLATE_GRAPH_TYPEDEFS(Graph);

  typedef typename GR::template NodeMap<BoolNodeMap*> BoolNodeMatrix;

  typedef typename std::vector<Node> NodeVector;
  typedef typename std::set<Node> NodeSet;
  typedef typename NodeSet::const_iterator NodeSetIt;
  typedef typename std::vector<NodeSet> NodeSetVector;
  typedef typename NodeSetVector::const_iterator NodeSetVectorIt;

  typedef typename lemon::Dijkstra<GR, WeightArcMap> DijkstraType;
  typedef typename std::list<Node> NodeList;
  typedef typename NodeList::const_iterator NodeListIt;

  typedef typename lemon::Suurballe<GR, WeightArcMap> SuurballeType;
  typedef typename SuurballeType::Path PathType;
  typedef typename GR::template NodeMap<PathType> PathMap;
  typedef typename GR::template NodeMap<PathMap*> PathMatrix;

  typedef typename GR::template NodeMap<NodeSet> NodeSetMap;


public:
  MwcsAnalyze(const MwcsGraphType& mwcsGraph);
  MwcsAnalyze(const Graph& g,
              const WeightNodeMap& weight,
              const IntNodeMap& comp);
  ~MwcsAnalyze();

  void analyze(bool eqClasses = false);
  void analyzeS(int k);
  void analyzeNegHubs();

  int getEqClassesCount() const { return static_cast<int>(_eqClasses.size()); }
  const NodeSetVector& getEqClasses() const { return _eqClasses; }
  int getEqClass(Node n) const { return _eqClassMap[n]; }
  void print(std::ostream& out) const
  {
    int eqClassIdx = 0;
    for (NodeSetVectorIt nodeSetIt = _eqClasses.begin(); nodeSetIt != _eqClasses.end(); nodeSetIt++, eqClassIdx++)
    {
      if (nodeSetIt->size() == 1) continue;
      out << "// Equivalence class " << eqClassIdx << " size " << nodeSetIt->size() << std::endl;
      out << "// ";
      for (NodeSetIt nodeIt = nodeSetIt->begin(); nodeIt != nodeSetIt->end(); nodeIt++)
      {
        out << _g.id(*nodeIt) << " ";
      }
      out << std::endl;
    }
  }
  
  void print(const PathType& p, std::ostream& out) const
  {
    NodeList path;
    bool f = true;
    for (typename PathType::ArcIt a(p); a != lemon::INVALID; ++a)
    {
      if (f)
      {
        path.push_back(_g.source(a));
        f = false;
      }
      path.push_back(_g.target(a));
    }
    
    double cumWeight = 0;
    double maxCumWeight = -std::numeric_limits<double>::max();
    
    bool first = true;
    for (NodeListIt nodeIt = path.begin(); nodeIt != path.end(); nodeIt++)
    {
      cumWeight += _weight[*nodeIt];
      if (cumWeight > maxCumWeight)
        maxCumWeight = cumWeight;
      
      if (!first)
      {
        std::cout << " -> ";
      }
      else
      {
        first = false;
      }
      
      out << _g.id(*nodeIt) << " (" << _weight[*nodeIt] << ")";
    }
    out << " : " << cumWeight << "/" << maxCumWeight << std::endl;
  }

  bool ok(Node i, Node j) const
  {
    if (_ok[i] != NULL)
      return (*_ok[i])[j];
    else
      return false;
  }

  const PathType& path(Node i, Node j) const
  {
    assert(_path[i] != NULL);
    return (*_path[i])[j];
  }

  void printNegHubs(const MwcsGraphType& mwcsGraph, std::ostream& out) const
  {
    for (NodeIt v(_g); v != lemon::INVALID; ++v)
    {
      if (_beneficial[v])
      {
        out << "// NegHub: " << mwcsGraph.getLabel(v) << "\t"
            << _weight[v] << "\t"
            << _benefit[v] << std::endl;

        const NodeSet& posNeighbors = _posNeighbors[v];
        for (NodeSetIt nodeIt = posNeighbors.begin();
             nodeIt != posNeighbors.end(); nodeIt++)
        {
          out << mwcsGraph.getLabel(*nodeIt) << "\t"
              << _weight[*nodeIt] << std::endl;
        }
        out << std::endl;
      }
    }
  }

private:
  const Graph& _g;
  const WeightNodeMap& _weight;
  const IntNodeMap& _comp;
  WeightArcMap _arcWeight;
  NodeSetVector _eqClasses;
  IntNodeMap _eqClassMap;
  BoolNodeMatrix _ok;
  PathMatrix _path;

  DoubleNodeMap _benefit;
  NodeSetMap _posNeighbors;
  BoolNodeMap _beneficial;
  int _nBenificialNegHubs;
  NodeSet _beneficialNegHubs;
  NodeVector _rouletteWheel;

  double initArcWeights();
  bool isPathOK(const PathType& p) const;

public:
  const NodeSet& getBeneficialNegHubs() const
  {
    return _beneficialNegHubs;
  }

  double getBenefit(Node negHub) const
  {
    return _benefit[negHub];
  }

  int getNumberOfBeneficialNegHubs() const
  {
    return _nBenificialNegHubs;
  }

  const NodeSet& getPosNeighbors(Node negHub) const
  {
    return _posNeighbors[negHub];
  }

  const NodeVector& getRouletteWheel() const
  {
    return _rouletteWheel;
  }
};

template<typename GR, typename WGHT>
inline MwcsAnalyze<GR, WGHT>::~MwcsAnalyze()
{
  for (NodeIt n(_g); n != lemon::INVALID; ++n)
  {
    delete _ok[n];
    delete _path[n];
  }
}

template<typename GR, typename WGHT>
inline MwcsAnalyze<GR, WGHT>::MwcsAnalyze(const MwcsGraphType& mwcsGraph)
  : _g(mwcsGraph.getGraph())
  , _weight(mwcsGraph.getScores())
  , _comp(mwcsGraph.getComponentMap())
  , _arcWeight(_g)
  , _eqClasses()
  , _eqClassMap(_g, -1)
  , _ok(_g, NULL)
  , _path(_g, NULL)
  , _benefit(_g)
  , _posNeighbors(_g)
  , _beneficial(_g, false)
  , _nBenificialNegHubs(0)
  , _beneficialNegHubs()
  , _rouletteWheel()
{
}

template<typename GR, typename WGHT>
inline MwcsAnalyze<GR, WGHT>::MwcsAnalyze(const Graph& g,
                                          const WeightNodeMap& weight,
                                          const IntNodeMap& comp)
  : _g(g)
  , _weight(weight)
  , _comp(comp)
  , _arcWeight(_g)
  , _eqClasses()
  , _eqClassMap(_g, -1)
  , _ok(_g, NULL)
  , _path(_g, NULL)
  , _benefit(_g)
  , _posNeighbors(_g)
  , _beneficial(_g, false)
  , _nBenificialNegHubs(0)
  , _beneficialNegHubs()
  , _rouletteWheel()
{
}

template<typename GR, typename WGHT>
inline bool MwcsAnalyze<GR, WGHT>::isPathOK(const PathType& p) const
{
  NodeList path;
  bool f = true;
  for (typename PathType::ArcIt a(p); a != lemon::INVALID; ++a)
  {
    if (f)
    {
      path.push_back(_g.source(a));
      f = false;
    }
    path.push_back(_g.target(a));
  }

  //if (path.size() && path.front() != i)
  //  return false;

  //bool first = true;

  double cumWeight = 0;
  double maxCumWeight = -std::numeric_limits<double>::max();

  for (NodeListIt nodeIt = path.begin(); nodeIt != path.end(); nodeIt++)
  {
    cumWeight += _weight[*nodeIt];
    if (cumWeight > maxCumWeight)
      maxCumWeight = cumWeight;

    //if (!first)
    //{
    //  std::cout << " -> ";
    //}
    //else
    //{
    //  first = false;
    //}

    //std::cout << _g.id(*nodeIt);
  }
  //std::cout << " : " << cumWeight << "/" << maxCumWeight << std::endl;

  return cumWeight == maxCumWeight && cumWeight >= 0;
}

template<typename GR, typename WGHT>
inline void MwcsAnalyze<GR, WGHT>::analyzeS(int k)
{
  initArcWeights();
  SuurballeType sb(_g, _arcWeight);

  for (NodeIt i(_g); i != lemon::INVALID; ++i)
  {
    if (_weight[i] <= 0)
      continue;

    int comp_i = _comp[i];
    delete _ok[i];
    _ok[i] = new BoolNodeMap(_g, false);

    // compute single source shortest path from i
    sb.init(i);

    for (NodeIt j(_g); j != lemon::INVALID; ++j)
    {
      if (i == j)
        continue;

      if (_weight[j] <= 0 || comp_i != _comp[j])
        continue;

      int n = sb.start(j, k);
      for (int idx = 0; idx < n; idx++)
      {
        PathType path = sb.path(idx);
        if (isPathOK(path))
        {
          (*_ok[i])[j] = true;
          break;
        }
      }
    }
  }
}

template<typename GR, typename WGHT>
inline void MwcsAnalyze<GR, WGHT>::analyze(bool eqClasses)
{
  initArcWeights();
  lemon::Dijkstra<GR, WeightArcMap> dijkstra(_g, _arcWeight);

  int count = 0;
  for (NodeIt i(_g); i != lemon::INVALID; ++i)
  {
    if (_weight[i] <= 0)
      continue;

    int comp_i = _comp[i];
    _ok[i] = new BoolNodeMap(_g, false);
    _path[i] = new PathMap(_g);

    // compute single source shortest path from i
    int localCount = 0;
    do
    {
      dijkstra.init();
      dijkstra.run(i);

      for (NodeIt j(_g); j != lemon::INVALID; ++j)
      {
        if (i == j)
          continue;

        if (_weight[j] <= 0 || comp_i != _comp[j])
          continue;

        if (dijkstra.reached(j))
        {
          PathType p = dijkstra.path(j);
          if (isPathOK(p))
          {
            (*_ok[i])[j] = true;
            (*_path[i])[j] = p;
            count++;
          }
        }
      }
    } while (localCount != 0);
  }

  if (g_verbosity >= VERBOSE_ESSENTIAL)
    std::cout << "// Identified " << count << " dependent node pairs" << std::endl;

  // determine all equivalence classes
  if (eqClasses)
  {
    lemon::mapFill(_g, _eqClassMap, -1);
    int eqClassIdx = 0;
    _eqClasses.clear();

    for (NodeIt i(_g); i != lemon::INVALID; ++i)
    {
      if (_weight[i] <= 0)
        continue;

      for (NodeIt j = i; j != lemon::INVALID; ++j)
      {
        if (_weight[j] <= 0)
          continue;

        if (i != j)
        {
          bool ok_i_j = (*_ok[i])[j];
          bool ok_j_i = (*_ok[j])[i];
          
          if (_g.id(i) == 3522 && _g.id(j) == 3437)
          {
            print(path(i, j), std::cerr);
            print(path(j, i), std::cerr);
          }

          if (ok_i_j && ok_j_i)
          {
            print(path(i, j), std::cerr);
            print(path(j, i), std::cerr);
            
            if (_eqClassMap[i] == -1 && _eqClassMap[j] == -1)
            {
              // add to new eqClass
              _eqClasses.push_back(NodeSet());
              _eqClassMap[i] = _eqClassMap[j] = eqClassIdx;
              _eqClasses[eqClassIdx].insert(i);
              _eqClasses[eqClassIdx].insert(j);

              eqClassIdx++;
            }
            else if (_eqClassMap[i] == -1)
            {
              _eqClassMap[i] = _eqClassMap[j];
              _eqClasses[_eqClassMap[i]].insert(i);
            }
            else if (_eqClassMap[j] == -1)
            {
              _eqClassMap[j] = _eqClassMap[i];
              _eqClasses[_eqClassMap[j]].insert(j);
            }
            else
            {
//              assert(_eqClassMap[i] == _eqClassMap[j]);
              // merge j into i
              for (NodeSetIt it = _eqClasses[_eqClassMap[j]].begin(), it_end = _eqClasses[_eqClassMap[j]].end();
                   it != it_end; ++it)
              {
                _eqClassMap[*it] = _eqClassMap[i];
              }
              _eqClasses[_eqClassMap[i]].insert(_eqClasses[_eqClassMap[j]].begin(), _eqClasses[_eqClassMap[j]].end());
            }
            std::cerr << _eqClassMap[i] << std::endl << std::endl;
          }
        }
      }

      // wtf is this?
//      if (_eqClassMap[i] == -1)
//      {
//        _eqClassMap[i] = eqClassIdx;
//        _eqClasses.push_back(NodeSet());
//        _eqClasses[eqClassIdx].insert(i);
//        eqClassIdx++;
//      }
    }
  }
}

template<typename GR, typename WGHT>
inline void MwcsAnalyze<GR, WGHT>::analyzeNegHubs()
{
  // 1. start by identifying all negative hubs
  _nBenificialNegHubs = 0;
  _beneficialNegHubs.clear();
  for (NodeIt v(_g); v != lemon::INVALID; ++v)
  {
    _beneficial[v] = false;
    _benefit[v] = 0;

    NodeSet& posNeighbors = _posNeighbors[v];
    posNeighbors.clear();

    if (_weight[v] >= 0) continue;

    double maxPosWeight = 0;
    _benefit[v] = _weight[v];
    for (IncEdgeIt e(_g, v); e != lemon::INVALID; ++e)
    {
      Node w = _g.oppositeNode(v, e);
      double weight_w = _weight[w];
      if (weight_w >= 0)
      {
        posNeighbors.insert(w);
        _benefit[v] += weight_w;
        if (weight_w > maxPosWeight)
        {
          maxPosWeight = weight_w;
        }
      }
    }

    if (_benefit[v] >= maxPosWeight)
    {
      _beneficial[v] = true;
      _nBenificialNegHubs++;
      _beneficialNegHubs.insert(v);
    }
  }

  // 2. construct roulette wheel
  _rouletteWheel.clear();

  for (NodeSetIt nodeIt = _beneficialNegHubs.begin();
       nodeIt != _beneficialNegHubs.end(); nodeIt++)
  {
    int b = round(_benefit[*nodeIt]);
    for (int i = 0; i < b; i++) _rouletteWheel.push_back(*nodeIt);
  }
}

template<typename GR, typename WGHT>
inline double MwcsAnalyze<GR, WGHT>::initArcWeights()
{
  double offset = 0;

  // first determine the offset
  for (NodeIt n(_g); n != lemon::INVALID; ++n)
  {
    if (_weight[n] < offset)
    {
      offset = _weight[n];
    }
  }

  // now compute arc weights
  for (ArcIt a(_g); a != lemon::INVALID; ++a)
  {
    Node v = _g.target(a);
    _arcWeight[a] = _weight[v] - offset;
  }

  return -offset;
}

} // namespace mwcs
} // namespace nina

#endif // MWCSANALYZE_H
