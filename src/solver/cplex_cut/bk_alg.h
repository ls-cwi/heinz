/*
 *  bk_alg.h
 *
 *   Created on: 12-may-2013
 *       Author: M. El-Kebir
 */
// BIG FAT WARNING: reverse arcs not supported
// (uncomment ArcLookUp to add support)

#ifndef BK_ALG_H
#define BK_ALG_H

#include <lemon/core.h>
#include <limits>
#include <ostream>
#include <maxflow-v3.01/graph.h>

namespace nina {

template<typename DGR>
class BkFlowAlg
{
public:
  typedef DGR Digraph;

  TEMPLATE_DIGRAPH_TYPEDEFS(Digraph);
  typedef DoubleArcMap CapacityMap;

  BkFlowAlg(const Digraph& g,
            const CapacityMap& cap)
    : _g(g)
    , _cap(_g)
    , _source(lemon::INVALID)
    , _target(lemon::INVALID)
    , _bkNode(_g, -1)
    , _bkArc(_g, NULL)
    , _pBK(NULL)
    , _flow()
  {
    lemon::mapCopy(_g, cap, _cap);
    init();
  }

  virtual ~BkFlowAlg()
  {
    delete _pBK;
  }

  void incCap(Arc a, double c);
  void setCap(const CapacityMap& cap);
  double resCap(Arc a) const;
  double revResCap(Arc a) const;
  double cap(Arc a) const;
  double run(bool reuse = false);
  double flow(Arc a) const;
  bool cut(Arc a) const;
  bool cut(Node v) const;
  double maxFlow() const;
  Node getSource() const { return _source; }
  Node getTarget() const { return _target; }
  void setSource(Node source);
  void setTarget(Node target);

  void printFlow(std::ostream& out, bool cutOnly = false) const;
  void printCut(std::ostream& out) const;

private:
  typedef bk::Graph<double, double, double> BkGraphType;
  typedef BkGraphType::arc BkArc;
  typedef lemon::ArcLookUp<Digraph> ArcLookUpType;
  typedef typename Digraph::template ArcMap<BkArc*> ToBkArcMap;
  typedef IntNodeMap ToBkNodeMap;

  const Digraph& _g;
  CapacityMap _cap;
  Node _source;
  Node _target;
  ToBkNodeMap _bkNode;
  ToBkArcMap _bkArc;

  BkGraphType* _pBK;
  double _flow;

  void init();
};

template<typename DGR>
void BkFlowAlg<DGR>::init()
{
  const int n = lemon::countNodes(_g);
  const int m = lemon::countArcs(_g);
  //const ArcLookUpType lookUp(_g);

  delete _pBK;
  _pBK = new BkGraphType(n, m);
  _pBK->add_node(n);

  // copy nodes
  int i = 0;
  for (NodeIt v(_g); v != lemon::INVALID; ++v)
  {
    _bkNode[v] = i++;
  }

  // copy arcs
  for (ArcIt a(_g); a != lemon::INVALID; ++a)
  {
    Node v = _g.source(a);
    Node w = _g.target(a);

    // check if reverse arc also occurs
    // NO SELF LOOPS!
    //Arc revA = lookUp(w, v);
    //if (revA != lemon::INVALID && _g.id(v) < _g.id(w))
    //{
    //  BkArc* pBkArc = _pBK->add_edge(_bkNode[v], _bkNode[w], _cap[a], _cap[revA]);
    //  _bkArc[a] = pBkArc;
    //  _bkArc[revA] = ++pBkArc;
    //}
    //else
    {
      _bkArc[a] = _pBK->add_edge(_bkNode[v], _bkNode[w], _cap[a], 0);
    }
  }
  //const double max = std::numeric_limits<double>::max();
  //_pBK->add_tweights(_bkNode[_source], max, 0);
  //_pBK->add_tweights(_bkNode[_target], 0, max);
}

template<typename DGR>
double BkFlowAlg<DGR>::run(bool reuse)
{
  _flow = _pBK->maxflow(reuse);
  return _flow;
}

template<typename DGR>
double BkFlowAlg<DGR>::maxFlow() const
{
  return _flow;
}

template<typename DGR>
bool BkFlowAlg<DGR>::cut(Arc a) const
{
  Node v = _g.source(a);
  Node w = _g.target(a);
  return cut(v) && !cut(w);
}

template<typename DGR>
bool BkFlowAlg<DGR>::cut(Node v) const
{
  return _pBK->what_segment(_bkNode[v]) == BkGraphType::SOURCE;
}

template<typename DGR>
double BkFlowAlg<DGR>::flow(Arc a) const
{
  const BkArc* pBkArc = _bkArc[a];
  return _cap[a] - pBkArc->r_cap;
}

template<typename DGR>
void BkFlowAlg<DGR>::incCap(Arc a, double c)
{
  BkArc* pBkArc = _bkArc[a];
  pBkArc->r_cap += c;
  _cap[a] += c;

  _pBK->mark_node(_bkNode[_g.source(a)]);
}

template<typename DGR>
void BkFlowAlg<DGR>::setCap(const CapacityMap& cap)
{
  //const ArcLookUpType lookUp(_g);

  lemon::mapCopy(_g, cap, _cap);
  for (ArcIt a(_g); a != lemon::INVALID; ++a)
  {
    double cap_a = _cap[a];
    BkArc* pBkArc = _bkArc[a];

    //if (lookUp(_g.target(a), _g.source(a)) == lemon::INVALID)
    {
      //assert(pBkArc->r_cap == cap_a);
      pBkArc->r_cap = cap_a;
      //assert(pBkArc->sister->r_cap == 0);
      pBkArc->sister->r_cap = 0;
    }
  }

  for (NodeIt v(_g); v != lemon::INVALID; ++v)
  {
    //assert(_pBK->get_trcap(_bkNode[v]) == 0);
    _pBK->set_trcap(_bkNode[v], 0);
  }

  const double max = std::numeric_limits<double>::max();
  if (_source != lemon::INVALID)
    _pBK->add_tweights(_bkNode[_source], max, 0);
  if (_target != lemon::INVALID)
    _pBK->add_tweights(_bkNode[_target], 0, max);
}

template<typename DGR>
double BkFlowAlg<DGR>::resCap(Arc a) const
{
  BkArc* pBkArc = _bkArc[a];
  return pBkArc->r_cap;
}

template<typename DGR>
double BkFlowAlg<DGR>::revResCap(Arc a) const
{
  BkArc* pBkArc = _bkArc[a];
  return pBkArc->sister->r_cap;
}

template<typename DGR>
double BkFlowAlg<DGR>::cap(Arc a) const
{
  return _cap[a];
}

template<typename DGR>
void BkFlowAlg<DGR>::setSource(Node source)
{
  if (_source != lemon::INVALID)
    _pBK->set_trcap(_bkNode[_source], 0);

  const double max = std::numeric_limits<double>::max();
  _source = source;
  _pBK->add_tweights(_bkNode[_source], max, 0);
}

template<typename DGR>
void BkFlowAlg<DGR>::setTarget(Node target)
{
  if (_target != lemon::INVALID)
  {
    _pBK->set_trcap(_bkNode[_target], 0);
  }

  const double max = std::numeric_limits<double>::max();
  _target = target;
  _pBK->add_tweights(_bkNode[_target], 0, max);
}

template<typename DGR>
void BkFlowAlg<DGR>::printFlow(std::ostream& out, bool cutOnly) const
{
  for (ArcIt a(_g); a != lemon::INVALID; ++a)
  {
    Node v = _g.source(a);
    Node w = _g.target(a);
    if (cut(a))
    {
      out << _g.id(a) << " : ";
      out << _g.id(v) << " -> " << _g.id(w)
          << " * " << flow(a) << "/" << _cap[a]
          << std::endl;
    }
    else if (!cutOnly)
    {
      out << _g.id(a) << " : ";
      out << _g.id(v) << " -> " << _g.id(w)
          << "   " << flow(a) << "/" << _cap[a]
          << std::endl;
    }
  }
}

template<typename DGR>
void BkFlowAlg<DGR>::printCut(std::ostream& out) const
{
  for (NodeIt v(_g); v != lemon::INVALID; ++v)
  {
    out << _g.id(v) << (cut(v) ? ": source" : ": target") << std::endl;
  }
}

} // namespace nina

#endif // BK_ALG_H
