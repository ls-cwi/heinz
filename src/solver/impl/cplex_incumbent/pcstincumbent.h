/*
 * pcstincumbent.h
 *
 *  Created on: 25-nov-2014
 *      Author: M. El-Kebir
 */

#ifndef PCSTINCUMBENT_H
#define PCSTINCUMBENT_H

#include <ilcplex/ilocplex.h>
#include <ilcplex/ilocplexi.h>
#include <ilconcert/ilothread.h>
#include <limits>
#include <set>
#include "utils.h"
#include "mwcsgraph.h"

namespace nina {
namespace mwcs {

template<typename GR,
         typename NWGHT = typename GR::template NodeMap<double>,
         typename NLBL = typename GR::template NodeMap<std::string>,
         typename EWGHT = typename GR::template EdgeMap<double> >
class PcstIncumbent : public IloCplex::IncumbentCallbackI
{
public:
  typedef GR Graph;
  typedef NWGHT WeightNodeMap;
  typedef NLBL LabelNodeMap;
  typedef EWGHT WeightEdgeMap;
  
  typedef MwcsGraph<Graph, WeightNodeMap, LabelNodeMap, WeightEdgeMap> MwcsGraphType;

  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  
  typedef std::set<Node> NodeSet;
  typedef typename NodeSet::const_iterator NodeSetIt;
  
public:
  PcstIncumbent(IloEnv env,
                const MwcsGraphType& mwcsGraph,
                IloBoolVarArray x,
                IloFastMutex* pMutex)
    : IloCplex::IncumbentCallbackI(env)
    , _mwcsGraph(mwcsGraph)
    , _x(x)
    , _pMutex(pMutex)
  {
  }
  
  PcstIncumbent(const PcstIncumbent& other)
    : IloCplex::IncumbentCallbackI(other._env)
    , _mwcsGraph(other._mwcsGraph)
    , _x(other._x)
    , _pMutex(other._pMutex)
  {
  }
  
  static double reEvaluate(const MwcsGraphType& mwcsGraph,
                           const NodeSet& solution);
  
protected:
  const MwcsGraphType& _mwcsGraph;
  IloBoolVarArray _x;
  IloFastMutex* _pMutex;
  
  virtual void main()
  {
    lock();
    if (getObjValue() > _highestObj)
    {
      _highestObj = getObjValue();
      *g_pOut << "Solution " << g_timer.realTime() << " " << getObjValue() << std::endl;
    }
    unlock();
  }
  
  virtual IloCplex::CallbackI* duplicateCallback() const
  {
    return (new (_env) PcstIncumbent(*this));
  }
  
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
  
  static double _highestObj;
};

template<typename GR, typename NWGHT, typename NLBL, typename EWGHT>
double PcstIncumbent<GR, NWGHT, NLBL, EWGHT>::_highestObj = -std::numeric_limits<double>::max();
  
} // namespace mwcs
} // namespace nina

#endif // INCUMBENT_H
