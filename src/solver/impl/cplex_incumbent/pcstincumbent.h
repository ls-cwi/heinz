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
#include <lemon/tolerance.h>
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
                double pT,
                IloFastMutex* pMutex)
    : IloCplex::IncumbentCallbackI(env)
    , _pT(pT)
    , _pMutex(pMutex)
  {
  }
  
  PcstIncumbent(const PcstIncumbent& other)
    : IloCplex::IncumbentCallbackI(other._env)
    , _pT(other._pT)
    , _pMutex(other._pMutex)
  {
  }

protected:
  const double _pT;
  IloFastMutex* _pMutex;
  
  virtual void main()
  {
    lock();
    if (getObjValue() > _highestObj)
    {
      _highestObj = getObjValue();
      *g_pOut << "Solution " << g_timer.realTime() << " " << -1 * _highestObj + _pT << std::endl;
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
