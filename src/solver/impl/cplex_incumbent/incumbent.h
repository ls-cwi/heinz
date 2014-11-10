/*
 * incumbent.h
 *
 *  Created on: 10-nov-2014
 *      Author: M. El-Kebir
 */

#ifndef INCUMBENT_H
#define INCUMBENT_H

#include <ilcplex/ilocplex.h>
#include <ilcplex/ilocplexi.h>
#include <ilconcert/ilothread.h>
#include "utils.h"

namespace nina {
namespace mwcs {

class Incumbent : public IloCplex::IncumbentCallbackI
{
public:
  Incumbent(IloEnv env,
            IloFastMutex* pMutex)
    : IloCplex::IncumbentCallbackI(env)
    , _pMutex(pMutex)
  {
  }
  
  Incumbent(const Incumbent& other)
    : IloCplex::IncumbentCallbackI(other._env)
    , _pMutex(other._pMutex)
  {
  }
  
protected:
  IloFastMutex* _pMutex;
  
  virtual void main()
  {
    *g_pOut << "Solution " << g_timer.realTime() << " " << getObjValue() << std::endl;
  }
  
  virtual IloCplex::CallbackI* duplicateCallback() const
  {
    return (new (_env) Incumbent(*this));
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
};
  
} // namespace mwcs
} // namespace nina

#endif // INCUMBENT_H
