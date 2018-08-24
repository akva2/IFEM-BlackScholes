// $Id$
//==============================================================================
//!
//! \file BlackScholes.h
//!
//! \date Aug 24 2018
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand implementations for Black-Scholes.
//!
//==============================================================================

#include "BlackScholes.h"
#include "FiniteElement.h"
#include "TimeDomain.h"
#include "ElmMats.h"
#include "AnaSol.h"
#include "Function.h"
#include "Vec3Oper.h"
#include "Utilities.h"


BlackScholes::BlackScholes(unsigned short int n,
                           TimeIntegration::Method method)
  : IntegrandBase(n), bdf(TimeIntegration::Order(method))
{
  primsol.resize(1+bdf.getActualOrder());
}


bool BlackScholes::evalInt (LocalIntegral& elmInt,
                            const FiniteElement& fe,
                            const TimeDomain& time,
                            const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  double V = 0;
  for (int t = 0; t< bdf.getOrder(); ++t) {
    double val = fe.N.dot(elMat.vec[t+1]);
    V += val*bdf[1+t]/time.dt;
  }

  WeakOps::Mass(elMat.A.front(), fe, bdf[0]/time.dt  - riskFree);
  WeakOps::Laplacian(elMat.A.front(), fe, 0.5*sigma*sigma*(X[0]*X[0]));
  Vec3 a;
  a[0] = riskFree*X[0];
  WeakOps::Advection(elMat.A.front(), fe, a);

  WeakOps::Source(elMat.b.front(), fe, -V);

  return true;
}


std::string BlackScholes::getField1Name (size_t, const char* prefix) const
{
  if (!prefix) return "Value";

  return prefix + std::string(" Value");
}
