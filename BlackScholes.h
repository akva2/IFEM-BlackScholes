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

#ifndef _BLACKSCHOLES_H
#define _BLACKSCHOLES_H

#include "BDF.h"
#include "IntegrandBase.h"
#include "TimeIntUtils.h"
#include "EqualOrderOperators.h"
#include <array>
#include <memory>


/*!
  \brief Class representing the integrand of Black-Scholes.

  \details Time stepping is done using BDF1/BDF2
*/

class BlackScholes : public IntegrandBase
{
public:
  using WeakOps = EqualOrderOperators::Weak;     //!< Convenience renaming

  //! \brief The default constructor initializes all pointers to zero.
  //! \param[in] n Number of spatial dimensions
  //! \param[in] method The time integration method to use
  //! \param[in] itg_type The integrand type to use
  //! \param[in] useALE If \e true, use ALE formulation
  BlackScholes(unsigned short int n = 3,
               TimeIntegration::Method method = TimeIntegration::BE);

  //! \brief Empty destructor.
  virtual ~BlackScholes() = default;

  using IntegrandBase::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] time Parameters for time-dependent simulations
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
               const TimeDomain& time, const Vec3& X) const override;

  //! \brief Advances the time stepping scheme.
  void advanceStep() { bdf.advanceStep(); }

  //! \brief Returns the name of the primary solution field.
  //! \param[in] prefix Name prefix
  std::string getField1Name(size_t, const char* prefix = 0) const override;

  //! \brief Set the risk-free rate.
  void setRiskFree(double r) { riskFree = r; }
  //! \brief 
  void setSigma(double r) { sigma= r; }

protected:
  TimeIntegration::BDF bdf; //!< BDF helper class

  double riskFree = 0.1; //!< Rate of risk-free returns
  double sigma = 0.1; //!< Standard deviation of returns
};

#endif
