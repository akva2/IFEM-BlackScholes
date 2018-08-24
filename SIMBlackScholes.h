// $Id$
//==============================================================================
//!
//! \file SIMBlackScholes.h
//!
//! \date Aug 24 2018
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Simulation driver for Isogeometric FE analysis of Black-Scholes.
//!
//==============================================================================

#ifndef _SIM_BLACKSCHOLES_H_
#define _SIM_BLACKSCHOLES_H_

#include "IFEM.h"
#include "BlackScholes.h"
#include "Utilities.h"
#include "TimeStep.h"
#include "tinyxml.h"
#include "DataExporter.h"


/*!
  \brief Driver class for isogeometric FE analysis of Black-Scholes problems.
*/

template<class Dim> class SIMBlackScholes : public Dim
{
public:
  //! \brief Default constructor.
  SIMBlackScholes(TimeIntegration::Method method) : Dim(1), bs(Dim::dimension, method)
  {
    Dim::myProblem = &bs;
  }

  //! \brief Destructor.
  virtual ~SIMBlackScholes()
  {
    Dim::myProblem = nullptr;
    Dim::myInts.clear();
  }

  //! \brief Returns the name of this simulator (for use in the HDF5 export).
  std::string getName() const override { return "BlackScholes"; }

  //! \brief Register fields for data export.
  void registerFields(DataExporter& exporter)
  {
    int results = DataExporter::PRIMARY;

    if (!Dim::opt.pSolOnly)
      results |= DataExporter::SECONDARY;

    exporter.registerField("u", "primary", DataExporter::SIM, results);
    exporter.setFieldValue("u", this, &sol.front());
  }

  //! \brief Opens a new VTF-file and writes the model geometry to it.
  //! \param[in] fileName File name used to construct the VTF-file name from
  //! \param[out] geoBlk Running geometry block counter
  //! \param[out] nBlock Running result block counter
  bool saveModel(char* fileName, int& geoBlk, int& nBlock)
  {
    if (Dim::opt.format < 0)
      return true;

    nBlock = 0;
    return this->writeGlvG(geoBlk,fileName);
  }

  //! \brief Saves the converged results to VTF file of a given time step.
  //! \param[in] nBlock Running VTF block counter
  bool saveStep(const TimeStep& tp, int& nBlock)
  {
    if (Dim::opt.format < 0)
      return true;

    // Write solution fields
    int iDump = 1 + tp.step/Dim::opt.saveInc;
    if (!this->writeGlvS(sol.front(),iDump,nBlock,tp.time.t))
      return false;

    return this->writeGlvStep(iDump,tp.time.t);
  }

  //! \brief Initializes for time-dependent simulation.
  bool init()
  {
    if (!this->setMode(SIM::DYNAMIC))
      return false;

    sol.resize(bs.getNoSolutions());

    for (auto&  v : sol)
      v.resize(this->getNoDOFs());

    return true;
  }

  //! \brief Computes the solution for the current time step.
  bool solveStep(TimeStep& tp)
  {
    if (Dim::msgLevel >= 0)
      IFEM::cout <<"\n  step = "<< tp.step <<"  time = "<< tp.time.t << std::endl;

    Vector dummy;
    this->updateDirichlet(tp.time.t, &dummy);

    if (!this->assembleSystem(tp.time,sol))
      return false;

    if (!this->solveSystem(sol.front(), Dim::msgLevel-1,"Value       "))
      return false;

    return true;
  }

  //! \brief Advances the time step one step forward.
  bool advanceStep(TimeStep&)
  {
    // Update vectors between time steps
    const int nNusols = sol.size();
    for (int n = nNusols-1; n > 0; n--)
      sol[n] = sol[n-1];

    bs.advanceStep();

    return true;
  }

  using Dim::parse;
  //! \brief Parses a data section from an XML element.
  bool parse(const TiXmlElement* elem) override
  {
    if (strcasecmp(elem->Value(),"blackscholes"))
      return this->Dim::parse(elem);

    double r;
    if (utl::getAttribute(elem,"r",r))
      bs.setRiskFree(r);
    if (utl::getAttribute(elem,"sigma",r))
      bs.setSigma(r);

    return true;
  }

private:
  BlackScholes bs; //!< Black-Scholes integrand
  Vectors sol;     //!< Internal solution vectors
};

#endif
