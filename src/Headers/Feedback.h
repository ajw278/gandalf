//=================================================================================================
//  Feedback.h
//  ...
//
//  This file is part of GANDALF :
//  Graphical Astrophysics code for N-body Dynamics And Lagrangian Fluids
//  https://github.com/gandalfcode/gandalf
//  Contact : gandalfcode@gmail.com
//
//  Copyright (C) 2013  D. A. Hubber, G. Rosotti
//
//  GANDALF is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 2 of the License, or
//  (at your option) any later version.
//
//  GANDALF is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  General Public License (http://www.gnu.org/licenses) for more details.
//=================================================================================================


#ifndef _FEEDBACK_H_
#define _FEEDBACK_H_


#include <string>
#include "Precision.h"
#include "Constants.h"
#include "InlineFuncs.h"
#include "Hydrodynamics.h"
using namespace std;



//=================================================================================================
//  Class Feedback
/// \brief   ...
/// \details ...
/// \author  ...
/// \date    18/05/2015
//=================================================================================================
template <int ndim>
class Feedback
{
 public:

  Feedback() {};
  virtual ~Feedback() {};

  //virtual int AddNewParticles(FLOAT, FLOAT, Hydrodynamics<ndim> *) = 0;
  virtual int AddWindMassFlux(FLOAT, FLOAT, Hydrodynamics<ndim> *) = 0;

};



//=================================================================================================
//  Class NullFeedback
/// \brief   ...
/// \details ...
/// \author  ...
/// \date    18/05/2015
//=================================================================================================
template <int ndim>
class NullFeedback : public Feedback<ndim>
{
 public:

  NullFeedback() {};
  virtual ~NullFeedback() {};

  //virtual int AddNewParticles(FLOAT, FLOAT, Hydrodynamics<ndim> *) {return 0;}
  virtual int AddWindMassFlux(FLOAT, FLOAT, Hydrodynamics<ndim> *) {};

};



//=================================================================================================
//  Class HotWindFeedback
/// \brief   ...
/// \details ...
/// \author  ...
/// \date    18/05/2015
//=================================================================================================
template <int ndim, template<int> class ParticleType>
class HotWindFeedback : public Feedback<ndim>
{
 public:

  HotWindFeedback();
  virtual ~HotWindFeedback() {};

  //virtual int AddNewParticles(FLOAT, FLOAT, Hydrodynamics<ndim> *);
  virtual int AddWindMassFlux(FLOAT, FLOAT, Hydrodynamics<ndim> *);

};

#endif
