//=================================================================================================
//  LocalIsotherm.cpp
//  Contains functions for a locally isothermal EOS.
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
 
 
 
#include <math.h>
#include "Debug.h"
#include "Constants.h"
#include "EnergyEquation.h"
#include "EOS.h"
#include "Hydrodynamics.h"
#include "Sinks.h"
#include "SphIntegration.h"

using namespace std;
 
 
 
//=================================================================================================
//  LocalIsotherm::LocalIsotherm()
/// LocalIsotherm class constructor
//=================================================================================================
template <int ndim, template <int> class ParticleType>
LocalIsotherm<ndim, ParticleType>::LocalIsotherm(FLOAT temp0aux, FLOAT gamma_aux, FLOAT templaw_aux, FLOAT mu_bar_aux, SimUnits *units, Nbody<ndim> *nbody_aux):
  EnergyEquation<ndim>(1),
  gammam1(gamma_aux-1.0),
  templaw(templaw_aux),
  temp0(temp0aux/units->temp.outscale),
  mu_bar(mu_bar_aux),
  nbody(nbody_aux)
{
}
 
//=================================================================================================
//  LocalIsotherm::EndTimestep
/// Record all important thermal quantities at the end of the step for start of the new timestep.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void LocalIsotherm<ndim,ParticleType>::EndTimestep
 (const int n,                         ///< [in] Integer time in block time struct
  const int Npart,                     ///< [in] Number of particles
  const FLOAT t,                       ///< [in] Current simulation time
  const FLOAT timestep,                ///< [in] Base timestep value
  Particle<ndim>* part_gen)            ///< [inout] Pointer to SPH particle array
{
  int dn;                              // Integer time since beginning of step
  int i;                               // Particle counter                        // .

  //FLOAT dt_therm;                      // ..
  ParticleType<ndim>* partdata = static_cast<ParticleType<ndim>* > (part_gen);
 
  timing->StartTimingSection("LOCAL_ISO");
 
  //-----------------------------------------------------------------------------------------------
//NOT USING 'temp' AS A VARIABLE - WHAT DOES THIS MEAN FOR a) FUNCTIONALITY OF ENERGY INTEGRATION AND b) THE #pragma omp CALL (deleted temp from private for time being)
#pragma omp parallel for default(none) private(dn,i) shared(partdata)
  for (i=0; i<Npart; i++) {
    ParticleType<ndim> &part = partdata[i];
    if (part.flags.is_dead()) continue;
    dn = n - part.nlast;
 
    if (dn == part.nstep) { 
//Is this assignment right? Seems to be cf. 'LocalIsotherm.cpp'. When are u and u0 used? Both are assigned in Radws.
      part.u = Radialu(part.r);
      part.u0 = part.u;
    }
 
  }
  //-----------------------------------------------------------------------------------------------
  timing->EndTimingSection("LOCAL_ISO");
 
  return;
}
 
 
//=================================================================================================
//  LocalIsotherm::Radialu
/// RadialTemp returns temperature due to closest star particle
//=================================================================================================
template <int ndim, template <int> class ParticleType>
FLOAT LocalIsotherm<ndim,ParticleType>::Radialu
 (const FLOAT *partpos)                      
{
  int i=0;
  int j=0;
  FLOAT stardistmin=1e30;
  FLOAT stardist=0.0; 
  NbodyParticle<ndim>** star = nbody->nbodydata;

  for(i=0; i<nbody->Nnbody; i++){
     for(j=0; j<ndim; j++){
       stardist += (partpos[j]-star[i]->r[j])*(partpos[j]-star[i]->r[j]);
     }
     if (stardist<stardistmin) stardistmin= stardist;
     stardist=0.0;
  }
  stardistmin = sqrt(stardistmin) ;

  return temp0*pow(stardistmin,-templaw)/gammam1/mu_bar;
}
 
 
template class LocalIsotherm<1, GradhSphParticle>;
template class LocalIsotherm<2, GradhSphParticle>;
template class LocalIsotherm<3, GradhSphParticle>;
template class LocalIsotherm<1, SM2012SphParticle>;
template class LocalIsotherm<2, SM2012SphParticle>;
template class LocalIsotherm<3, SM2012SphParticle>;
