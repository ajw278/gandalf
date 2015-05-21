//=================================================================================================
//  Feedback.cpp
//  Contains all functions for managing the tree for hydrodynamical particles.
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


#include <assert.h>
#include <iostream>
#include <string>
#include "Feedback.h"
#include "chealpix.h"
#if defined _OPENMP
#include <omp.h>
#endif
using namespace std;



//=================================================================================================
//  HotWindFeedback::HotWindFeedback
/// ...
//=================================================================================================
template <int ndim, template<int> class ParticleType>
HotWindFeedback<ndim,ParticleType>::HotWindFeedback()
{
  ifirstshell = -1;
  ilastshell = -1;
  tlastshell = 0.0;
}



//=================================================================================================
//  HotWindFeedback::AddNewParticles
/// ...
//=================================================================================================
/*template <int ndim, template<int> class ParticleType>
int HotWindFeedback<ndim,ParticleType>::AddNewParticles
 (FLOAT t,
  FLOAT timestep,
  Hydrodynamics<ndim> *hydro)
{

  return 0;
}*/



//=================================================================================================
//  HotWindFeedback::AddMassFlux
/// ...
//=================================================================================================
template <int ndim, template<int> class ParticleType>
int HotWindFeedback<ndim,ParticleType>::AddWindMassFlux
 (FLOAT t,
  FLOAT timestep,
  Hydrodynamics<ndim> *hydro)
{
  FLOAT rsource[ndim];
  ParticleType<ndim> *partdata = static_cast<ParticleType<ndim>* > (hydro->GetParticleArray());

  int i,k;
  int level = 0;
  //int Nwind = 12;
  long nside = (long) (pow(2,level));
  int Nwind = (int) (12*nside*nside);


  cout << "WINDS!!" << endl;


  int Nnew = 0;
  FLOAT mwind = 0.0001;
  FLOAT mdot = 0.01;
  FLOAT vwind = 10.0;
  FLOAT vshell;
  FLOAT rshell = 0.01;
  FLOAT mu_ion = 0.5;
  FLOAT kb = 1.0;
  FLOAT gamma = 1.6666666666;
  double rnew[ndim];
  FLOAT dt;
  for (k=0; k<ndim; k++) rsource[k] = 0.0;

  FLOAT dt_wind = (FLOAT) Nwind*mwind/mdot;
  vshell = rshell/dt_wind;
  FLOAT Twind = 3.0*mu_ion*vwind*vwind/16.0/kb;
  FLOAT uwind = Twind / (gamma - 1.0);


  //assert(dt_wind < timestep);


  // First, add mass to existing particles
  if (ifirstshell != -1) {

    if (t > tlastshell + dt_wind) dt = t - tlastshell + dt_wind;
    else dt = timestep;

    for (i=ifirstshell; i<=ilastshell; i++) {
      for (k=0; k<ndim; k++) partdata[i].r0[k] = partdata[i].r[k];
      for (k=0; k<ndim; k++) partdata[i].v0[k] = partdata[i].v[k];
      for (k=0; k<ndim; k++) partdata[i].r[k] = rsource[k] + partdata[i].v[k]*(t - tlastshell);
      partdata[i].m += mdot*dt/(FLOAT) Nwind;
      if (t > tlastshell + dt_wind) partdata[i].itype = gas;
    }

  }


  // Add new particles if we need to
  //-----------------------------------------------------------------------------------------------
  if (t > tlastshell) {

    ifirstshell = hydro->Nhydro;
    ilastshell = hydro->Nhydro + Nwind - 1;

    for (long ipix=0; ipix<Nwind; ipix++) {
      int inew = hydro->Nhydro + ipix;
      pix2vec_nest(nside, ipix, rnew);
      for (k=0; k<ndim; k++) partdata[inew].v[k] = vshell*rnew[k];
      for (k=0; k<ndim; k++) partdata[inew].r[k] = rsource[k] + partdata[inew].v[k]*(t - tlastshell);

      ParticleType<ndim> &part = partdata[inew];
      partdata[inew].m     = mdot*(t - tlastshell)/(FLOAT) Nwind;
      partdata[inew].u     = uwind;
      partdata[inew].itype = wind;
      partdata[inew].h     = 0.5*rshell;
      partdata[inew].rho   = part.m/(4.0*pi*pow(0.5*rshell,3)/3.0);

      cout << "New particle : " << part.m << "    " << part.u << "    " << vshell << endl;
      cout << "rnew         : " << rnew[0] << "   " << rnew[1] << "   " << rnew[2] << endl;

    }


    tlastshell = tlastshell + dt_wind;
    hydro->Nhydro += Nwind;
    Nnew = Nwind;

    cout << "Created new wind particles : " << Nwind << "   " << hydro->Nhydro << endl;

  }


  cout << "t : " << t << "   " << tlastshell << "   " << dt_wind << endl;


  return Nnew;
}



template class HotWindFeedback<1, GradhSphParticle>;
template class HotWindFeedback<2, GradhSphParticle>;
template class HotWindFeedback<3, GradhSphParticle>;

template class HotWindFeedback<1, MeshlessFVParticle>;
template class HotWindFeedback<2, MeshlessFVParticle>;
template class HotWindFeedback<3, MeshlessFVParticle>;
