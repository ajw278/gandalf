//=================================================================================================
//  MfvCommon.cpp
//  Contains all functions for calculating Meshless Finite-Volume Hydrodynamics quantities.
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


#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <math.h>
#include "Precision.h"
#include "MeshlessFV.h"
#include "Particle.h"
#include "Parameters.h"
#include "SmoothingKernel.h"
#include "EOS.h"
#include "Debug.h"
#include "Exception.h"
#include "InlineFuncs.h"
using namespace std;



//=================================================================================================
//  MfvCommon::MfvCommon
/// MfvCommon class constructor.  Calls main SPH class constructor and also
/// sets additional kernel-related quantities
//=================================================================================================
template <int ndim, template<int> class kernelclass>
MfvCommon<ndim, kernelclass>::MfvCommon
  (int hydro_forces_aux, int self_gravity_aux, FLOAT _accel_mult, FLOAT _courant_mult,
   FLOAT _h_fac, FLOAT h_converge_aux, FLOAT gamma_aux, string gas_eos_aux, string KernelName,
   int size_part, SimUnits &units, Parameters *params):
  MeshlessFV<ndim>(hydro_forces_aux, self_gravity_aux, _accel_mult, _courant_mult, _h_fac,
                   h_converge_aux, gamma_aux, gas_eos_aux, KernelName, size_part, units, params),
  kern(kernelclass<ndim>(KernelName))
{
  this->kernp      = &kern;
  this->kernfac    = (FLOAT) 1.0;
  this->kernfacsqd = (FLOAT) 1.0;
  this->kernrange  = this->kernp->kernrange;
}



//=================================================================================================
//  MfvCommon::~MfvCommon
/// MfvCommon class destructor
//=================================================================================================
template <int ndim, template<int> class kernelclass>
MfvCommon<ndim, kernelclass>::~MfvCommon()
{
}



//=================================================================================================
//  MfvCommon::ComputeH
/// Compute the value of the smoothing length of particle 'i' by iterating the relation :
/// h = h_fac*(m/rho)^(1/ndim).
/// Uses the previous value of h as a starting guess and then uses either a Newton-Rhapson solver,
/// or fixed-point iteration, to converge on the correct value of h.  The maximum tolerance used
/// for deciding whether the iteration has converged is given by the 'h_converge' parameter.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
int MfvCommon<ndim, kernelclass>::ComputeH
 (const int i,                         ///< [in] id of particle
  const int Nneib,                     ///< [in] No. of potential neighbours
  const FLOAT hmax,                    ///< [in] Max. h permitted by neib list
  FLOAT *m,                            ///< [in] Array of neib. masses
  FLOAT *mu,                           ///< [in] Array of m*u (not needed here)
  FLOAT *drsqd,                        ///< [in] Array of neib. distances squared
  FLOAT *gpot,                         ///< [in] Array of neib. grav. potentials
  MeshlessFVParticle<ndim> &part,      ///< [inout] Particle i data
  Nbody<ndim> *nbody)                  ///< [in] Main N-body object
{
  int j;                               // Neighbour id
  int iteration = 0;                   // h-rho iteration counter
  int iteration_max = 30;              // Max. no of iterations
  FLOAT h_lower_bound = (FLOAT) 0.0;   // Lower bound on h
  FLOAT h_upper_bound = hmax;          // Upper bound on h
  FLOAT invhsqd;                       // (1 / h)^2
  FLOAT ssqd;                          // Kernel parameter squared, (r/h)^2


  // If there are sink particles present, check if the particle is inside one.
  // If so, then adjust the iteration bounds and ensure they are valid (i.e. hmax is large enough)
  if (part.sinkid != -1) {
    h_lower_bound = hmin_sink;
    //h_lower_bound = nbody->stardata[part.sinkid].h;  //hmin_sink;
    if (hmax < hmin_sink) return -1;
  }

  // Some basic sanity-checking in case of invalid input into routine
  assert(Nneib > 0);
  assert(hmax > (FLOAT) 0.0);
  assert(!part.flags.is_dead());
  assert(part.m > (FLOAT) 0.0);


  // Main smoothing length iteration loop
  //===============================================================================================
  do {

    // Initialise all variables for this value of h
    iteration++;
    part.ndens    = (FLOAT) 0.0;
    part.invomega = (FLOAT) 0.0;
    part.zeta     = (FLOAT) 0.0;
    part.invh     = (FLOAT) 1.0/part.h;
    part.hfactor  = pow(part.invh,ndim);
    invhsqd       = part.invh*part.invh;

    // Loop over all nearest neighbours in list to calculate density, omega and zeta.
    //---------------------------------------------------------------------------------------------
    for (j=0; j<Nneib; j++) {
      ssqd           = drsqd[j]*invhsqd;
      part.ndens    += kern.w0_s2(ssqd);
      part.invomega += part.invh*kern.womega_s2(ssqd);
      part.zeta     += m[j]*kern.wzeta_s2(ssqd);
    }
    //---------------------------------------------------------------------------------------------

    part.ndens    *= part.hfactor;
    part.invomega *= part.hfactor;
    part.zeta     *= invhsqd;
    part.volume    = (FLOAT) 1.0/part.ndens;
    part.rho       = part.m*part.ndens;
    if (part.rho > (FLOAT) 0.0) part.invrho = (FLOAT) 1.0/part.rho;

    // If h changes below some fixed tolerance, exit iteration loop
    if (part.rho > (FLOAT) 0.0 && part.h > h_lower_bound &&
        fabs(part.h - h_fac*pow(part.volume,MeshlessFV<ndim>::invndim)) < h_converge) break;


    // Use fixed-point iteration, i.e. h_new = h_fac*(m/rho_old)^(1/ndim), for now.  If this does
    // not converge in a reasonable number of iterations (iteration_max), then assume something is
    // wrong and switch to a bisection method, which should be guaranteed to converge, albeit much
    // more slowly.  (N.B. will implement Newton-Raphson soon)
    //---------------------------------------------------------------------------------------------
    if (iteration < iteration_max) {
      part.h = h_fac*pow(part.volume,MeshlessFV<ndim>::invndim);
    }
    else if (iteration == iteration_max) {
      part.h = (FLOAT) 0.5*(h_lower_bound + h_upper_bound);
    }
    else if (iteration < 5*iteration_max) {
      if (part.ndens < small_number || part.ndens*pow(part.h,ndim) > pow(h_fac,ndim)) {
        h_upper_bound = part.h;
      }
      else {
        h_lower_bound = part.h;
      }
      part.h = (FLOAT) 0.5*(h_lower_bound + h_upper_bound);
    }
    else {
      cout << "H ITERATION : " << iteration << "    h : " << part.h
           << "   rho : " << part.rho << "   h_upper " << h_upper_bound << "    hmax :  " << hmax
           << "   h_lower : " << h_lower_bound << "    " << part.hfactor << "    m : " << part.m
           << "     " << part.m*part.hfactor*kern.w0(0.0) << "    " << Nneib << endl;
      string message = "Problem with convergence of h-rho iteration";
      ExceptionHandler::getIstance().raise(message);
    }

    // If the smoothing length is too large for the neighbour list, exit routine and flag neighbour
    // list error in order to generate a larger neighbour list (not properly implemented yet).
    if (part.h > hmax) return 0;

  } while (part.h > h_lower_bound && part.h < h_upper_bound);
  //===============================================================================================


  // Compute other terms once number density and smoothing length are known
  part.h         = max(h_fac*powf(part.volume, (FLOAT) MeshlessFV<ndim>::invndim), h_lower_bound);
  part.invh      = (FLOAT) 1.0/part.h;
  part.hfactor   = pow(part.invh, ndim+1);
  part.hrangesqd = kernfacsqd*kern.kernrangesqd*part.h*part.h;
  part.div_v     = (FLOAT) 0.0;
  part.invomega  = (FLOAT) 1.0 + (FLOAT) MeshlessFV<ndim>::invndim*part.h*part.invomega/part.ndens;
  part.invomega  = (FLOAT) 1.0/part.invomega;
  part.zeta      = -(FLOAT) MeshlessFV<ndim>::invndim*part.h*part.zeta*part.invomega/part.ndens;

  // Calculate the minimum neighbour potential (used later to identify new sinks)
  if (create_sinks == 1) {
    part.flags.set_flag(potmin);
    for (j=0; j<Nneib; j++) {
      if (gpot[j] > (FLOAT) 1.000000001*part.gpot &&
          drsqd[j]*invhsqd < kern.kernrangesqd) part.flags.unset_flag(potmin);
    }
  }

  // Set important thermal variables here
  this->ComputeThermalProperties(part);
  this->UpdatePrimitiveVector(part);
  //this->UpdateArrayVariables(part);


  // If h is invalid (i.e. larger than maximum h), then return error code (0)
  if (part.h <= hmax) return 1;
  else return -1;
}



//=================================================================================================
//  MfvCommon::ComputeDerivatives
/// Compute Psi factors required for computing derivatives needed in Meshess FV equations.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void MfvCommon<ndim, kernelclass>::ComputePsiFactors
 (const int i,                                 ///< [in] id of particle
  const int Nneib,                             ///< [in] No. of neins in neibpart array
  int *neiblist,                               ///< [in] id of gather neibs in neibpart
  FLOAT *drmag,                                ///< [in] Distances of gather neighbours
  FLOAT *invdrmag,                             ///< [in] Inverse distances of gather neibs
  FLOAT *dr,                                   ///< [in] Position vector of gather neibs
  MeshlessFVParticle<ndim> &part,              ///< [inout] Particle i data
  MeshlessFVParticle<ndim> *neibpart)          ///< [inout] Neighbour particle data
{
  int j;                                       // Neighbour list id
  int jj;                                      // Aux. neighbour loop counter
  int k;                                       // Dimension counter
  FLOAT draux[ndim];                           // Relative position vector
  FLOAT drsqd;                                 // Distance squared
  FLOAT E[ndim][ndim];                         // E-matrix for computing normalised B-matrix
  const FLOAT invhsqd = part.invh*part.invh;   // Local copy of 1/h^2

  // Zero all matrices
  for (k=0; k<ndim; k++) {
    for (int kk=0; kk<ndim; kk++) {
      E[k][kk] = (FLOAT) 0.0;
      part.B[k][kk] = (FLOAT) 0.0;
    }
  }


  // Loop over all potential neighbours in the list
  //-----------------------------------------------------------------------------------------------
  for (jj=0; jj<Nneib; jj++) {
    j = neiblist[jj];

    for (k=0; k<ndim; k++) draux[k] = neibpart[j].r[k] - part.r[k];
    drsqd = DotProduct(draux, draux, ndim);

    for (k=0; k<ndim; k++) {
      for (int kk=0; kk<ndim; kk++) {
        E[k][kk] += draux[k]*draux[kk]*part.hfactor*kern.w0_s2(drsqd*invhsqd)/part.ndens;
      }
    }
  }
  //-----------------------------------------------------------------------------------------------


  // Invert the matrix (depending on dimensionality)
  if (ndim == 1) {
    part.B[0][0] = (FLOAT) 1.0/E[0][0];
  }
  else if (ndim == 2) {
    const FLOAT invdet = (FLOAT) 1.0/(E[0][0]*E[1][1] - E[0][1]*E[1][0]);
    part.B[0][0] = invdet*E[1][1];
    part.B[0][1] = -(FLOAT) 1.0*invdet*E[0][1];
    part.B[1][0] = -(FLOAT) 1.0*invdet*E[1][0];
    part.B[1][1] = invdet*E[0][0];
  }
  else if (ndim == 3) {
    const FLOAT invdet = (FLOAT) 1.0/(E[0][0]*(E[1][1]*E[2][2] - E[2][1]*E[1][2]) -
                                      E[0][1]*(E[1][0]*E[2][2] - E[1][2]*E[2][0]) +
                                      E[0][2]*(E[1][0]*E[2][1] - E[1][1]*E[2][0]));
    part.B[0][0] = (E[1][1]*E[2][2] - E[2][1]*E[1][2])*invdet;
    part.B[0][1] = (E[0][2]*E[2][1] - E[0][1]*E[2][2])*invdet;
    part.B[0][2] = (E[0][1]*E[1][2] - E[0][2]*E[1][1])*invdet;
    part.B[1][0] = (E[1][2]*E[2][0] - E[1][0]*E[2][2])*invdet;
    part.B[1][1] = (E[0][0]*E[2][2] - E[0][2]*E[2][0])*invdet;
    part.B[1][2] = (E[1][0]*E[0][2] - E[0][0]*E[1][2])*invdet;
    part.B[2][0] = (E[1][0]*E[2][1] - E[2][0]*E[1][1])*invdet;
    part.B[2][1] = (E[2][0]*E[0][1] - E[0][0]*E[2][1])*invdet;
    part.B[2][2] = (E[0][0]*E[1][1] - E[1][0]*E[0][1])*invdet;
  }


  return;
}



//=================================================================================================
//  MfvCommon::ComputeDerivatives
/// Compute derivatives required for Meshless FV equations.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void MfvCommon<ndim, kernelclass>::ComputeGradients
 (const int i,                         ///< [in] id of particle
  const int Nneib,                     ///< [in] No. of neins in neibpart array
  int *neiblist,                       ///< [in] id of gather neibs in neibpart
  FLOAT *drmag,                        ///< [in] Distances of gather neighbours
  FLOAT *invdrmag,                     ///< [in] Inverse distances of gather neibs
  FLOAT *dr,                           ///< [in] Position vector of gather neibs
  MeshlessFVParticle<ndim> &part,      ///< [inout] Particle i data
  MeshlessFVParticle<ndim> *neibpart)  ///< [inout] Neighbour particle data
{
  int j;                               // Neighbour list id
  int jj;                              // Aux. neighbour counter
  int k;                               // Dimension counter
  int var;                             // Particle state vector summation variable
  FLOAT draux[ndim];                   // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT dv[ndim];                      // Relative velocity vector
  FLOAT dvdr;                          // Dot product of dv and dr
  FLOAT psitilda[ndim];                // Normalised gradient psi factor
  const FLOAT invhsqd = part.invh*part.invh;   // Local copy of 1/h^2


  // Initialise/zero all variables to be updated in this routine
  part.vsig_max = (FLOAT) 0.0;
  for (k=0; k<ndim; k++) part.vreg[k] = (FLOAT) 0.0;
  for (k=0; k<ndim; k++) {
    for (var=0; var<nvar; var++) {
      part.grad[var][k] = (FLOAT) 0.0;
    }
  }

  // Loop over all potential neighbours in the list
  //-----------------------------------------------------------------------------------------------
  for (jj=0; jj<Nneib; jj++) {
    j = neiblist[jj];

    for (k=0; k<ndim; k++) draux[k] = neibpart[j].r[k] - part.r[k];
    for (k=0; k<ndim; k++) dv[k] = neibpart[j].v[k] - part.v[k];
    dvdr = DotProduct(dv, draux, ndim);
    drsqd = DotProduct(draux, draux, ndim);

    // Calculate psitilda values
    for (k=0; k<ndim; k++) psitilda[k] = (FLOAT) 0.0;
    for (k=0; k<ndim; k++) {
      for (int kk=0; kk<ndim; kk++) {
        psitilda[k] += part.B[k][kk]*draux[kk]*part.hfactor*kern.w0_s2(drsqd*invhsqd)/part.ndens;
      }
    }

    // Calculate contribution to gradient from neighbour
    for (var=0; var<nvar; var++) {
      for (k=0; k<ndim; k++) {
        part.grad[var][k] += (neibpart[j].Wprim[var] - part.Wprim[var])*psitilda[k];
      }
    }

    // Calculate maximum signal velocity
    part.vsig_max = max(part.vsig_max, part.sound + neibpart[j].sound -
                                       min((FLOAT) 0.0, dvdr/(sqrtf(drsqd) + small_number)));

    for (k=0; k<ndim; k++) part.vreg[k] -= draux[k]*kern.w0_s2(drsqd*invhsqd);

  }
  //-----------------------------------------------------------------------------------------------

  for (k=0; k<ndim; k++) part.vreg[k] *= part.invh*part.sound;  //pow(part.invh, ndim);



  // Apply the slope limiter
  //-----------------------------------------------------------------------------------------------
  limiter->CellLimiter(part, neibpart, neiblist, Nneib) ;


  assert(part.vsig_max >= part.sound);

  return;
}



//=================================================================================================
//  MfvCommon::CopyDataToGhosts
/// Copy any newly calculated data from original SPH particles to ghosts.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void MfvCommon<ndim, kernelclass>::CopyDataToGhosts
 (DomainBox<ndim> &simbox,
  MeshlessFVParticle<ndim> *partdata)      ///< [inout] Neighbour particle data
{
  int i;                                   // Particle id
  int iorig;                               // Original (real) particle id
  int itype;                               // Ghost particle type
  int j;                                   // Ghost particle counter

  debug2("[MfvCommon::CopyDataToGhosts]");


  //-----------------------------------------------------------------------------------------------
//#pragma omp parallel for default(none) private(i,iorig,itype,j) shared(simbox,sph,partdata)
  for (j=0; j<this->NPeriodicGhost; j++) {
    i = this->Nhydro + j;
    iorig = partdata[i].iorig;
    itype = partdata[i].flags.get();
    assert(itype != none) ;

    partdata[i]        = partdata[iorig];
    partdata[i].iorig  = iorig;
    partdata[i].flags  = type_flag(itype);
    partdata[i].flags.unset_flag(active);


    // Modify ghost position based on ghost type
    // Ghosts of ghosts refer only to their previous ghosts not the base cell, so
    // only update one direction.
    if (ndim > 2) {
      if (itype & z_periodic_lhs) {
        partdata[i].r[2] += simbox.boxsize[2];
        continue ;
      }
      else if (itype & z_periodic_rhs) {
    	partdata[i].r[2] -= simbox.boxsize[2];
        continue ;
      }
      else if (itype & z_mirror_lhs) {
    	partdata[i].reflect(2, simbox.boxmin[2]) ;
        continue ;
      }
      else if (itype & z_mirror_rhs) {
      	partdata[i].reflect(2, simbox.boxmax[2]) ;
        continue ;
      }

    }
    if (ndim > 1) {
      if (itype & y_periodic_lhs) {
    	partdata[i].r[1] += simbox.boxsize[1];
    	continue ;
      }
      else if (itype & y_periodic_rhs) {
    	partdata[i].r[1] -= simbox.boxsize[1];
    	continue ;
      }
      else if (itype & y_mirror_lhs) {
      	partdata[i].reflect(1, simbox.boxmin[1]) ;
    	continue ;
      }
      else if (itype & y_mirror_rhs) {
    	partdata[i].reflect(1, simbox.boxmax[1]) ;
        continue ;
      }
    }

    if (itype & x_periodic_lhs) {
      partdata[i].r[0] += simbox.boxsize[0];
      continue ;
    }
    else if (itype & x_periodic_rhs) {
      partdata[i].r[0] -= simbox.boxsize[0];
      continue ;
    }
    else if (itype & x_mirror_lhs) {
      partdata[i].reflect(0, simbox.boxmin[0]) ;
      continue ;
    }
    else if (itype & x_mirror_rhs) {
      partdata[i].reflect(0, simbox.boxmax[0]) ;
      continue ;
    }
  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  MfvCommon::ComputeSmoothedGravForces
/// Compute SPH neighbour force pairs for
/// (i) All neighbour interactions of particle i with i.d. j > i,
/// (ii) Active neighbour interactions of particle j with i.d. j > i
/// (iii) All inactive neighbour interactions of particle i with i.d. j < i.
/// This ensures that particle-particle pair interactions are computed once only for efficiency.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void MfvCommon<ndim, kernelclass>::ComputeSmoothedGravForces
 (const int i,                         ///< [in] id of particle
  const int Nneib,                     ///< [in] No. of neins in neibpart array
  int *neiblist,                       ///< [in] id of gather neibs in neibpart
  MeshlessFVParticle<ndim> &parti,     ///< [inout] Particle i data
  MeshlessFVParticle<ndim> *neibpart)  ///< [inout] Neighbour particle data
  //MeshlessFVParticle<ndim> &part,      ///< [inout] Particle i data
  //MeshlessFVParticle<ndim> *neib_gen)  ///< [inout] Neighbour particle data
{
  int j;                               // Neighbour list id
  int jj;                              // Aux. neighbour counter
  int k;                               // Dimension counter
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drmag;                         // Distance
  FLOAT invdrmag;                      // 1 / distance
  FLOAT gaux;                          // Aux. grav. potential variable
  FLOAT paux;                          // Aux. pressure force variable
  //MeshlessFVParticle<ndim>& parti = static_cast<MeshlessFVParticle<ndim>& > (part);
  //MeshlessFVParticle<ndim>* neibpart = static_cast<MeshlessFVParticle<ndim>* > (neib_gen);


  // Loop over all potential neighbours in the list
  //-----------------------------------------------------------------------------------------------
  for (jj=0; jj<Nneib; jj++) {
    j = neiblist[jj];
    assert(!neibpart[j].flags.is_dead());

    for (k=0; k<ndim; k++) dr[k] = neibpart[j].r[k] - parti.r[k];
    drmag = sqrt(DotProduct(dr, dr, ndim) + small_number);
    invdrmag = (FLOAT) 1.0/drmag;
    for (k=0; k<ndim; k++) dr[k] *= invdrmag;

    // Main SPH gravity terms
    //---------------------------------------------------------------------------------------------
    paux = (FLOAT) 0.5*(parti.invh*parti.invh*kern.wgrav(drmag*parti.invh) +
                        parti.zeta*parti.hfactor*kern.w1(drmag*parti.invh) +
                        neibpart[j].invh*neibpart[j].invh*kern.wgrav(drmag*neibpart[j].invh) +
                        neibpart[j].zeta*neibpart[j].hfactor*kern.w1(drmag*neibpart[j].invh));
    gaux = (FLOAT) 0.5*(parti.invh*kern.wpot(drmag*parti.invh) +
                        neibpart[j].invh*kern.wpot(drmag*neibpart[j].invh));

    // Add total hydro contribution to acceleration for particle i
    for (k=0; k<ndim; k++) parti.agrav[k] += neibpart[j].m*dr[k]*paux;
    parti.gpot += neibpart[j].m*gaux;

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  MfvCommon::ComputeDirectGravForces
/// Compute the contribution to the total gravitational force of particle 'i'
/// due to 'Nneib' neighbouring particles in the list 'neiblist'.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void MfvCommon<ndim, kernelclass>::ComputeDirectGravForces
 (const int i,                         ///< id of particle
  const int Ndirect,                   ///< No. of nearby 'gather' neighbours
  int *directlist,                     ///< id of gather neighbour in neibpart
  MeshlessFVParticle<ndim> &part,      ///< Particle i data
  MeshlessFVParticle<ndim> *neib_gen)  ///< Neighbour particle data
{
  int j;                               // Neighbour particle id
  int jj;                              // Aux. neighbour loop counter
  int k;                               // Dimension counter
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT invdrmag;                      // 1 / distance
  FLOAT invdr3;                        // 1 / dist^3
  MeshlessFVParticle<ndim>& parti = static_cast<MeshlessFVParticle<ndim>& > (part);
  MeshlessFVParticle<ndim>* neibdata = static_cast<MeshlessFVParticle<ndim>* > (neib_gen);


  // Loop over all neighbouring particles in list
  //-----------------------------------------------------------------------------------------------
  for (jj=0; jj<Ndirect; jj++) {
    j = directlist[jj];
    assert(!neibdata[j].flags.is_dead());

    for (k=0; k<ndim; k++) dr[k] = neibdata[j].r[k] - parti.r[k];
    drsqd    = DotProduct(dr,dr,ndim) + small_number;
    invdrmag = (FLOAT) 1.0/sqrt(drsqd);
    invdr3   = invdrmag*invdrmag*invdrmag;

    // Add contribution to current particle
    for (k=0; k<ndim; k++) parti.agrav[k] += neibdata[j].m*dr[k]*invdr3;
    parti.gpot += neibdata[j].m*invdrmag;

    assert(drsqd >= parti.hrangesqd && drsqd >= neibdata[j].hrangesqd);

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  MfvCommon::ComputeStarGravForces
/// Computes contribution of gravitational force and potential due to stars.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void MfvCommon<ndim, kernelclass>::ComputeStarGravForces
 (const int N,                         ///< [in] No. of stars
  NbodyParticle<ndim> **nbodydata,     ///< [in] Array of star pointers
  MeshlessFVParticle<ndim> &part)      ///< [inout] SPH particle reference
{
  int j;                               // Star counter
  int k;                               // Dimension counter
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drmag;                         // Distance
  FLOAT drsqd;                         // Distance squared
  //FLOAT drdt;                          // Rate of change of relative distance
  //FLOAT dv[ndim];                      // Relative velocity vector
  FLOAT invdrmag;                      // 1 / drmag
  FLOAT invhmean;                      // 1 / hmean
  FLOAT ms;                            // Star mass
  FLOAT paux;                          // Aux. force variable
  MeshlessFVParticle<ndim>& parti = static_cast<MeshlessFVParticle<ndim>& > (part);

  // Loop over all stars and add each contribution
  //-----------------------------------------------------------------------------------------------
  for (j=0; j<N; j++) {

    //if (fixed_sink_mass) ms = msink_fixed;
    //else
    ms = nbodydata[j]->m;

    for (k=0; k<ndim; k++) dr[k] = nbodydata[j]->r[k] - parti.r[k];
    //for (k=0; k<ndim; k++) dv[k] = nbodydata[j]->v[k] - parti.v[k];
    drsqd    = DotProduct(dr,dr,ndim) + small_number;
    drmag    = sqrt(drsqd);
    invdrmag = (FLOAT) 1.0/drmag;
    invhmean = (FLOAT) 2.0/(parti.h + nbodydata[j]->h);
    //drdt     = DotProduct(dv, dr, ndim)*invdrmag;
    paux     = ms*invhmean*invhmean*kern.wgrav(drmag*invhmean)*invdrmag;

    // Add total hydro contribution to acceleration for particle i
    for (k=0; k<ndim; k++) parti.agrav[k] += paux*dr[k];
    //for (k=0; k<ndim; k++) parti.adot[k] += paux*dv[k] - (FLOAT) 3.0*paux*drdt*invdrmag*dr[k] +
    //  (FLOAT) 2.0*twopi*ms*drdt*kern.w0(drmag*invhmean)*powf(invhmean,ndim)*invdrmag*dr[k];
    parti.gpot += ms*invhmean*kern.wpot(drmag*invhmean);

    assert(drmag > (FLOAT) 0.0);
    assert(drmag*invhmean > (FLOAT) 0.0);

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



template class MfvCommon<1, M4Kernel>;
template class MfvCommon<2, M4Kernel>;
template class MfvCommon<3, M4Kernel>;
template class MfvCommon<1, QuinticKernel>;
template class MfvCommon<2, QuinticKernel>;
template class MfvCommon<3, QuinticKernel>;
template class MfvCommon<1, TabulatedKernel>;
template class MfvCommon<2, TabulatedKernel>;
template class MfvCommon<3, TabulatedKernel>;
