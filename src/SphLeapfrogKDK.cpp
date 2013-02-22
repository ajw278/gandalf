// ============================================================================
// SphLeapfrogKDK.cpp
// ============================================================================


#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include "Dimensions.h"
#include "Sph.h"
#include "SphKernel.h"
#include "SphIntegration.h"
#include "SphParticle.h"
#include "EOS.h"
#include "Debug.h"
using namespace std;



// ============================================================================
// SphLeapfrogKDK::SphLeapfrogKDK()
// ============================================================================
SphLFKDK::SphLFKDK(int ndimaux, int vdimaux, 
		   double accel_mult_aux, double courant_mult_aux) :
  SphIntegration(ndimaux, vdimaux, accel_mult_aux, courant_mult_aux)
{
}



// ============================================================================
// SphLeapfrogKDK::~SphLeapfrog()
// ============================================================================
SphLFKDK::~SphLFKDK()
{
}



// ============================================================================
// SphLeapfrogKDK::AdvanceParticles
// ============================================================================
void SphLFKDK::AdvanceParticles(int n, int level_step, 
				int Nsph, SphParticle *sph, double dt)
{
  int i;
  int k;
  int nstep;

  debug2("[SphLFKDK::AdvanceParticles]");

  for (i=0; i<Nsph; i++) {
    nstep = pow(2,level_step - sph[i].level);
    for (k=0; k<ndim; k++) sph[i].r[k] +=
      sph[i].v[k]*dt + 0.5*sph[i].a[k]*dt*dt;
    for (k=0; k<vdim; k++) sph[i].v[k] += sph[i].a[k]*dt;
    if (n%nstep == 0) sph[i].active = true;
  }

  return;
}
 


// ============================================================================
// SphLeapfrogKDK::CorrectionTerms
// ============================================================================
void SphLFKDK::CorrectionTerms(int n, int level_step, 
			       int Nsph, SphParticle *sph, double dt)
{
  int i;
  int k;
  int nstep;

  debug2("[SphLFKDK::CorrectionTerms]");

  for (i=0; i<Nsph; i++) {
    nstep = pow(2,level_step - sph[i].level);
    if (n%nstep == 0)
      for (k=0; k<ndim; k++) 
	sph[i].v[k] += 0.5*(sph[i].a[k] - sph[i].a0[k])*dt*(DOUBLE) nstep;
  }

  return;
}



// ============================================================================
// SphLeapfrogKDK::EndTimestep
// ============================================================================
void SphLFKDK::EndTimestep(int n, int level_step, int Nsph, SphParticle *sph)
{
  int i,k,nstep;

  debug2("[SphLFKDK::EndTimestep]");

  for (i=0; i<Nsph; i++) {
    nstep = pow(2,level_step - sph[i].level);
    if (n%nstep == 0) {
      for (k=0; k<vdim; k++) sph[i].a0[k] = sph[i].a[k];
      sph[i].active = false;
    }
  }

  return;
}