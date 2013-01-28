// ============================================================================
// SphIntegration.cpp
// ============================================================================

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include "Sph.h"
#include "SphKernel.h"
#include "SphIntegration.h"
#include "SphParticle.h"
#include "EOS.h"
#include "Debug.h"
using namespace std;



// ============================================================================
// SphIntegration::SphIntegration
// ============================================================================
SphIntegration::SphIntegration()
{
}



// ============================================================================
// SphIntegration::~SphIntegration
// ============================================================================
SphIntegration::~SphIntegration()
{
}



// ============================================================================
// ..
// ============================================================================
double SphIntegration::Timestep(SphParticle &part, EOS *eos)
{
  double timestep;
  double amag;

  //Courant condition
  timestep = 0.1*part.h/(eos->SoundSpeed(part) + );

  //Acceleration condition
  amag = sqrt(part.a[0]*part.a[0] + part.a[1]*part.a[1] + part.a[2]*part.a[2]);
  timestep = min(timestep, accel_mult*sqrt(part.h/(amag + small_number_dp));

  //TODO: implement energy condition (once we will be solving the energy equation)

  return timestep;
}