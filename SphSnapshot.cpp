// ============================================================================
// SNAPSHOT.CPP
// ============================================================================


#include <cstdio>
#include <iostream>
#include "SphSnapshot.h"
#include "Sph.h"
#include "SphParticle.h"
using namespace std;


// ============================================================================
// SNAPSHOT::SNAPSHOT
// ============================================================================
SphSnapshot::SphSnapshot()
{
  allocated = false;
  nallocated = 0;
  ndim = 3;
  Nsph = 0;
  t = 0.0f;
}



// ============================================================================
// SNAPSHOT::~SNAPSHOT
// ============================================================================
SphSnapshot::~SphSnapshot()
{
}



// ============================================================================
// SNAPSHOT::ALLOCATEBUFFERMEMORY
// ============================================================================
void SphSnapshot::AllocateBufferMemory(void)
{
  if (ndim == 1) {
    x = new float[Nsph];
    vx = new float[Nsph];
    ax = new float[Nsph];
  }
  else if (ndim == 2) {
    x = new float[Nsph];
    y = new float[Nsph];
    vx = new float[Nsph];
    vy = new float[Nsph];
    ax = new float[Nsph];
    ay = new float[Nsph];
  }
  else if (ndim == 3) {
    x = new float[Nsph];
    y = new float[Nsph];
    z = new float[Nsph];
    vx = new float[Nsph];
    vy = new float[Nsph];
    vz = new float[Nsph];
    ax = new float[Nsph];
    ay = new float[Nsph];
    az = new float[Nsph];
  }
  
  m = new float[Nsph];
  h = new float[Nsph];
  rho = new float[Nsph];
  u = new float[Nsph];

  allocated = true;
  nallocated = 3*ndim + 4;

  return;
}



// ============================================================================
// SNAPSHOT::DEALLOCATEBUFFERMEMORY
// ============================================================================
void SphSnapshot::DeallocateBufferMemory(void)
{
  if (ndim == 1) {
    delete[] x;
    delete[] vx;
    delete[] ax;
  }
  else if (ndim == 2) {
    delete[] x;
    delete[] y;
    delete[] vx;
    delete[] vy;
    delete[] ax;
    delete[] ay;
  }
  else if (ndim == 3) {
    delete[] x;
    delete[] y;
    delete[] z;
    delete[] vx;
    delete[] vy;
    delete[] vz;
    delete[] ax;
    delete[] ay;
    delete[] az;
  }
  
  delete[] m;
  delete[] h;
  delete[] rho;
  delete[] u;

  allocated = false;
  nallocated = 0;

  return;
}



// ============================================================================
// SNAPSHOT::CALCULATEMEMORYUSAGE
// ============================================================================
int SphSnapshot::CalculateMemoryUsage(void)
{
  return Nsph*nallocated*sizeof(float);
}



// ============================================================================
// SNAPSHOT::COPYDATAFROMSIMULATION
// ============================================================================
void SphSnapshot::CopyDataFromSimulation(int ndimaux, int Nsphaux, 
					 SphParticle *sphaux)
{
  Nsph = Nsphaux;
  ndim = ndimaux;

  AllocateBufferMemory();

  for (int i=0; i<Nsph; i++) {

    if (ndim == 1) {
      x[i] = sphaux[i].r[0];
      vx[i] = sphaux[i].v[0];
      ax[i] = sphaux[i].a[0];
    }
    else if (ndim == 2) {
      x[i] = sphaux[i].r[0];
      y[i] = sphaux[i].r[1];
      vx[i] = sphaux[i].v[0];
      vy[i] = sphaux[i].v[1];
      ax[i] = sphaux[i].a[0];
      ay[i] = sphaux[i].a[1];
    }
    else if (ndim == 3) {
      x[i] = sphaux[i].r[0];
      y[i] = sphaux[i].r[1];
      z[i] = sphaux[i].r[2];
      vx[i] = sphaux[i].v[0];
      vy[i] = sphaux[i].v[1];
      vz[i] = sphaux[i].v[2];
      ax[i] = sphaux[i].a[0];
      ay[i] = sphaux[i].a[1];
      az[i] = sphaux[i].a[2];
    }

    m[i] = sphaux[i].m;
    h[i] = sphaux[i].h;
    rho[i] = sphaux[i].rho;
    u[i] = sphaux[i].u;

  }

  return;
}



// ============================================================================
// SphSnapshot::GETARRAY
// ============================================================================
void SphSnapshot::ExtractArray(string name, float** out_array, 
			       int* size_array)
{
  if (name == "x") *out_array = x;
  if (name == "y") *out_array = y;
  if (name == "z") *out_array = z;
  if (name == "vx") *out_array = vx;
  if (name == "vy") *out_array = vy;
  if (name == "vz") *out_array = vz;
  if (name == "ax") *out_array = ax;
  if (name == "ay") *out_array = ay;
  if (name == "az") *out_array = az;
  if (name == "m") *out_array = m;
  if (name == "h") *out_array = h;
  if (name == "rho") *out_array = rho;
  if (name == "u") *out_array = u;

  *size_array = Nsph;

  return;
}