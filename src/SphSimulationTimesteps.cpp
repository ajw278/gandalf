//=============================================================================
//  SphSimulationTimesteps.cpp
//  Contains all functions for controlling the SPH simulation timestep structure
//=============================================================================


#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <math.h>
#include <cstdio>
#include <cstring>
#include "Precision.h"
#include "Exception.h"
#include "SphSimulation.h"
#include "Parameters.h"
#include "InlineFuncs.h"
#include "Debug.h"
using namespace std;



//=============================================================================
//  SphSimulation::ComputeGlobalTimestep
//  Computes global timestep for SPH simulation.
//=============================================================================
template <int ndim>
void SimulationDim<ndim>::ComputeGlobalTimestep(void)
{
  int i;                            // Particle counter
  DOUBLE dt = big_number_dp;        // Particle timestep
  DOUBLE dt_min = big_number_dp;    // Local copy of minimum timestep

  debug2("[SphSimulation::ComputeGlobalTimestep]");

  // --------------------------------------------------------------------------
  if (n == nresync) {

    n = 0;
    level_max = 0;
    level_step = level_max + integration_step - 1;
    nresync = integration_step;

    // Find minimum timestep from all SPH particles
    // ------------------------------------------------------------------------
#pragma omp parallel default(shared) private(i,dt)
    {
#pragma omp for
      for (i=0; i<sph->Nsph; i++)
        dt = min(dt,sphint->Timestep(sph->sphdata[i],sph->hydro_forces));

      // If integrating energy equation, include energy timestep
      if (simparams->stringparams["gas_eos"] == "energy_eqn") {
#pragma omp for
        for (i=0; i<sph->Nsph; i++)
          dt = min(dt,uint->Timestep(sph->sphdata[i]));
      }

#pragma omp critical
      if (dt < dt_min) dt_min = dt;
    }
    // ------------------------------------------------------------------------


    // Now compute minimum timestep due to stars/systems
    for (i=0; i<nbody->Nstar; i++)
      dt_min = min(dt_min,nbody->Timestep(nbody->nbodydata[i]));

    
    // Set all particles to same timestep
    timestep = dt_min;
    for (i=0; i<sph->Nsph; i++) {
      sph->sphdata[i].level = 0;
      sph->sphdata[i].nstep = pow(2,level_step - sph->sphdata[i].level);
      sph->sphdata[i].dt = timestep;
    }
    for (i=0; i<nbody->Nnbody; i++) {
      nbody->nbodydata[i]->level = 0;
      nbody->nbodydata[i]->nstep = pow(2,level_step - nbody->nbodydata[i]->level);
      nbody->nbodydata[i]->dt = timestep;
    }

  }
  // --------------------------------------------------------------------------

  return;
}



// ============================================================================
// SphSimulation::ComputeBlockTimesteps
// ..
// ============================================================================
template <int ndim>
void SimulationDim<ndim>::ComputeBlockTimesteps(void)
{
  int i;                                    // Particle counter
  int istep;                                // ??
  int level;                                // Particle timestep level
  int last_level;                           // Previous timestep level
  int level_max_old;                        // Old level_max
  int level_max_sph = 0;                    // level_max for SPH particles only
  int nstep;                                // ??
  DOUBLE dt;                                // Aux. timestep variable
  DOUBLE dt_max_sph = 0.0;                  // Maximum SPH particle timestep

  debug2("[SphSimulation::ComputeBlockTimesteps]");

  timestep = big_number;

  // Synchronise all timesteps and reconstruct block timestep structure.
  // ==========================================================================
  if (n == nresync) {

    n = 0;

    // Find minimum timestep from all SPH particles
    for (i=0; i<sph->Nsph; i++) {
      dt = sphint->Timestep(sph->sphdata[i],sph->hydro_forces);
      if (dt < timestep) timestep = dt;
      if (dt > dt_max_sph) dt_max_sph = dt;
      sph->sphdata[i].dt = dt;
    }
    
    // If integrating energy equation, include energy timestep
    if (sph->gas_eos == "energy_eqn") {
      for (i=0; i<sph->Nsph; i++) {
	dt = uint->Timestep(sph->sphdata[i]);
	if (dt < timestep) timestep = dt;
	sph->sphdata[i].dt = min(sph->sphdata[i].dt,dt);
      }
    }

    // Calculate new block timestep levels
    level_max = Nlevels - 1;
    level_step = level_max + integration_step - 1;
    dt_max = timestep*powf(2.0,level_max);
    nresync = pow(2,level_step);
    timestep = dt_max / (DOUBLE) nresync;
    
    // Calculate the maximum level occupied by all SPH particles
    level_max_sph = max((int) (invlogetwo*log(dt_max/dt_max_sph)) + 1, 0);

    // If enforcing a single SPH timestep, set it here.  Otherwise, populate 
    // the timestep levels with SPH particles.
    if (sph_single_timestep == 1)
      for (i=0; i<sph->Nsph; i++) sph->sphdata[i].level = level_max_sph;
    else {
      for (i=0; i<sph->Nsph; i++) {
	level = min((int) (invlogetwo*log(dt_max/dt)) + 1, level_max);
	level = max(level,0);
	sph->sphdata[i].level = level;
	sph->sphdata[i].dt = 
	  (DOUBLE) pow(2,level_step - sph->sphdata[i].level)*timestep;
      }
    }

  }

  // If not resynchronising, check if any SPH particles need to move up 
  // or down timestep levels
  // ==========================================================================
  else {

    level_max_old = level_max;
    level_max = 0;

    // Find all SPH particles at the beginning of a new timestep
    // ------------------------------------------------------------------------
    for (i=0; i<sph->Nsph; i++) {
      last_level = sph->sphdata[i].level;

      nstep = pow(2,level_step - last_level);
      istep = pow(2,level_step - last_level + 1);

      // Skip particles that are not at end of step
      if (n%nstep == 0) {
	last_level = sph->sphdata[i].level;
	dt = sphint->Timestep(sph->sphdata[i],sph->hydro_forces);
	if (sph->gas_eos == "energy_eqn") 
	  dt = min(dt,uint->Timestep(sph->sphdata[i]));
	sph->sphdata[i].dt = dt;
	level = min((int) (invlogetwo*log(dt_max/dt)) + 1, level_max);
	level = max(level,0);

	// Move up one level (if levels are correctly synchronised) or 
	// down several levels if required
	if (level < last_level && last_level > 1 && n%istep == 0) 
	  sph->sphdata[i].level--;
	else if (level > last_level) {
	  sph->sphdata[i].level = level;
	}
      }

      // Find maximum level of all SPH particles
      level_max_sph = max(level_max_sph,sph->sphdata[i].level);
      level_max = max(level_max,sph->sphdata[i].level);
    }
    // ------------------------------------------------------------------------
      

    // Set fixed SPH timestep level here in case maximum has changed
    if (sph_single_timestep == 1)
      for (i=0; i<sph->Nsph; i++) sph->sphdata[i].level = level_max_sph;


    // Update all timestep variables if we have removed or added any levels
    // ------------------------------------------------------------------------
    if (level_max != level_max_old) {

      // Increase maximum timestep level if correctly synchronised
      istep = pow(2,level_step - level_max_old + 1);
      if (level_max <= level_max_old - 1 && level_max_old > 1 && n%istep == 0)
	level_max = level_max_old - 1;
      else if (level_max == level_max_old)
	level_max = level_max_old;
      level_step = level_max + integration_step - 1;

      // Adjust integer time if levels added or removed
      if (level_max > level_max_old)
	n *= pow(2,level_max - level_max_old);
      else if (level_max < level_max_old)
	n /= pow(2,level_max_old - level_max);

    }
    // ------------------------------------------------------------------------

    nresync = pow(2,level_step);
    timestep = dt_max / (DOUBLE) nresync;

    for (i=0; i<sph->Nsph; i++) sph->sphdata[i].dt = 
      (DOUBLE) pow(2,level_step - sph->sphdata[i].level)*timestep;

  }
  // ==========================================================================

#if defined(VERIFY_ALL)
  VerifyBlockTimesteps();
#endif

  return;
}



#if defined(VERIFY_ALL)
// ============================================================================
// SphSimulation::VerifyBlockTimesteps
// Perform simple sanity check on block timesteps to check if 
// values are consistent.
// ============================================================================
template <int ndim>
void SimulationDim<ndim>::VerifyBlockTimesteps(void)
{
  debug2("[SphSimulation::VerifyBlockTimesteps]");

  // Check integer timestep variables are valid
  if (n < 0 || n > nresync) {
    cout << "Invalid integer timestep value : " 
	 << n << "   " << nresync << endl;
    exit(0);
  }

  // Check all particles occupy valid timestep levels
  for (int i=0; i<sph->Nsph; i++) {
    if (sph->sphdata[i].level > level_max && sph->sphdata[i].level < 0) {
      cout << "Invalid SPH timestep level : " << i << "   "
	   << sph->sphdata[i].level << "   " << level_max << endl;
      exit(0);
    }
  }

  return;
}
#endif
