//=================================================================================================
//  ForcelessSimulation.cpp
//  Contains all main functions controlling force free simulation work-flow.
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


#include <iostream>
#include <iomanip>
#include <ostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <time.h>
#include <cstdio>
#include <cstring>
#include "Precision.h"
#include "CodeTiming.h"
#include "Exception.h"
#include "Debug.h"
#include "InlineFuncs.h"
#include "Simulation.h"
#include "Parameters.h"
#include "Nbody.h"
#include "RandomNumber.h"
#include "Sph.h"
#include "RiemannSolver.h"
#include "Ghosts.h"
#include "Sinks.h"
using namespace std;



// Create template class instances of the main SphSimulation object for
// each dimension used (1, 2 and 3)
template class ForcelessSimulation<1>;
template class ForcelessSimulation<2>;
template class ForcelessSimulation<3>;



//=================================================================================================
//  ForcelessSimulation::ProcessParameters
/// Process all the options chosen in the parameters file, setting various
/// simulation variables and creating important simulation objects.
//=================================================================================================
template <int ndim>
void ForcelessSimulation<ndim>::ProcessParameters(void)
{
  // Local references to parameter variables for brevity
  map<string, int> &intparams = simparams->intparams;
  map<string, double> &floatparams = simparams->floatparams;
  map<string, string> &stringparams = simparams->stringparams;
  string sim = stringparams["sim"];
  string gas_eos = stringparams["gas_eos"];
  string gas_radiation = stringparams["radiation"];

  debug2("[ForcelessSimulation::ProcessParameters]");


  // Sanity check for valid dimensionality
  if (ndim < 1 || ndim > 3) {
    std::ostringstream message;
    message << "Invalid dimensionality chosen : ndim = " << ndim;
    ExceptionHandler::getIstance().raise(message.str());
  }

  // Set-up random number generator object
  //-----------------------------------------------------------------------------------------------
  if (stringparams["rand_algorithm"] == "xorshift") {
    randnumb = new XorshiftRand(intparams["randseed"]);
  }
  else if (stringparams["rand_algorithm"] == "none") {
    randnumb = new DefaultSystemRand(intparams["randseed"]);
  }
  else {
    string message = "Unrecognised parameter : rand_algorithm= " +
      stringparams["rand_algorithm"];
    ExceptionHandler::getIstance().raise(message);
  }


  // Set-up all output units for scaling parameters
  simunits.SetupUnits(simparams);


  // Boundary condition variables
  //-----------------------------------------------------------------------------------------------
  simbox.boundary_lhs[0] = setBoundaryType(stringparams["boundary_lhs[0]"]);
  simbox.boundary_rhs[0] = setBoundaryType(stringparams["boundary_rhs[0]"]);
  simbox.boxmin[0] = floatparams["boxmin[0]"]/simunits.r.outscale;
  simbox.boxmax[0] = floatparams["boxmax[0]"]/simunits.r.outscale;

  if (ndim > 1) {
    simbox.boundary_lhs[1] = setBoundaryType(stringparams["boundary_lhs[1]"]);
    simbox.boundary_rhs[1] = setBoundaryType(stringparams["boundary_rhs[1]"]);
    simbox.boxmin[1] = floatparams["boxmin[1]"]/simunits.r.outscale;
    simbox.boxmax[1] = floatparams["boxmax[1]"]/simunits.r.outscale;
  }

  if (ndim == 3) {
    simbox.boundary_lhs[2] = setBoundaryType(stringparams["boundary_lhs[2]"]);
    simbox.boundary_rhs[2] = setBoundaryType(stringparams["boundary_rhs[2]"]);
    simbox.boxmin[2] = floatparams["boxmin[2]"]/simunits.r.outscale;
    simbox.boxmax[2] = floatparams["boxmax[2]"]/simunits.r.outscale;
  }

  for (int k=0; k<ndim; k++) {
    simbox.boxsize[k] = simbox.boxmax[k] - simbox.boxmin[k];
    simbox.boxhalf[k] = 0.5*simbox.boxsize[k];
  }

  // Process all N-body parameters and set-up main N-body objects
  this->ProcessNbodyParameters();

  // Set-up main SPH objects depending on which SPH algorithm we are using
  ProcessSphParameters();
  hydro = sph;


  // Set external potential field object and set pointers to object
  if (stringparams["external_potential"] == "none") {
    extpot = new NullPotential<ndim>();
  }
  else if (stringparams["external_potential"] == "vertical") {
    extpot = new VerticalPotential<ndim>
      (intparams["kgrav"], floatparams["avert"], simbox.boxmin[intparams["kgrav"]]);
  }
  else if (stringparams["external_potential"] == "plummer") {
    extpot = new PlummerPotential<ndim>(floatparams["mplummer"], floatparams["rplummer"]);
  }
  else {
    string message = "Unrecognised parameter : external_potential = "
      + simparams->stringparams["external_potential"];
    ExceptionHandler::getIstance().raise(message);
  }

  sph->extpot = extpot;
  nbody->extpot = extpot;


  // Create Ewald periodic gravity object
  periodicBoundaries = IsAnyBoundaryPeriodic(simbox);
  if (periodicBoundaries && intparams["self_gravity"] == 1) {
    ewaldGravity = true;
    ewald = new Ewald<ndim>
      (simbox, intparams["gr_bhewaldseriesn"], intparams["in"], intparams["nEwaldGrid"],
       floatparams["ewald_mult"], floatparams["ixmin"], floatparams["ixmax"],
       floatparams["EFratio"], timing);
    simbox.PeriodicGravity = true ;
  }
  else{
    simbox.PeriodicGravity = false ;
    if (IsAnyBoundaryReflecting(simbox) && intparams["self_gravity"]){
      ExceptionHandler::getIstance().raise("Error: Reflecting boundaries and self-gravity is not "
    		                               "supported") ;
    }
  }


  // Set all other SPH parameter variables
  sph->Nhydromax       = intparams["Nhydromax"];
  sph->create_sinks    = intparams["create_sinks"];
  sph->fixed_sink_mass = intparams["fixed_sink_mass"];
  sph->msink_fixed     = floatparams["m1"];
  sph->alpha_visc_min  = floatparams["alpha_visc_min"];


  // Set important variables for N-body objects
  nbody->Nstarmax       = intparams["Nstarmax"];
  nbody_single_timestep = intparams["nbody_single_timestep"];
  nbodytree.gpehard     = floatparams["gpehard"];
  nbodytree.gpesoft     = floatparams["gpesoft"];
  //nbody->perturbers     = intparams["perturbers"];
  //if (intparams["sub_systems"] == 1) subsystem->perturbers = intparams["perturbers"];



  // Sink particles
  //-----------------------------------------------------------------------------------------------
  sinks = new Sinks<ndim>(sphneib);
  sink_particles             = intparams["sink_particles"];
  sinks->sink_particles      = intparams["sink_particles"];
  sinks->create_sinks        = intparams["create_sinks"];
  sinks->smooth_accretion    = intparams["smooth_accretion"];
  sinks->alpha_ss            = floatparams["alpha_ss"];
  sinks->smooth_accrete_frac = floatparams["smooth_accrete_frac"];
  sinks->smooth_accrete_dt   = floatparams["smooth_accrete_dt"];
  sinks->sink_radius_mode    = stringparams["sink_radius_mode"];
  sinks->rho_sink            = floatparams["rho_sink"];
  sinks->rho_sink            /= (simunits.rho.outscale*simunits.rho.outcgs);

  if (sinks->sink_radius_mode == "fixed") {
    sinks->sink_radius = floatparams["sink_radius"]/simunits.r.outscale;
  }
  else {
    sinks->sink_radius = floatparams["sink_radius"];
  }

  // Sanity-check for various sink particle values
  if (intparams["sink_particles"] == 1 &&
      (stringparams["nbody"] != "lfkdk" && stringparams["nbody"] != "lfdkd")) {
    string message = "Invalid parameter : nbody must use lfkdk or lfdkd when "
      "using accreting sink particles";
    ExceptionHandler::getIstance().raise(message);
  }

#if defined MPI_PARALLEL
  sinks->SetMpiControl(mpicontrol);
#endif

  // Set other important simulation variables
  dt_litesnap         = floatparams["dt_litesnap"]/simunits.t.outscale;
  dt_python           = floatparams["dt_python"];
  dt_snap             = floatparams["dt_snap"]/simunits.t.outscale;
  extra_sink_output   = intparams["extra_sink_output"];
  level_diff_max      = intparams["level_diff_max"];
  litesnap            = intparams["litesnap"];
  Nlevels             = intparams["Nlevels"];
  ndiagstep           = intparams["ndiagstep"];
  noutputstep         = intparams["noutputstep"];
  nradstep            = intparams["nradstep"];
  nrestartstep        = intparams["nrestartstep"];
  ntreebuildstep      = intparams["ntreebuildstep"];
  ntreestockstep      = intparams["ntreestockstep"];
  Nstepsmax           = intparams["Nstepsmax"];
  out_file_form       = stringparams["out_file_form"];
  pruning_level_min   = intparams["pruning_level_min"];
  pruning_level_max   = intparams["pruning_level_max"];
  run_id              = stringparams["run_id"];
  sph_single_timestep = intparams["sph_single_timestep"];
  tmax_wallclock      = floatparams["tmax_wallclock"];
  tend                = floatparams["tend"]/simunits.t.outscale;
  tlitesnapnext       = floatparams["tlitesnapfirst"]/simunits.t.outscale;
  tsnapnext           = floatparams["tsnapfirst"]/simunits.t.outscale;
  tdyn_mult			  = floatparams["tdyn_mult"];


  // Set pointers to timing object
  nbody->timing   = timing;
  if (sim == "sph" || sim == "gradhsph" || sim == "sm2012sph" || sim=="forceless") {
    sinks->timing    = timing;
    sphint->timing  = timing;
    sphneib->timing = timing;
    uint->timing    = timing;
    radiation->timing = timing;
  }


  // Flag that we've processed all parameters already
  ParametersProcessed = true;


  return;
}



//=================================================================================================
//  ForcelessSimulation::PostInitialConditionsSetup
/// Call routines for calculating all initial SPH and N-body quantities
/// once initial conditions have been set-up.
//=================================================================================================
template <int ndim>
void ForcelessSimulation<ndim>::PostInitialConditionsSetup(void)
{
  int i;                               // Particle counter
  int k;                               // Dimension counter
  FLOAT adot[ndim];                    // Dummy adot array
  SphParticle<ndim> *partdata;         // Pointer to main SPH data array
  debug2("[ForcelessSimulation::PostInitialConditionsSetup]");

  // Set iorig
  if (rank == 0) {
    for (i=0; i<sph->Nhydro; i++) sph->GetSphParticlePointer(i).iorig = i;
  }

  // Copy information from the stars to the sinks
  if (simparams->intparams["sink_particles"]==1) {
    sinks->Nsink = nbody->Nstar;
    sinks->AllocateMemory(sinks->Nsink);
    for (int i=0; i<sinks->Nsink; i++) {
      sinks->sink[i].star   = &(nbody->stardata[i]);
      sinks->sink[i].istar  = i;
      sinks->sink[i].radius = hydro->kernp->kernrange*nbody->stardata[i].h;
      //sinks->sink[i].radius = simparams->floatparams["sink_radius"];
    }
  }

  // Perform initial MPI decomposition
  //-----------------------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
  mpicontrol->CreateInitialDomainDecomposition(sph, nbody, simparams, simbox,
                                               this->initial_h_provided);
  this->AllocateParticleMemory();
#endif

  // Set pointer to SPH particle data
  partdata = sph->GetSphParticleArray();

  // Set time variables here (for now)
  nresync = 0;
  n = 0;

  // Set initial smoothing lengths and create initial ghost particles
  //-----------------------------------------------------------------------------------------------
  sph->Nghost = 0;
  sph->Nghostmax = sph->Nhydromax - sph->Nhydro;
  sph->Ntot = sph->Nhydro;
  for (i=0; i<sph->Nhydro; i++) sph->GetSphParticlePointer(i).flags.set_flag(active);

  // Set initial artificial viscosity alpha values
  if (simparams->stringparams["time_dependent_avisc"] == "none") {
    for (i=0; i<sph->Nhydro; i++) sph->GetSphParticlePointer(i).alpha = sph->alpha_visc;
  }
  else {
    for (i=0; i<sph->Nhydro; i++) sph->GetSphParticlePointer(i).alpha = sph->alpha_visc_min;
  }

  // Compute instantaneous mean mass (used for smooth sink accretion)
  sph->mmean = (FLOAT) 0.0;
  for (i=0; i<sph->Nhydro; i++) sph->mmean += sph->GetSphParticlePointer(i).m;
  sph->mmean /= (FLOAT) sph->Nhydro;

  // Compute minimum smoothing length of sinks
  sph->hmin_sink = big_number;
  for (i=0; i<sinks->Nsink; i++) {
    sph->hmin_sink = min(sph->hmin_sink, (FLOAT) sinks->sink[i].star->h);
  }

  // If the smoothing lengths have not been provided beforehand, then
  // calculate the initial values here
  //sphneib->neibcheck = false;
  if (!this->initial_h_provided) {
    sph->InitialSmoothingLengthGuess();
    sphneib->BuildTree(rebuild_tree, 0, ntreebuildstep, ntreestockstep,
                       sph->Ntot, sph->Nhydromax, timestep, partdata, sph);
    sphneib->UpdateAllSphProperties(sph->Nhydro, sph->Ntot, partdata, sph, nbody);
  }
  else {
    sphneib->BuildTree(rebuild_tree, 0, ntreebuildstep, ntreestockstep,
                       sph->Ntot, sph->Nhydromax, timestep, partdata, sph);
  }

#ifdef MPI_PARALLEL
  mpicontrol->UpdateAllBoundingBoxes(sph->Nhydro, sph, sph->kernp);
#endif

  // Search ghost particles
  sphneib->SearchBoundaryGhostParticles((FLOAT) 0.0, simbox, sph);
  sphneib->BuildGhostTree(true, 0, ntreebuildstep, ntreestockstep, sph->Ntot,
                          sph->Nhydromax, timestep, partdata, sph);
#ifdef MPI_PARALLEL
  mpicontrol->UpdateAllBoundingBoxes(sph->Nhydro+sph->NPeriodicGhost, sph, sph->kernp);
  for (int i=0; i<sph->Nhydro+sph->NPeriodicGhost; i++) {
    SphParticle<ndim>& parti = sph->GetSphParticlePointer(i);
    parti.hrangesqd = sph->kernfacsqd*sph->kernp->kernrangesqd*parti.h*parti.h;
  }
  MpiGhosts->SearchGhostParticles((FLOAT) 0.0, simbox, sph);
  sphneib->BuildMpiGhostTree(true, 0, ntreebuildstep, ntreestockstep, sph->Ntot,
                             sph->Nhydromax, timestep, partdata, sph);
#endif

  // Zero accelerations
  for (i=0; i<sph->Nhydro; i++) sph->GetSphParticlePointer(i).flags.set_flag(active);

  // Update neighbour tree
  rebuild_tree = true;
  sphneib->BuildTree(rebuild_tree, 0, ntreebuildstep, ntreestockstep,
                     sph->Ntot, sph->Nhydromax, timestep, partdata, sph);
  level_step = 1;


  // For Eigenvalue MAC, need non-zero values
  for (i=0; i<sph->Nhydro; i++) sph->GetSphParticlePointer(i).gpot = big_number;

  // Calculate all SPH properties
  sphneib->UpdateAllSphProperties(sph->Nhydro, sph->Ntot, partdata, sph, nbody);

#ifdef MPI_PARALLEL
  mpicontrol->UpdateAllBoundingBoxes(sph->Nhydro, sph, sph->kernp);
#endif

  // Regularise particle positions (if selected in parameters file)
  if (simparams->intparams["regularise_particle_ics"] == 1) {
    RegulariseParticleDistribution(simparams->intparams["Nreg"]);
  }

  // Search ghost particles
  sphneib->SearchBoundaryGhostParticles((FLOAT) 0.0, simbox, sph);
  sphneib->BuildGhostTree(true, 0, ntreebuildstep, ntreestockstep, sph->Ntot,
                          sph->Nhydromax, timestep, partdata, sph);
#ifdef MPI_PARALLEL
  mpicontrol->UpdateAllBoundingBoxes(sph->Nhydro + sph->NPeriodicGhost, sph, sph->kernp);
  MpiGhosts->SearchGhostParticles((FLOAT) 0.0, simbox, sph);
  sphneib->BuildMpiGhostTree(true, 0, ntreebuildstep, ntreestockstep, sph->Ntot,
                             sph->Nhydromax, timestep, partdata, sph);
#endif

  // Update neighbour tree
  rebuild_tree = true;
  sphneib->BuildTree(rebuild_tree, 0, ntreebuildstep, ntreestockstep,
                     sph->Ntot, sph->Nhydromax, timestep, partdata, sph);
  sphneib->SearchBoundaryGhostParticles((FLOAT) 0.0, simbox, sph);
  sphneib->BuildGhostTree(true, 0, ntreebuildstep, ntreestockstep, sph->Ntot,
                          sph->Nhydromax, timestep, partdata, sph);
  //sphneib->neibcheck = true;

    // Communicate pruned trees for MPI
#ifdef MPI_PARALLEL
  sphneib->BuildPrunedTree(rank, sph->Nhydromax, simbox, mpicontrol->mpinode, partdata);
#endif

  // Compute all initial N-body terms
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<nbody->Nstar; i++) {
    for (k=0; k<ndim; k++) nbody->stardata[i].a[k]     = (FLOAT) 0.0;
    for (k=0; k<ndim; k++) nbody->stardata[i].adot[k]  = (FLOAT) 0.0;
    for (k=0; k<ndim; k++) nbody->stardata[i].a2dot[k] = (FLOAT) 0.0;
    for (k=0; k<ndim; k++) nbody->stardata[i].a3dot[k] = (FLOAT) 0.0;
    nbody->stardata[i].gpot   = (FLOAT) 0.0;
    nbody->stardata[i].gpe    = (FLOAT) 0.0;
    nbody->stardata[i].tlast  = t;
    nbody->stardata[i].active = true;
    nbody->stardata[i].level  = 0;
    nbody->stardata[i].nstep  = 0;
    nbody->stardata[i].nlast  = 0;
    nbody->nbodydata[i]       = &(nbody->stardata[i]);
  }
  nbody->Nnbody = nbody->Nstar;


  // Read-in stellar properties table here
  nbody->LoadStellarPropertiesTable(&simunits);
  nbody->UpdateStellarProperties();

  // Compute all initial SPH force terms
  //-----------------------------------------------------------------------------------------------

  // Zero accelerations (here for now)
  for (i=0; i<sph->Ntot; i++) {
    SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
    part.tlast     = t;
    part.flags.unset_flag(active);
    part.level     = 0;
    part.levelneib = 0;
    part.nstep     = 0;
    part.nlast     = 0;
    part.dalphadt  = (FLOAT) 0.0;
    part.div_v     = (FLOAT) 0.0;
    part.dudt      = (FLOAT) 0.0;
    part.gpot      = (FLOAT) 0.0;
    part.mu_bar    = (FLOAT) simparams->floatparams["mu_bar"];
    for (k=0; k<ndim; k++) part.a[k] = (FLOAT) 0.0;
    for (k=0; k<ndim; k++) part.agrav[k] = (FLOAT) 0.0;
  }
  for (i=0; i<sph->Nhydro; i++) sph->GetSphParticlePointer(i).flags.set_flag(active);

  // Copy all other data from real SPH particles to ghosts
  LocalGhosts->CopyHydroDataToGhosts(simbox, sph);

  sphneib->BuildTree(true, 0, ntreebuildstep, ntreestockstep, sph->Ntot,
                     sph->Nhydromax, timestep, partdata, sph);
  sphneib->SearchBoundaryGhostParticles((FLOAT) 0.0, simbox, sph);
  sphneib->BuildGhostTree(true, 0, ntreebuildstep, ntreestockstep, sph->Ntot,
                          sph->Nhydromax, timestep, partdata, sph);
#ifdef MPI_PARALLEL
  mpicontrol->UpdateAllBoundingBoxes(sph->Nhydro + sph->NPeriodicGhost, sph, sph->kernp);
  MpiGhosts->SearchGhostParticles((FLOAT) 0.0, simbox, sph);
  sphneib->BuildMpiGhostTree(true, 0, ntreebuildstep, ntreestockstep, sph->Ntot,
                             sph->Nhydromax, timestep, partdata, sph);
#endif

  // Calculate all SPH properties
  sphneib->UpdateAllSphProperties(sph->Nhydro, sph->Ntot, partdata, sph, nbody);


#ifdef MPI_PARALLEL
  if (sph->self_gravity == 1) {
    sphneib->UpdateGravityExportList(rank, sph->Nhydro, sph->Ntot, partdata, sph, nbody, simbox);
  }
  else {
    sphneib->UpdateHydroExportList(rank, sph->Nhydro, sph->Ntot, partdata, sph, nbody, simbox);
  }

  mpicontrol->ExportParticlesBeforeForceLoop(sph);
#endif


  for (i=0; i<sph->Nhydro; i++) {
    SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
    part.ionfrac = (FLOAT) 0.9999999;
  }
  // Update the radiation field
  for (int jj=0; jj<10; jj++) {
    radiation->UpdateRadiationField(sph->Nhydro, nbody->Nnbody, sinks->Nsink,
                                    partdata, nbody->nbodydata, sinks->sink);
  }


  // Update thermal properties (if radiation field has altered them)
  /*for (i=0; i<sph->Nhydro; i++) {
    SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
    sph->ComputeThermalProperties(part);
  }*/


  // Calculate SPH gravity and hydro forces, depending on which are activated
 /* if (sph->hydro_forces == 1 && sph->self_gravity == 1) {
    sphneib->UpdateAllSphForces(sph->Nhydro, sph->Ntot, partdata, sph, nbody, simbox, ewald);
  }
  else if (sph->self_gravity == 1) {
    sphneib->UpdateAllSphGravForces(sph->Nhydro, sph->Ntot, partdata, sph, nbody, simbox, ewald);
  }
  else if (sph->hydro_forces == 1) {
    sphneib->UpdateAllSphHydroForces(sph->Nhydro, sph->Ntot, partdata, sph, nbody, simbox);
  }
  else{
    ExceptionHandler::getIstance().raise("Error: No forces included in simulation");
  }*/
  sphneib->UpdateAllSphForces(sph->Nhydro, sph->Ntot, partdata, sph, nbody, simbox, ewald);

#if defined MPI_PARALLEL
  mpicontrol->GetExportedParticlesAccelerations(sph);
#endif

  // Add external potential for all active SPH particles
  for (i=0; i<sph->Nhydro; i++) {
    SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
    sph->extpot->AddExternalPotential(part.r, part.v, part.a, adot, part.gpot);
  }

/*  // Compute the dust forces if present.
  if (sphdust != NULL){
    // Copy properties from original particles to ghost particles
    LocalGhosts->CopyHydroDataToGhosts(simbox, sph);
#ifdef MPI_PARALLEL
    MpiGhosts->CopyHydroDataToGhosts(simbox, sph);
#endif
    sphdust->UpdateAllDragForces(sph->Nhydro, sph->Ntot, partdata) ;
  }*/

  // Set initial accelerations
  for (i=0; i<sph->Nhydro; i++) {
      SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
      for (k=0; k<ndim; k++) part.r0[k] = part.r[k];
      for (k=0; k<ndim; k++) part.v0[k] = part.v[k];
      for (k=0; k<ndim; k++) part.a0[k] = part.a[k];
      part.flags.unset_flag(active);
    }

  LocalGhosts->CopyHydroDataToGhosts(simbox,sph);
#ifdef MPI_PARALLEL
  MpiGhosts->CopyHydroDataToGhosts(simbox,sph);
#endif


  // Compute initial N-body forces
  //-----------------------------------------------------------------------------------------------
  if (sph->self_gravity == 1 && sph->Nhydro > 0) {
    sphneib->UpdateAllStarGasForces(sph->Nhydro, sph->Ntot, partdata, sph, nbody);
#if defined MPI_PARALLEL
    // We need to sum up the contributions from the different domains
    mpicontrol->ComputeTotalStarGasForces(nbody);
#endif
  }

  if (nbody->nbody_softening == 1) {
    nbody->CalculateDirectSmoothedGravForces(nbody->Nnbody, nbody->nbodydata);
  }
  else {
    nbody->CalculateDirectGravForces(nbody->Nnbody, nbody->nbodydata);
  }
  nbody->CalculateAllStartupQuantities(nbody->Nnbody, nbody->nbodydata);

  for (i=0; i<nbody->Nnbody; i++) {
    if (nbody->nbodydata[i]->active) {
      nbody->extpot->AddExternalPotential(nbody->nbodydata[i]->r, nbody->nbodydata[i]->v,
                                          nbody->nbodydata[i]->a, nbody->nbodydata[i]->adot,
                                          nbody->nbodydata[i]->gpot);
    }
  }

  // Set particle values for initial step (e.g. r0, v0, a0, u0)
  uint->EndTimestep(n, sph->Nhydro, t, timestep, partdata);
  sphint->EndTimestep(n, sph->Nhydro, t, timestep, partdata);
  nbody->EndTimestep(n, nbody->Nstar, t, timestep, nbody->nbodydata);

  this->CalculateDiagnostics();
  this->diag0 = this->diag;
  this->setup = true;

  return;
}

//=================================================================================================
//  ForcelessSimulation::MainLoop
/// Main  forcless simulation integration loop.
//=================================================================================================
template <int ndim>
void ForcelessSimulation<ndim>::MainLoop(void)
{
  //int activecount = 0;                 // Flag if we need to recompute particles
  int i;                               // Particle loop counter
  int it;                              // Time-symmetric iteration counter
  int k;                               // Dimension counter
  FLOAT adot[ndim];                    // Dummy adot variable
  FLOAT tghost;                        // Approx. ghost particle lifetime
  SphParticle<ndim> *partdata = sph->GetSphParticleArray();

  debug2("[ForcelessSimulation::MainLoop]");

  // Compute timesteps for all particles
  if (Nlevels == 1) this->ComputeGlobalTimestep();
  else this->ComputeBlockTimesteps();

  // Advance time variables
  n = n + 1;
  Nsteps = Nsteps + 1;
  t = t + timestep;
  if (n == nresync) Nblocksteps = Nblocksteps + 1;
  if (n%integration_step == 0) Nfullsteps = Nfullsteps + 1;

  // Advance SPH and N-body particles' positions and velocities
  uint->EnergyIntegration(n, sph->Nhydro, (FLOAT) t, (FLOAT) timestep, partdata);
  sphint->AdvanceParticles(n, sph->Nhydro, (FLOAT) t, (FLOAT) timestep, partdata);
  nbody->AdvanceParticles(n, nbody->Nnbody, t, timestep, nbody->nbodydata);

  // Check all boundary conditions
  // (DAVID : Move this function to sphint and create an analagous one
  //  for N-body.  Also, only check this on tree-build steps)
  sphint->CheckBoundaries(simbox,sph);


  // Perform the load-balancing step for MPI simulations.  First update the pruned trees on all
  // processors, then compute the new load-balanced MPI domains and finally transfer the
  // particles to the new domains.
  //-----------------------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
  if (Nsteps%ntreebuildstep == 0 || rebuild_tree) {
    sphneib->BuildPrunedTree(rank, sph->Nhydromax, simbox, mpicontrol->mpinode, partdata);
    mpicontrol->UpdateAllBoundingBoxes(sph->Nhydro, sph, sph->kernp);
    mpicontrol->LoadBalancing(sph, nbody);
    sphneib->InitialiseCellWorkCounters();
  }
#endif


  // Rebuild or update local neighbour and gravity tree
  sphneib->BuildTree(rebuild_tree, Nsteps, ntreebuildstep, ntreestockstep,
                     sph->Ntot, sph->Nhydromax, timestep, partdata, sph);

  // Search for new ghost particles and create on local processor
  //-----------------------------------------------------------------------------------------------
  if (Nsteps%ntreebuildstep == 0 || rebuild_tree) {
    tghost = timestep*(FLOAT)(ntreebuildstep - 1);
    sphneib->SearchBoundaryGhostParticles(tghost, simbox, sph);
    sphneib->BuildGhostTree(rebuild_tree, Nsteps, ntreebuildstep, ntreestockstep,
                            sph->Ntot, sph->Nhydromax, timestep, partdata, sph);

  // Re-build and communicate the new pruned trees (since the trees will necessarily change
  // once there has been communication of particles to new domains)
#ifdef MPI_PARALLEL
    sphneib->BuildPrunedTree(rank, sph->Nhydromax, simbox, mpicontrol->mpinode, partdata);
    mpicontrol->UpdateAllBoundingBoxes(sph->Nhydro + sph->NPeriodicGhost, sph, sph->kernp);
    MpiGhosts->SearchGhostParticles(tghost, simbox, sph);
    sphneib->BuildMpiGhostTree(rebuild_tree, Nsteps, ntreebuildstep, ntreestockstep,
                               sph->Ntot, sph->Nhydromax, timestep, partdata, sph);
#endif
  }
  // Otherwise copy properties from original particles to ghost particles
  else {
    LocalGhosts->CopyHydroDataToGhosts(simbox, sph);
#ifdef MPI_PARALLEL
    MpiGhosts->CopyHydroDataToGhosts(simbox, sph);
#endif
    sphneib->BuildGhostTree(rebuild_tree, Nsteps, ntreebuildstep, ntreestockstep,
                            sph->Ntot, sph->Nhydromax, timestep, partdata, sph);
  }
  // Iterate if we need to immediately change SPH particle timesteps
  // (e.g. due to feedback, or sudden change in neighbour timesteps)
  //-----------------------------------------------------------------------------------------------

    // Update cells containing active particles
    sphneib->UpdateActiveParticleCounters(partdata, sph);

      // Zero accelerations (here for now)
      for (i=0; i<sph->Nhydro; i++) {
        SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
        if (part.flags.check_flag(active)) {
          part.levelneib = 0;
          part.dalphadt  = (FLOAT) 0.0;
          part.div_v     = (FLOAT) 0.0;
          part.dudt      = (FLOAT) 0.0;
          part.gpot      = (FLOAT) 0.0;
          for (k=0; k<ndim; k++) part.a[k] = (FLOAT) 0.0;
          for (k=0; k<ndim; k++) part.agrav[k] = (FLOAT) 0.0;
          for (k=0; k<ndim; k++) part.a_dust[k] = (FLOAT) 0.0;
        }
      }

      // Calculate all SPH properties
      if (t >= tsnapnext)
        sphneib->UpdateAllSphProperties(sph->Nhydro, sph->Ntot, partdata, sph, nbody);
      // Calculate SPH gravity and hydro forces, depending on which are activated
      // Call ForcelessSimulation version, not GradhSph one??
      sphneib->UpdateAllSphForces(sph->Nhydro, sph->Ntot, partdata, sph, nbody, simbox, ewald);

      // Add external potential for all active SPH particles
      for (i=0; i<sph->Nhydro; i++) {
        SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
        if (part.flags.check_flag(active)) {
          sph->extpot->AddExternalPotential(part.r, part.v, part.a, adot, part.gpot);
        }
      }

      // Checking if acceleration or other values are invalid
      for (i=0; i<sph->Nhydro; i++) {
        SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
        if (part.flags.check_flag(active)) {
          for (k=0; k<ndim; k++) assert(part.r[k] == part.r[k]);
          for (k=0; k<ndim; k++) assert(part.v[k] == part.v[k]);
          for (k=0; k<ndim; k++) assert(part.a[k] == part.a[k]);
          assert(part.gpot == part.gpot);
        }
      }

    // Zero all active flags once accelerations have been computed
    for (i=0; i<sph->Nhydro; i++) {
      SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
      part.flags.unset_flag(active);
    }

  //-----------------------------------------------------------------------------------------------


  // Iterate for P(EC)^n schemes for N-body particles
  //-----------------------------------------------------------------------------------------------
  for (it=0; it<nbody->Npec; it++) {

    // Zero all acceleration terms
    for (i=0; i<nbody->Nnbody; i++) {
      if (nbody->nbodydata[i]->active) {
        for (k=0; k<ndim; k++) nbody->nbodydata[i]->a[k]     = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) nbody->nbodydata[i]->adot[k]  = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) nbody->nbodydata[i]->a2dot[k] = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) nbody->nbodydata[i]->a3dot[k] = (FLOAT) 0.0;
        nbody->nbodydata[i]->gpot = (FLOAT) 0.0;
        nbody->nbodydata[i]->gpe = (FLOAT) 0.0;
      }
    }
    if (sink_particles == 1) {
      for (i=0; i<sinks->Nsink; i++) {
        if (sinks->sink[i].star->active) {
          for (k=0; k<ndim; k++) sinks->sink[i].fhydro[k] = (FLOAT) 0.0;
        }
      }
    }

    // Calculate forces, force derivatives etc.., for active stars/systems
    if (nbody->nbody_softening == 1) {
      nbody->CalculateDirectSmoothedGravForces(nbody->Nnbody,nbody->nbodydata);
    }
    else {
      nbody->CalculateDirectGravForces(nbody->Nnbody,nbody->nbodydata);
    }

    for (i=0; i<nbody->Nnbody; i++) {
      if (nbody->nbodydata[i]->active) {
        nbody->extpot->AddExternalPotential(nbody->nbodydata[i]->r, nbody->nbodydata[i]->v,
                                            nbody->nbodydata[i]->a, nbody->nbodydata[i]->adot,
                                            nbody->nbodydata[i]->gpot);
      }
    }

    // Calculate correction step for all stars at end of step, except the
    // final iteration (since correction is computed in EndStep also).
    //if (it < nbody->Npec - 1)
    nbody->CorrectionTerms(n,nbody->Nnbody,t,timestep,nbody->nbodydata);

  }
  //-----------------------------------------------------------------------------------------------


  rebuild_tree = false;

  // End-step terms for all SPH particles
  if (sph->Nhydro > 0) {
    uint->EndTimestep(n, sph->Nhydro, (FLOAT) t, (FLOAT) timestep, partdata);
    sphint->EndTimestep(n, sph->Nhydro, (FLOAT) t, (FLOAT) timestep, partdata);
  }

  // End-step terms for all star particles
  if (nbody->Nstar > 0) nbody->EndTimestep(n, nbody->Nnbody, t, timestep, nbody->nbodydata);

  // Search for new sink particles (if activated) and accrete to existing sinks
  if (sink_particles == 1) {
    if (sinks->create_sinks == 1 && (rebuild_tree || Nfullsteps%ntreebuildstep == 0)) {
      sinks->SearchForNewSinkParticles(n, t ,sph, nbody);
    }
    if (sinks->Nsink > 0) {
      sph->mmean = (FLOAT) 0.0;
      for (i=0; i<sph->Nhydro; i++) sph->mmean += sph->GetSphParticlePointer(i).m;
      sph->mmean /= (FLOAT) sph->Nhydro;
      sph->hmin_sink = big_number;
      for (i=0; i<sinks->Nsink; i++) {
        sph->hmin_sink = min(sph->hmin_sink, (FLOAT) sinks->sink[i].star->h);
      }

      sinks->AccreteMassToSinks(n, timestep, partdata, sph, nbody);
      nbody->UpdateStellarProperties();
      if (extra_sink_output) WriteExtraSinkOutput();
    }
    // If we will output a snapshot (regular or for restarts), then delete all accreted particles
    if ((t >= tsnapnext && sinks->Nsink > 0) || n == nresync || kill_simulation ||
         timing->WallClockTime() - timing->tstart_wall > (FLOAT) 0.99*tmax_wallclock) {
      sph->DeleteDeadParticles();
      rebuild_tree = true;
    }
  }

  return;
}



//=================================================================================================
//  ForcelessSimulation::ComputeGlobalTimestep
/// Computes global timestep for SPH simulation.  Calculates the minimum
/// timestep for all SPH and N-body particles in the simulation.
//=================================================================================================
template <int ndim>
void ForcelessSimulation<ndim>::ComputeGlobalTimestep(void)
{
  int i;                               // Particle counter
  DOUBLE dt;                           // Particle timestep
  DOUBLE dt_min = big_number_dp;       // Local copy of minimum timestep
  DOUBLE dt_nbody;                     // Aux. minimum N-body timestep
  DOUBLE dt_sph;                       // Aux. minimum SPH timestep

  debug2("[ForcelessSimulation::ComputeGlobalTimestep]");
  timing->StartTimingSection("GLOBAL_TIMESTEPS");


  // Only update timestep when all particles are synced at end of last step.
  //-----------------------------------------------------------------------------------------------
  if (n == nresync) {

    n            = 0;
    level_max    = 0;
    level_step   = level_max + integration_step - 1;
    nresync      = integration_step;
    dt_min_nbody = big_number_dp;
    dt_min_hydro   = big_number_dp;

    // Find minimum timestep from all SPH particles
    //---------------------------------------------------------------------------------------------
#pragma omp parallel default(none) private(i,dt,dt_nbody,dt_sph) shared(dt_min)
    {
      dt       = big_number_dp;
      dt_nbody = big_number_dp;
      dt_sph   = big_number_dp;

#pragma omp for
      for (i=0; i<sph->Nhydro; i++) {
        SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
        part.level     = 0;
        part.levelneib = 0;
        part.nstep     = pow(2,level_step - part.level);
        part.nlast     = n;
        part.tlast     = t;
        part.dt        = Timestep(tdyn_mult, part);
        dt             = min(dt,part.dt);
        dt_sph         = min(dt_sph,part.dt);
      }

      // Now compute minimum timestep due to stars/systems
#pragma omp for
      for (i=0; i<nbody->Nnbody; i++) {
        nbody->nbodydata[i]->level = 0;
        nbody->nbodydata[i]->nstep = pow(2,level_step - nbody->nbodydata[i]->level);
        nbody->nbodydata[i]->nlast = n;
        nbody->nbodydata[i]->tlast = t;
        nbody->nbodydata[i]->dt    = nbody->Timestep(nbody->nbodydata[i]);
        dt       = min(dt,nbody->nbodydata[i]->dt);
        dt_nbody = min(dt_nbody,nbody->nbodydata[i]->dt);
      }

#pragma omp critical
      {
        if (dt < dt_min) dt_min = dt;
        if (dt_sph < dt_min_hydro) dt_min_hydro = dt_sph;
        if (dt_nbody < dt_min_nbody) dt_min_nbody = dt_nbody;
      }

    }
    //---------------------------------------------------------------------------------------------

#ifdef MPI_PARALLEL
    dt = dt_min;
    MPI_Allreduce(&dt, &dt_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif
    timestep = dt_min;

    // Set minimum timestep for all SPH and N-body particles
    for (i=0; i<sph->Nhydro; i++) sph->GetSphParticlePointer(i).dt = timestep;
    for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->dt = timestep;

  }
  //-----------------------------------------------------------------------------------------------

  timing->EndTimingSection("GLOBAL_TIMESTEPS");

  return;
}



//=================================================================================================
//  ForcelessSimulation::ComputeBlockTimesteps
/// Compute timesteps for all particles using hierarchical block timesteps.
//=================================================================================================
template <int ndim>
void ForcelessSimulation<ndim>::ComputeBlockTimesteps(void)
{
  int i;                                     // Particle counter
  int istep;                                 // Aux. variable for changing steps
  int last_level;                            // Previous timestep level
  int level;                                 // Particle timestep level
  int level_max_aux;                         // Aux. maximum level variable
  int level_max_nbody = 0;                   // level_max for star particles only
  int level_max_old;                         // Old level_max
  int level_max_sph = 0;                     // level_max for SPH particles only
  int level_min_sph = 9999999;               // level_min for SPH particles
  int level_nbody;                           // local thread var. for N-body level
  int level_sph;                             // local thread var. for SPH level
  int nfactor;                               // Increase/decrease factor of n
  int nstep;                                 // Particle integer step-size
  DOUBLE dt;                                 // Aux. timestep variable
  DOUBLE dt_min = big_number_dp;             // Minimum timestep
  DOUBLE dt_min_aux;                         // Aux. minimum timestep variable
  DOUBLE dt_nbody;                           // Aux. minimum N-body timestep
  DOUBLE dt_sph;                             // Aux. minimum SPH timestep

  debug2("[ForcelessSimulation::ComputeBlockTimesteps]");
  timing->StartTimingSection("BLOCK_TIMESTEPS");


  dt_min_nbody = big_number_dp;
  dt_min_hydro = big_number_dp;


  // Synchronise all timesteps and reconstruct block timestep structure.
  //===============================================================================================
  if (n == nresync) {

    n = 0;
    timestep = big_number_dp;

#pragma omp parallel default(none) private(dt,dt_min_aux,dt_nbody,dt_sph,i)
    {
      // Initialise all timestep and min/max variables
      dt_min_aux = big_number_dp;
      dt_sph     = big_number_dp;
      dt_nbody   = big_number_dp;

      // Find minimum timestep from all SPH particles
#pragma omp for
      for (i=0; i<sph->Nhydro; i++) {
        SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
        if (part.flags.is_dead()) continue;
        dt         = Timestep(tdyn_mult, part);
        dt_min_aux = min(dt_min_aux,dt);
        dt_sph     = min(dt_sph,dt);
        part.dt    = dt;
      }

      // Now compute minimum timestep due to stars/systems
#pragma omp for
      for (i=0; i<nbody->Nnbody; i++) {
        dt         = nbody->Timestep(nbody->nbodydata[i]);
        dt_min_aux = min(dt_min_aux,dt);
        dt_nbody   = min(dt_nbody,dt);
        nbody->nbodydata[i]->dt = dt;
      }

#pragma omp critical
      {
        timestep     = min(timestep,dt_min_aux);
        dt_min_hydro = min(dt_min_hydro,dt_sph);
        dt_min_nbody = min(dt_min_nbody,dt_nbody);
      }
#pragma omp barrier
    }


    // For MPI, determine the global minimum timestep over all processors
#ifdef MPI_PARALLEL
    dt = timestep;
    MPI_Allreduce(&dt, &timestep, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    dt = dt_min_hydro;
    MPI_Allreduce(&dt, &dt_min_hydro, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    dt = dt_min_nbody;
    MPI_Allreduce(&dt, &dt_min_nbody, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif


    // Calculate new block timestep levels
    level_max  = Nlevels - 1;
    level_step = level_max + integration_step - 1;
    dt_max     = timestep*powf(2.0, level_max);

    // Calculate the maximum level occupied by all SPH particles
    level_max_sph   = min(ComputeTimestepLevel(dt_min_hydro, dt_max), level_max);
    level_max_nbody = min(ComputeTimestepLevel(dt_min_nbody, dt_max), level_max);

    // Populate timestep levels with N-body particles.
    // Ensures that N-body particles occupy levels lower than all SPH particles
    for (i=0; i<nbody->Nnbody; i++) {
      dt = nbody->nbodydata[i]->dt;
      level = min(ComputeTimestepLevel(dt, dt_max), level_max);
      nbody->nbodydata[i]->level = max(level, level_max_sph);
      nbody->nbodydata[i]->nlast = n;
      nbody->nbodydata[i]->nstep = pow(2, level_step - nbody->nbodydata[i]->level);
      nbody->nbodydata[i]->tlast = t;
    }

    // If particles are sink neighbours, set to same timesteps as sinks
    for (i=0; i<sph->Nhydro; i++) {
      SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
      if (part.sinkid != -1) {
        if (sinks->sink[part.sinkid].star->level - part.level > level_diff_max) {
          part.level     = sinks->sink[part.sinkid].star->level - level_diff_max;
          part.levelneib = sinks->sink[part.sinkid].star->level;
          level_max_sph  = max(level_max_sph, part.level);
        }
      }
    }

    // If enforcing a single SPH timestep, set it here.
    // Otherwise, populate the timestep levels with SPH particles.
    if (sph_single_timestep == 1) {
      for (i=0; i<sph->Nhydro; i++) {
        SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
        if (part.flags.is_dead()) continue;
        part.level     = level_max_sph;
        part.levelneib = level_max_sph;
        part.nlast     = n;
        part.tlast     = t;
        part.nstep     = pow(2, level_step - part.level);
      }
      level_min_sph = level_max_sph;
    }
    else {
      for (i=0; i<sph->Nhydro; i++) {
        SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
        if (part.flags.is_dead()) continue;
        dt             = part.dt;
        level          = min(ComputeTimestepLevel(dt, dt_max), level_max);
        part.level     = level;
        part.levelneib = level;
        part.nlast     = n;
        part.tlast     = t;
        part.nstep     = pow(2, level_step - part.level);
        level_min_sph  = min(level_min_sph, part.level);
      }
    }


    nresync = pow(2, level_step);
    assert(nresync > 0);
    timestep = dt_max / (DOUBLE) nresync;

  }
  // If not resynchronising, check if any SPH/N-body particles need to move
  // up or down timestep levels.
  //===============================================================================================
  else {

    level_max_old   = level_max;
    level_max       = 0;
    level_max_nbody = 0;
    level_max_sph   = 0;


#pragma omp parallel default(shared) private(dt,dt_nbody,dt_sph,i)\
     private(istep,last_level,level,level_max_aux,level_nbody,level_sph,nstep,nfactor)
    //shared(dt_min,imin,level_max_nbody,level_max_sph,level_min_sph)
    //shared(cout)
    {
      dt_sph        = big_number_dp;
      dt_nbody      = big_number_dp;
      level_max_aux = 0;
      level_nbody   = 0;
      level_sph     = 0;


      // Find all SPH particles at the beginning of a new timestep
      //-------------------------------------------------------------------------------------------
#pragma omp for
      for (i=0; i<sph->Nhydro; i++) {
        SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
        if (part.flags.is_dead()) continue;

        // SPH particles whose timestep has been artificially reduced by Saitoh & Makino scheme.
        if (part.nlast == n && part.nstep != pow(2,level_step - part.level)) {
          dt             = Timestep(tdyn_mult, part);
          level          = max(ComputeTimestepLevel(dt, dt_max), part.levelneib - level_diff_max);
          part.level     = max(part.level, level);
          part.levelneib = part.level;
          part.dt        = dt;
          part.nlast     = n;
          part.tlast     = t;
          part.nstep     = pow(2, level_step - part.level);
        }
        // SPH particles that have naturally reached the end of their step
        else if (part.nlast == n) {
          nstep      = part.nstep;
          last_level = part.level;

          // Compute new timestep value and level number
          dt    = Timestep(tdyn_mult, part);
          level = max(ComputeTimestepLevel(dt, dt_max), part.levelneib - level_diff_max);

          // Move up one level (if levels are correctly synchronised) or
          // down several levels if required
          if (level < last_level && last_level > 1 && n%(2*nstep) == 0) {
            part.level = last_level - 1;
          }
          else if (level > last_level) {
            part.level = level;
          }
          else {
            part.level = last_level;
          }

          part.levelneib = level;
          part.dt        = dt;
          part.nlast     = n;
          part.tlast     = t;
          part.nstep     = pow(2, level_step - part.level);
        }

        // Find maximum level of all SPH particles
        level_sph     = max(level_sph, part.level);
        level_max_aux = max(level_max_aux, part.level);

        dt_sph = min(dt_sph, part.dt);
      }
      //-------------------------------------------------------------------------------------------


#pragma omp critical
      {
        dt_min        = min(dt_min, dt_sph);
        dt_min_hydro  = min(dt_min_hydro, dt_sph);
        level_max     = max(level_max, level_max_aux);
        level_max_sph = max(level_max_sph, level_sph);
      }
#pragma omp barrier

#if defined MPI_PARALLEL
#pragma omp master
      {
        level = level_max_sph;
        MPI_Allreduce(&level, &level_max_sph, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      }
#pragma omp barrier
#endif

      // Now find all N-body particles at the beginning of a new timestep
      //-------------------------------------------------------------------------------------------
#pragma omp for
      for (i=0; i<nbody->Nnbody; i++) {

        // Skip particles that are not at end of step
        if (nbody->nbodydata[i]->nlast == n) {
          nstep = nbody->nbodydata[i]->nstep;
          last_level = nbody->nbodydata[i]->level;

          // Compute new timestep value and level number
          dt    = nbody->Timestep(nbody->nbodydata[i]);
          level = max(ComputeTimestepLevel(dt, dt_max), level_max_sph);

          // Move up one level (if levels are correctly synchronised) or
          // down several levels if required
          if (level < last_level && level > level_max_sph && last_level > 1 && n%(2*nstep) == 0) {
            nbody->nbodydata[i]->level = last_level - 1;
          }
          else if (level > last_level) {
            nbody->nbodydata[i]->level = level;
          }
          else {
            nbody->nbodydata[i]->level = last_level;
          }

          nbody->nbodydata[i]->dt    = dt;
          nbody->nbodydata[i]->nlast = n;
          nbody->nbodydata[i]->nstep = pow(2, level_step - nbody->nbodydata[i]->level);
          nbody->nbodydata[i]->tlast = t;
        }

        // Find maximum level of all N-body particles
        level_nbody   = max(level_nbody, nbody->nbodydata[i]->level);
        level_max_aux = max(level_max_aux, nbody->nbodydata[i]->level);
        dt_nbody      = min(dt_nbody, nbody->nbodydata[i]->dt);
      }
      //-------------------------------------------------------------------------------------------

      #pragma omp critical
      {
        dt_min          = min(dt_min, dt_nbody);
        dt_min_nbody    = min(dt_min_nbody, dt_nbody);
        level_max       = max(level_max, level_max_aux);
        level_max_nbody = max(level_max_nbody, level_nbody);
      }
      #pragma omp barrier


      // Correct timestep levels for any particles that have entered a sink
      //-------------------------------------------------------------------------------------------
  #pragma omp for
      for (i=0; i<sph->Nhydro; i++) {
        SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
        if (part.flags.is_dead() || part.nlast != n) continue;
        if (part.sinkid != -1) {
          if (sinks->sink[part.sinkid].star->level - part.level > level_diff_max) {
            part.level = sinks->sink[part.sinkid].star->level - level_diff_max;
          }
        }
      }

    }

    // For MPI, find the global maximum timestep levels for each processor
#ifdef MPI_PARALLEL
    level = level_max;
    MPI_Allreduce(&level, &level_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    level = level_max_nbody;
    MPI_Allreduce(&level, &level_max_nbody, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    assert(level_max_sph >= 0);
#endif

    assert(!(isnan(dt_min)) && !(isinf(dt_min)));
    assert(!(isnan(dt_max)) && !(isinf(dt_max)));

    // Set fixed SPH timestep level here in case maximum has changed
    if (sph_single_timestep == 1) {
      for (i=0; i<sph->Nhydro; i++) {
        SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
        if (part.flags.is_dead()) continue;
        if (part.nlast == n) part.level = level_max_sph;
      }
    }

    // Update all timestep variables if we have removed or added any levels
    /*if (level_max != level_max_old) {

      // Increase maximum timestep level if correctly synchronised
      istep = pow(2, level_step - level_max_old + 1);
      if (level_max <= level_max_old - 1 && level_max_old > 1 && n%istep == 0) {
        level_max = level_max_old - 1;
      }
      else if (level_max < level_max_old) {
        level_max = level_max_old;
      }
    }*/

    istep = pow(2, level_step - level_max_old + 1);

    // Adjust integer time if levels are added or removed
    if (level_max > level_max_old) {
      nfactor = pow(2, level_max - level_max_old);
      n *= nfactor;
      for (i=0; i<sph->Nhydro; i++) {
        SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
        if (part.flags.is_dead()) continue;
        part.nstep *= nfactor;
        part.nlast *= nfactor;
      }
      for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nstep *= nfactor;
      for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nlast *= nfactor;
    }
    else if (level_max <= level_max_old - 1 && level_max_old > 1 && n%istep == 0) {
      level_max = level_max_old - 1;

      nfactor = pow(2, level_max_old - level_max);
      assert(n%nfactor == 0);
      n /= nfactor;
      for (i=0; i<sph->Nhydro; i++) {
        SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
        if (part.flags.is_dead()) continue;
        part.nlast /= nfactor;
        part.nstep /= nfactor;
      }
      for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nlast /= nfactor;
      for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nstep /= nfactor;
    }
    else {
      level_max = level_max_old;
    }

    level_step = level_max + integration_step - 1;
    nresync    = pow(2, level_step);
    timestep   = dt_max / (DOUBLE) nresync;

    // Update values of nstep for both SPH and star particles
    for (i=0; i<sph->Nhydro; i++) {
      SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
      if (part.flags.is_dead()) continue;
      if (part.nlast == n) part.nstep = pow(2, level_step - part.level);
    }
    for (i=0; i<nbody->Nnbody; i++) {
      if (nbody->nbodydata[i]->nlast == n) {
        nbody->nbodydata[i]->nstep = pow(2, level_step - nbody->nbodydata[i]->level);
      }
    }

    assert(level_max >= level_max_old - 1);

  }
  //===============================================================================================


  // Various asserts for debugging
  assert(timestep >= 0.0 && !(isinf(timestep)) && !(isnan(timestep)));
  assert(dt_max > 0.0 && !(isinf(dt_max)) && !(isnan(dt_max)));
  assert(level_step == level_max + integration_step - 1);
  assert(level_max_sph <= level_max);
  assert(level_max_nbody <= level_max);
  assert(n <= nresync);
  for (i=0; i<sph->Nhydro; i++) {
    SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
    if (part.flags.is_dead()) continue;
    assert(part.level <= level_max);
    assert(part.nlast <= n);
    assert(part.tlast <= t);
    assert(part.nstep == pow(2,level_step - part.level));
    assert(part.nlast != n || n%part.nstep == 0);
  }
  for (i=0; i<nbody->Nnbody; i++) {
    assert(nbody->nbodydata[i]->level <= level_max);
    assert(nbody->nbodydata[i]->nlast <= n);
    assert(nbody->nbodydata[i]->nstep == pow(2,level_step - nbody->nbodydata[i]->level));
    assert(nbody->nbodydata[i]->nlast != n || n%nbody->nbodydata[i]->nstep == 0);
    assert(nbody->nbodydata[i]->level >= level_max_sph);
    assert(nbody->nbodydata[i]->tlast <= t);
  }
  if (timestep <= 0.0) {
    cout << "Timestep fallen to zero : " << timestep << "    dtmax: " << dt_max
         << "    nresync " << nresync << endl;
    ExceptionHandler::getIstance().raise("Error : timestep fallen to zero");
  }

  timing->EndTimingSection("BLOCK_TIMESTEPS");

  return;


  // Some validations
  //-----------------------------------------------------------------------------------------------
  /*int *ninlevel;
  int Nactive=0;
  ninlevel = new int[level_max+1];
  //SphParticle<ndim>& part = sph->GetSphParticlePointer(imin);
  cout << "-----------------------------------------------------" << endl;
  cout << "Checking timesteps : " << level_max << "   " << level_max_sph << "    "
       << level_max_nbody << "    " << level_step << "   " << level_max_old << endl;
  cout << "n : " << n << endl;
  cout << "dt_min_hydro : " << dt_min_hydro << "    dt_min_nbody : " << dt_min_nbody
       << "    timestep : " << timestep << endl;
  cout << "hmin_sink : " << sph->hmin_sink << endl;
  //cout << "imin : " << imin << "    " << part.dt << "     " << part.m/sph->mmean << "    "
  //     << part.h << "    " << "    " << part.sound << "     " << part.div_v << "     "
  //     << part.h/(part.sound + part.h*fabs(part.div_v)) << "     "
  //     << sqrt(part.h/sqrt(DotProduct(part.a,part.a,ndim)))
  //     << endl;
  for (int l=0; l<=level_max; l++) ninlevel[l] = 0;
  for (i=0; i<sph->Nhydro; i++) {
    SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
    if (part.active) Nactive++;
    ninlevel[part.level]++;
  }
  cout << "No. of active SPH particles : " << Nactive << endl;
  cout << "SPH level occupancy" << endl;
  for (int l=0; l<=level_max; l++) cout << "level : " << l << "     N : " << ninlevel[l] << endl;

  for (int l=0; l<=level_max; l++) ninlevel[l] = 0;
  for (i=0; i<nbody->Nstar; i++) ninlevel[nbody->nbodydata[i]->level]++;
  cout << "N-body level occupancy" << endl;
  for (int l=0; l<=level_max; l++) cout << "level : " << l << "     N : " << ninlevel[l] << endl;

  delete[] ninlevel;

  if (timestep <= 0.0) {
    cout << "Timestep fallen to zero : " << timestep << endl;
    ExceptionHandler::getIstance().raise("Error : Timestep fallen to zero");
  }
  cin >> i;*/


  return;

}



//=================================================================================================
//  ForcelessSimulation::WriteExtraSinkOutput
/// For any simulations loaded into memory via a snapshot file, all particle
/// variables are converted into dimensionless code units here.
//=================================================================================================
template <int ndim>
void ForcelessSimulation<ndim>::WriteExtraSinkOutput(void)
{
  int k;                               // ..
  int s;                               // ..
  string filename;                     // Output snapshot filename
  string nostring;                     // String of number of snapshots
  stringstream ss;                     // Stream object for preparing filename
  ofstream outfile;                    // Stream of restart file


  //-----------------------------------------------------------------------------------------------
  for (s=0; s<sinks->Nsink; s++) {

    SinkParticle<ndim> &sink = sinks->sink[s];
    nostring = "";
    ss << setfill('0') << setw(5) << s;
    nostring = ss.str();
    filename = run_id + ".sink." + nostring;
    ss.str(std::string());

    outfile.open(filename.c_str(), std::ofstream::app);
    outfile << t << "    ";
    outfile << Nsteps << "    ";
    for (k=0; k<ndim; k++) outfile << sink.star->r[k] << "    ";
    for (k=0; k<ndim; k++) outfile << sink.star->v[k] << "    ";
    for (k=0; k<ndim; k++) outfile << sink.star->a[k] << "    ";
    for (k=0; k<3; k++) outfile << sink.angmom[k] << "    ";
    outfile << sink.star->m  << "    ";
    outfile << sink.menc     << "    ";
    outfile << sink.mmax     << "    ";
    outfile << sink.macctot  << "    ";
    outfile << sink.dmdt     << "    ";
    outfile << sink.ketot    << "    ";
    outfile << sink.gpetot   << "    ";
    outfile << sink.rotketot << "    ";
    outfile << sink.utot     << "    ";
    outfile << sink.taccrete << "    ";
    outfile << sink.trad     << "    ";
    outfile << sink.trot     << "    ";
    outfile << sink.tvisc    << "    ";
    outfile << endl;
    outfile.close();

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  ForcelessSimulation::RegulariseParticleDistribution
/// ...
//=================================================================================================
template <int ndim>
void ForcelessSimulation<ndim>::RegulariseParticleDistribution
 (const int Nreg)                                  ///< [in] No. of regularisation steps
{
  const FLOAT alphaReg = 0.2;                      // Particle displacement magnitude
  FLOAT *rreg = new FLOAT[ndim*sph->Nhydromax];    // Arrya of particle positions
  SphParticle<ndim> *partdata = sph->GetSphParticleArray();


  //===============================================================================================
  for (int ireg=0; ireg<Nreg; ireg++) {

    // Buid/re-build tree, create ghosts and update particle properties
    rebuild_tree = true;
    sphneib->BuildTree(rebuild_tree, 0, ntreebuildstep, ntreestockstep,
                       sph->Ntot, sph->Nhydromax, timestep, partdata, sph);
    sphneib->SearchBoundaryGhostParticles((FLOAT) 0.0, simbox, sph);
    sphneib->BuildGhostTree(true, 0, ntreebuildstep, ntreestockstep, sph->Ntot,
                            sph->Nhydromax, timestep, partdata, sph);
    sphneib->UpdateAllSphProperties(sph->Nhydro, sph->Ntot, partdata, sph, nbody);


    //=============================================================================================
#pragma omp parallel default(none) shared(partdata, rreg)
    {
      int k;
      FLOAT dr[ndim];
      FLOAT drsqd;
      int *neiblist = new int[sph->Nhydromax];


      //-------------------------------------------------------------------------------------------
#pragma omp for
      for (int i=0; i<sph->Nhydro; i++) {
        SphParticle<ndim> &part = sph->GetSphParticlePointer(i);
        FLOAT invhsqd = part.invh*part.invh;
        for (k=0; k<ndim; k++) rreg[ndim*i + k] = (FLOAT) 0.0;

        // Find list of gather neighbours
        int Nneib = sphneib->GetGatherNeighbourList(part.r, sph->kernrange*part.h, partdata,
                                                    sph->Ntot, sph->Nhydromax, neiblist);

        // Loop over all neighbours and calculate position correction for regularisation
        //-----------------------------------------------------------------------------------------
        for (int jj=0; jj<Nneib; jj++) {
          int j = neiblist[jj];
          SphParticle<ndim> &neibpart = sph->GetSphParticlePointer(j);

          for (k=0; k<ndim; k++) dr[k] = neibpart.r[k] - part.r[k];
          drsqd = DotProduct(dr, dr, ndim);
          if (drsqd >= part.hrangesqd) continue;
          for (k=0; k<ndim; k++) rreg[ndim*i + k] += dr[k]*sph->kernp->w0_s2(drsqd*invhsqd);

        }
        //-----------------------------------------------------------------------------------------

      }
      //-------------------------------------------------------------------------------------------


      // Apply all regularisation corrections to particle positions
      //-------------------------------------------------------------------------------------------
#pragma omp for
      for (int i=0; i<sph->Nhydro; i++) {
        SphParticle<ndim> &part = sph->GetSphParticlePointer(i);
        for (k=0; k<ndim; k++) part.r[k] -= alphaReg*rreg[ndim*i + k];
      }

      delete[] neiblist;

    }
    //=============================================================================================

    // Check that new positions don't fall outside the domain box
    sphint->CheckBoundaries(simbox, sph);

  }
  //================================================================================================

  delete[] rreg;

  return;
}

//=================================================================================================
//  ForcelessSimulation::Timestep
/// Default timestep size for forceless SPH particles
//  Dynamical time scale from orbit of sink particles
//=================================================================================================
template <int ndim>
DOUBLE ForcelessSimulation<ndim>::Timestep
 (const FLOAT tdyn_mult_aux, SphParticle<ndim> &part)            ///< [inout] Reference to SPH particle
{
  int j;                               // Star counter
  int k;                               // Dimension counter
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drmag;                         // Distance
  FLOAT drsqd;                         // Distance squared
  FLOAT ms;                            // Star mass
  FLOAT stardistmin ;				   // Minimum distance to star
  FLOAT mstar;						   // Minimum distance star mass
  DOUBLE timestep;                     // Minimum value of particle timesteps

  mstar =0.0;
  ms = 0.0;
  drmag = 0.0;
  stardistmin = 1e100 ; //std::numeric_limits<double>::max() ;
  timestep =0.0;

  NbodyParticle<ndim>** star = nbody->nbodydata;
  GradhSphParticle<ndim>& parti = static_cast<GradhSphParticle<ndim>& > (part);


  // Loop over all stars and add each contribution
    //-----------------------------------------------------------------------------------------------
    for (j=0; j<nbody->Nnbody; j++) {

      if (sph->fixed_sink_mass){
    	  ms = sph->msink_fixed;
      }
      else {
    	  ms = star[j]->m;
      }

      for (k=0; k<ndim; k++) dr[k] = star[j]->r[k] - parti.r[k];
      drsqd    = DotProduct(dr,dr,ndim) + small_number;
      drmag    = sqrt(drsqd);
      if (drmag<stardistmin) {
    	  stardistmin= drmag;
    	  mstar = ms;
      }
    }

   timestep = tdyn_mult_aux*sqrt(stardistmin*stardistmin*stardistmin/mstar);

  if (timestep > 1.0e20) {
    cout << "Timestep problem : " << timestep  << endl;
    ExceptionHandler::getIstance().raise("Error : Hydro timestep too large");
  }

  assert(!(isnan(timestep)) && !(isinf(timestep)));
  return timestep;
}


//=================================================================================================
//  ForcelessSimulation::ProcessSphParameters
/// Process all the options chosen in the parameters file, setting various
/// simulation variables and creating important simulation objects.
//=================================================================================================
template <int ndim>
void ForcelessSimulation<ndim>::ProcessSphParameters(void)
{
  aviscenum avisc = noav;              // Artificial viscosity enum
  acondenum acond = noac;              // Artificial conductivity enum
  eosenum eos_type = noeos;            // Gas EOS enum
  tdaviscenum tdavisc = notdav;        // Time-dependent viscosity enum

  // Local references to parameter variables for brevity
  map<string, int> &intparams = simparams->intparams;
  map<string, double> &floatparams = simparams->floatparams;
  map<string, string> &stringparams = simparams->stringparams;
  string KernelName = stringparams["kernel"];
  string gas_eos = stringparams["gas_eos"];
  string gas_radiation = stringparams["radiation"];

  debug2("[ForcelessSimulation::ProcessSphParameters]");


  // Set the enum for artificial viscosity
  if (stringparams["avisc"] == "none") {
    avisc = noav;
    tdavisc = notdav;
  }
  else if (stringparams["avisc"] == "mon97" && stringparams["time_dependent_avisc"] == "mm97") {
    avisc = mon97mm97;
    tdavisc = mm97;
  }
  else if (stringparams["avisc"] == "mon97" && stringparams["time_dependent_avisc"] == "cd2010") {
    avisc = mon97cd2010;
    tdavisc = cd2010;
  }
  else if (stringparams["avisc"] == "mon97") {
    avisc = mon97;
    tdavisc = notdav;
  }
  else {
    string message = "Unrecognised parameter : avisc = " + simparams->stringparams["avisc"] +
      "   or time_dependent_avisc : " + simparams->stringparams["time_dependent_avisc"];
    ExceptionHandler::getIstance().raise(message);
  }

  // Set the enum for artificial conductivity
  if (stringparams["acond"] == "none") {
    acond = noac;
  }
  else if (stringparams["acond"] == "wadsley2008") {
    acond = wadsley2008;
  }
  else if (stringparams["acond"] == "price2008") {
    acond = price2008;
  }
  else {
    string message = "Unrecognised parameter : acond = " + simparams->stringparams["acond"];
    ExceptionHandler::getIstance().raise(message);
  }

  // Set gas EOS values
  if (stringparams["gas_eos"] == "isothermal") {
    eos_type = isothermal;
  }
  else if (stringparams["gas_eos"] == "barotropic") {
    eos_type = barotropic;
  }
  else if (stringparams["gas_eos"] == "barotropic2") {
    eos_type = barotropic2;
  }
  else if (stringparams["gas_eos"] == "energy_eqn") {
    eos_type = energy_eqn;
  }
  else if (stringparams["gas_eos"] == "constant_temp") {
    eos_type = constant_temp;
  }
  else if (stringparams["gas_eos"] == "rad_ws" || stringparams["gas_eos"] == "radws") {
    eos_type = radws;
  }
  else {
    string message = "Unrecognised eos parameter : gas_eos = " + simparams->stringparams["gas_eos"];
    ExceptionHandler::getIstance().raise(message);
  }


  // Create 'grad-h' SPH object depending on choice of kernel
  //===============================================================================================
  if (intparams["tabulated_kernel"] == 1) {
    sph = new Forceless<ndim, TabulatedKernel>
      (intparams["hydro_forces"], intparams["self_gravity"], floatparams["alpha_visc"],
       floatparams["beta_visc"], floatparams["h_fac"], floatparams["h_converge"],
       avisc, acond, tdavisc, stringparams["gas_eos"], KernelName, simunits, simparams);
  }
  else if (intparams["tabulated_kernel"] == 0) {
    // Depending on the kernel, instantiate a different GradSph object
    if (KernelName == "m4") {
      sph = new Forceless<ndim, M4Kernel>
        (intparams["hydro_forces"], intparams["self_gravity"], floatparams["alpha_visc"],
         floatparams["beta_visc"], floatparams["h_fac"], floatparams["h_converge"],
         avisc, acond, tdavisc, stringparams["gas_eos"], KernelName, simunits, simparams);
    }
    else if (KernelName == "quintic") {
      sph = new Forceless<ndim, QuinticKernel>
        (intparams["hydro_forces"], intparams["self_gravity"], floatparams["alpha_visc"],
         floatparams["beta_visc"], floatparams["h_fac"], floatparams["h_converge"],
         avisc, acond, tdavisc, stringparams["gas_eos"], KernelName, simunits, simparams);
    }
    else if (KernelName == "gaussian") {
      sph = new Forceless<ndim, GaussianKernel>
        (intparams["hydro_forces"], intparams["self_gravity"], floatparams["alpha_visc"],
         floatparams["beta_visc"], floatparams["h_fac"], floatparams["h_converge"],
         avisc, acond, tdavisc, stringparams["gas_eos"], KernelName, simunits, simparams);
    }
    else {
      string message = "Unrecognised parameter : kernel = " + simparams->stringparams["kernel"];
      ExceptionHandler::getIstance().raise(message);
    }
  }
  else {
    string message = "Invalid option for the tabulated_kernel parameter: " +
      stringparams["tabulated_kernel"];
    ExceptionHandler::getIstance().raise(message);
  }


  // Create SPH particle integration object
  //-----------------------------------------------------------------------------------------------
  if (stringparams["sph_integration"] == "lfkdk") {
    sphint = new SphLeapfrogKDK<ndim, GradhSphParticle>
      (floatparams["accel_mult"], floatparams["courant_mult"],
       floatparams["energy_mult"], eos_type, tdavisc);
  }
  else if (stringparams["sph_integration"] == "lfdkd") {
    sphint = new SphLeapfrogDKD<ndim, GradhSphParticle>
      (floatparams["accel_mult"], floatparams["courant_mult"],
       floatparams["energy_mult"], eos_type, tdavisc);
    integration_step = max(integration_step,2);
  }
  else {
    string message = "Unrecognised parameter : sph_integration = "
      + simparams->stringparams["sph_integration"];
    ExceptionHandler::getIstance().raise(message);
  }


  // Energy integration object
  //-----------------------------------------------------------------------------------------------
 /* if (stringparams["energy_integration"] == "Radws" ||
      stringparams["energy_integration"] == "radws"||
      stringparams["energy_integration"] == "rad_ws") {
    uint = new EnergyRadws<ndim, GradhSphParticle>
      (floatparams["energy_mult"], stringparams["radws_table"],
       floatparams["temp_ambient"], &simunits, sph->eos);
  }
  else if (stringparams["energy_integration"] == "lociso") {
     uint = new LocalIsotherm<ndim,GradhSphParticle>(floatparams["temp0"], floatparams["gamma_eos"],floatparams["templaw"], floatparams["mu_bar"], &simunits, nbody);
  }
  else if (stringparams["energy_integration"] == "null" ||
           stringparams["energy_integration"] == "none") {
    uint = new NullEnergy<ndim>(floatparams["energy_mult"]);
  }
  else {
    string message = "Unrecognised parameter : energy_integration = "
      + simparams->stringparams["energy_integration"];
    ExceptionHandler::getIstance().raise(message);
  }*/

  //For forceless - always null energy
  uint = new NullEnergy<ndim>(floatparams["energy_mult"]);

  //-----------------------------------------------------------------------------------------------
#if defined MPI_PARALLEL
  if (stringparams["mpi_decomposition"] == "kdtree") {
    mpicontrol = new MpiKDTreeDecomposition<ndim,GradhSphParticle>();
  }
  else {
    string message = "Unrecognised parameter : mpi_decomposition = "
      + simparams->stringparams["mpi_decomposition"];
    ExceptionHandler::getIstance().raise(message);
  }
  mpicontrol->timing = timing;
  rank = mpicontrol->rank;
  Nmpi = mpicontrol->Nmpi;
#endif



  // Create neighbour searching object based on chosen method in params file
  //-----------------------------------------------------------------------------------------------
  // Here I do a horrible hack to get at the underlying tree, needed for the dust.
  // TreeBase<ndim> * t = NULL, * gt = NULL, *mpit = NULL ;
/*
  if (stringparams["neib_search"] == "bruteforce") {
    sphneib = new GradhSphBruteForce<ndim,GradhSphParticle>
     (sph->kernp->kernrange, &simbox, sph->kernp, timing);
  }
  else if (stringparams["neib_search"] == "kdtree") {
    sphneib = new GradhSphKDTree<ndim,GradhSphParticle,KDTreeCell>
     (intparams["Nleafmax"], Nmpi, intparams["pruning_level_min"], intparams["pruning_level_max"],
      floatparams["thetamaxsqd"], sph->kernp->kernrange, floatparams["macerror"],
      stringparams["gravity_mac"], stringparams["multipole"], &simbox, sph->kernp, timing, sph->types);

   // typedef GradhSphKDTree<ndim,GradhSphParticle,KDTreeCell> TreeType ;
    //TreeType *pTree = reinterpret_cast<TreeType*>(sphneib) ;
    //t = pTree->tree ; gt = pTree->ghosttree ;
//#ifdef MPI_PARALLEL
  //  mpit = pTree->mpighosttree ;
//#endif
  }
  else if (stringparams["neib_search"] == "octtree" && gas_radiation == "treeray" && ndim == 3) {
    sphneib = new GradhSphOctTree<ndim,GradhSphParticle,TreeRayCell>
     (intparams["Nleafmax"], Nmpi, intparams["pruning_level_min"], intparams["pruning_level_max"],
      floatparams["thetamaxsqd"], sph->kernp->kernrange, floatparams["macerror"],
      stringparams["gravity_mac"], stringparams["multipole"], &simbox, sph->kernp, timing, sph->types);

    //typedef GradhSphOctTree<ndim,GradhSphParticle,TreeRayCell> TreeType ;
   // TreeType *pTree = reinterpret_cast<TreeType*>(sphneib) ;
   // t = pTree->tree ; gt = pTree->ghosttree ;
//#ifdef MPI_PARALLEL
 //   mpit = pTree->mpighosttree ;
//#endif
  }
  else if (stringparams["neib_search"] == "octtree") {
    sphneib = new GradhSphOctTree<ndim,GradhSphParticle,OctTreeCell>
     (intparams["Nleafmax"], Nmpi, intparams["pruning_level_min"], intparams["pruning_level_max"],
      floatparams["thetamaxsqd"], sph->kernp->kernrange, floatparams["macerror"],
      stringparams["gravity_mac"], stringparams["multipole"], &simbox, sph->kernp, timing, sph->types);

    //typedef GradhSphOctTree<ndim,GradhSphParticle,OctTreeCell> TreeType ;
    //TreeType *pTree = reinterpret_cast<TreeType*>(sphneib) ;
   // t = pTree->tree ; gt = pTree->ghosttree ;
//#ifdef MPI_PARALLEL
 //   mpit = pTree->mpighosttree ;
//#endif
  }
  else {
    string message = "Unrecognised parameter : neib_search = "
      + simparams->stringparams["neib_search"];
    ExceptionHandler::getIstance().raise(message);
  }
  //sphneib->kernp = sph->kernp;
  sphneib->kernfac = sph->kernfac;
  //sphneib->kernrange = sph->kernp->kernrange;


#if defined MPI_PARALLEL
  mpicontrol->SetNeibSearch(sphneib);
#endif*/

  //Forceless Tree option - currently take GradhSph tree
  sphneib = new ForcelessTree<ndim,GradhSphParticle,KDTreeCell>
  (intparams["Nleafmax"], Nmpi, intparams["pruning_level_min"], intparams["pruning_level_max"],
   floatparams["thetamaxsqd"], sph->kernp->kernrange, floatparams["macerror"],
   stringparams["gravity_mac"], stringparams["multipole"], &simbox, sph->kernp, timing, sph->types);

  // Radiation transport object
  //-----------------------------------------------------------------------------------------------
 /* if (gas_radiation == "treeray" && ndim == 3) {
    radiation = new TreeRay<ndim,1,GradhSphParticle,TreeRayCell>
      (Nmpi, intparams["on_the_spot"], intparams["nside"], intparams["ilNR"], intparams["ilNTheta"],
       intparams["ilNPhi"], intparams["ilNNS"], intparams["ilFinePix"], floatparams["maxDist"],
       floatparams["rayRadRes"], floatparams["relErr"], stringparams["errControl"],
       simbox, &simunits, simparams, sphneib);
  }
  else if (gas_radiation == "ionisation") {
    radiation = new MultipleSourceIonisation<ndim,GradhSphParticle>
      (sphneib, floatparams["mu_bar"], floatparams["mu_ion"], floatparams["temp0"],
       floatparams["temp_ion"], floatparams["Ndotmin"], floatparams["gamma_eos"],
       pow(simunits.r.outscale*simunits.r.outcgs, 3.)/
       pow(simunits.m.outscale*simunits.m.outcgs, 2.),
       simunits.temp.outscale, pow(simunits.r.outscale*simunits.r.outcgs,-4)*
       pow(simunits.t.outscale*simunits.t.outcgs,+2)/simunits.m.outscale*simunits.m.outcgs);
  }
  else if (gas_radiation == "monoionisation") {
    radiation = new MonochromaticIonisationMonteCarlo<ndim,1,GradhSphParticle,MonoIonTreeCell>
      (intparams["Nleafmax"], intparams["Nraditerations"], intparams["Nradlevels"],
       floatparams["Nphotonratio"], floatparams["temp_ion"], floatparams["arecomb"],
       floatparams["NLyC"], stringparams["rand_algorithm"], &simunits, sph->eos);
  }
  else if (gas_radiation == "none") {
    radiation = new NullRadiation<ndim>();
  }
  else {
    string message = "Unrecognised parameter : radiation = " + gas_radiation;
    ExceptionHandler::getIstance().raise(message);
  }*/

  //Forceless - radiation always null
  radiation = new NullRadiation<ndim>();

  // Create ghost particle object
  //-----------------------------------------------------------------------------------------------
  if (IsAnyBoundarySpecial(simbox)) {
    LocalGhosts = new PeriodicGhostsSpecific<ndim,GradhSphParticle >();
  }
  else {
    LocalGhosts = new NullGhosts<ndim>();
  }
#ifdef MPI_PARALLEL
  MpiGhosts = new MpiGhostsSpecific<ndim, GradhSphParticle>(mpicontrol);
#endif


  // Depending on the dimensionality, calculate expected neighbour number
  //-----------------------------------------------------------------------------------------------
  if (ndim == 1) {
    sph->Ngather = (int) (2.0*sph->kernp->kernrange*sph->h_fac);
  }
  else if (ndim == 2) {
    sph->Ngather = (int) (pi*pow(sph->kernp->kernrange*sph->h_fac, 2));
  }
  else if (ndim == 3) {
    sph->Ngather = (int) (4.0*pi*pow(sph->kernp->kernrange*sph->h_fac, 3)/3.0);
  }

  // Setup the dust force object
  //-----------------------------------------------------------------------------------------------
  //sphdust = DustFactory<ndim, GradhSphParticle>::ProcessParameters(simparams, timing, sph->types, t, gt, mpit) ;

  return;
}

//=================================================================================================
//  ForcelessSimulation::SmoothParticleQuantity
/// ...
//=================================================================================================
template <int ndim>
void ForcelessSimulation<ndim>::SmoothParticleQuantity
 (const int Npart,
  FLOAT *values)
{

  return;
}
