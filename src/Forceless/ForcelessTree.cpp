//=================================================================================================
//  ForcelessTree.cpp
//  Contains all functions for building, stocking and walking for the
//  binary KD tree for SPH particles.
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


#include <cstdlib>
#include <cassert>
#include <iostream>
#include <string>
#include <math.h>
#include "Precision.h"
#include "Exception.h"
#include "SphNeighbourSearch.h"
#include "Sph.h"
#include "Parameters.h"
#include "InlineFuncs.h"
#include "Particle.h"
#include "Debug.h"
#if defined _OPENMP
#include <omp.h>
#endif
using namespace std;



//=================================================================================================
//  ForcelessTree::ForcelessTree
/// ForcelessTree constructor.  Initialises various variables.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
ForcelessTree<ndim,ParticleType,TreeCell>::ForcelessTree
 (int _Nleafmax, int _Nmpi, int _pruning_level_min, int _pruning_level_max, FLOAT _thetamaxsqd,
  FLOAT _kernrange, FLOAT _macerror, string _gravity_mac, string _multipole,
  DomainBox<ndim>* _box, SmoothingKernel<ndim>* _kern, CodeTiming* _timing):
 NeighbourSearch<ndim>(_kernrange, _box, _kern, _timing),
 SphTree<ndim,ParticleType,TreeCell>
  (_Nleafmax, _Nmpi, _pruning_level_min, _pruning_level_max, _thetamaxsqd,
   _kernrange, _macerror, _gravity_mac, _multipole, _box, _kern, _timing)
{
}



//=================================================================================================
//  ForcelessTree::~ForcelessTree
/// ForcelessTree destructor.  Deallocates tree memory upon object destruction.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
ForcelessTree<ndim,ParticleType,TreeCell>::~ForcelessTree()
{
  if (tree->allocated_tree) {
    this->DeallocateMemory();
    tree->DeallocateTreeMemory();
  }
}



//=================================================================================================
//  ForcelessKDTree::ForcelessKDTree
/// ForcelessKDTree constructor.  Initialises various variables and creates tree objects.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
ForcelessKDTree<ndim,ParticleType,TreeCell>::ForcelessKDTree
 (int _Nleafmax, int _Nmpi, int _pruning_level_min, int _pruning_level_max, FLOAT _thetamaxsqd,
  FLOAT _kernrange, FLOAT _macerror, string _gravity_mac, string _multipole,
  DomainBox<ndim>* _box, SmoothingKernel<ndim>* _kern, CodeTiming* _timing, ParticleTypeRegister& types):
 NeighbourSearch<ndim>(_kernrange, _box, _kern, _timing),
 ForcelessTree<ndim,ParticleType,TreeCell>
  (_Nleafmax, _Nmpi, _pruning_level_min, _pruning_level_max, _thetamaxsqd,
   _kernrange, _macerror, _gravity_mac, _multipole, _box, _kern, _timing)
{
  // Set-up main tree object
  tree = new KDTree<ndim,ParticleType,TreeCell>(_Nleafmax, _thetamaxsqd, _kernrange,
                                                _macerror, _gravity_mac, _multipole, *_box, types);

  // Set-up ghost-particle tree object
  ghosttree = new KDTree<ndim,ParticleType,TreeCell>(_Nleafmax, _thetamaxsqd, _kernrange,
                                                     _macerror, _gravity_mac, _multipole, *_box, types);

#ifdef MPI_PARALLEL
  // Set-up ghost-particle tree object
  mpighosttree = new KDTree<ndim,ParticleType,TreeCell>(_Nleafmax, _thetamaxsqd, _kernrange,
                                                        _macerror, _gravity_mac, _multipole, *_box, types);

  // Set-up multiple pruned trees, one for each MPI process
  KDTree<ndim,ParticleType,TreeCell>** prunedtree_derived = new KDTree<ndim,ParticleType,TreeCell>*[Nmpi];
  prunedtree = (Tree<ndim,ParticleType,TreeCell> **) prunedtree_derived;
  KDTree<ndim,ParticleType,TreeCell>** sendprunedtree_derived = new KDTree<ndim,ParticleType,TreeCell>*[Nmpi];
  sendprunedtree = (Tree<ndim,ParticleType,TreeCell> **) sendprunedtree_derived;

  for (int i=0; i<Nmpi; i++) {
    prunedtree[i] = new KDTree<ndim,ParticleType,TreeCell>
     (_Nleafmax, _thetamaxsqd, _kernrange, _macerror, _gravity_mac, _multipole, *_box, types);
  }
  for (int i=0; i<Nmpi; i++) {
    sendprunedtree[i] = new KDTree<ndim,ParticleType,TreeCell>
     (_Nleafmax, _thetamaxsqd, _kernrange, _macerror, _gravity_mac, _multipole, *_box, types);
  }
#endif
}


//=================================================================================================
//  ForcelessTree::UpdateAllSphProperties
/// Update all gather SPH properties (e.g. rho, div_v) for all active particles in domain.
/// Loops over all cells containing active particles, performs a tree walk for all particles in
/// the cell, and then calls SPH class routine to compute properties from neighbours.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void ForcelessTree<ndim,ParticleType,TreeCell>::UpdateAllSphProperties
 (int Nhydro,                              ///< [in] No. of SPH particles
  int Ntot,                                ///< [in] No. of SPH + ghost particles
  SphParticle<ndim> *sph_gen,              ///< [inout] Pointer to SPH ptcl array
  Sph<ndim> *sph,                          ///< [in] Pointer to SPH object
  Nbody<ndim> *nbody)                      ///< [in] Pointer to N-body object
{
  int cactive;                             // No. of active tree cells
  TreeCell<ndim> *celllist;                // List of active tree cells
  ParticleType<ndim> *sphdata = static_cast<ParticleType<ndim>* > (sph_gen);
#ifdef MPI_PARALLEL
  int Nactivetot = 0;                      // Total number of active particles
  double twork = timing->WallClockTime();  // Start time (for load balancing)
#endif

  debug2("[ForcelessTree::UpdateAllSphProperties]");
  timing->StartTimingSection("SPH_FORCELESS_PROPERTIES");


  // Find list of all cells that contain active particles
  celllist = new TreeCell<ndim>[tree->gtot];
  cactive = tree->ComputeActiveCellList(celllist);
  assert(cactive <= tree->gtot);

  // If there are no active cells, return to main loop
  if (cactive == 0) {
    delete[] celllist;
    timing->EndTimingSection("SPH_FORCELESS_PROPERTIES");
    return;
  }


  // Set-up all OMP threads
  //===============================================================================================
#pragma omp parallel default(none) shared(cactive,celllist,cout,nbody,sph,sphdata,Ntot)
  {
#if defined _OPENMP
    const int ithread = omp_get_thread_num();
#else
    const int ithread = 0;
#endif
    int celldone;                              // Flag if cell is done
    int cc;                                    // Aux. cell counter
    int i;                                     // Particle id
    int j;                                     // Aux. particle counter
    int jj;                                    // Aux. particle counter
    int k;                                     // Dimension counter
    int Nactive;                               // No. of active particles in cell
    int Ngather;                               // No. of gather neighbours
    int Nneib;                                 // No. of neighbours from tree-walk
    int okflag;                                // Flag if particle is done
    FLOAT draux[ndim];                         // Aux. relative position vector var
    FLOAT drsqdaux;                            // Distance squared
    FLOAT hrangesqd;                           // Kernel extent
    FLOAT hmax;                                // Maximum smoothing length
    FLOAT rp[ndim];                            // Local copy of particle position
    FLOAT *mu2 = 0;                            // Trimmed array (dummy for grad-h)
    int Nneibmax = Nneibmaxbuf[ithread];       // Local copy of neighbour buffer size
    int* activelist = activelistbuf[ithread];  // Local array of active particle ids
    int* neiblist = new int[Nneibmax];         // Local array of neighbour particle ids
    int* ptype    = new int[Nneibmax];         // Local array of particle types
    FLOAT* gpot   = new FLOAT[Nneibmax];       // Local array of particle potentials
    FLOAT* gpot2  = new FLOAT[Nneibmax];       // Local reduced array of neighbour potentials
    FLOAT* drsqd  = new FLOAT[Nneibmax];       // Local array of distances (squared)
    FLOAT* m      = new FLOAT[Nneibmax];       // Local array of particle masses
    FLOAT* m2     = new FLOAT[Nneibmax];       // Local reduced array of neighbour masses
    FLOAT* r      = new FLOAT[Nneibmax*ndim];  // Local array of particle positions
    ParticleType<ndim>* activepart = activepartbuf[ithread];   // Local array of active particles


    // Loop over all active cells
    //=============================================================================================
#pragma omp for schedule(guided)
    for (cc=0; cc<cactive; cc++) {
      TreeCell<ndim>& cell = celllist[cc];
      celldone = 1;
      hmax = cell.hmax;

      // If hmax is too small so the neighbour lists are invalid, make hmax
      // larger and then recompute for the current active cell.
      //-------------------------------------------------------------------------------------------
      do {
        hmax = (FLOAT) 1.05*hmax;
        celldone = 1;

        // Find list of active particles in current cell
        Nactive = tree->ComputeActiveParticleList(cell, sphdata, activelist);
        for (j=0; j<Nactive; j++) activepart[j] = sphdata[activelist[j]];

        // Compute neighbour list for cell from particles on all trees
        Nneib = 0;
        Nneib = tree->ComputeGatherNeighbourList(cell,sphdata,hmax,Nneibmax,Nneib,neiblist);
        Nneib = ghosttree->ComputeGatherNeighbourList(cell,sphdata,hmax,Nneibmax,Nneib,neiblist);
#ifdef MPI_PARALLEL
        Nneib = mpighosttree->ComputeGatherNeighbourList(cell,sphdata,hmax,Nneibmax,Nneib,neiblist);
#endif

        // If there are too many neighbours so the buffers are filled,
        // reallocate the arrays and recompute the neighbour lists.
        while (Nneib == -1) {
          delete[] r;
          delete[] m2;
          delete[] m;
          delete[] drsqd;
          delete[] gpot2;
          delete[] gpot;
          delete[] ptype;
          delete[] neiblist;
          Nneibmax = 2*Nneibmax;
          neiblist = new int[Nneibmax];
          ptype    = new int[Nneibmax];
          gpot     = new FLOAT[Nneibmax];
          gpot2    = new FLOAT[Nneibmax];
          drsqd    = new FLOAT[Nneibmax];
          m        = new FLOAT[Nneibmax];
          m2       = new FLOAT[Nneibmax];
          r        = new FLOAT[Nneibmax*ndim];
          Nneib = 0;
          Nneib = tree->ComputeGatherNeighbourList(cell,sphdata,hmax,Nneibmax,Nneib,neiblist);
          Nneib = ghosttree->ComputeGatherNeighbourList(cell,sphdata,hmax,Nneibmax,Nneib,neiblist);
#ifdef MPI_PARALLEL
          Nneib = mpighosttree->ComputeGatherNeighbourList(cell,sphdata,hmax,
                                                           Nneibmax,Nneib,neiblist);
#endif
        };

        // Make local copies of important neib information (mass and position)
        for (jj=0; jj<Nneib; jj++) {
          j         = neiblist[jj];
          gpot[jj]  = sphdata[j].gpot;
          m[jj]     = sphdata[j].m;
          ptype[jj] = sphdata[j].ptype;
          for (k=0; k<ndim; k++) r[ndim*jj + k] = sphdata[j].r[k];
        }

        // Loop over all active particles in the cell
        //-----------------------------------------------------------------------------------------
        for (j=0; j<Nactive; j++) {
          i = activelist[j];
          for (k=0; k<ndim; k++) rp[k] = activepart[j].r[k];

          // Set gather range as current h multiplied by some tolerance factor
          hrangesqd = kernrangesqd*hmax*hmax;
          Ngather = 0;

          // Compute distance (squared) to all
          //---------------------------------------------------------------------------------------
          for (jj=0; jj<Nneib; jj++) {

            // Only include particles of appropriate types in density calculation
            if (!sph->types[activepart[j].ptype].hmask[ptype[jj]]) continue ;

            for (k=0; k<ndim; k++) draux[k] = r[ndim*jj + k] - rp[k];
            drsqdaux = DotProduct(draux,draux,ndim) + small_number;

            // Record distance squared and masses for all potential gather neighbours
            if (drsqdaux <= hrangesqd) {
              gpot[Ngather]  = gpot[jj];
              drsqd[Ngather] = drsqdaux;
              m2[Ngather]    = m[jj];
              Ngather++;
            }

          }
          //---------------------------------------------------------------------------------------

          // Validate that gather neighbour list is correct
#if defined(VERIFY_ALL)
          if (neibcheck) this->CheckValidNeighbourList(i, Ntot, Nneib, neiblist, sphdata, "gather");
#endif

          // Compute smoothing length and other gather properties for ptcl i
          okflag = sph->ComputeH(i, Ngather, hmax, m2, mu2, drsqd, gpot, activepart[j], nbody);

          // If h-computation is invalid, then break from loop and recompute
          // larger neighbour lists
          if (okflag == 0) {
            celldone = 0;
            break;
          }

        }
        //-----------------------------------------------------------------------------------------

      } while (celldone == 0);
      //-------------------------------------------------------------------------------------------

      // Once cell is finished, copy all active particles back to main memory
      for (j=0; j<Nactive; j++) sphdata[activelist[j]] = activepart[j];


    }
    //=============================================================================================


    // Free-up all memory
    delete[] r;
    delete[] m2;
    delete[] m;
    delete[] drsqd;
    delete[] gpot2;
    delete[] gpot;
    delete[] ptype;
    delete[] neiblist;

  }
  //===============================================================================================

  // Compute time spent in routine and in each cell for load balancing
#ifdef MPI_PARALLEL
  twork = timing->WallClockTime() - twork;
  for (int cc=0; cc<cactive; cc++) Nactivetot += celllist[cc].Nactive;
  for (int cc=0; cc<cactive; cc++) {
    int c = celllist[cc].id;
    tree->celldata[c].worktot += twork*(DOUBLE) tree->celldata[c].Nactive / (DOUBLE) Nactivetot;
  }
#ifdef OUTPUT_ALL
  cout << "Time computing smoothing lengths : " << twork << "     Nactivetot : " << Nactivetot << endl;
#endif
#endif

  delete[] celllist;

  // Update tree smoothing length values here
  tree->UpdateHmaxValues(tree->celldata[0],sphdata);

  timing->EndTimingSection("SPH_FORCELESS_PROPERTIES");

  return;
}


//=================================================================================================
// ForcelessTree::UpdateAllSphForces
/// Compute all forces on active SPH particles (hydro + gravity) for periodic boundary conditions.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void ForcelessTree<ndim,ParticleType,TreeCell>::UpdateAllSphForces
 (int Nhydro,                          ///< [in] No. of SPH particles
  int Ntot,                            ///< [in] No. of SPH + ghost particles
  SphParticle<ndim> *sph_gen,          ///< [inout] Pointer to SPH ptcl array
  Sph<ndim> *sph,                      ///< [in] Pointer to SPH object
  Nbody<ndim> *nbody,                  ///< [in] Pointer to N-body object
  DomainBox<ndim> &simbox,             ///< [in] Simulation domain box
  Ewald<ndim> *ewald)                  ///< [in] Ewald gravity object pointer
{
  int cactive;                         // No. of active cells
  TreeCell<ndim> *celllist;            // List of active cells
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  timing->StartTimingSection("FORCELESS_ALL_FORCES");


  // Find list of all cells that contain active particles
#if defined (MPI_PARALLEL)
  celllist = new TreeCell<ndim>[tree->Ncellmax];
#else
  celllist = new TreeCell<ndim>[tree->gtot];
#endif
  cactive = tree->ComputeActiveCellList(celllist);

  // If there are no active cells, return to main loop
  if (cactive == 0) {
    delete[] celllist;
    timing->EndTimingSection("FORCELESS_ALL_FORCES");
    return;
  }


  // Set-up all OMP threads
  //===============================================================================================
  #pragma omp parallel default(none) shared(celllist,cactive,ewald,nbody,simbox,sph,sphdata,cout)
  {
#if defined _OPENMP
    const int ithread = omp_get_thread_num();
#else
    const int ithread = 0;
#endif
    int cc;                                      // Aux. cell counter
    int i;                                       // Particle id
    int j;                                       // Aux. particle counter
    int k;                                       // Dimension counter
    int Nactive;                                 // ..
    //int Nhydroaux;                               // ..
    //FLOAT aperiodic[ndim];                       // ..
    int *activelist  = activelistbuf[ithread];   // ..
    ParticleType<ndim>* activepart = activepartbuf[ithread];   // ..

    // Loop over all active cells
    //=============================================================================================
#pragma omp for schedule(guided)
    for (cc=0; cc<cactive; cc++) {
      TreeCell<ndim> &cell = celllist[cc];

      // Find list of active particles in current cell
      Nactive = tree->ComputeActiveParticleList(cell, sphdata, activelist);

      // Make local copies of active particles
      for (j=0; j<Nactive; j++) activepart[j] = sphdata[activelist[j]];

      // Zero/initialise all summation variables for active particles
      for (j=0; j<Nactive; j++) {
        activepart[j].div_v     = (FLOAT) 0.0;
        activepart[j].dudt      = (FLOAT) 0.0;
        activepart[j].levelneib = 0;
        activepart[j].gpot      = 0;
        for (k=0; k<ndim; k++) activepart[j].a[k]     = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) activepart[j].agrav[k] = (FLOAT) 0.0;
      }


      // Compute all star forces for active particles
      for (j=0; j<Nactive; j++) {
        sph->ComputeStarGravForces(nbody->Nnbody, nbody->nbodydata, activepart[j]);
      }

      // Add all active particles contributions to main array
      for (j=0; j<Nactive; j++) {
        i = activelist[j];
        for (k=0; k<ndim; k++) sphdata[i].a[k]     += activepart[j].a[k] + activepart[j].agrav[k];
        for (k=0; k<ndim; k++) sphdata[i].agrav[k] += activepart[j].agrav[k];
        sphdata[i].gpot  += activepart[j].gpot;
      }
    }
  }
  //===============================================================================================

  delete[] celllist;

  timing->EndTimingSection("FORCELESS_ALL_FORCES");

  return;
}


template class ForcelessTree<1,GradhSphParticle,KDTreeCell>;
template class ForcelessTree<2,GradhSphParticle,KDTreeCell>;
template class ForcelessTree<3,GradhSphParticle,KDTreeCell>;

template class ForcelessKDTree<1,GradhSphParticle,KDTreeCell>;
template class ForcelessKDTree<2,GradhSphParticle,KDTreeCell>;
template class ForcelessKDTree<3,GradhSphParticle,KDTreeCell>;
