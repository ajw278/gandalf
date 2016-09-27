//=================================================================================================
//  SphNeighbourSearch.h
//  Header file containing class definitions for all SPH neighbour searching algorithms.
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


#ifndef _SPH_NEIGHBOUR_SEARCH_H_
#define _SPH_NEIGHBOUR_SEARCH_H_


#include <assert.h>
#include <iostream>
#include <string>
#include <vector>
#include "Precision.h"
#include "Constants.h"
#include "CodeTiming.h"
#include "GhostNeighbours.hpp"
#include "InlineFuncs.h"
#include "Nbody.h"
#include "NeighbourSearch.h"
#include "SmoothingKernel.h"
#include "Particle.h"
#include "Sph.h"
#include "DomainBox.h"
#include "Ewald.h"
#include "Parameters.h"
#include "KDTree.h"
#include "OctTree.h"
#include "BruteForceTree.h"
#if defined MPI_PARALLEL
#include "MpiExport.h"
#include "MpiNode.h"
#endif
using namespace std;


// Forward declaration of Sph to break circular dependency
template <int ndim>
class Sph;

// Forward declare MpiNode to break circular dependency
#if defined MPI_PARALLEL
template <int ndim>
class MpiNode;
#endif



//=================================================================================================
//  Class SphNeighbourSearch
/// \brief   SphNeighbourSearch class definition.
/// \details Class for creating the SPH neighbour search data structure, and for computing local
///          neighbour lists and calling SPH functions (e.g. computing h, SPH forces, etc..).
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim>
class SphNeighbourSearch : public virtual NeighbourSearch<ndim>
{
/*#if defined MPI_PARALLEL
protected:
  vector<int> ids_active_particles;
#endif*/
 public:

  using NeighbourSearch<ndim>::neibcheck;
  using NeighbourSearch<ndim>::box;
  using NeighbourSearch<ndim>::timing;
  using NeighbourSearch<ndim>::kernp;
  using NeighbourSearch<ndim>::kernfac;
  using NeighbourSearch<ndim>::kernrange;
  using NeighbourSearch<ndim>::kernrangesqd;


  //-----------------------------------------------------------------------------------------------
  SphNeighbourSearch(FLOAT kernrangeaux, DomainBox<ndim> *boxaux,
                     SmoothingKernel<ndim> *kernaux, CodeTiming *timingaux) :
    NeighbourSearch<ndim>(kernrangeaux, boxaux, kernaux, timingaux) {};
  virtual ~SphNeighbourSearch() {};


  //-----------------------------------------------------------------------------------------------
  virtual void UpdateAllSphProperties(int, int, SphParticle<ndim> *,
                                      Sph<ndim> *, Nbody<ndim> *) = 0;
  virtual void UpdateAllSphHydroForces(int, int, SphParticle<ndim> *, Sph<ndim> *,
                                       Nbody<ndim> *, DomainBox<ndim> &) =0 ;
  virtual void UpdateAllSphForces(int, int, SphParticle<ndim> *, Sph<ndim> *,
                                  Nbody<ndim> *, DomainBox<ndim> &, Ewald<ndim> *) = 0;
  virtual void UpdateAllSphGravForces(int, int, SphParticle<ndim> *, Sph<ndim> *,
                                      Nbody<ndim> *, DomainBox<ndim> &, Ewald<ndim> *) = 0;

  //virtual void UpdateAllStarGasForces(int, int, SphParticle<ndim> *,
  //                                    Sph<ndim> *, Nbody<ndim> *) = 0;

};


//=================================================================================================
//  Class SphTree
/// \brief   Class containing tree for efficient SPH neighbour searching and gravity calculations.
/// \details Class containing tree for efficient SPH neighbour searching and gravity calculations.
/// \author  D. A. Hubber
/// \date    08/01/2014
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
class SphTree : public SphNeighbourSearch<ndim>, public HydroTree<ndim,ParticleType,TreeCell>
{
#if defined MPI_PARALLEL
  vector<vector<int> > ids_sent_particles;
protected:
  using NeighbourSearch<ndim>::ids_active_particles;
  using NeighbourSearch<ndim>::N_imported_part_per_proc;
#endif
 public:

  using HydroTree<ndim,ParticleType,TreeCell>::activelistbuf;
  using HydroTree<ndim,ParticleType,TreeCell>::activepartbuf;
  using HydroTree<ndim,ParticleType,TreeCell>::allocated_buffer;
  using HydroTree<ndim,ParticleType,TreeCell>::box;
  using HydroTree<ndim,ParticleType,TreeCell>::cellbuf;
  using HydroTree<ndim,ParticleType,TreeCell>::gravity_mac;
  using HydroTree<ndim,ParticleType,TreeCell>::kernp;
  using HydroTree<ndim,ParticleType,TreeCell>::kernrange;
  using HydroTree<ndim,ParticleType,TreeCell>::kernrangesqd;
  using HydroTree<ndim,ParticleType,TreeCell>::levelneibbuf;
  using HydroTree<ndim,ParticleType,TreeCell>::multipole;
  using HydroTree<ndim,ParticleType,TreeCell>::neibcheck;
  using HydroTree<ndim,ParticleType,TreeCell>::neibpartbuf;
  using HydroTree<ndim,ParticleType,TreeCell>::Ngravcellmaxbuf;
  using HydroTree<ndim,ParticleType,TreeCell>::Nleafmax;
  using HydroTree<ndim,ParticleType,TreeCell>::Nneibmaxbuf;
  using HydroTree<ndim,ParticleType,TreeCell>::Nthreads;
  using HydroTree<ndim,ParticleType,TreeCell>::Ntot;
  using HydroTree<ndim,ParticleType,TreeCell>::Ntotmax;
  using HydroTree<ndim,ParticleType,TreeCell>::Ntotmaxold;
  using HydroTree<ndim,ParticleType,TreeCell>::Ntotold;
  using HydroTree<ndim,ParticleType,TreeCell>::timing;
  using HydroTree<ndim,ParticleType,TreeCell>::tree;
  using HydroTree<ndim,ParticleType,TreeCell>::ghosttree;
#ifdef MPI_PARALLEL
  using HydroTree<ndim,ParticleType,TreeCell>::mpighosttree;
  using HydroTree<ndim,ParticleType,TreeCell>::Nmpi;
  using HydroTree<ndim,ParticleType,TreeCell>::prunedtree;
  using HydroTree<ndim,ParticleType,TreeCell>::sendprunedtree;
#endif


  //-----------------------------------------------------------------------------------------------
  SphTree(int _Nleafmax, int _Nmpi, int _pruning_level_min, int _pruning_level_max,
          FLOAT _thetamaxsqd, FLOAT _kernrange, FLOAT _macerror,
          string _gravity_mac, string _multipole,
          DomainBox<ndim> *_box, SmoothingKernel<ndim> *_kern, CodeTiming *_timing) :
    NeighbourSearch<ndim>(_kernrange, _box, _kern, _timing),
    SphNeighbourSearch<ndim>(_kernrange, _box, _kern, _timing),
    HydroTree<ndim,ParticleType,TreeCell>(_Nleafmax, _Nmpi, _pruning_level_min, _pruning_level_max,
                                          _thetamaxsqd, _kernrange, _macerror, _gravity_mac,
                                          _multipole, _box, _kern, _timing) {};
  virtual ~SphTree() {};


  //-----------------------------------------------------------------------------------------------
  void UpdateAllSphProperties(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *) = 0;
  void UpdateAllSphHydroForces(int, int, SphParticle<ndim> *, Sph<ndim> *,
                              Nbody<ndim> *, DomainBox<ndim> &) =0;
  void UpdateAllSphForces(int, int, SphParticle<ndim> *, Sph<ndim> *,
                          Nbody<ndim> *, DomainBox<ndim> &, Ewald<ndim> *) =0;
  void UpdateAllSphGravForces(int, int, SphParticle<ndim> *, Sph<ndim> *,
                              Nbody<ndim> *, DomainBox<ndim> &, Ewald<ndim> *) =0;
  //void UpdateAllStarGasForces(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *) =0;

};


//=================================================================================================
//  Class GradhSphTree
/// \brief   Class containing tree for computing grad-h SPH summation and force loops.
/// \details Class containing tree for computing grad-h SPH summation and force loops.
/// \author  D. A. Hubber
/// \date    08/01/2014
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
class GradhSphTree : public SphTree<ndim,ParticleType,TreeCell>
{
 public:

  using SphTree<ndim,ParticleType,TreeCell>::activelistbuf;
  using SphTree<ndim,ParticleType,TreeCell>::activepartbuf;
  using SphTree<ndim,ParticleType,TreeCell>::allocated_buffer;
  using SphTree<ndim,ParticleType,TreeCell>::box;
  using SphTree<ndim,ParticleType,TreeCell>::cellbuf;
  using SphTree<ndim,ParticleType,TreeCell>::gravity_mac;
  using SphTree<ndim,ParticleType,TreeCell>::kernp;
  using SphTree<ndim,ParticleType,TreeCell>::kernrange;
  using SphTree<ndim,ParticleType,TreeCell>::kernrangesqd;
  using SphTree<ndim,ParticleType,TreeCell>::levelneibbuf;
  using SphTree<ndim,ParticleType,TreeCell>::multipole;
  using SphTree<ndim,ParticleType,TreeCell>::neibcheck;
  using SphTree<ndim,ParticleType,TreeCell>::neibpartbuf;
  using SphTree<ndim,ParticleType,TreeCell>::Ngravcellmaxbuf;
  using SphTree<ndim,ParticleType,TreeCell>::Nleafmax;
  using SphTree<ndim,ParticleType,TreeCell>::Nneibmaxbuf;
  using SphTree<ndim,ParticleType,TreeCell>::Nthreads;
  using SphTree<ndim,ParticleType,TreeCell>::Ntot;
  using SphTree<ndim,ParticleType,TreeCell>::Ntotmax;
  using SphTree<ndim,ParticleType,TreeCell>::Ntotmaxold;
  using SphTree<ndim,ParticleType,TreeCell>::Ntotold;
  using SphTree<ndim,ParticleType,TreeCell>::timing;
  using SphTree<ndim,ParticleType,TreeCell>::tree;
  using SphTree<ndim,ParticleType,TreeCell>::ghosttree;
#ifdef MPI_PARALLEL
  using SphTree<ndim,ParticleType,TreeCell>::mpighosttree;
  using SphTree<ndim,ParticleType,TreeCell>::Nmpi;
  using SphTree<ndim,ParticleType,TreeCell>::prunedtree;
  using SphTree<ndim,ParticleType,TreeCell>::sendprunedtree;
#endif


  //-----------------------------------------------------------------------------------------------
  GradhSphTree(int, int, int, int, FLOAT, FLOAT, FLOAT, string, string,
               DomainBox<ndim> *, SmoothingKernel<ndim> *, CodeTiming *,ParticleTypeRegister& types);
  virtual ~GradhSphTree();


  //-----------------------------------------------------------------------------------------------
  void UpdateAllSphProperties(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphHydroForces(int, int, SphParticle<ndim> *, Sph<ndim> *,
                               Nbody<ndim> *, DomainBox<ndim> &);
  void UpdateAllSphForces(int, int, SphParticle<ndim> *, Sph<ndim> *,
                          Nbody<ndim> *, DomainBox<ndim> &, Ewald<ndim> *);
  void UpdateAllSphGravForces(int, int, SphParticle<ndim> *, Sph<ndim> *,
                              Nbody<ndim> *, DomainBox<ndim> &, Ewald<ndim> *);

};

//=================================================================================================
//  Class ForcelessTree
/// \brief   Class containing tree for computing grad-h SPH summation and force loops.
/// \details Class containing tree for computing grad-h SPH summation and force loops.
/// \author  D. A. Hubber
/// \date    08/01/2014
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
class ForcelessTree : public SphTree<ndim,ParticleType,TreeCell>
{
 public:

  using SphTree<ndim,ParticleType,TreeCell>::activelistbuf;
  using SphTree<ndim,ParticleType,TreeCell>::activepartbuf;
  using SphTree<ndim,ParticleType,TreeCell>::allocated_buffer;
  using SphTree<ndim,ParticleType,TreeCell>::box;
  using SphTree<ndim,ParticleType,TreeCell>::cellbuf;
  using SphTree<ndim,ParticleType,TreeCell>::gravity_mac;
  using SphTree<ndim,ParticleType,TreeCell>::kernp;
  using SphTree<ndim,ParticleType,TreeCell>::kernrange;
  using SphTree<ndim,ParticleType,TreeCell>::kernrangesqd;
  using SphTree<ndim,ParticleType,TreeCell>::levelneibbuf;
  using SphTree<ndim,ParticleType,TreeCell>::multipole;
  using SphTree<ndim,ParticleType,TreeCell>::neibcheck;
  using SphTree<ndim,ParticleType,TreeCell>::neibpartbuf;
  using SphTree<ndim,ParticleType,TreeCell>::Ngravcellmaxbuf;
  using SphTree<ndim,ParticleType,TreeCell>::Nleafmax;
  using SphTree<ndim,ParticleType,TreeCell>::Nneibmaxbuf;
  using SphTree<ndim,ParticleType,TreeCell>::Nthreads;
  using SphTree<ndim,ParticleType,TreeCell>::Ntot;
  using SphTree<ndim,ParticleType,TreeCell>::Ntotmax;
  using SphTree<ndim,ParticleType,TreeCell>::Ntotmaxold;
  using SphTree<ndim,ParticleType,TreeCell>::Ntotold;
  using SphTree<ndim,ParticleType,TreeCell>::timing;
  using SphTree<ndim,ParticleType,TreeCell>::tree;
  using SphTree<ndim,ParticleType,TreeCell>::ghosttree;
#ifdef MPI_PARALLEL
  using SphTree<ndim,ParticleType,TreeCell>::mpighosttree;
  using SphTree<ndim,ParticleType,TreeCell>::Nmpi;
  using SphTree<ndim,ParticleType,TreeCell>::prunedtree;
  using SphTree<ndim,ParticleType,TreeCell>::sendprunedtree;
#endif


  //-----------------------------------------------------------------------------------------------
  ForcelessTree(int, int, int, int, FLOAT, FLOAT, FLOAT, string, string,
               DomainBox<ndim> *, SmoothingKernel<ndim> *, CodeTiming *, ParticleTypeRegister&);
  virtual ~ForcelessTree();

  void RaiseError(void){
	  ExceptionHandler::getIstance().raise("Error: update forces called in forceless simulation");
  };

  //-----------------------------------------------------------------------------------------------
  void UpdateAllSphProperties(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphForces(int, int, SphParticle<ndim> *, Sph<ndim> *,
                            Nbody<ndim> *, DomainBox<ndim> &, Ewald<ndim> *);

  void UpdateAllSphHydroForces(int, int, SphParticle<ndim> *, Sph<ndim> *,
                                 Nbody<ndim> *, DomainBox<ndim> &){RaiseError();};
  void UpdateAllSphGravForces(int, int, SphParticle<ndim> *, Sph<ndim> *,
                                Nbody<ndim> *, DomainBox<ndim> &, Ewald<ndim> *){RaiseError();};

  SphParticle<ndim> *part_gen;

};


//=================================================================================================
//  Class SM2012SphTree
/// \brief   Class containing tree for computing Saitoh & Makino (2012) SPH force loops.
/// \details Class containing tree for computing Saitoh & Makino (2012) SPH force loops.
/// \author  D. A. Hubber
/// \date    08/01/2014
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
class SM2012SphTree: public SphTree<ndim,ParticleType,TreeCell>
{
 public:

  using SphTree<ndim,ParticleType,TreeCell>::activelistbuf;
  using SphTree<ndim,ParticleType,TreeCell>::activepartbuf;
  using SphTree<ndim,ParticleType,TreeCell>::allocated_buffer;
  using SphTree<ndim,ParticleType,TreeCell>::box;
  using SphTree<ndim,ParticleType,TreeCell>::cellbuf;
  using SphTree<ndim,ParticleType,TreeCell>::gravity_mac;
  using SphTree<ndim,ParticleType,TreeCell>::kernp;
  using SphTree<ndim,ParticleType,TreeCell>::kernrange;
  using SphTree<ndim,ParticleType,TreeCell>::kernrangesqd;
  using SphTree<ndim,ParticleType,TreeCell>::levelneibbuf;
  using SphTree<ndim,ParticleType,TreeCell>::multipole;
  using SphTree<ndim,ParticleType,TreeCell>::neibcheck;
  using SphTree<ndim,ParticleType,TreeCell>::neibpartbuf;
  using SphTree<ndim,ParticleType,TreeCell>::Ngravcellmaxbuf;
  using SphTree<ndim,ParticleType,TreeCell>::Nleafmax;
  using SphTree<ndim,ParticleType,TreeCell>::Nneibmaxbuf;
  using SphTree<ndim,ParticleType,TreeCell>::Ntot;
  using SphTree<ndim,ParticleType,TreeCell>::Ntotmax;
  using SphTree<ndim,ParticleType,TreeCell>::Ntotmaxold;
  using SphTree<ndim,ParticleType,TreeCell>::Ntotold;
  using SphTree<ndim,ParticleType,TreeCell>::timing;
  using SphTree<ndim,ParticleType,TreeCell>::tree;
  using SphTree<ndim,ParticleType,TreeCell>::ghosttree;
#ifdef MPI_PARALLEL
  using SphTree<ndim,ParticleType,TreeCell>::mpighosttree;
  using SphTree<ndim,ParticleType,TreeCell>::Nmpi;
  using SphTree<ndim,ParticleType,TreeCell>::prunedtree;
  using SphTree<ndim,ParticleType,TreeCell>::sendprunedtree;
#endif


  //-----------------------------------------------------------------------------------------------
  SM2012SphTree(int, int, int, int, FLOAT, FLOAT, FLOAT, string, string,
                DomainBox<ndim> *, SmoothingKernel<ndim> *, CodeTiming *,
                ParticleTypeRegister&);


  //-----------------------------------------------------------------------------------------------
  void UpdateAllSphProperties(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *) {};
  void UpdateAllSphHydroForces(int, int, SphParticle<ndim> *, Sph<ndim> *,
                               Nbody<ndim> *, DomainBox<ndim> &) {};
  void UpdateAllSphForces(int, int, SphParticle<ndim> *, Sph<ndim> *,
                          Nbody<ndim> *, DomainBox<ndim> &, Ewald<ndim> *) {};
  void UpdateAllSphGravForces(int, int, SphParticle<ndim> *, Sph<ndim> *,
                              Nbody<ndim> *, DomainBox<ndim> &, Ewald<ndim> *) {};

};

#endif
