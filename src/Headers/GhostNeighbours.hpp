//=================================================================================================
//  GhostNeighbours.hpp
//  Contains the definitions of a class used for constructing periodic and mirror ghosts
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

#ifndef _GHOST_NEIGHBOURS_H_
#define _GHOST_NEIGHBOURS_H_

#include "DomainBox.h"
#include "Precision.h"
#include "Particle.h"
#include "Tree.h"

//=================================================================================================
/// \brief  Make the ghost particles based upon the boundary conditions
/// \author R. A. Booth
/// \date   26/10/2015
/// \return The number of neighbours found
//=================================================================================================
template<int ndim>
class GhostNeighbourFinder
{
private:

  const DomainBox<ndim>& _domain ;
  Box<ndim> _cell ;
  FLOAT _centre[ndim] ;
  bool _mirror_bound[ndim][2] ;
  bool _need_mirror[ndim][2] ;
  bool _periodic_bound[ndim] ;
  bool _any_periodic, _any_mirror, _need_mirrors, _any_special;

  void 	SetBoundaryFlags() {
    for (int k=0; k <ndim; k++) {
      _periodic_bound[k] = _domain.boundary_lhs[k] == periodicBoundary ;
      if (_periodic_bound[k])
    	  assert(_domain.boundary_rhs[k] == periodicBoundary) ;

      _mirror_bound[k][0] = _domain.boundary_lhs[k] == mirrorBoundary ;
      _mirror_bound[k][1] = _domain.boundary_rhs[k] == mirrorBoundary ;

      _any_periodic |= _periodic_bound[k] ;
      _any_mirror   |= (_mirror_bound[k][0] | _mirror_bound[k][1]) ;

      _need_mirror[k][0] = _need_mirror[k][1] = false ;
    }
  }

public:

	GhostNeighbourFinder(const DomainBox<ndim>& simbox)
	: _domain(simbox),
	  _any_periodic(false), _any_mirror(false), _need_mirrors(false), _any_special(false)
	{
	  SetBoundaryFlags() ;
	}

	GhostNeighbourFinder(const DomainBox<ndim>& simbox, const TreeCellBase<ndim>& cell)
	: _domain(simbox),
	  _any_periodic(false), _any_mirror(false), _need_mirrors(false), _any_special(false)
	{
	  SetBoundaryFlags() ;
	  SetTargetCell(cell) ;
	}


	//=================================================================================================
	/// \brief  Set a target cell to compute mirrors for.
	/// \author R. A. Booth
	/// \date   27/10/2015
	/// \return A boolean saying whether the boxes overlap
	//=================================================================================================
	void SetTargetCell(const TreeCellBase<ndim>& cell)
	{
	  _need_mirrors = false ;
	  for (int k=0; k < ndim; k++){
		_centre[k] = cell.rcell[k] ;
		_cell.boxmin[k] = cell.hboxmin[k] ;
		_cell.boxmax[k] = cell.hboxmax[k] ;

		if (_any_mirror){
		  // Compute whether we need a mirror in the gather case
		  _need_mirror[k][0] = (_mirror_bound[k][0] & (_cell.boxmin[k] < _domain.boxmin[k]));
		  _need_mirror[k][1] = (_mirror_bound[k][1] & (_cell.boxmax[k] > _domain.boxmax[k]));
		  _need_mirrors     |= (_need_mirror[k][0] | _need_mirror[k][1]) ;
	    }
	  }
	  _any_special = _need_mirrors | _any_periodic ;
	}

	//=================================================================================================
	/// \brief  Find the smallest separation vector dr[] in a periodic box
	/// \author D. A. Hubber, G. Rosotti
	/// \date   12/11/2013
	/// \return A boolean saying whether the boxes overlap
	//=================================================================================================
	type_flag NearestPeriodicVector(FLOAT dr[ndim]) const
	{
	  type_flag bound ;
   	  if (_any_periodic)
		{
   		  for (int k=0; k<ndim; k++) {
   		    if (_periodic_bound[k]) {
   		      if (dr[k] > _domain.boxhalf[k]) {
   		        dr[k] -=_domain.boxsize[k];
   		        bound.set_flag(periodic_bound_flags[k][1]) ;
   		      }
   		      else if (dr[k] < -_domain.boxhalf[k]) {
   		        dr[k] += _domain.boxsize[k];
   		        bound.set_flag(periodic_bound_flags[k][0]) ;
   		      }
   		    }
   		  }
		}
   	  return bound ;
	}

	//=================================================================================================
	/// \brief  Find 'periodic' correction vector.
	/// \author D. A. Hubber, G. Rosotti
	/// \date   12/11/2013
	/// \return A boolean saying whether the boxes overlap
	//=================================================================================================
	type_flag PeriodicDistanceCorrection(const FLOAT dr[ndim], FLOAT dr_corr[ndim]) const
	{
	  type_flag bound ;
	  if (_any_periodic)
		{
	   	  for (int k=0; k<ndim; k++) {
	   	    if (_periodic_bound[k]) {
	   	      if (dr[k] > _domain.boxhalf[k]) {
	   	        dr_corr[k] =- _domain.boxsize[k];
	   	        bound.set_flag(periodic_bound_flags[k][1]) ;
	   	      }
	   	      else if (dr[k] < -_domain.boxhalf[k]) {
	   	        dr_corr[k] = _domain.boxsize[k];
	   	        bound.set_flag(periodic_bound_flags[k][0]) ;
	   	      }
	   	      else
	   	    	dr_corr[k] = 0 ;
	   	    }
	   	  }
		}
	  return bound ;
	}

	//=================================================================================================
	/// \brief  Find the maximum number of mirrors required
	/// \author R. A. Booth
	/// \date   28/10/2015
	/// \return An integer specifying the max number of neighbours (>=1)
	//=================================================================================================
	int MaxNumGhosts() const
	{
	  int NumGhosts = 1;
		if (_any_mirror){
		  for (int k=0; k < ndim; k++){
			NumGhosts *= 1 + _mirror_bound[k][0] + _mirror_bound[k][1] ;
		}
	  }
	  return NumGhosts ;
	}

	//=================================================================================================
	/// \brief Construct the centres and reflection signs of a cell. This list will include the
	///        original cell or the periodic neighbour of the original cell.
	/// \author R. A. Booth
	/// \date   27/10/2015
	/// \return The number of neighbours found
	//===============================================================================================
	template<template <int> class ParticleType>
	int ConstructGhostsScatterGather(const ParticleType<ndim>& p, ParticleType<ndim>* ngbs) const
	{
	  // First find the nearest periodic mirror
	  ngbs[0] = p ;
	  if (_any_periodic)
		_MakePeriodicGhost(ngbs[0]) ;


	  // Number of Ghost cells
	  int Nghost = 1 ;

	  if (_any_mirror)
		Nghost = _MakeReflectedScatterGatherGhosts(ngbs) ;

	  return Nghost ;
	}

    //=================================================================================================
	/// \brief Construct the centres and reflection signs of a cell. This list will include the
	///        original cell or the periodic neighbour of the original cell.
	/// \author R. A. Booth
	/// \date   27/10/2015
	/// \return The number of neighbours found
	//===============================================================================================
	template<template <int> class ParticleType>
	int ConstructGhostsGather(const ParticleType<ndim>& p, ParticleType<ndim>* ngbs) const
	{
	  // First find the nearest periodic mirror
	  ngbs[0] = p ;
	  if (_any_periodic)
		_MakePeriodicGhost(ngbs[0]) ;


	  // Number of Ghost cells
	  int Nghost = 1 ;
	  // Now recursively reflect the cells
	  if (_need_mirrors)
		Nghost = _MakeReflectedGhostsGather(ngbs) ;

	  return Nghost ;
	}

	template< template<int> class TreeCell, template<int> class ParticleType>
	void CorrectGhostParticlePosition(const TreeCell<ndim>& cell, const FLOAT * r, const int * sign,
			                          ParticleType<ndim>& p) const
	{
		if (_any_special){
		  for (int k=0; k < ndim; k++){
		    FLOAT dr_k = p.r[k] - cell.rcell[k] ;
		    p.r[k] = r[k] + dr_k * sign[k] ;
		    p.v[k] *= sign[k] ;
		  }
		}
	}

private:
    //=================================================================================================
	//  _MakePeriodicGhost
	/// \brief Do the actual construction of the nearest periodic ghost.
	/// \author R. A. Booth
	/// \date   27/10/2015
	/// \return The number of neighbours found
	//===============================================================================================
	template<template <int> class ParticleType>
	void _MakePeriodicGhost(ParticleType<ndim>& p) const {
	  FLOAT dr[ndim] ;

	  for (int k=0; k <ndim; k++)
	    dr[k] = p.r[k] - _centre[k] ;

	  type_flag bound_flag =  NearestPeriodicVector(dr) ;
	  for (int k=0; k <ndim; k++)
	    p.r[k] = _centre[k] + dr[k];

	  p.flags.set_flag(bound_flag.get()) ;

	}

    //=================================================================================================
	//  _MakeReflectedGhostsGather
	/// \brief Do the actual construction of the mirror ghosts within the gather range. Assumes the
	/// first particle is already saved in ngbs.
	/// \author R. A. Booth
	/// \date   27/10/2015
	/// \return The number of neighbours found
	//===============================================================================================
	template<template <int> class ParticleType>
	int _MakeReflectedGhostsGather(ParticleType<ndim>* ngbs) const {
	  int nc = 1 ;
	  // Loop over the possible directions for reflections
	  for (int k = 0; k < ndim; k++){
		// Save the current number of images
		int Nghost = nc ;
		// Do reflections on the left edge
		if (_need_mirror[k][0]){
		  for (int n=0; n < Nghost; n++){
			double rk = 2*_domain.boxmin[k] - ngbs[n].r[k] ;
			if (rk > _cell.boxmin[k]) {
			  ngbs[nc] = ngbs[n] ;
			  ngbs[nc].reflect(k, _domain.boxmin[k]) ;
			  ngbs[nc].flags.set_flag(mirror_bound_flags[k][0]) ;
			  nc++ ;
			}
		  }
		}
		// Do reflections on the right edge
		if (_need_mirror[k][1]){
		  for (int n=0; n < Nghost; n++){
			double rk = 2*_domain.boxmax[k] - ngbs[n].r[k] ;
			if (rk < _cell.boxmax[k]) {
			  ngbs[nc] = ngbs[n] ;
			  ngbs[nc].reflect(k, _domain.boxmax[k]) ;
			  ngbs[nc].flags.set_flag(mirror_bound_flags[k][1]) ;
			  nc++ ;
			}
		  }
		}
	  }
	  return nc ;
	}

	//=================================================================================================
	//  _MakeReflectedGhostsScatterGather
	/// \brief Do the actual construction of the mirror ghosts within the scatter-gather range. Assumes
	/// the  first particle is already saved in ngbs.
	/// \author R. A. Booth
	/// \date   27/10/2015
	/// \return The number of neighbours found
	//===============================================================================================
	template<template <int> class ParticleType>
	int _MakeReflectedScatterGatherGhosts(ParticleType<ndim>* ngbs) const {
	  int nc = 1 ;
	  double h2 = ngbs[0].hrangesqd ;
	  // Loop over the possible directions for reflections
	  for (int k = 0; k < ndim; k++){
		// Save the current number of images
		int Nghost = nc ;

		// Do reflections on the left edge
		if (_mirror_bound[k][0]){
		  FLOAT dx = 2*_domain.boxmin[k] - ngbs[0].r[k] - _cell.boxmin[k] ;
		  if (dx*dx < h2){
			for (int n=0; n < Nghost; n++){
			  ngbs[nc] = ngbs[n] ;
			  ngbs[nc].reflect(k, _domain.boxmin[k]) ;
			  ngbs[nc].flags.set_flag(mirror_bound_flags[k][0]) ;
			  nc++;
			}
		  }
		}
		// Do reflections on the right edge
		if (_mirror_bound[k][1]){
		  FLOAT dx = 2*_domain.boxmax[k] - ngbs[0].r[k] - _cell.boxmax[k] ;
		  if (dx*dx < h2){
			for (int n=0; n < Nghost; n++){
			 ngbs[nc] = ngbs[n] ;
			 ngbs[nc].reflect(k, _domain.boxmax[k]) ;
			 ngbs[nc].flags.set_flag(mirror_bound_flags[k][1]) ;
			 nc++ ;
			}
		  }
		}
	  }
	  return nc ;
	}
};


//=================================================================================================
//  ParticleCellProxy
/// Proxy class that we can use to create the Ghosts for the brute force search
//=================================================================================================
template<int ndim>
struct ParticleCellProxy
: public TreeCellBase<ndim>
{
	using TreeCellBase<ndim>::hboxmin ;
	using TreeCellBase<ndim>::hboxmax ;
	using TreeCellBase<ndim>::rcell ;


	ParticleCellProxy(const Particle<ndim>& p)
	{
		FLOAT dr =  sqrt(p.hrangesqd) ;
		for (int k=0; k<ndim;k++)
		{
			rcell[k] = p.r[k] ;

			hboxmax[k] = rcell[k] + dr ;
			hboxmin[k] = rcell[k] - dr ;

		}
	}
} ;


//=================================================================================================
//  DomainCellProxy
/// Proxy class that we can use to create the Ghosts for the brute force search
//=================================================================================================
template<int ndim>
struct DomainCellProxy
: public TreeCellBase<ndim>
{
	using TreeCellBase<ndim>::hboxmin ;
	using TreeCellBase<ndim>::hboxmax ;
	using TreeCellBase<ndim>::rcell ;

	DomainCellProxy(const DomainBox<ndim>& box)
	{
		for (int k=0; k<ndim;k++)
		{
			rcell[k] = 0.5 * (box.boxmax[k] + box.boxmin[k]) ;
			FLOAT dr = 0.6 * (box.boxmax[k] - box.boxmin[k]) ;

			hboxmax[k] = rcell[k] + dr ;
			hboxmax[k] = rcell[k] - dr ;
		}
	}
} ;


#endif//_GHOST_NEIGHBOURS_H_
