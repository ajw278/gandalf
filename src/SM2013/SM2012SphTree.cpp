//=================================================================================================
//  SM2012SphTree.cpp
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
//  SM2012SphTree::SM2012SphTree
/// SM2012SphTree constructor.  Initialises various variables.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
SM2012SphTree<ndim,ParticleType,TreeCell>::SM2012SphTree
 (int _Nleafmax, int _Nmpi, int _pruning_level_min, int _pruning_level_max, FLOAT _thetamaxsqd,
  FLOAT _kernrange, FLOAT _macerror, string _gravity_mac, string _multipole,
  DomainBox<ndim>* _box, SmoothingKernel<ndim>* _kern, CodeTiming* _timing, ParticleTypeRegister& types):
 NeighbourSearch<ndim>(_kernrange, _box, _kern, _timing),
 SphTree<ndim,ParticleType,TreeCell>
  (_Nleafmax, _Nmpi, _pruning_level_min, _pruning_level_max, _thetamaxsqd,
   _kernrange, _macerror, _gravity_mac, _multipole, _box, _kern, _timing)
{
  // Set-up main tree object
  tree = new_tree<ndim,ParticleType,TreeCell>(_Nleafmax, _thetamaxsqd, _kernrange,
   		  	  	  	  	  	  	  	  	  	  _macerror, _gravity_mac, _multipole, *_box, types);

  // Set-up ghost-particle tree object
  ghosttree = new_tree<ndim,ParticleType,TreeCell>(_Nleafmax, _thetamaxsqd, _kernrange,
		  _macerror, _gravity_mac, _multipole, *_box, types);

#ifdef MPI_PARALLEL
  // Set-up ghost-particle tree object
  mpighosttree = new_tree<ndim,ParticleType,TreeCell>(_Nleafmax, _thetamaxsqd, _kernrange,
		  _macerror, _gravity_mac, _multipole, *_box, types);

  // Set-up multiple pruned trees, one for each MPI process
  prunedtree = new Tree<ndim,ParticleType,TreeCell>*[Nmpi] ;
  // new_tree_array<ndim,ParticleType,TreeCell>(Nmpi);
  sendprunedtree =  new Tree<ndim,ParticleType,TreeCell>*[Nmpi] ;
  //new_tree_array<ndim,ParticleType,TreeCell>(Nmpi);

  for (int i=0; i<Nmpi; i++) {
	prunedtree[i] = new_tree<ndim,ParticleType,TreeCell>(_Nleafmax, _thetamaxsqd, _kernrange,
   	  	  	  	  	  	  	  	  	  	  	  	  	  	 _macerror, _gravity_mac, _multipole, *_box, types);
  }
  for (int i=0; i<Nmpi; i++) {
	sendprunedtree[i] = new_tree<ndim,ParticleType,TreeCell>(_Nleafmax, _thetamaxsqd, _kernrange,
       														 _macerror, _gravity_mac, _multipole, *_box, types);
  }
#endif

}



template class SM2012SphTree<1,SM2012SphParticle,BruteForceTreeCell>;
template class SM2012SphTree<2,SM2012SphParticle,BruteForceTreeCell>;
template class SM2012SphTree<3,SM2012SphParticle,BruteForceTreeCell>;
template class SM2012SphTree<1,SM2012SphParticle,KDTreeCell>;
template class SM2012SphTree<2,SM2012SphParticle,KDTreeCell>;
template class SM2012SphTree<3,SM2012SphParticle,KDTreeCell>;
template class SM2012SphTree<1,SM2012SphParticle,OctTreeCell>;
template class SM2012SphTree<2,SM2012SphParticle,OctTreeCell>;
template class SM2012SphTree<3,SM2012SphParticle,OctTreeCell>;

