//=================================================================================================
//  Feedback.cpp
//  Contains all functions for managing the tree for hydrodynamical particles.
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
#include <string>
#include "Feedback.h"
#if defined _OPENMP
#include <omp.h>
#endif
using namespace std;



//=============================================================================
//  HotWindFeedback::HotWindFeedback
/// ...
//=============================================================================
template <int ndim, template<int> class ParticleType>
HotWindFeedback<ndim,ParticleType>::HotWindFeedback()
{
}



//=============================================================================
//  HotWindFeedback::AddNewParticles
/// ...
//=============================================================================
template <int ndim, template<int> class ParticleType>
int HotWindFeedback<ndim,ParticleType>::AddNewParticles()
{
  return 0;
}



template class HotWindFeedback<1, GradhSphParticle>;
template class HotWindFeedback<2, GradhSphParticle>;
template class HotWindFeedback<3, GradhSphParticle>;

template class HotWindFeedback<1, MeshlessFVParticle>;
template class HotWindFeedback<2, MeshlessFVParticle>;
template class HotWindFeedback<3, MeshlessFVParticle>;
