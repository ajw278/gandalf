//=================================================================================================
//  SlopeLimiter.h
//  Contains all routines for ...
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


#ifndef _SLOPE_LIMITER_H_
#define _SLOPE_LIMITER_H_


#include <assert.h>
#include <iostream>
#include <math.h>
#include "Exception.h"
#include "InlineFuncs.h"
#include "Precision.h"
using namespace std;



//=================================================================================================
//  Class SlopeLimiter
/// \brief   Parent class for all slope limiters
/// \details ...
/// \author  D. A. Hubber
/// \date    23/03/2015
//=================================================================================================
template <int ndim, template<int> class ParticleType>
class SlopeLimiter
{
 public:

  SlopeLimiter() {};
  virtual ~SlopeLimiter() {};


  virtual void CellLimiter(ParticleType<ndim> &parti,
                           const ParticleType<ndim>* neibpart, const int* neiblist, int Nneib) {} ;

  virtual void ComputeLimitedSlopes(ParticleType<ndim> &parti, ParticleType<ndim> &partj,
                                    FLOAT draux[ndim], FLOAT gradW[ndim+2][ndim], FLOAT dW[ndim+2]) {
    for (int var=0; var<ndim+2; var++) {
      dW[var] = DotProduct(parti.grad[var], draux, ndim);
      for (int k=0; k<ndim; k++) gradW[var][k] = parti.grad[var][k];
    }
  }

};



//=================================================================================================
//  Class ZeroSlopeLimiter
/// \brief   Imposes 1st-order (Godunov method) by zeroing slopes.
/// \details Imposes 1st-order (Godunov method) by zeroing slopes.
/// \author  D. A. Hubber
/// \date    28/08/2015
//=================================================================================================
template <int ndim, template<int> class ParticleType>
class ZeroSlopeLimiter : public SlopeLimiter<ndim,ParticleType>
{
 public:

  ZeroSlopeLimiter() {};
  ~ZeroSlopeLimiter() {};

  //===============================================================================================
  void ComputeLimitedSlopes(ParticleType<ndim> &parti, ParticleType<ndim> &partj,
                            FLOAT draux[ndim], FLOAT gradW[ndim+2][ndim], FLOAT dW[ndim+2])
  {
    for (int var=0; var<ndim+2; var++) {
      dW[var] = (FLOAT) 0.0;
      for (int k=0; k<ndim; k++) gradW[var][k] = (FLOAT) 0.0;
    }
  }

};



//=================================================================================================
//  Class NullLimiter
/// \brief   Null slope limiter.  Extrapolates variables fully without limiting their values.
/// \details Null slope limiter.  Extrapolates variables fully without limiting their values.
/// \author  D. A. Hubber
/// \date    23/03/2015
//=================================================================================================
template <int ndim, template<int> class ParticleType>
class NullLimiter : public SlopeLimiter<ndim,ParticleType>
{
 public:

  NullLimiter() {};
  ~NullLimiter() {};
};


//=================================================================================================
//  Class TVDScalarLimiter
/// \brief   TVD Scalar limiter.
/// \details The standard generalisation of the minmod or monotized central slope limiters to
///          multiple dimensions. For each neighbour, the slope is reconstructed to either the
///          neighbouring particle (minmod) or to the face (monotized central slopes). The slope is
///          is then limited to ensure the reconstruction lies between the values for the particle
///          and it's neighbour. This repeated for each primitive variable and each neighbour.
///          The TESS slope limiter (Hess & Springel, 2011) is equivalent to this limiter.
/// \author  R. A. Booth
/// \date    14/09/2016
//=================================================================================================
template <int ndim, template<int> class ParticleType>
class TVDScalarLimiter : public SlopeLimiter<ndim,ParticleType>
{
 public:

  TVDScalarLimiter(bool limit_at_edge=true)
    : _edge_limit(limit_at_edge)
  {};

  //===============================================================================================
  void CellLimiter(ParticleType<ndim> &parti,
                   const ParticleType<ndim>* neibpart, const int* neiblist, int Nneib) {
    double dr[ndim] ;
    double alpha[ndim+2] ;
    int j, jj, k, var ;

    for (var=0; var<ndim+2; var++) alpha[var] = 1.0 ;

    for (jj=0; jj<Nneib; jj++) {
      j = neiblist[jj];

      for (k=0; k<ndim; k++) dr[k] = neibpart[j].r[k] - parti.r[k];

      // Calculate min and max values of primitives for slope limiters
      for (var=0; var<ndim+2; var++) {
        // Reconstruct slope to cell edge or neighbouring particle.
        double dW = DotProduct(parti.grad[var], dr, ndim) ;
        if (_edge_limit) dW *= 0.51 ;

        double dWcell = neibpart[j].Wprim[var] - parti.Wprim[var] ;

        // Limit the reconstructive value to be between the two cell values
        if (dW != 0)
          alpha[var] = min(alpha[var], max(0., min(1., dWcell / dW))) ;
      }
    }

    for (var=0; var<ndim+2; var++)
      for (int k=0; k<ndim; k++)
        parti.grad[var][k] *= alpha[var] ;
  }

private:
  bool _edge_limit ;
};

//=================================================================================================
//  Class ScalarLimiter.
/// \brief   Scalar limiter with relaxed constraints. The limiter merely assures that the maximum
///          and minimum reconstructed values lie within the maximum / minimum values of all of the
///          neighbours. Slopes can either be reconstructed to particle locations (minmod style) or
///          or cell faces (mon-cen style, preferred).
///          The Balsara (2004) limiter with psi = 0.5 or 1 (minmod or mon-cen style) and the
///          original Gaburov & Nitadori (2011) limiter are of this form.
/// \details
/// \author  R. A. Booth
/// \date    14/09/2016
//=================================================================================================
template <int ndim, template<int> class ParticleType>
class ScalarLimiter: public SlopeLimiter<ndim,ParticleType>
{
 public:

  ScalarLimiter(bool limit_at_edge=true)
    : _edge_limit(limit_at_edge)
  {};

  //===============================================================================================
  void CellLimiter(ParticleType<ndim> &parti,
                   const ParticleType<ndim>* neibpart, const int* neiblist, int Nneib) {
    double dr[ndim], drmax ;
    double Wmax[ndim+2], Wmin[ndim+2] ;
    int j, jj, k, var ;

    // Get the max / min values over the neighbours
    for (var=0; var<ndim+2; var++) {
      Wmax[var] = Wmin[var] = parti.Wprim[var] ;
    }

    // Compute the maximum particle volume and max/min values in
    // that volume.
    drmax = 0 ;
    for (jj=0; jj<Nneib; jj++) {
      j = neiblist[jj];

      for (k=0; k<ndim; k++) dr[k] = neibpart[j].r[k] - parti.r[k];

      for (var=0; var<ndim+2; var++) {
        Wmax[var] = max(Wmax[var], neibpart[j].Wprim[var]) ;
        Wmin[var] = min(Wmin[var], neibpart[j].Wprim[var]) ;
        drmax = max(drmax, sqrt(DotProduct(dr,dr, ndim))) ;
      }
    }

    // TODO:
    //   Factor of 2 here is a hack assuming m4 kernel
    drmax = max(drmax, 2 * parti.h) ;
    if (_edge_limit) drmax *= 0.51 ;

    for (var=0; var<ndim+2; var++) {
      double gradW = sqrt(DotProduct(parti.grad[var], parti.grad[var], ndim)) ;
      double dW = drmax * gradW ;

      double dWmax = Wmax[var] - parti.Wprim[var] ;
      double dWmin = parti.Wprim[var] - Wmin[var] ;

      double alpha = 1 ;
      if (dW != 0)
        alpha = max(0., min(1., min(dWmax/dW, dWmin/dW))) ;

      for (int k=0; k<ndim; k++)
        parti.grad[var][k] *= alpha ;
    }
  }

private:
  bool _edge_limit ;
};

//=================================================================================================
//  Class Springel2009Limiter.
/// \brief   Limiter used in the AREPO paper. The limiter is similar to the ScalarLimiter, except it
///          only limits the gradients based on values that are actually reconstructed rather than
///          all possible values that could be reconstructed, making it slightly less diffusive.
/// \details
/// \author  R. A. Booth, D. A. Hubber
/// \date    14/09/2016
//=================================================================================================
template <int ndim, template<int> class ParticleType>
class Springel2009Limiter: public SlopeLimiter<ndim,ParticleType>
{
 public:

  Springel2009Limiter(bool limit_at_edge=true)
    : _edge_limit(limit_at_edge)
  {};

  //===============================================================================================
  void CellLimiter(ParticleType<ndim> &parti,
                   const ParticleType<ndim>* neibpart, const int* neiblist, int Nneib) {
    double dr[ndim] ;
    double dWmax[ndim+2], dWmin[ndim+2] ;
    double alpha[ndim+2] ;
    int j, jj, k, var ;

    // Get the max / min values over the neighbours
    for (var=0; var<ndim+2; var++) {
      alpha[var] = 1.0 ;
      dWmax[var] = dWmin[var] = parti.Wprim[var] ;
    }

    // Compute the maximum particle volume and max/min values in
    // that volume.
    for (jj=0; jj<Nneib; jj++) {
      j = neiblist[jj];
      for (var=0; var<ndim+2; var++) {
        dWmax[var] = max(dWmax[var], neibpart[j].Wprim[var]) ;
        dWmin[var] = min(dWmin[var], neibpart[j].Wprim[var]) ;
      }
    }
    for (var=0; var<ndim+2; var++) {
      dWmax[var] -= parti.Wprim[var] ;
      dWmin[var] -= parti.Wprim[var] ;
    }


    for (jj=0; jj<Nneib; jj++) {
      j = neiblist[jj];

      for (k=0; k<ndim; k++) dr[k] = neibpart[j].r[k] - parti.r[k];

      for (var=0; var<ndim+2; var++) {
        double dW = DotProduct(parti.grad[var], dr, ndim) ;

        if (_edge_limit) dW *= 0.51 ;

        if (dW > 0)
          alpha[var] = min(alpha[var], dWmax[var] / dW) ;
        else if (dW < 0)
          alpha[var] = min(alpha[var], dWmin[var] / dW) ;
      }
    }

    for (var=0; var<ndim+2; var++)
     for (k =0; k < ndim; k++)
       parti.grad[var][k] *= alpha[var] ;
  }

private:
  bool _edge_limit ;
};


//=================================================================================================
//  Class GizmoLimiter
/// \brief   Slope limiter from the GIZMO paper.
/// \details This limiter consists of first applying the non-TVD ScalarLimiter to limit the
///          gradient for a given particle. Subsequently a per-face limiter is applied that ensures
///          each of the 1D Riemann problems are (approximately) TVD, even if the full problem is
//           not.
/// \author  D. A. Hubber
/// \date    23/03/2015
//=================================================================================================
template <int ndim, template<int> class ParticleType>
class GizmoLimiter : public ScalarLimiter<ndim,ParticleType>
{
 public:

  GizmoLimiter() {};
  ~GizmoLimiter() {};

  //===============================================================================================
  void ComputeLimitedSlopes(ParticleType<ndim> &parti, ParticleType<ndim> &partj,
                            FLOAT draux[ndim], FLOAT gradW[ndim+2][ndim], FLOAT dW[ndim+2])
  {
    int var;
    FLOAT alpha = (FLOAT) 0.0;
    FLOAT dr[ndim];
    FLOAT phiminus;
    FLOAT phiplus;
    FLOAT phimid;
    const FLOAT psi1 = (FLOAT) 0.5;
    const FLOAT psi2 = (FLOAT) 0.375;

    //---------------------------------------------------------------------------------------------
    for (var=0; var<ndim+2; var++) {

      dW[var] = DotProduct(parti.grad[var], draux, ndim) ;
      for (int k=0; k<ndim; k++) gradW[var][k] = parti.grad[var][k];


      for (int k=0; k<ndim; k++) dr[k] = partj.r[k] - parti.r[k];
      FLOAT drmag = sqrt(DotProduct(dr, dr, ndim));

      const FLOAT delta1 = psi1*fabs(parti.Wprim[var] - partj.Wprim[var]);
      const FLOAT delta2 = psi2*fabs(parti.Wprim[var] - partj.Wprim[var]);
      const FLOAT phimin = min(parti.Wprim[var], partj.Wprim[var]);
      const FLOAT phimax = max(parti.Wprim[var], partj.Wprim[var]);
      const FLOAT phibar = parti.Wprim[var] + (partj.Wprim[var] - parti.Wprim[var])*
        sqrt(DotProduct(draux, draux, ndim))/drmag;
      const FLOAT phimid0 = parti.Wprim[var] + dW[var];

      if (sgn(phimin - delta1) == sgn(phimin)) {
        phiminus = phimin - delta1;
      }
      else {
        phiminus = phimin / ((FLOAT) 1.0 + delta1/fabs(phimin));
      }

      if (sgn(phimax + delta1) == sgn(phimax)) {
        phiplus = phimax + delta1;
      }
      else {
        phiplus = phimax / ((FLOAT) 1.0 + delta1/fabs(phimax));
      }

      if (parti.Wprim[var] < partj.Wprim[var]) {
        phimid = max(phiminus, min(phibar + delta2, phimid0));
      }
      else if (parti.Wprim[var] > partj.Wprim[var]) {
        phimid = min(phiplus, max(phibar - delta2, phimid0));
      }
      else {
        phimid = parti.Wprim[var];
      }

      FLOAT drsqd = DotProduct(draux, draux, ndim);
      dW[var] = phimid - parti.Wprim[var];
    }
    //---------------------------------------------------------------------------------------------
  }

};

#endif
