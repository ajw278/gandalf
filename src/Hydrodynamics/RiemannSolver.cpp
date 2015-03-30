//=================================================================================================
//  RiemannSolver.cpp
//  Contains all available Riemann solver functions.
//=================================================================================================


#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <map>
#include <string>
#include "Constants.h"
#include "RiemannSolver.h"
using namespace std;


static const int Niterationmax = 100;            // ..
static const FLOAT tolerance = 1.0e-6;           // Iteration tolerance


//=================================================================================================
//  RiemannSolver::RiemannSolver
/// Constructor for RiemannSolver class.
//=================================================================================================
template <int ndim>
RiemannSolver<ndim>::RiemannSolver(FLOAT gamma_aux):
  gamma(gamma_aux),
  invgamma(1.0/gamma_aux),
  g1(0.5*(gamma_aux - 1.0)/gamma_aux),
  g2(0.5*(gamma_aux + 1.0)/gamma_aux),
  g3(2.0*gamma_aux/(gamma_aux - 1.0)),
  g4(2.0/(gamma_aux - 1.0)),
  g5(2.0/(gamma_aux + 1.0)),
  g6((gamma_aux - 1.0)/(gamma_aux + 1.0)),
  g7(0.5*(gamma_aux - 1.0)),
  g8(gamma_aux - 1.0),
  g9(1.0/(gamma_aux - 1.0))
{
}



//=================================================================================================
//  RiemannSolver::ComputeRotationMatrices
/// ...
//=================================================================================================
template <int ndim>
void RiemannSolver<ndim>::ComputeRotationMatrices
 (const FLOAT runit[ndim],
  FLOAT rotMat[ndim][ndim],
  FLOAT invRotMat[ndim][ndim])
{
  if (ndim == 1) {
    rotMat[0][0]    = runit[0];
    invRotMat[0][0] = runit[0];
  }
  else if (ndim == 2) {
    FLOAT theta = atan2(runit[1], runit[0]);
    rotMat[0][0] = cos(theta);
    rotMat[0][1] = -sin(theta);
    rotMat[1][0] = sin(theta);
    rotMat[1][1] = cos(theta);
    invRotMat[0][0] = cos(theta);
    invRotMat[0][1] = sin(theta);
    invRotMat[1][0] = -sin(theta);
    invRotMat[1][1] = cos(theta);
  }

  return;
}



//=================================================================================================
//  RiemannSolver::RotateVector
/// ...
//=================================================================================================
template <int ndim>
void RiemannSolver<ndim>::RotateVector
 (FLOAT rotMat[ndim][ndim],
  FLOAT vec[ndim])
{
  if (ndim == 1) {
    //cout << "Before : " << vec[0] << "    " << rotMat[0][0] << endl;
    vec[0] = rotMat[0][0]*vec[0];
    //cout << "After  : " << vec[0] << endl;
  }
  else if (ndim == 2) {
    FLOAT oldVec[ndim];
    for (int k=0; k<ndim; k++) oldVec[k] = vec[k];
    vec[0] = rotMat[0][0]*oldVec[0] + rotMat[0][1]*oldVec[1];
    vec[1] = rotMat[1][0]*oldVec[0] + rotMat[1][1]*oldVec[1];
  }

  return;
}



//=================================================================================================
//  ExactRiemannSolver::SolveRiemannProblem
/// Exact Riemann solver, based on approach outlined by Toro (1999).
//=================================================================================================
template <int ndim>
void ExactRiemannSolver<ndim>::ComputeStarRegion
 (const FLOAT pl,                      ///< LHS pressure
  const FLOAT pr,                      ///< RHS pressure
  const FLOAT dl,                      ///< LHS density
  const FLOAT dr,                      ///< RHS density
  const FLOAT cl,                      ///< LHS sound speed
  const FLOAT cr,                      ///< RHS sound speed
  const FLOAT ul,                      ///< LHS velocity
  const FLOAT ur,                      ///< RHS velocity
  FLOAT &pstar,                        ///< Intermediate pressure state
  FLOAT &ustar)                        ///< Velocity of intermediate state
{
  const FLOAT Al = g5/dl;
  const FLOAT Ar = g5/dr;
  const FLOAT Bl = pl*g6;
  const FLOAT Br = pr*g6;

  int iteration = 0;                   // No. of iterations (should be configured in parameter file)
  FLOAT cup;                           // Primitive variable ???
  FLOAT fl;                            // Left-state variable
  FLOAT flprime;                       // Left-state variable gradient
  FLOAT fr;                            // Right-state variable
  FLOAT frprime;                       // Right-state variable gradient
  FLOAT gel;                           // Two-shock variable
  FLOAT ger;                           // Two-shock variable
  FLOAT quser;                         // First guess of solution
  FLOAT pmin;                          // Max. pressure
  FLOAT pmax;                          // Min. pressure
  FLOAT pold;                          // Old intermediate pressure
  FLOAT ppv;                           // Primitive Variable pressure
  // Two-rarefaction RS variables
  FLOAT pq;                            // ..
  FLOAT um;                            // ..
  FLOAT ptl;                           // ..
  FLOAT ptr;                           // ..

  //std::cout << "All values : press : " << pl << "   " << pr << "     vel : " << ul << "    " << ur
  //<< "    press : " << pl << "    " << pr << std::endl;

  //// Iteration variables

  // Make guess of solution for first iteration
  quser = 2.0;

  // Compute guess from Primitive Variable RS
  cup  = 0.25*(dl + dr)*(cl + cr);
  ppv  = 0.5*(pl + pr) + 0.5*(ul - ur)*cup;
  ppv  = max(0.0,ppv);
  pmin = min(pl,pr);
  pmax = max(pl,pr);

  // Select guess from PVRS
  if (pmax/pmin <= quser && pmin <= ppv && ppv <= pmax) {
    pstar = ppv;
  }
  // Select two-rarefaction RS
  else if (ppv < pmin) {
    pq    = pow(pl/pr,g1);
    um    = (pq*ul/cl + ur/cr + g4*(pq - 1.0))/(pq/cl + 1.0/cr);
    ptl   = 1.0 + g7*(ul - um)/cl;
    ptr   = 1.0 + g7*(um - ur)/cr;
    pstar = 0.5*(pl*powf(ptl,g3) + powf(pr*ptr,g3));
  }
  // Select two-shock RS with PVRS as estimate
  else {
    gel   = sqrt((g5/dl)/(g6*pl + ppv));
    ger   = sqrt((g5/dr)/(g6*pr + ppv));
    pstar = (gel*pl + ger*pr - (ur - ul))/(gel + ger);
  }


  // Main iteration loop
  //-----------------------------------------------------------------------------------------------
  do {
    iteration++;

    // Calculate contribution to f and fprime for LHS
    if (pstar > pl) {
      fl = (pstar - pl)*sqrt(Al/(pstar + Bl));
      flprime = sqrt(Al/(pstar + Bl))*(1.0 - 0.5*(pstar - pl)/(pstar + Bl));
    }
    else {
      fl = g4*cl*(pow(pstar/pl,g1) - 1.0);
      flprime = pow(pstar/pl,-g2)/(dl*cl);
    }

    // Calculate contribution to f and fprime for RHS
    if (pstar > pr) {
      fr = (pstar - pr)*sqrt(Ar/(pstar + Br));
      frprime = sqrt(Ar/(pstar + Br))*(1.0 - 0.5*(pstar - pr)/(pstar + Br));
    }
    else {
      fr = g4*cr*(pow(pstar/pr,g1) - 1.0);
      frprime = pow(pstar/pr,-g2)/(dr*cr);
    }

    // Perform Newton-Raphson iteration
    pold = pstar;
    pstar = pstar - (fl + fr + ur - ul)/(flprime + frprime);

    // Check if convergence has been achieved
    if (pstar < small_number) pstar = small_number; //this could also be smaller than zero
    else if (2.0*fabs(pstar - pold)/(pstar + pold) < tolerance) break;

    // Check the star variables have not become NaNs
    if (pstar != pstar || ustar != ustar) {
      std::cout << "Checking Riemann values : " << pstar << "   "
                << ustar << "   iteration : " << iteration << std::endl;
      std::cout << "rho : " << dl << "   " << dr << "     vel : " << ul << "    "
                << ur << "    press : " << pl << "    " << pr << std::endl;
      exit(0);
    }


  } while (iteration < Niterationmax);
  //-----------------------------------------------------------------------------------------------

  // Compute velocity of star region
  ustar = 0.5*(ul + ur) + 0.5*(fr - fl);

  return;
}



//=================================================================================================
//  ExactRiemannSolver::SampleExactSolution
//=================================================================================================
template <int ndim>
void ExactRiemannSolver<ndim>::SampleExactSolution
 (const FLOAT pstar,                   ///< ..
  const FLOAT ustar,                   ///< ..
  const FLOAT s,                       ///< ..
  const FLOAT pl,                      ///< ..
  const FLOAT pr,                      ///< ..
  const FLOAT dl,                      ///< ..
  const FLOAT dr,                      ///< ..
  const FLOAT cl,                      ///< ..
  const FLOAT cr,                      ///< ..
  const FLOAT ul,                      ///< ..
  const FLOAT ur,                      ///< ..
  FLOAT &p,                            ///< ..
  FLOAT &d,                            ///< ..
  FLOAT &u)                            ///< ..
{

  FLOAT cm;                            // ..
  FLOAT sl,sr;                         // ..
  FLOAT sh, st;                        // Rarefaction head and tail speed
  FLOAT c;                             // Parameter to determine fan position


  // If sampling point lies to the left of the contact discontinuity
  //-----------------------------------------------------------------------------------------------
  if (s <= ustar) {

    // Left rarefaction
    //---------------------------------------------------------------------------------------------
    if (pstar <= pl) {
      sh = ul - cl;

      // Sampled point is left data state
      if (s <= sh) {
        d = dl;
        u = ul;
        p = pl;
      }
      else {
        cm = cl*powf(pstar/pl,g1);
        st = ustar - cm;

        // Sample point is star left state
        if (s > st) {
          d = dl*powf(pstar/pl,invgamma);
          u = ustar;
          p = pstar;
        }
        // Sampled point is inside left fan
        else {
          u = g5*(cl + g7*ul + s);
          c = g5*(cl + g7*(ul - s));
          d = dl*powf(c/cl,g4);
          p = pl*powf(c/cl,g3);
        }

      }

    }

    // left shock
    //---------------------------------------------------------------------------------------------
    else {
      FLOAT pml = pstar/pl;
      sl = ul - cl*sqrt(g2*pml + g1);

      // Sampled point is left data state
      if (s <= sl) {
        d = dl;
        u = ul;
        p = pl;
      }
      // Sampled point is star left state
      else {
        d = dl*(pml + g6)/(pml*g6 + 1.0);
        u = ustar;
        p = pstar;
      }
    }
    //---------------------------------------------------------------------------------------------

  }
  // Sampling point lies to the right of the contact discontinuity
  //-----------------------------------------------------------------------------------------------
  else {

    // Right shock
    //---------------------------------------------------------------------------------------------
    if (pstar >= pr) {

      FLOAT pmr = pstar/pr;
      sr = ur + cr*sqrt(g2*pmr + g1);

      // Sampled point is right data state
      if (s >= sr) {
        d = dr;
        u = ur;
        p = pr;
      }
      // Sampled point is star right state
      else {
        d = dr*(pmr + g6)/(pmr*g6 + 1.0);
        u = ustar;
        p = pstar;
      }

    }
    // Right rarefaction
    //---------------------------------------------------------------------------------------------
    else {

      sh = ur + cr;

      // Sampled point is left data state
      if (s >= sh) {
        d = dr;
        u = ur;
        p = pr;
      }
      else {
        cm = cr*powf(pstar/pr,g1);
        st = ustar + cm;

        // Sample point is star right state
        if (s <= st) {
          d = dr*powf(pstar/pr,invgamma);
          u = ustar;
          p = pstar;
        }
        // Sampled point is inside left fan
        else {
          u = g5*(-cr + g7*ur + s);
          c = g5*(cr - g7*(ur - s));
          d = dr*powf(c/cr,g4);
          p = pr*powf(c/cr,g3);
        }

      }

    }

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  ExactRiemannSolver::ComputeFluxes
/// Exact Riemann solver, based on approach outlined by Toro (1999).
//=================================================================================================
template <int ndim>
void ExactRiemannSolver<ndim>::ComputeFluxes
 (const FLOAT Wleft[nvar],             ///< LHS primitive state
  const FLOAT Wright[nvar],            ///< RHS primitive state
  const FLOAT runit[ndim],             ///< ..
  FLOAT vface[ndim],                   ///< ..
  FLOAT flux[nvar][ndim])              ///< flux vector
{
  const FLOAT cl = sqrt(gamma*Wleft[ipress]/Wleft[irho]);    // ..
  const FLOAT cr = sqrt(gamma*Wright[ipress]/Wright[irho]);  // ..

  int k,kv;
  FLOAT pstar;                         // Pressure in star region
  FLOAT ustar;                         // Velocity in star region
  FLOAT etot;                          // Total specific energy
  FLOAT p,d,u;                         // Primitive variables at s=0 from Riemann solver
  FLOAT rotMat[ndim][ndim];
  FLOAT invRotMat[ndim][ndim];
  FLOAT uleft[ndim];
  FLOAT uright[ndim];
  FLOAT Wface[ndim];

  // Compute rotation matrices
  this->ComputeRotationMatrices(runit, rotMat, invRotMat);

  for (k=0; k<ndim; k++) uleft[k] = Wleft[k];
  for (k=0; k<ndim; k++) uright[k] = Wright[k];

  this->RotateVector(invRotMat, uleft);
  this->RotateVector(invRotMat, uright);

  //cout << "ul : " << Wleft[0] << "   " << uleft[0] << "    ur : " << Wright[0] << "   " << uright[0] << endl;
  //cout << "Rotation matrix : " << rotMat[0][0] << "   " << invRotMat[0][0] << endl;

  //const FLOAT ul = DotProduct(Wleft, runit, ndim);
  //const FLOAT ur = DotProduct(Wright, runit, ndim);

  // Compute p and u values at interface (in star region)
  ComputeStarRegion(Wleft[ipress], Wright[ipress], Wleft[irho], Wright[irho],
                    cl, cr, uleft[0], uright[0], pstar, ustar);
  SampleExactSolution(pstar, ustar, 0.0, Wleft[ipress], Wright[ipress], Wleft[irho],
                      Wright[irho], cl, cr, uleft[0], uright[0], p, d, u);

  for (kv=0; kv<ndim; kv++) Wface[kv] = 0.0;
  Wface[irho]   = d;
  Wface[ivx]    = u;
  Wface[ipress] = p;

  // Rotate velocity to original frame
  this->RotateVector(rotMat, Wface);
  //cout << "ul : " << Wleft[0] << "   " << Wright[0] << "   " << Wface[0] << endl;

  // Compute fluxes
  /*etot           = p/(gamma - 1.0) + 0.5*d*u*u;
  flux[irho][0]  = d*u;
  flux[ivx][0]   = p + d*u*u;
  flux[ietot][0] = u*(etot + p);*/

  // Compute fluxes in moving frame
  FLOAT ekin = 0.0;
  for (kv=0; kv<ndim; kv++) ekin += Wface[kv]*Wface[kv];
  for (k=0; k<ndim; k++) {
    for (kv=0; kv<ndim; kv++) flux[kv][k] = Wface[irho]*Wface[k]*Wface[kv];
    flux[k][k]     = Wface[irho]*Wface[k]*Wface[k] + Wface[ipress];
    flux[irho][k]  = Wface[irho]*(Wface[k]);
    flux[ietot][k] = (Wface[ipress]/(gamma - 1.0) + 0.5*Wface[irho]*ekin)*(Wface[k])
      + Wface[ipress]*Wface[k];
  }

  // Add corrections for transforming back to original lab frame
  /*for (k=0; k<ndim; k++) {
    flux[ietot][k] += (FLOAT) 0.5*DotProduct(vface, vface, ndim)*flux[irho][k] + vface[k]*flux[ivx][k];
    flux[ivx][k] += vface[k]*flux[irho][k];
  }*/

  // Add corrections for transforming back to original lab frame
  for (k=0; k<ndim; k++) {
    flux[ietot][k] += (FLOAT) 0.5*DotProduct(vface, vface, ndim)*flux[irho][k] + DotProduct(vface, flux[k], ndim);
  }
  for (k=0; k<ndim; k++) {
    for (kv=0; kv<ndim; kv++) flux[kv][k] += vface[kv]*flux[irho][k];
  }


  /*FLOAT ekin = 0.0;
  for (kv=0; kv<ndim; kv++) ekin += Wface[kv]*Wface[kv];
  for (k=0; k<ndim; k++) {
    for (kv=0; kv<ndim; kv++) flux[kv][k] = Wface[irho]*(Wface[k]*Wface[kv] - vface[k]*Wface[kv]);
    flux[k][k]     = Wface[irho]*(Wface[k]*Wface[k] - Wface[k]*vface[k]) + Wface[ipress];
    flux[irho][k]  = Wface[irho]*(Wface[k] - vface[k]);
    flux[ietot][k] = (Wface[ipress]/(gamma - 1.0) + 0.5*Wface[irho]*ekin)*(Wface[k] - vface[k])
      + Wface[ipress]*Wface[k];
  }*/

  /*if (ndim == 1) {
    etot           = Wface[ipress]/(gamma - 1.0) + 0.5*Wface[irho]*Wface[ivx]*Wface[ivx];
    flux[irho][0]  = Wface[irho]*Wface[ivx];
    flux[ivx][0]   = Wface[ipress] + Wface[irho]*Wface[ivx]*Wface[ivx];
    flux[ietot][0] = Wface[ivx]*(etot + Wface[ipress]);
  }*/


  /*for (int var=0; var<nvar; var++) {
    this->RotateVector(rotMat, flux[var]);
  }*/


  return;
}



//=================================================================================================
//  HllcRiemannSolver::SolveRiemannProblem
/// HLLc Riemann solver for pstar and ustar.
//=================================================================================================
template <int ndim>
void HllcRiemannSolver<ndim>::ComputeFluxes
 (const FLOAT Wleft[nvar],               ///< LHS primitive state
  const FLOAT Wright[nvar],              ///< RHS primitive state
  const FLOAT runit[ndim],               ///< ..
  FLOAT vface[ndim],                     ///< ..
  FLOAT flux[nvar][ndim])                ///< flux vector
{
  //// State variables
  const FLOAT dl = Wleft[irho];    //qleft[0];
  const FLOAT dr = Wright[irho];   //qright[0];
  const FLOAT ul = DotProduct(Wleft, runit, ndim);
  const FLOAT ur = DotProduct(Wright, runit, ndim);
  //const FLOAT ul = Wleft[k];       //qleft[1];
  //const FLOAT ur = Wright[k];      //qright[1];
  const FLOAT pl = Wleft[ipress];  //qleft[2];
  const FLOAT pr = Wright[ipress]; //qright[2];
  const FLOAT cl = sqrt(gamma*pl/dl);
  const FLOAT cr = sqrt(gamma*pr/dr);

  FLOAT ekin;                            // Kinetic energy
  FLOAT etotl;                           // Left total energy
  FLOAT etotr;                           // Right total energy
  FLOAT Sl;                              // LHS wave speed estimator
  FLOAT Sr;                              // RHS wave speed estimator
  FLOAT pstar;                           // Pressure in star region
  FLOAT ustar;                           // Velocity in star region
  FLOAT ddu;
  FLOAT dstarl;                          // Left star region density
  FLOAT etotstarl;                       // Left star region total energy
  FLOAT dstarr;                          // Right star region density
  FLOAT etotstarr;                       // Right star region total energy
  FLOAT df;                              // Flux density
  FLOAT uf;                              // Flux velocity
  FLOAT pf;                              // Flux pressure
  FLOAT etotf;                           // Flux total energy
  int k;

  // Compute total energy
  ekin = 0.5*dl*ul*ul;
#if NDIM>1
  ekin = ekin + 0.5*dl*qleft[4]*qleft[4];
#endif
#if NDIM>2
  ekin = ekin + 0.5*dl*qleft[5]*qleft[5];
#endif
  etotl = pl/(gamma - 1.0) + ekin;

  ekin = 0.5*dr*ur*ur;
#if NDIM>1
  ekin = ekin + 0.5*dr*qright[4]*qright[4];
#endif
#if NDIM>2
  ekin = ekin + 0.5*dr*qright[5]*qright[5];
#endif
  etotr = pr/(gamma - 1.0) + ekin;

  // Compute HLL wave speed
  Sl = min(ul,ur) - max(cl,cr);
  Sr = max(ul,ur) + max(cl,cr);

  ddu = dl*(Sl - ul) - dr*(Sr - ur);

  // Compute intermediate ('star') velocity and pressure
  ustar = ((pr - pl) + dl*ul*(Sl - ul) - dr*ur*(Sr - ur))/ddu;
  pstar = (dl*(Sl - ul)*pr - dr*(Sr - ur)*pr + dl*dr*(Sr - ur)*(Sl - ul)*(ul - ur))
          /ddu;

  // Left star region variables
  dstarl = dl*(Sl - ul)/(Sl - ustar);
  etotstarl = dstarl*(etotl/dl + (ustar - ul)*(ustar + pl/(dl*(Sl - ul))));

  // Right star region variables
  dstarr = dr*(Sr - ur)/(Sr - ustar);
  etotstarr = dstarr*(etotr/dr + (ustar - ur)*(ustar + pr/(dr*(Sr - ur))));

  // Sample solution at x/t=0
  if( Sl > 0.0 ) {
    df = dl;
    uf = ul;
    pf = pl;
    etotf = etotl;
  }
  else if( ustar > 0.0 ) {
    df = dstarl;
    uf = ustar;
    pf = pstar;
    etotf = etotstarl;
  }
  else if( Sr > 0.0 ) {
    df = dstarr;
    uf = ustar;
    pf = pstar;
    etotf = etotstarr;
  }
  else {
    df = dr;
    uf = ur;
    pf = pr;
    etotf = etotr;
  }

  // Compute the conservative fluxes
  /*flux[irho]  = df*uf;
  flux[ivx]   = df*uf*uf + pf;
  flux[ietot] = (etotf + pf)*uf;*/
/*#if NDIM>1
  do(ivar=3; ivar<nvar; ivar++) {
    if( flux[0] > 0.0d0) {
      flux[ivar] = flux[0]*qleft[ivar];
    }
    else {
      flux[ivar] = flux[0]*qright[ivar];
    }
  }
#endif
*/

  return;
}



template class RiemannSolver<1>;
template class RiemannSolver<2>;
template class RiemannSolver<3>;
template class ExactRiemannSolver<1>;
template class ExactRiemannSolver<2>;
template class ExactRiemannSolver<3>;
template class HllcRiemannSolver<1>;
template class HllcRiemannSolver<2>;
template class HllcRiemannSolver<3>;
