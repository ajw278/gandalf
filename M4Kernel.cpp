// ============================================================================
// M4Kernel.Cpp
// ============================================================================


#include <cmath>
#include <iostream>
#include "Constants.h"
#include "Dimensions.h"
#include "SphKernel.h"
using namespace std;



// ============================================================================
// M4Kernel::M4Kernel
// ============================================================================
M4Kernel::M4Kernel(int ndimaux)
{
  kernrange = 2.0;
  invkernrange = 0.5;
  kernrangesqd = 4.0;
  if (ndimaux == 1) kernnorm = twothirds;
  else if (ndimaux == 2) kernnorm = invpi*10.0/7.0;
  else if (ndimaux == 3) kernnorm = invpi;
#if !defined(FIXED_DIMENSIONS)
  ndimpr = (float) ndimaux;
#endif
}



// ============================================================================
// M4Kernel::~M4Kernel
// ============================================================================
M4Kernel::~M4Kernel()
{
}



// ============================================================================
// M4Kernel::Setup
// ============================================================================
//void M4Kernel::Setup(int ndim)
//{
//  return;
//}



// ============================================================================
// M4Kernel::w0
// ============================================================================
float M4Kernel::w0(float s)
{
  if (s < 1.0f)
    return kernnorm*(1.0 - 1.5*s*s + 0.75*s*s*s);
  else if (s < 2.0)
    return 0.25f*kernnorm*powf(2.0 - s,3.0);
  else
    return 0.0f;
}



// ============================================================================
// M4Kernel::w1
// ============================================================================
float M4Kernel::w1(float s)
{
  if (s < 1.0f)
    return kernnorm*(-3.0*s + 2.25*s*s);
  else if (s < 2.0)
    return -0.75*kernnorm*(2.0 - s)*(2.0 - s);
  else
    return 0.0;
}



// ============================================================================
// M4Kernel::womega
// ============================================================================
float M4Kernel::womega(float s)
{
  if (s < 1.0)
    return kernnorm*(-ndimpr + 1.5*(ndimpr + 2.0)*s*s - 
		     0.75*(ndimpr + 3.0)*pow(s,3));
  else if (s < 2.0)
    return kernnorm*(-2.0*ndimpr + 3.0*(ndimpr + 1.0)*s - 1.50*
		     (ndimpr + 2.0)*s*s + 0.25*(ndimpr + 3.0)*pow(s,3));
  else
    return 0.0;
}



// ============================================================================
// KERNEL::W1_TC
// ============================================================================
/*float Kernel::w1_m4_tc(float s)
{
  if (s < 1.0f)
    return kernnorm*(-3.0f*s + 2.25f*s*s);
  else if (s < 2.0f)
    return -0.75f*kernnorm*(2.0f - s)*(2.0f - s);
  else
    return 0.0f;
    }



// ============================================================================
// KERNEL::WO_QUINTIC
// ============================================================================
float Kernel::w0_quintic(float s)
{
  if (s < 1.0f)
    return kernnorm*(66.0f - 60.0f*s*s + 30.0f*powf(s,4) - 10.0f*powf(s,5));
  else if (s < 2.0f)
    return kernnorm*(51.0f+ 75.0f*s - 210.0f*s*s + 150.0f*pow(s,3) -
		     45.0f*powf(s,4) + 5.0f*powf(s,5));
  else if (s < 3.0f)
    return kernnorm*(243.0f - 405.0f*s + 270.0f*s*s - 
		     90.0f*powf(s,3) + 15.0f*powf(s,4) - powf(s,5));
  else
    return 0.0f;
}



// ============================================================================
// KERNEL::W1_QUINTIC
// ============================================================================
float Kernel::w1_quintic(float s) { return 0.0f; }


*/