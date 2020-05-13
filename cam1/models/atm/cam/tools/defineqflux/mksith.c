#include <math.h>

#include "defineqflux.h"

void mksith (double siccnt[PLAT][PLON], 
	     double sicthk[PLAT][PLON], 
	     double flat[PLAT],
	     int nlon[PLAT],
	     int n)
{
  int i, j;

  double r;
  double x;
  double y;
  double zfac;
  double theta;

  const double sicvar[12] = {.942, 1.000, 1.058, 1.124, 1.161, 1.175,
                             1.058, 0.931, 0.883, 0.880, 0.876, 0.912};
  const double sicmin = 0.25;  // minimum sea-ice thickness for initial dataset
  const double r1  = 3.5;
  const double phi = (120./360.)*2.*M_PI;
  const double psi = (  7./360.)*2.*M_PI;
  const double nx  = cos (phi) * sin (psi);
  const double ny  = sin (phi) * sin (psi);
  const double nz  = cos (psi);

  double sicmax;              // maximum sea-ice thickness for northern hemisphere
  
  // Northern hemisphere

  for (j = PLAT/2; j < PLAT; j++) {
    for (i = 0; i < nlon[j]; i++) {
      sicthk[j][i] = 0.;

      if (flat[j] >= +70.0)
	sicmax = 3.0;
      else if (flat[j] >= 50.)
	sicmax = sicmin + (3.0 - sicmin)*(flat[j]-50.)/20.;
      else 
	sicmax = sicmin;

      if (siccnt[j][i] > 0.) {
	if (siccnt[j][i] > 0.9)
	  sicthk[j][i] = sicmax;
	else if (siccnt[j][i] > 0.1)
	  sicthk[j][i] = sicmin + (siccnt[j][i] - 0.10)*(sicmax - sicmin)/0.80;
	else
	  sicthk[j][i] = sicmin;
      }

      // next, modify thickness based on seasonal variation:

      sicthk[j][i] = sicthk[j][i] * sicvar[n];
    }

    // make modification for arctic pack ice, poleward of 70N,
    // so that the largest thicknesses are just north west of Greenland:

    for (i = 0; i < nlon[j]; i++) {
      if (flat[j] >= 70.0) {
	r     = (90. - flat[j]) / r1;
	theta =  i*2.*M_PI / nlon[j];
	x     =  r * cos(theta);
	y     =  r * sin(theta);
	zfac  =  1. + (-(x*nx + y*ny)/nz);

	if (flat[j] < +78.0)
	  zfac =  1. + (-(x*nx + y*ny) / nz) * exp((flat[j] - 78.0) / r1);

	sicthk[j][i] = sicthk[j][i]*zfac;
      }
    }
  }

  // Southern hemisphere
  
  for (j = 0; j < PLAT/2; j++) {
    for (i = 0; i < nlon[j]; i++)
      if (siccnt[j][i] > 0.)
	sicthk[j][i] = siccnt[j][i] < 0.25 ? 0.25 + siccnt[j][i] : 0.50;
  }

  return;
}

      
