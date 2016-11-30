#define _USE_MATH_DEFINES
 
#include <stdio.h>
#include <cmath>  
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>

#include "logging.h"

#define MAXIT 10000 // Maximum number of iterations to run
#define TOL 1e-7 // Tolerance for summation difference
#define INTTOL 1e-7 // Tolerance for integration error

struct gap_params{
  int n; // Dimensionless Matsubara index (n+0.5)
  double d; // Dimensionless gap ($\delta$)
  double T; // Temperature
  gsl_spline* jacSpline; // Spline interpolation of the jacobian
  gsl_interp_accel* jacSplineAcc;
  gsl_spline* gapSpline; // Spline interpolation of the gap
  gsl_interp_accel* gapSplineAcc;
  gsl_spline* vkSpline; // Spline interpolation of the velocity
  gsl_interp_accel* vkSplineAcc; 
};

struct dos_params{
  double omega;
  double Delta;
  gsl_spline* jacSpline; // Spline interpolation of the jacobian
  gsl_interp_accel* jacSplineAcc;
  gsl_spline* gapSpline; // Spline interpolation of the gap
  gsl_interp_accel* gapSplineAcc;
};

double AIntegrand(double phi, void* p){
  struct gap_params *params = (struct gap_params*)p;
  double jac, gap, d;
  int n = params->n;
  d = params->d;
  jac = gsl_spline_eval(params->jacSpline, phi, params->jacSplineAcc);
  gap = gsl_spline_eval(params->gapSpline, phi, params->gapSplineAcc);

  return jac*gap*gap/sqrt(d*d*gap*gap + (n+0.5)*(n+0.5));
}

double rhoSfIntegrand(double phi, void* p){
  struct gap_params *params = (struct gap_params*)p;
  double jac, gap, vk, d, T;
  int n = params->n;
  d = params->d;
  T = params->T;
  jac = gsl_spline_eval(params->jacSpline, phi, params->jacSplineAcc);
  gap = gsl_spline_eval(params->gapSpline, phi, params->gapSplineAcc);
  vk = gsl_spline_eval(params->vkSpline, phi, params->vkSplineAcc);

  return 2*M_PI*T*jac*vk*vk*d*d*gap*gap/pow(d*d*gap*gap + 4.0*M_PI*M_PI*T*T*(n+0.5)*(n+0.5), 1.5);
}

double rhoSfNormIntegrand(double phi, void* p){
  struct gap_params *params = (struct gap_params*)p;
  double jac = gsl_spline_eval(params->jacSpline, phi, params->jacSplineAcc);
  double vk = gsl_spline_eval(params->vkSpline, phi, params->vkSplineAcc);
  double gap = gsl_spline_eval(params->gapSpline, phi, params->gapSplineAcc);
  return jac*vk*vk;
}

double ANormIntegrand(double phi, void* p){
  struct gap_params *params = (struct gap_params*)p;
  double jac = gsl_spline_eval(params->jacSpline, phi, params->jacSplineAcc);
  double gap = gsl_spline_eval(params->gapSpline, phi, params->gapSplineAcc);
  return jac*gap*gap;
}

double dosIntegrand(double phi, void* p){
  struct dos_params *params = (struct dos_params*)p;
  double omega = params->omega;
  double d = params->Delta;
  double jac = gsl_spline_eval(params->jacSpline, phi, params->jacSplineAcc);
  double gap = gsl_spline_eval(params->gapSpline, phi, params->gapSplineAcc);

  return jac*omega/sqrt(omega*omega + d*d*gap*gap);
}

double dosNormIntegrand(double phi, void* p){
  struct dos_params *params = (struct dos_params*)p;
  return gsl_spline_eval(params->jacSpline, phi, params->jacSplineAcc);
}

extern "C" void sd_omega_tilde(double* delta, double* omega, double* tp, double* phi0, double* phi1, double* jacAngles, double* jac, int* nJac, double* gapAngles, double* gap, int *nGap, double* result){
  // Setup logging
  openLog("librhosf.log");
  debug_print("%s| Begin: sd_omega_tilde\n", getTime());
  debug_print("%s| Settings: INTTOL: %g\n", getTime(), INTTOL);
  gsl_error_handler_t* old_handler = gsl_set_error_handler(&error_handler);

  debug_print("%s| Setting up integration workspace\n", getTime());
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

  debug_print("%s| Initializing jacobian interpolation\n", getTime());
  const gsl_interp_type *type = gsl_interp_cspline;
  gsl_interp_accel *jacSplineAcc = gsl_interp_accel_alloc();
  gsl_interp_accel *gapSplineAcc = gsl_interp_accel_alloc();
  gsl_spline *jacSpline = gsl_spline_alloc(type, *nJac);
  gsl_spline *gapSpline = gsl_spline_alloc(type, *nGap);

  debug_print("%s| Initializing jacobian interpolation\n", getTime());
  gsl_spline_init(jacSpline, jacAngles, jac, *nJac);
  debug_print("%s| Initializing gap interpolation\n", getTime());
  gsl_spline_init(gapSpline, gapAngles, gap, *nGap);

  struct dos_params params;
  params.jacSpline = jacSpline;
  params.jacSplineAcc = jacSplineAcc;
  params.gapSpline = gapSpline;
  params.gapSplineAcc = gapSplineAcc;

  gsl_function F;
  F.params = &params;

  double dos, error, norm;
  debug_print("%s| Integrating dosNormIntegrand\n", getTime());
  F.function = &dosNormIntegrand;
  gsl_integration_qag(&F, *phi0, *phi1, 0, INTTOL, 1000, 2, w, &norm, &error);
  
  F.function = &dosIntegrand;
  double omega_t, omega_t_new;
  omega_t = *omega;
  int n;
  for(n = 0; n < MAXIT; n++){
    params.omega = omega_t;
    params.Delta = *delta;
    gsl_integration_qag(&F, *phi0, *phi1, 0, INTTOL, 1000, 2, w, &dos, &error);
    omega_t_new = *omega + (*tp)*dos/norm;
    if(fabs(omega_t_new - omega_t) < TOL)
      break;
    omega_t = omega_t_new;
  }
  if(n < MAXIT)
    debug_print("%s| omega_tilde(delta=%g) ended on step %d, with result: %g\n", getTime(), *delta, n, omega_t_new);
  else
    // Use fprint here so that the warning is always logged
    fprintf(fpLog, "%s| WARNING: omega_tilde(delta=%g) ended on step %d=MAXIT, with result: %g. The result may not be valid to the prescribed tolerance: %g\n", getTime(), *delta, n, omega_t_new, TOL);

  *result = omega_t_new;

  debug_print("%s| Freeing memory for sd_omega_tilde\n", getTime());
  gsl_integration_workspace_free(w);
  gsl_spline_free(jacSpline);
  gsl_spline_free(gapSpline);
  gsl_interp_accel_free(jacSplineAcc);
  gsl_interp_accel_free(gapSplineAcc);
  debug_print("%s| End: sd_omega_tilde\n", getTime());
  closeLog();
}

extern "C" void sd_rhoSf(double* d, double* T, int *nd, double* phi0, double* phi1, double* jacAngles, double* jac, int *nJac, double* gapAngles, double* gap, int *nGap, double* vkAngles, double* vk, int* nVk, double* result){
  double norm, error, res, gapRes, newRes;
  // Setup logging
  openLog("librhosf.log");
  debug_print("%s| Begin: sd_rhoSf\n", getTime());
  debug_print("%s| Settings: MAXIT: %d, TOL: %g, INTTOL: %g\n", getTime(), MAXIT, TOL, INTTOL);
  gsl_error_handler_t* old_handler = gsl_set_error_handler(&error_handler);

  debug_print("%s| Setting up integration workspace\n", getTime());
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

  debug_print("%s| Initializing jacobian interpolation\n", getTime());
  const gsl_interp_type *type = gsl_interp_cspline;
  gsl_interp_accel *jacSplineAcc = gsl_interp_accel_alloc();
  gsl_interp_accel *gapSplineAcc = gsl_interp_accel_alloc();
  gsl_interp_accel *vkSplineAcc = gsl_interp_accel_alloc();
  gsl_spline *jacSpline = gsl_spline_alloc(type, *nJac);
  gsl_spline *gapSpline = gsl_spline_alloc(type, *nGap);
  gsl_spline *vkSpline = gsl_spline_alloc(type, *nVk);

  debug_print("%s| Initializing jacobian interpolation\n", getTime());
  gsl_spline_init(jacSpline, jacAngles, jac, *nJac);
  debug_print("%s| Initializing gap interpolation\n", getTime());
  gsl_spline_init(gapSpline, gapAngles, gap, *nGap);
  debug_print("%s| Initializing vk interpolation\n", getTime());
  gsl_spline_init(vkSpline, vkAngles, vk, *nVk);

  struct gap_params params;
  params.jacSpline = jacSpline;
  params.jacSplineAcc = jacSplineAcc;
  params.gapSpline = gapSpline;
  params.gapSplineAcc = gapSplineAcc;
  params.vkSpline = vkSpline;
  params.vkSplineAcc = vkSplineAcc;

  gsl_function F;
  F.params = &params;

  debug_print("%s| Integrating rhoSfNormIntegrand\n", getTime());
  F.function = &rhoSfNormIntegrand;
  gsl_integration_qag(&F, *phi0, *phi1, 0, INTTOL, 1000, 2, w, &norm, &error);
  
  F.function = &rhoSfIntegrand;
  int n;
  for(int j = 0; j < *nd; j++){
    debug_print("%s| Calculating rho_sf(d=%g)\n", getTime(), d[j]);
    gapRes = 0.0;
    for(n = 0; n < MAXIT; n++){
      debug_print("%s| Calculating Matsubara frequency (n=%d+0.5)\n", getTime(), n);
      params.n = n;
      params.d = d[j];
      params.T = T[j];
      gsl_integration_qag(&F, *phi0, *phi1, 0, INTTOL, 1000, 2, w, &res, &error);
      newRes = gapRes + res/norm;
      if(fabs(newRes - gapRes) < TOL)
        break;
      gapRes = newRes;
    }
    if(n < MAXIT)
      debug_print("%s| rhoSf(d=%g) ended on step %d, with result: %g\n", getTime(), d[j], n, gapRes);
    else
      // Use fprint here so that the warning is always logged
      fprintf(fpLog, "%s| WARNING: rhoSf(d=%g) ended on step %d=MAXIT, with result: %g. The result may not be valid to the prescribed tolerance: %g\n", getTime(), d[j], n, gapRes, TOL);
    result[j] = gapRes;
  }

  debug_print("%s| Freeing memory for sd_rhoSf\n", getTime());
  gsl_integration_workspace_free(w);
  gsl_spline_free(jacSpline);
  gsl_spline_free(gapSpline);
  gsl_spline_free(vkSpline);
  gsl_interp_accel_free(jacSplineAcc);
  gsl_interp_accel_free(gapSplineAcc);
  gsl_interp_accel_free(vkSplineAcc);
  debug_print("%s| End: sd_rhoSf\n", getTime());
  closeLog();
}

extern "C" void sd_A(double* d, int *nd, double* phi0, double* phi1, double* jacAngles, double* jac, int *nJac, double* gapAngles, double* gap, int *nGap, double* result){

  double error, norm, res, gapRes, newRes;
  // Setup logging
  openLog("librhosf.log");
  debug_print("%s| Begin: sd_A\n", getTime());
  debug_print("%s| Settings: MAXIT: %d, TOL: %g, INTTOL: %g\n", getTime(), MAXIT, TOL, INTTOL);
  gsl_error_handler_t* old_handler = gsl_set_error_handler(&error_handler);

  debug_print("%s| Setting up integration workspace\n", getTime());
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

  debug_print("%s| Creating interpolation objects for jacobian and gap angular dependence\n", getTime());
  const gsl_interp_type *type = gsl_interp_cspline;
  gsl_interp_accel *jacSplineAcc = gsl_interp_accel_alloc();
  gsl_interp_accel *gapSplineAcc = gsl_interp_accel_alloc();
  gsl_spline *jacSpline = gsl_spline_alloc(type, *nJac);
  gsl_spline *gapSpline = gsl_spline_alloc(type, *nGap);

  debug_print("%s| Initializing jacobian interpolation\n", getTime());
  gsl_spline_init(jacSpline, jacAngles, jac, *nJac);
  debug_print("%s| Initializing gap interpolation\n", getTime());
  gsl_spline_init(gapSpline, gapAngles, gap, *nGap);

  struct gap_params params;
  params.jacSpline = jacSpline;
  params.jacSplineAcc = jacSplineAcc;
  params.gapSpline = gapSpline;
  params.gapSplineAcc = gapSplineAcc;

  gsl_function F;
  F.params = &params;

  debug_print("%s| Integrating ANormIntegrand\n", getTime());
  F.function = &ANormIntegrand;
  gsl_integration_qag(&F, *phi0, *phi1, 0, INTTOL, 1000, 2, w, &norm, &error);
  
  F.function = &AIntegrand;
  int n;
  for(int j = 0; j < *nd; j++){
    debug_print("%s| Calculating A(d=%g)\n", getTime(), d[j]);
    gapRes = 0.0;
    for(n = 0; n < MAXIT; n++){
      debug_print("%s| Calculating Matsubara frequency (n=%d+0.5)\n", getTime(), n);
      params.n = n;
      params.d = d[j];
      gsl_integration_qag(&F, *phi0, *phi1, 0, INTTOL, 1000, 2, w, &res, &error);
      newRes = gapRes + 1.0/(n+0.5) - res/norm;
      if(fabs(newRes - gapRes) < TOL){
        break;
      }
      gapRes = newRes;
    }
    if(n < MAXIT)
      debug_print("%s| A(d=%g) ended on step %d, with result: %g\n", getTime(), d[j], n, gapRes);
    else
      // Use fprint here so that the warning is always logged
      fprintf(fpLog, "%s| WARNING: A(d=%g) ended on step %d=MAXIT, with result: %g. The result may not be valid to the prescribed tolerance: %g\n", getTime(), d[j], n, gapRes, TOL);
    result[j] = gapRes;
  }

  debug_print("%s| Freeing memory for sd_A\n", getTime());
  gsl_integration_workspace_free(w);
  gsl_spline_free(jacSpline);
  gsl_spline_free(gapSpline);
  gsl_interp_accel_free(jacSplineAcc);
  gsl_interp_accel_free(gapSplineAcc);
  gsl_set_error_handler(old_handler);
  debug_print("%s| End: sd_A\n", getTime());
  closeLog();
}
