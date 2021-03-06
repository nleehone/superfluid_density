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
  double omega_t;
  double d; // Dimensionless gap ($\delta$)
  double T; // Temperature
  double Delta;
  double Gamma; 
  double c;
  //double P;
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
  double jac, gap, d, P;
  double omega_t = params->omega_t;
  //P = params->P;
  int n = params->n;
  d = params->d;
  jac = gsl_spline_eval(params->jacSpline, phi, params->jacSplineAcc);
  gap = gsl_spline_eval(params->gapSpline, phi, params->gapSplineAcc);

  return jac*gap*gap/sqrt(d*d*gap*gap + omega_t*omega_t);
}

double rhoSfIntegrand(double phi, void* p){
  struct gap_params *params = (struct gap_params*)p;
  double jac, gap, vk, d, T;
  int n = params->n;
  double omega_t = params->omega_t;
  d = params->d;
  T = params->T;
  jac = gsl_spline_eval(params->jacSpline, phi, params->jacSplineAcc);
  gap = gsl_spline_eval(params->gapSpline, phi, params->gapSplineAcc);
  vk = gsl_spline_eval(params->vkSpline, phi, params->vkSplineAcc);

  return jac*vk*vk*d*d*gap*gap/pow(d*d*gap*gap + omega_t*omega_t, 1.5);
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
  struct gap_params *params = (struct gap_params*)p;
  double omega = params->omega_t;
  double d = params->Delta;
  double jac = gsl_spline_eval(params->jacSpline, phi, params->jacSplineAcc);
  double gap = gsl_spline_eval(params->gapSpline, phi, params->gapSplineAcc);

  return jac*omega/sqrt(omega*omega + d*d*gap*gap);
}

double dosNormIntegrand(double phi, void* p){
  struct gap_params *params = (struct gap_params*)p;
  double jac = gsl_spline_eval(params->jacSpline, phi, params->jacSplineAcc);
  return jac;
}

double omega_tilde2(double omega, double Delta, double Gamma1, double Gamma2, double phi0, double phi1, struct gap_params* params, gsl_integration_workspace* w){
  int n;
  gsl_function F;
  F.params = params;
  double dos, error, norm;
  F.function = &dosNormIntegrand;
  gsl_integration_qag(&F, phi0, phi1, 0, INTTOL, 1000, 2, w, &norm, &error);
  F.function = &dosIntegrand;
  double omega_t = omega;
  double omega_t_new;
  for(n = 0; n < MAXIT; n++){
    params->omega_t = omega_t;
    params->Delta = Delta;
    gsl_integration_qag(&F, phi0, phi1, 0, INTTOL, 1000, 2, w, &dos, &error);
    dos /= norm;
    omega_t_new = omega + M_PI*Gamma1*dos + M_PI*Gamma2/dos; 
    if(fabs(omega_t_new - omega_t) < TOL)
      break;
    omega_t = omega_t_new;
  }
  return omega_t_new;
}

double omega_tilde(double omega, double Delta, double Gamma, double c, double phi0, double phi1, struct gap_params* params, gsl_integration_workspace* w){
  int n;
  gsl_function F;
  F.params = params;
  double dos, error, norm;
  F.function = &dosNormIntegrand;
  gsl_integration_qag(&F, phi0, phi1, 0, INTTOL, 1000, 2, w, &norm, &error);
  F.function = &dosIntegrand;
  double omega_t = omega;
  double omega_t_new;
  for(n = 0; n < MAXIT; n++){
    params->omega_t = omega_t;
    params->Delta = Delta;
    gsl_integration_qag(&F, phi0, phi1, 0, INTTOL, 1000, 2, w, &dos, &error);
    dos /= norm;
    if(c >= 0.0)
      omega_t_new = omega + M_PI*Gamma*dos/(c*c + dos*dos);
    else
      omega_t_new = omega + M_PI*Gamma*dos;
    if(fabs(omega_t_new - omega_t) < TOL)
      break;
    omega_t = omega_t_new;
  }
  return omega_t_new;
}

extern "C" void sd_rhoSf(double* d, double* T, int *nd, double* gamma, double* c, double* phi0, double* phi1, double* angles, double* jac, double* gap, double* vk, int* nAngles, double* result){
  double norm, error, res, gapRes, newRes;
  // Setup logging
  openLog("librhosf.log");
  gsl_error_handler_t* old_handler = gsl_set_error_handler(&error_handler);

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

  const gsl_interp_type *type = gsl_interp_cspline;
  gsl_interp_accel *jacSplineAcc = gsl_interp_accel_alloc();
  gsl_interp_accel *gapSplineAcc = gsl_interp_accel_alloc();
  gsl_interp_accel *vkSplineAcc = gsl_interp_accel_alloc();
  gsl_spline *jacSpline = gsl_spline_alloc(type, *nAngles);
  gsl_spline *gapSpline = gsl_spline_alloc(type, *nAngles);
  gsl_spline *vkSpline = gsl_spline_alloc(type, *nAngles);

  gsl_spline_init(jacSpline, angles, jac, *nAngles);
  gsl_spline_init(gapSpline, angles, gap, *nAngles);
  gsl_spline_init(vkSpline, angles, vk, *nAngles);

  struct gap_params params;
  params.jacSpline = jacSpline;
  params.jacSplineAcc = jacSplineAcc;
  params.gapSpline = gapSpline;
  params.gapSplineAcc = gapSplineAcc;
  params.vkSpline = vkSpline;
  params.vkSplineAcc = vkSplineAcc;
  //params.P = *P;
  params.c = *c;

  gsl_function F;
  F.params = &params;

  F.function = &rhoSfNormIntegrand;
  gsl_integration_qag(&F, *phi0, *phi1, 0, INTTOL, 1000, 2, w, &norm, &error);
  
  F.function = &rhoSfIntegrand;
  int n;
  for(int j = 0; j < *nd; j++){
    gapRes = 0.0;
    for(n = 0; n < MAXIT; n++){
      params.Gamma = *gamma/(2.0*M_PI*T[j]);
      params.omega_t = (n+0.5);
      if(*gamma > 0)
        params.omega_t = omega_tilde(params.omega_t, d[j], *gamma/(2.0*M_PI*T[j]), *c, *phi0, *phi1, &params, w);
      else if(*gamma < 0)
        params.omega_t = omega_tilde2(params.omega_t, d[j], -*gamma/(2.0*M_PI*T[j]), *c, *phi0, *phi1, &params, w);
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

  gsl_integration_workspace_free(w);
  gsl_spline_free(jacSpline);
  gsl_spline_free(gapSpline);
  gsl_spline_free(vkSpline);
  gsl_interp_accel_free(jacSplineAcc);
  gsl_interp_accel_free(gapSplineAcc);
  gsl_interp_accel_free(vkSplineAcc);
  closeLog();
}

extern "C" void sd_A(double* d, int *nd, double* Gamma, double* t, double* c, double* phi0, double* phi1, double* angles, double* jac, double* gap, int *nAngles, double* result){

  double error, norm, res, gapRes, newRes;
  // Setup logging
  openLog("librhosf.log");
  gsl_error_handler_t* old_handler = gsl_set_error_handler(&error_handler);

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

  const gsl_interp_type *type = gsl_interp_cspline;
  gsl_interp_accel *jacSplineAcc = gsl_interp_accel_alloc();
  gsl_interp_accel *gapSplineAcc = gsl_interp_accel_alloc();
  gsl_spline *jacSpline = gsl_spline_alloc(type, *nAngles);
  gsl_spline *gapSpline = gsl_spline_alloc(type, *nAngles);

  gsl_spline_init(jacSpline, angles, jac, *nAngles);
  gsl_spline_init(gapSpline, angles, gap, *nAngles);

  struct gap_params params;
  params.jacSpline = jacSpline;
  params.jacSplineAcc = jacSplineAcc;
  params.gapSpline = gapSpline;
  params.gapSplineAcc = gapSplineAcc;
  params.Gamma = *Gamma;
  params.c = *c;

  gsl_function F;
  F.params = &params;

  F.function = &ANormIntegrand;
  gsl_integration_qag(&F, *phi0, *phi1, 0, INTTOL, 1000, 2, w, &norm, &error);
  
  F.function = &AIntegrand;
  int n;
  for(int j = 0; j < *nd; j++){
    gapRes = 0.0;
    for(n = 0; n < MAXIT; n++){
      params.omega_t = (n+0.5);
      params.omega_t = omega_tilde(params.omega_t, d[j], *Gamma/(2*M_PI*(*t)), *c, *phi0, *phi1, &params, w);
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

  gsl_integration_workspace_free(w);
  gsl_spline_free(jacSpline);
  gsl_spline_free(gapSpline);
  gsl_interp_accel_free(jacSplineAcc);
  gsl_interp_accel_free(gapSplineAcc);
  gsl_set_error_handler(old_handler);
  closeLog();
}
