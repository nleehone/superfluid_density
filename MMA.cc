// MMA.cc
// Mathematica interface code
#include "WolframLibrary.h"
#include "WolframCompileLibrary.h"

DLLEXPORT mint WolframLibrary_getVersion(){
  return WolframLibraryVersion;
}

DLLEXPORT int WolframLibrary_initialize(WolframLibraryData libData){
  return 0;
}

extern "C"{
  void sd_A_clean(mreal* d, mint* nd, mreal* phi0, mreal* phi1, mreal* jacAngles, mreal* jac, mint* nJac, mreal* gapAngles, mreal* gap, mint* nGap, mreal* result);
  void sd_A_born(mreal* d, mint* nd, mreal* tp, mreal* phi0, mreal* phi1, mreal* jacAngles, mreal* jac, mint* nJac, mreal* gapAngles, mreal* gap, mint* nGap, mreal* result);
  void sd_omega_tilde(mreal* delta, mreal* omega, mreal* tp, mreal* phi0, mreal* phi1, mreal* jacAngles, mreal* jac, mint* nJac, mreal* gapAngles, mreal* gap, mint* nGap, mreal* result);
  void sd_rhoSf(mreal* d, mreal* T, mint* nd, mreal* Gamma, mreal* c, mreal* phi0, mreal* phi1, mreal* jacAngles, mreal* jac, mint* nJac, mreal* gapAngles, mreal* gap, mint* nGap, mreal* vkAngles, mreal* vk, mint* nVk, mreal* result);
}

EXTERN_C DLLEXPORT int rho_sf(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res){
  mreal phi0, phi1, Gamma, c;
  MTensor jac, jacAngles, gap, gapAngles, deltas, vk, vkAngles, result, temperatures;

  deltas = MArgument_getMTensor(Args[0]);
  temperatures = MArgument_getMTensor(Args[1]);
  Gamma = MArgument_getReal(Args[2]);
  c = MArgument_getReal(Args[3]);
  phi0 = MArgument_getReal(Args[4]);
  phi1 = MArgument_getReal(Args[5]);
  jacAngles = MArgument_getMTensor(Args[6]);
  jac = MArgument_getMTensor(Args[7]);
  gapAngles = MArgument_getMTensor(Args[8]);
  gap = MArgument_getMTensor(Args[9]);
  vkAngles = MArgument_getMTensor(Args[10]);
  vk = MArgument_getMTensor(Args[11]);

  libData->MTensor_new(MType_Real, 1, MTensor_getDimensionsMacro(deltas), &result);

  sd_rhoSf(
      MTensor_getRealDataMacro(deltas), 
      MTensor_getRealDataMacro(temperatures), 
      MTensor_getDimensionsMacro(deltas), 
      &Gamma,
      &c,
      &phi0, 
      &phi1, 
      MTensor_getRealDataMacro(jacAngles), 
      MTensor_getRealDataMacro(jac), 
      MTensor_getDimensionsMacro(jac), 
      MTensor_getRealDataMacro(gapAngles), 
      MTensor_getRealDataMacro(gap), 
      MTensor_getDimensionsMacro(gap), 
      MTensor_getRealDataMacro(vkAngles), 
      MTensor_getRealDataMacro(vk), 
      MTensor_getDimensionsMacro(vk), 
      MTensor_getRealDataMacro(result)
      );

  MArgument_setMTensor(Res, result);
  return LIBRARY_NO_ERROR;
}

EXTERN_C DLLEXPORT int A_born(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res){
  mreal phi0, phi1, tp;
  MTensor jac, jacAngles, gap, gapAngles, deltas, result;

  deltas = MArgument_getMTensor(Args[0]);
  tp = MArgument_getReal(Args[1]);
  phi0 = MArgument_getReal(Args[2]);
  phi1 = MArgument_getReal(Args[3]);
  jacAngles = MArgument_getMTensor(Args[4]);
  jac = MArgument_getMTensor(Args[5]);
  gapAngles = MArgument_getMTensor(Args[6]);
  gap = MArgument_getMTensor(Args[7]);

  libData->MTensor_new(MType_Real, 1, MTensor_getDimensionsMacro(deltas), &result);

  sd_A_born(MTensor_getRealDataMacro(deltas), MTensor_getDimensionsMacro(deltas), &tp, &phi0, &phi1, MTensor_getRealDataMacro(jacAngles), MTensor_getRealDataMacro(jac), MTensor_getDimensionsMacro(jac), MTensor_getRealDataMacro(gapAngles), MTensor_getRealDataMacro(gap), MTensor_getDimensionsMacro(gap), MTensor_getRealDataMacro(result));

  MArgument_setMTensor(Res, result);
  return LIBRARY_NO_ERROR;
}

EXTERN_C DLLEXPORT int A_clean(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res){
  mreal phi0, phi1;
  MTensor jac, jacAngles, gap, gapAngles, deltas, result;

  deltas = MArgument_getMTensor(Args[0]);
  phi0 = MArgument_getReal(Args[1]);
  phi1 = MArgument_getReal(Args[2]);
  jacAngles = MArgument_getMTensor(Args[3]);
  jac = MArgument_getMTensor(Args[4]);
  gapAngles = MArgument_getMTensor(Args[5]);
  gap = MArgument_getMTensor(Args[6]);

  libData->MTensor_new(MType_Real, 1, MTensor_getDimensionsMacro(deltas), &result);

  sd_A_clean(MTensor_getRealDataMacro(deltas), MTensor_getDimensionsMacro(deltas), &phi0, &phi1, MTensor_getRealDataMacro(jacAngles), MTensor_getRealDataMacro(jac), MTensor_getDimensionsMacro(jac), MTensor_getRealDataMacro(gapAngles), MTensor_getRealDataMacro(gap), MTensor_getDimensionsMacro(gap), MTensor_getRealDataMacro(result));

  MArgument_setMTensor(Res, result);
  return LIBRARY_NO_ERROR;
}

EXTERN_C DLLEXPORT int omega_tilde(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res){
  mreal phi0, phi1, delta, omega, tp, result;
  MTensor jac, jacAngles, gap, gapAngles;

  delta = MArgument_getReal(Args[0]);
  omega = MArgument_getReal(Args[1]);
  tp = MArgument_getReal(Args[2]);
  phi0 = MArgument_getReal(Args[3]);
  phi1 = MArgument_getReal(Args[4]);
  jacAngles = MArgument_getMTensor(Args[5]);
  jac = MArgument_getMTensor(Args[6]);
  gapAngles = MArgument_getMTensor(Args[7]);
  gap = MArgument_getMTensor(Args[8]);

  sd_omega_tilde(&delta, &omega, &tp, &phi0, &phi1, MTensor_getRealDataMacro(jacAngles), MTensor_getRealDataMacro(jac), MTensor_getDimensionsMacro(jac), MTensor_getRealDataMacro(gapAngles), MTensor_getRealDataMacro(gap), MTensor_getDimensionsMacro(gap), &result);

  MArgument_setReal(Res, result);
  return LIBRARY_NO_ERROR;
}
