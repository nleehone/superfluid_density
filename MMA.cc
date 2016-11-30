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
  void sd_A(mreal* d, mint* nd, mreal* phi0, mreal* phi1, mreal* jacAngles, mreal* jac, mint* nJac, mreal* gapAngles, mreal* gap, mint* nGap, mreal* result);
  void sd_rhoSf(mreal* d, mint* nd, mreal* phi0, mreal* phi1, mreal* jacAngles, mreal* jac, mint* nJac, mreal* gapAngles, mreal* gap, mint* nGap, mreal* vkAngles, mreal* vk, mint* nVk, mreal* result);
}

EXTERN_C DLLEXPORT int rho_sf(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res){
  mreal phi0, phi1;
  MTensor jac, jacAngles, gap, gapAngles, deltas, vk, vkAngles, result;

  deltas = MArgument_getMTensor(Args[0]);
  phi0 = MArgument_getReal(Args[1]);
  phi1 = MArgument_getReal(Args[2]);
  jacAngles = MArgument_getMTensor(Args[3]);
  jac = MArgument_getMTensor(Args[4]);
  gapAngles = MArgument_getMTensor(Args[5]);
  gap = MArgument_getMTensor(Args[6]);
  vkAngles = MArgument_getMTensor(Args[7]);
  vk = MArgument_getMTensor(Args[8]);

  libData->MTensor_new(MType_Real, 1, MTensor_getDimensionsMacro(deltas), &result);

  sd_rhoSf(MTensor_getRealDataMacro(deltas), MTensor_getDimensionsMacro(deltas), &phi0, &phi1, MTensor_getRealDataMacro(jacAngles), MTensor_getRealDataMacro(jac), MTensor_getDimensionsMacro(jac), MTensor_getRealDataMacro(gapAngles), MTensor_getRealDataMacro(gap), MTensor_getDimensionsMacro(gap), MTensor_getRealDataMacro(vkAngles), MTensor_getRealDataMacro(vk), MTensor_getDimensionsMacro(vk), MTensor_getRealDataMacro(result));

  MArgument_setMTensor(Res, result);
  return LIBRARY_NO_ERROR;
}

EXTERN_C DLLEXPORT int A(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res){
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

  sd_A(MTensor_getRealDataMacro(deltas), MTensor_getDimensionsMacro(deltas), &phi0, &phi1, MTensor_getRealDataMacro(jacAngles), MTensor_getRealDataMacro(jac), MTensor_getDimensionsMacro(jac), MTensor_getRealDataMacro(gapAngles), MTensor_getRealDataMacro(gap), MTensor_getDimensionsMacro(gap), MTensor_getRealDataMacro(result));

  MArgument_setMTensor(Res, result);
  return LIBRARY_NO_ERROR;
}

