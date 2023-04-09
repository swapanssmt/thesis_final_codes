
#include <string>
#define _USE_MATH_DEFINES

#include <cmath>
#include <limits>
#include <inttypes.h>
#include <time.h>

#include "mex.h"
#include "Array.hpp"
#include "ArrayMEX.hpp"
#include "MC2D.hpp"
#include "../versionstring.h"


#ifdef _MAKE_CTRL_C_POSSIBLE_
extern "C" bool utIsInterruptPending();
#endif

// Compiling (from MATLAB prompt):
//   mex MC2Dmex.cpp
//
// To compile with OpenMP (multithread) support (from MATLAB prompt):
//   mex -DUSE_OMP cpp/2d/MC2Dmex.cpp CFLAGS="$CFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp"
// Do not use OpenMP version if the MATLAB does not support the compiler used

time_t starting_time;

void finalchecks(int csum, int Nphoton) {
  if (csum != Nphoton)
  {
    mexPrintf("WARNING: RUN WAS ABORTED OR PARALLEL COMPUTING ENVIRONMENT IS NOT WORKING CORRECTLY. \n");
    // destroy progress bar
    mexEvalString("delete(mcwaitbar);");
  }
}

void finalchecks_destroy_bar(int csum, int Nphoton) {
   finalchecks(csum, Nphoton);
}


bool Progress_with_bar(double perc){
  //  printf("  %d %%\r", perc);
  mxArray *result;
  result=mexGetVariable("base", "abort_photonMC");
  if(result != NULL) {
    if(mxIsLogicalScalarTrue(result)) {
      mxDestroyArray(result);
      return false;
    }
  }

  #ifdef _MAKE_CTRL_C_POSSIBLE_
  if(utIsInterruptPending()) {
      mxDestroyArray(result);
      return false;
  }
  #endif

  time_t now;
  time(&now);
  double timedifference = difftime(now,starting_time);

  char matlabstring[5012];
  
  if(timedifference > 0) {
    
    double remainingtime = (100.0-perc)/(perc/timedifference);
    double hours = floor(remainingtime/(60*60));
    double minutes = floor((remainingtime - hours*60*60)/60);
    double seconds = (remainingtime - hours*60*60 - minutes*60);    
    
    sprintf(&matlabstring[0], "waitbar(%f,mcwaitbar,'%i hours %i minutes and %i seconds left');\n", perc / 100.0, (int) hours, (int) minutes, (int) ceil(seconds)); 
  //  mexPrintf("%s",matlabstring);
  } else {
     sprintf(&matlabstring[0],  "waitbar(0, mcwaitbar,'Estimating the time left');\n");    
  }

  mexEvalString(matlabstring);
  
  fflush(stdout);
  
  if(result != NULL) mxDestroyArray(result);
  
  return true;
}

bool Progress(double perc){
  mexPrintf("  %f %%\r", perc);

  return true;
}


void mexFunction(int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs)
{
  mexPrintf("                 ValoMC-2D\n");
  char infobuf[5012];
  version_string(infobuf);
  mexPrintf("%s",infobuf);

  if ((nrhs != 25) || ((nlhs != 14) && (nlhs != 15)))
  {
    mexPrintf("nrhs %i nlhs %i", nrhs, nlhs);
    mexErrMsgTxt("Syntax:\n [vsol, bsol, ebsol, R_vsol, R_bsol, R_ebsol, F_vsol, F_bsol, F_ebsol, F_R_vsol, F_R_bsol, F_R_ebsol,simulationtime, seed, [HN]] = MC2Dmex(H, HN, BH, r, BCType, BCIntensity, BCLightDirectionType, BCLNormal, BCn, mua_ex_sol, mua_ex_f, mua_em_sol, mus_ex, mus_em, g, n, f, phase0, Nphoton, Qyield_f, Tau_f, NBin2Dtheta, GaussianSigma, disablepbar, rndseed)\n");
  }

  // Parse input
  Array<int_fast64_t> H, HN, BH;
  Array<double> r, mua_ex_sol, mua_ex_f, mua_em_sol, mus_ex, mus_em, g, n;
  Array<char> BCType, BCLightDirectionType;
  Array<double> BCLNormal, BCn, f, Qyield_f, Tau_f, BCIntensity, phase0;
  Array<int_fast64_t> Nphoton;
  Array<int_fast64_t> NBin2Dtheta;
  Array<double> GaussianSigma;
  Array<int_fast64_t> disable_pbar;
  Array<uint_fast64_t> rndseed;

  Convert_mxArray(prhs[0], H);
  Convert_mxArray(prhs[1], HN);
  Convert_mxArray(prhs[2], BH);
  Convert_mxArray(prhs[3], r);
  Convert_mxArray(prhs[4], BCType);
  Convert_mxArray(prhs[5], BCIntensity);
  Convert_mxArray(prhs[6], BCLightDirectionType); // [AL]: New array for normal types
  Convert_mxArray(prhs[7], BCLNormal);
  Convert_mxArray(prhs[8], BCn);
  Convert_mxArray(prhs[9], mua_ex_sol);
  Convert_mxArray(prhs[10], mua_ex_f);
  Convert_mxArray(prhs[11], mua_em_sol);
  Convert_mxArray(prhs[12], mus_ex);
  Convert_mxArray(prhs[13], mus_em);
  Convert_mxArray(prhs[14], g);
  Convert_mxArray(prhs[15], n);
  Convert_mxArray(prhs[16], f);
  Convert_mxArray(prhs[17], phase0);
  Convert_mxArray(prhs[18], Nphoton);
  Convert_mxArray(prhs[19], Qyield_f);
  Convert_mxArray(prhs[20], Tau_f);
  Convert_mxArray(prhs[21], NBin2Dtheta);
  Convert_mxArray(prhs[22], GaussianSigma);
  Convert_mxArray(prhs[23], disable_pbar);
  Convert_mxArray(prhs[24], rndseed);

  // make negative phase0 positive by adding a multiple of 2*pi
  // Set parameters to MC

  MC2D MC;
  MC.H = H;
  MC.HN = HN;
  MC.BH = BH;
  MC.r = r;
  MC.BCType = BCType;
  MC.BCIntensity = BCIntensity; // [AL]
  MC.BCLightDirectionType = BCLightDirectionType; // [AL]
  MC.BCLNormal = BCLNormal;
  MC.BCn = BCn;
  MC.mua_ex_sol = mua_ex_sol;
  MC.mua_ex_f = mua_ex_f;
  MC.mua_em_sol = mua_em_sol;
  MC.mus_ex = mus_ex;
  MC.mus_em = mus_em;
  MC.g = g;
  MC.n = n;
  MC.f = f[0];
  MC.Nphoton = Nphoton[0];
  MC.Qyield_f = Qyield_f[0];
  MC.Tau_f = Tau_f[0];
  MC.NBin2Dtheta = NBin2Dtheta[0];
  MC.GaussianSigma = GaussianSigma;
  MC.phase0 = phase0[0];

  if(MC.phase0 < 0) {
    MC.phase0 += 2*M_PI*ceil(-MC.phase0 / (2*M_PI));
    mexPrintf("Transformed negative phase0 to positive %f\n", MC.phase0);
  }

  if(rndseed[1]) {
     MC.seed = rndseed[0];
  } else {
     MC.seed = (unsigned long) time(NULL);
  }
  // Initialize
  mexPrintf("Initializing MC2D...\n");
  try {
    MC.ErrorChecks();
    MC.Init();
  } catch(mcerror e) {
    std::string message = "Error in initializing MC2D: " + std::string(errorstring(e)) + "\n"; 
    mexErrMsgTxt(message.c_str());
    return;
  }
  // Create a wait bar
  
  time(&starting_time);

  // Compute
  if(disable_pbar[0] == 0) {
     mexPrintf("Computing... \n");
     mexEvalString("assignin('base','abort_photonMC', false);");
     mexEvalString("mcwaitbar = waitbar(0,'Please wait..', 'name', 'Running simulation', 'CreateCancelBtn','abort_photonMC=true;');");

     MC.MonteCarlo(Progress_with_bar, finalchecks_destroy_bar);
     mexPrintf("...done\n");
     printf("\n"); fflush(stdout);
  } else {
     mexPrintf("Computing... \n");
     MC.MonteCarlo(Progress, finalchecks);
     mexPrintf("...done\n");
     printf("\n"); fflush(stdout);
  }

  // Show lossage
  if(MC.loss) mexPrintf(" %ld photons lost during computation!\n", MC.loss);
  mexPrintf(" %ld fluorescence photons generated during computation!\n", MC.N_F_Photons);

  // Copy solution from MC to output
  Array<double> vsolr, vsoli, bsolr, bsoli;
  Array<double> dbsolr, dbsoli; // [AL]

  // ****************** MODIFY******************
  Array<double> R_vsolr, R_vsoli, R_bsolr, R_bsoli;
  Array<double> R_dbsolr, R_dbsoli; // [AL]
  //*****************************************
  //*********************** modify fluoroscence ***********************************
  Array<double> F_vsolr, F_vsoli, F_bsolr, F_bsoli;
  Array<double> F_dbsolr, F_dbsoli; // [AL]

  // ****************** MODIFY******************
  Array<double> F_R_vsolr, F_R_vsoli, F_R_bsolr, F_R_bsoli;
  Array<double> F_R_dbsolr, F_R_dbsoli; // [AL]
  //*****************************************
  //*******************************************************************************

  
  Convert_mxArray(&plhs[0], vsolr, vsoli, MC.ER.Nx, MC.ER.Ny);
  Convert_mxArray(&plhs[1], bsolr, bsoli, MC.EBR.Nx, MC.EBR.Ny);
  Convert_mxArray(&plhs[2], dbsolr, dbsoli, MC.DEBR.Nx, MC.DEBR.Ny);

  //*********************************modify************************
  Convert_mxArray(&plhs[3], R_vsolr, R_vsoli, MC.R_ER.Nx, MC.R_ER.Ny);
  Convert_mxArray(&plhs[4], R_bsolr, R_bsoli, MC.R_EBR.Nx, MC.R_EBR.Ny);
  Convert_mxArray(&plhs[5], R_dbsolr, R_dbsoli, MC.R_DEBR.Nx, MC.R_DEBR.Ny);
  //**************************************************************
  // ************************** Modify for fluoroscence *********************
  Convert_mxArray(&plhs[6], F_vsolr, F_vsoli, MC.F_ER.Nx, MC.F_ER.Ny);
  Convert_mxArray(&plhs[7], F_bsolr, F_bsoli, MC.F_EBR.Nx, MC.F_EBR.Ny);
  Convert_mxArray(&plhs[8], F_dbsolr, F_dbsoli, MC.F_DEBR.Nx, MC.F_DEBR.Ny);

  //*********************************modify************************
  Convert_mxArray(&plhs[9], F_R_vsolr, F_R_vsoli, MC.F_R_ER.Nx, MC.F_R_ER.Ny);
  Convert_mxArray(&plhs[10], F_R_bsolr, F_R_bsoli, MC.F_R_EBR.Nx, MC.F_R_EBR.Ny);
  Convert_mxArray(&plhs[11], F_R_dbsolr, F_R_dbsoli, MC.F_R_DEBR.Nx, MC.F_R_DEBR.Ny);
  // ***************************************************************************
  plhs[12]=mxCreateDoubleMatrix(1,1,mxREAL); // [AL]
  time_t now;
  time(&now);
  *mxGetPr(plhs[12])=(double) difftime(now,starting_time);

  const mwSize dims[] = {1,1};
  plhs[13] = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
  *((unsigned long*) mxGetData(plhs[13])) = MC.seed;

  long ii;
  for(ii = 0; ii < MC.ER.N; ii++){
    vsolr[ii] = MC.ER[ii];
    vsoli[ii] = MC.EI[ii];
  }
  for(ii = 0; ii < MC.EBR.N; ii++){
    bsolr[ii] = MC.EBR[ii];
    bsoli[ii] = MC.EBI[ii];
  }
  for(ii = 0; ii < MC.DEBR.N; ii++){
    dbsolr[ii] = MC.DEBR[ii];
    dbsoli[ii] = MC.DEBI[ii];
  }
  //***************************MODIFY*****************
  for(ii = 0; ii < MC.R_ER.N; ii++){
    R_vsolr[ii] = MC.R_ER[ii];
    R_vsoli[ii] = MC.R_EI[ii];
  }
  for(ii = 0; ii < MC.R_EBR.N; ii++){
    R_bsolr[ii] = MC.R_EBR[ii];
    R_bsoli[ii] = MC.R_EBI[ii];
  }
  for(ii = 0; ii < MC.R_DEBR.N; ii++){
    R_dbsolr[ii] = MC.R_DEBR[ii];
    R_dbsoli[ii] = MC.R_DEBI[ii];
  }

  //************************************************
  //*************** MODIFY FLUOROSCENCE ************************
  for(ii = 0; ii < MC.F_ER.N; ii++){
    F_vsolr[ii] = MC.F_ER[ii];
    F_vsoli[ii] = MC.F_EI[ii];
  }
  for(ii = 0; ii < MC.F_EBR.N; ii++){
    F_bsolr[ii] = MC.F_EBR[ii];
    F_bsoli[ii] = MC.F_EBI[ii];
  }
  for(ii = 0; ii < MC.F_DEBR.N; ii++){
    F_dbsolr[ii] = MC.F_DEBR[ii];
    F_dbsoli[ii] = MC.F_DEBI[ii];
  }
  //***************************MODIFY*****************
  for(ii = 0; ii < MC.F_R_ER.N; ii++){
    F_R_vsolr[ii] = MC.F_R_ER[ii];
    F_R_vsoli[ii] = MC.F_R_EI[ii];
  }
  for(ii = 0; ii < MC.F_R_EBR.N; ii++){
    F_R_bsolr[ii] = MC.F_R_EBR[ii];
    F_R_bsoli[ii] = MC.F_R_EBI[ii];
  }
  for(ii = 0; ii < MC.F_R_DEBR.N; ii++){
    F_R_dbsolr[ii] = MC.F_R_DEBR[ii];
    F_R_dbsoli[ii] = MC.F_R_DEBI[ii];
  }
  //****************************************************************

  // Copy topology neighbourhood
  if(nlhs == 15){
    Array<long> HNo;
    Convert_mxArray(&plhs[14], HNo, MC.HN.Nx, MC.HN.Ny);
    for(ii = 0; ii < MC.HN.N; ii++) HNo[ii] = MC.HN[ii]; 
  }

  if(disable_pbar[0] == 0) {
    mexEvalString("delete(mcwaitbar);");
  }
  mexPrintf("Done\n");
}
