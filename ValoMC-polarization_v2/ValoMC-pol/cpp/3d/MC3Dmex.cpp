
#include <string>
#define _USE_MATH_DEFINES
#define VALOMC_MEX
#include <cmath>
#include <limits>
#include <inttypes.h>
#include <string>
#include <vector>

#include "mex.h"
#include "Array.hpp"
#include "ArrayMEX.hpp"
#include "MC3D.hpp"
#include "../versionstring.h"

#include "matrix.h"

// Compiling (from MATLAB prompt):
//   mex MC3Dmex.cpp
//
// To compile with OpenMP (multithread) support (from MATLAB prompt):
//   mex -DUSE_OMP MC3Dmex.cpp CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
// Do not use OpenMP version if the MATLAB does not support the compiler used

time_t starting_time;


#ifdef _MAKE_CTRL_C_POSSIBLE_
extern "C" bool utIsInterruptPending();
#endif

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
  time_t now;
  time(&now);
  double timedifference = difftime(now,starting_time);
  
  #ifdef _MAKE_CTRL_C_POSSIBLE_
  if(utIsInterruptPending()) {
      mxDestroyArray(result);
      return false;
  }
  #endif

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
  mexPrintf("                 ValoMC-3D\n");
  char infobuf[5012];
  version_string(infobuf);
  mexPrintf("%s",infobuf);
  
  if ((nrhs != 26) || ((nlhs != 21) && (nlhs != 22)))
  {
    mexPrintf("nrhs %i nlhs %i", nrhs, nlhs);
    mexErrMsgTxt("Syntax:\n [I, Q, U, V, IB, QB, UB, VB, IR, QR, UR, VR, IT, QT, UT, VT, vsol, bsol, ebsol, simulationtime, rnseed, [HN]] = MC3Dmex(H, HN, BH, r, BCType, BCIntensity, BCLightDirectionType, BCLNormal, BCn, mua, mus, g, n, f, phase0, Nphoton, layer, s11, s12, s33, s43, S0,nangles,activate_pol, disablepbar, rnseed)\n");
  }
  mexPrintf("Initializing MC3D...\n");
  
  // Parse input
  Array<int_fast64_t> H, HN, BH, layer;
  Array<double> r, s11, s12, s33, s43, S0, mua, mus, g, n, phase0;
  Array<char> BCType, BCLightDirectionType;
  Array<double> BCLNormal, BCn, f, BCIntensity;
  Array<int_fast64_t> Nphoton, nangles, activate_pol;
  Array<double> GaussianSigma;
  Array<int_fast64_t> disable_pbar;
  Array<uint_fast64_t> rndseed;

  Convert_mxArray(prhs[0], H);
  Convert_mxArray(prhs[1], HN);
  Convert_mxArray(prhs[2], BH);
  Convert_mxArray(prhs[3], r);
  Convert_mxArray(prhs[4], BCType);
  Convert_mxArray(prhs[5], BCIntensity);    // [AL]: New array for light source intensity 
  Convert_mxArray(prhs[6], BCLightDirectionType); // [AL]: New array, determines if lightsource given relative to normal or not
  Convert_mxArray(prhs[7], BCLNormal);
  Convert_mxArray(prhs[8], BCn);
  Convert_mxArray(prhs[9], mua);
  Convert_mxArray(prhs[10], mus);
  Convert_mxArray(prhs[11], g);
  Convert_mxArray(prhs[12], n);
  Convert_mxArray(prhs[13], f);
  Convert_mxArray(prhs[14], phase0);
  Convert_mxArray(prhs[15], Nphoton);
  Convert_mxArray(prhs[16], layer);
  Convert_mxArray(prhs[17], s11);
  Convert_mxArray(prhs[18], s12);
  Convert_mxArray(prhs[19], s33);
  Convert_mxArray(prhs[20], s43);
  Convert_mxArray(prhs[21], S0);
  Convert_mxArray(prhs[22], nangles);
  Convert_mxArray(prhs[23], activate_pol);
  Convert_mxArray(prhs[24], disable_pbar);
  Convert_mxArray(prhs[25], rndseed);

//  Convert_mxArray(prhs[15], GaussianSigma); 

  // Set parameters to MC
  MC3D MC;
  MC.H = H;
  MC.HN = HN;
  MC.BH = BH;
  MC.r = r;
  MC.BCType = BCType;
  MC.BCIntensity = BCIntensity; // [AL]
  MC.BCLightDirectionType = BCLightDirectionType; // [AL]
  MC.BCLNormal = BCLNormal;
  MC.BCn = BCn;
  MC.mua = mua;
  MC.mus = mus;
  MC.g = g;
  MC.n = n;
  MC.f = f[0];
  MC.Nphoton = Nphoton[0];
  MC.layer =layer;
  MC.s11 = s11;
  MC.s12 = s12;
  MC.s33 = s33;
  MC.s43 = s43;
  MC.S0 = S0;
  MC.nangles = nangles[0];
  MC.activate_pol = activate_pol[0];
  MC.phase0 = phase0[0];
  //MC.GaussianSigma = GaussianSigma;
  //make negative phase0 positive

  if(MC.phase0 < 0) {
    MC.phase0 += 2*M_PI*ceil(-MC.phase0 / (2*M_PI));
  }
  if(rndseed[1]) {
     MC.seed = (unsigned long) rndseed[0];
  } else {
     MC.seed = (unsigned long) time(NULL);
  }
  // Initialize
  try {
    MC.ErrorChecks();
    MC.Init();
  } catch(mcerror e) {
    std::string message = "Error in initializing MC3D: " + std::string(errorstring(e)) + "\n"; 
    mexErrMsgTxt(message.c_str());
    return;
  }
  
  time(&starting_time);

  // Compute
  if(disable_pbar[0] == 0) {
     mexPrintf("Computing... \n");
    // Create a wait bar
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

  time_t now;

  // Show lossage
  if(MC.loss) mexPrintf(" %ld photons lost during computation!\n", MC.loss);

  // Copy solution from MC to output
  Array<double> out_I, out_Q, out_U, out_V, out_I_i, out_Q_i, out_U_i, out_V_i;
  Array<double> out_IB, out_QB, out_UB, out_VB, out_IB_i, out_QB_i, out_UB_i, out_VB_i;
  Array<double> out_IR, out_QR, out_UR, out_VR, out_IT, out_QT, out_UT, out_VT, out_IR_i, out_QR_i, out_UR_i, out_VR_i, out_IT_i, out_QT_i, out_UT_i, out_VT_i;
  Array<double> vsolr, vsoli, bsolr, bsoli;
  Array<double> dbsolr, dbsoli; // [AL]
  // start modify
  Convert_mxArray(&plhs[0], out_I, out_I_i, MC.I.Nx, MC.I.Ny);
  Convert_mxArray(&plhs[1], out_Q, out_Q_i, MC.Q.Nx, MC.Q.Ny);
  Convert_mxArray(&plhs[2], out_U, out_U_i, MC.U.Nx, MC.U.Ny);
  Convert_mxArray(&plhs[3], out_V, out_V_i, MC.V.Nx, MC.V.Ny);
  Convert_mxArray(&plhs[4], out_IB, out_IB_i, MC.IB.Nx, MC.IB.Ny);
  Convert_mxArray(&plhs[5], out_QB, out_QB_i, MC.QB.Nx, MC.QB.Ny);
  Convert_mxArray(&plhs[6], out_UB, out_UB_i, MC.UB.Nx, MC.UB.Ny);
  Convert_mxArray(&plhs[7], out_VB, out_VB_i, MC.VB.Nx, MC.VB.Ny);
  Convert_mxArray(&plhs[8], out_IR, out_IR_i, MC.IR.Nx, MC.IR.Ny);
  Convert_mxArray(&plhs[9], out_QR, out_QR_i, MC.QR.Nx, MC.QR.Ny);
  Convert_mxArray(&plhs[10], out_UR, out_UR_i, MC.UR.Nx, MC.UR.Ny);
  Convert_mxArray(&plhs[11], out_VR, out_VR_i, MC.VR.Nx, MC.VR.Ny);
  Convert_mxArray(&plhs[12], out_IT, out_IT_i, MC.IT.Nx, MC.IT.Ny);
  Convert_mxArray(&plhs[13], out_QT, out_QT_i, MC.QT.Nx, MC.QT.Ny);
  Convert_mxArray(&plhs[14],out_UT, out_UT_i, MC.UT.Nx, MC.UT.Ny);
  Convert_mxArray(&plhs[15],out_VT, out_VT_i, MC.VT.Nx, MC.VT.Ny);
  // end modify
  Convert_mxArray(&plhs[16], vsolr, vsoli, MC.ER.Nx, MC.ER.Ny);
  Convert_mxArray(&plhs[17], bsolr, bsoli, MC.EBR.Nx, MC.EBR.Ny);
  Convert_mxArray(&plhs[18], dbsolr, dbsoli, MC.DEBR.Nx, MC.DEBR.Ny);
  plhs[19]=mxCreateDoubleMatrix(1,1,mxREAL); // [AL]
  time(&now);

  *mxGetPr(plhs[19])=(double) difftime(now,starting_time);

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
  // MODIFY start
  for(ii = 0; ii < MC.I.N; ii++){
    out_I[ii] = MC.I[ii];
    out_I_i[ii] = MC.I_i[ii];
  }
  for(ii = 0; ii < MC.Q.N; ii++){
    out_Q[ii] = MC.Q[ii];
    out_Q_i[ii] = MC.Q_i[ii];
  }
  for(ii = 0; ii < MC.U.N; ii++){
    out_U[ii] = MC.U[ii];
    out_U_i[ii] = MC.U_i[ii];
  }
  for(ii = 0; ii < MC.V.N; ii++){
    out_V[ii] = MC.V[ii];
    out_V_i[ii] = MC.V_i[ii];
  }
  for(ii = 0; ii < MC.IB.N; ii++){
    out_IB[ii] = MC.IB[ii];
    out_IB_i[ii] = MC.IB_i[ii];
  }
  for(ii = 0; ii < MC.QB.N; ii++){
    out_QB[ii] = MC.QB[ii];
    out_QB_i[ii] = MC.QB_i[ii];
  }
  for(ii = 0; ii < MC.UB.N; ii++){
    out_UB[ii] = MC.UB[ii];
    out_UB_i[ii] = MC.UB_i[ii];
  }
  for(ii = 0; ii < MC.VB.N; ii++){
    out_VB[ii] = MC.VB[ii];
    out_VB_i[ii] = MC.VB_i[ii];
  }
  for(ii = 0; ii < MC.IR.N; ii++){
    out_IR[ii] = MC.IR[ii];
    out_IR_i[ii] = MC.IR_i[ii];
  }
  for(ii = 0; ii < MC.QR.N; ii++){
    out_QR[ii] = MC.QR[ii];
    out_QR_i[ii] = MC.QR_i[ii];
  }
  for(ii = 0; ii < MC.UR.N; ii++){
    out_UR[ii] = MC.UR[ii];
    out_UR_i[ii] = MC.UR_i[ii];
  }
  for(ii = 0; ii < MC.VR.N; ii++){
    out_VR[ii] = MC.VR[ii];
    out_VR_i[ii] = MC.VR_i[ii];
  }
  for(ii = 0; ii < MC.IT.N; ii++){
    out_IT[ii] = MC.IT[ii];
    out_IT_i[ii] = MC.IT_i[ii];
  }
  for(ii = 0; ii < MC.QT.N; ii++){
    out_QT[ii] = MC.QT[ii];
    out_QT_i[ii] = MC.QT_i[ii];
  }
  for(ii = 0; ii < MC.UT.N; ii++){
    out_UT[ii] = MC.UT[ii];
    out_UT_i[ii] = MC.UT_i[ii];
  }
  for(ii = 0; ii < MC.VT.N; ii++){
    out_VT[ii] = MC.VT[ii];
    out_VT_i[ii] = MC.VT_i[ii];
  }
  // modify end

  const mwSize dims[] = {1,1};
  plhs[20] = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
  *((unsigned long*) mxGetData(plhs[20])) = MC.seed;

  // Copy topology neighbourhood
  if(nlhs == 22){
    Array<long> HNo;
    Convert_mxArray(&plhs[21], HNo, MC.HN.Nx, MC.HN.Ny);
    for(ii = 0; ii < MC.HN.N; ii++) HNo[ii] = MC.HN[ii];
  }

  if(disable_pbar[0] == 0) {
    mexEvalString("delete(mcwaitbar);");
  }
  mexPrintf("Done\n");
}
