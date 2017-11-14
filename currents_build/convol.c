/* To be compile with : mex convol.c -lfftw3 */
/* You should have install fftw package first. See www.fftw.org */

#include <fftw3.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <mex.h>

/* gamma = convol(alpha,fft3k) is identical to ifftn(fftn(alpha).*fft3k) in matlab execpt that the size of fft3k is not the same to take into account its symmetry (since it's supposed to be the fft of a real matrix)*/

/* BE CAREFUL : fft3k should be a real matrix and the fft of real matrix. This is the case if fft3k = fftn(A) where A is real and symmetric.
/* alpha should be real, so is gamma*/

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]){

  const int *dim;
  int argu, Ndim, N, i, k, nx, ny, nz, nnx, n_iter;
  double *alpha, *gamma, *fft3k_R, *fft3k_I;
  double auxR;
  fftw_complex *aux;
  fftw_plan forward, backward;


  /* CHECK FOR PROPER INPUT ARGUMENTS */
  if(nrhs != 2){
    mexErrMsgTxt("Two inputs required.");
  }
  if(nlhs > 1){
    mexErrMsgTxt("Not more than one output required.");
  }
  /* the 0st input argument: alpha */
  Ndim = 4;
  argu = 0;
  if (mxGetNumberOfDimensions(prhs[argu]) != Ndim){
      mexErrMsgTxt("Alpha should be 4-dimensional array.");
  }
  dim = mxGetDimensions(prhs[argu]);
  alpha = mxGetPr(prhs[argu]);
  nx = dim[0];
  ny = dim[1];
  nz = dim[2];
  nnx = (int) (nx/2 + 1);
  N = nx*ny*nz;
  n_iter = dim[3];

   /* the first input argument: fft(noyau) */
  argu = 1;
  if (mxGetNumberOfDimensions(prhs[argu]) != (Ndim-1)){
      mexErrMsgTxt("noyau should be 3-dimensional array.");
  }
  /* fft3k should be of size nnx*ny*nz: keep only the (nx/2=1) rows of the matrix fftn(noyau)!*/
  if (mxGetDimensions(prhs[argu])[0] != nnx){
      mexErrMsgTxt("dimension of noyau along x-axis mismatch");
  }
  if (mxGetDimensions(prhs[argu])[1] != ny){
      mexErrMsgTxt("dimension of noyau along y-axis mismatch");
  }
  if (mxGetDimensions(prhs[argu])[2] != nz){
      mexErrMsgTxt("dimension of noyau along z-axis mismatch");
  }
  fft3k_R = mxGetPr(prhs[argu]);
/*  fft3k_I = mxGetPi(prhs[argu]);*/

  /* OUTPUT ARGUMENT */
  plhs[0] = mxCreateNumericArray(Ndim, dim, mxDOUBLE_CLASS, mxREAL);

  /********/
  /* MAIN */
  /********/
  aux = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) *nnx*ny*nz);

  for (i=0;i<n_iter;i++){

     forward = fftw_plan_dft_r2c_3d(nz,ny,nx,alpha + i*N, aux, FFTW_ESTIMATE);
     fftw_execute(forward);

 /*    for (k=0;k<(nnx*ny*nz);k++){
       auxR = aux[k][0];
       aux[k][0] = (fft3k_R[k]*auxR - fft3k_I[k]*aux[k][1])/N;
       aux[k][1] = (fft3k_R[k]*aux[k][1] + fft3k_I[k]*auxR)/N;
     }
*/

     /*ATTENTION : fft3k should be a real matrix !*/
     for (k=0;k<(nnx*ny*nz);k++){
       aux[k][0] *= fft3k_R[k]/N;
       aux[k][1] *= fft3k_R[k]/N;
     }

     gamma = mxGetPr(plhs[0])+i*N;
     backward = fftw_plan_dft_c2r_3d(nz,ny,nx, aux, gamma, FFTW_ESTIMATE);
     fftw_execute(backward);

  }

  fftw_destroy_plan(forward);
  fftw_destroy_plan(backward);
  fftw_cleanup();
  fftw_free(aux);
  return;

}
