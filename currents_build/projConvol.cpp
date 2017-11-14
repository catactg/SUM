/* To be compile with : mex convol.c -lfftw */
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
  int *outdim;
  int argu, Ndim, N, i, k, nx, ny, nz, nnx, Nnorm, c0X, c0Y, c0Z;
  double *alpha, *gamma, *fft3k_R, *taux, *cx, *grille_origine, *grille_long;
  double pas, pas3, rho000, rho100, rho010, rho001, rho110, rho101, rho011, rho111, deltaX, deltaY, deltaZ;
  fftw_complex *aux;
  fftw_plan forward, backward;


  /* CHECK FOR PROPER INPUT ARGUMENTS */
  if(nrhs != 6){
    mexErrMsgTxt("Two inputs required.");
  }
  if(nlhs > 1){
    mexErrMsgTxt("Not more than one output required.");
  }

  /* ------- the 0st input argument: taux ------- */
  argu = 0;
  if (mxGetNumberOfDimensions(prhs[argu]) != 2){
      mexErrMsgTxt("taux should be 2-dimensional array.");
  }
  dim = mxGetDimensions(prhs[argu]);
  taux = mxGetPr(prhs[argu]);
  Ndim = dim[0];
  Nnorm = dim[1];


   /* ------- the first input argument: cx ------- */
  argu = 1;
  if (mxGetNumberOfDimensions(prhs[argu]) != 2){
      mexErrMsgTxt("cx should be 2-dimensional array.");
  }
  /* cx should be of size 3*Npts*/
  if (mxGetM(prhs[argu]) != 3){
      mexErrMsgTxt("cx should have 3 rows");
  }
  if (mxGetN(prhs[argu]) != Nnorm){
      mexErrMsgTxt("cx should have the same number of columns as taux");
  }
  cx = mxGetPr(prhs[argu]);


  /* ------- the second input argument: grille.long ------- */
  argu = 2;
  if (mxGetN(prhs[argu]) != 3 || mxGetM(prhs[argu]) != 1){
      mexErrMsgTxt("grille.long should be 1x3 array.");
  }
  grille_long = mxGetPr(prhs[argu]);
  nx = (int) grille_long[0];
  ny = (int) grille_long[1];
  nz = (int) grille_long[2];
  nnx = (int) (nx/2 + 1);
  N = nx*ny*nz;

  /* ------- the third input argument: grille.pas ------- */
  argu = 3;
  if( !mxIsDouble(prhs[argu]) || mxIsComplex(prhs[argu]) ||
      mxGetN(prhs[argu])*mxGetM(prhs[argu])!=1 ) {
    mexErrMsgTxt("Input grille.pas must be a scalar.");
  }
  /*  get the scalar input pas */
  pas = mxGetScalar(prhs[argu]);
  if (pas <= 0)
	  mexErrMsgTxt("Input grille.pas must be a positive number.");
  pas3 = pas*pas*pas;

  /* ------- the fourth input argument: grille.origine ------- */
  argu = 4;
  if (mxGetM(prhs[argu]) != 3 || mxGetN(prhs[argu]) != 1){
      mexErrMsgTxt("grille.origine should be 3x1 array.");
  }
  grille_origine = mxGetPr(prhs[argu]);


   /* the fifth input argument: fft(noyau) */
  argu = 5;
  if (mxGetNumberOfDimensions(prhs[argu]) != 3){
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


  /* OUTPUT ARGUMENT */
  outdim = (int *) malloc(4*sizeof(int));
  outdim[0] = nx;
  outdim[1] = ny;
  outdim[2] = nz;
  outdim[3] = Ndim;
  plhs[0] = mxCreateNumericArray(4, outdim, mxDOUBLE_CLASS, mxREAL);


  /********/
  /* MAIN */
  /********/

  /* allocation */
  alpha = (double *) malloc(sizeof(double) * nx*ny*nz*Ndim);
  aux = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) *nnx*ny*nz);

  /* Tri-linear projection of the normals into the grid*/
  for (k=0;k<(nx*ny*nz*Ndim);k++){
      alpha[k] = 0;
  }

  for (k=0;k<Nnorm;k++){
        c0X = (int) ((cx[3*k]   - grille_origine[0])/pas);
        c0Y = (int) ((cx[1+3*k] - grille_origine[1])/pas);
        c0Z = (int) ((cx[2+3*k] - grille_origine[2])/pas);

        deltaX = cx[3*k]   - (grille_origine[0] + c0X*pas);
        deltaY = cx[1+3*k] - (grille_origine[1] + c0Y*pas);
        deltaZ = cx[2+3*k] - (grille_origine[2] + c0Z*pas);


        rho000 = (pas-deltaX) * (pas-deltaY)  * (pas-deltaZ)/pas3;
        rho100 = deltaX       * (pas-deltaY)  * (pas-deltaZ)/pas3;
        rho010 = (pas-deltaX) * deltaY        * (pas-deltaZ)/pas3;
        rho001 = (pas-deltaX) * (pas-deltaY)  * deltaZ      /pas3;
        rho110 = deltaX       * deltaY        * (pas-deltaZ)/pas3;
        rho011 = (pas-deltaX) * deltaY        * deltaZ      /pas3;
        rho101 = deltaX       * (pas-deltaY)  * deltaZ      /pas3;
        rho111 = deltaX       * deltaY        * deltaZ      /pas3;

        /*printf("%f %f %f %f %f %f %f %f\n",rho000,rho100,rho010,rho001,rho110,rho011,rho101,rho111);*/

        for (i=0;i<Ndim; i++){
          alpha[c0X   + c0Y*nx     + c0Z*nx*ny     + i*N] += rho000*taux[Ndim*k + i];
          alpha[c0X+1 + c0Y*nx     + c0Z*nx*ny     + i*N] += rho100*taux[Ndim*k + i];
          alpha[c0X   + (c0Y+1)*nx + c0Z*nx*ny     + i*N] += rho010*taux[Ndim*k + i];
          alpha[c0X   + c0Y*nx     + (c0Z+1)*nx*ny + i*N] += rho001*taux[Ndim*k + i];
          alpha[c0X+1 + (c0Y+1)*nx + c0Z*nx*ny     + i*N] += rho110*taux[Ndim*k + i];
          alpha[c0X   + (c0Y+1)*nx + (c0Z+1)*nx*ny + i*N] += rho011*taux[Ndim*k + i];
          alpha[c0X+1 + c0Y*nx     + (c0Z+1)*nx*ny + i*N] += rho101*taux[Ndim*k + i];
          alpha[c0X+1 + (c0Y+1)*nx + (c0Z+1)*nx*ny + i*N] += rho111*taux[Ndim*k + i];

          /*printf("%f %f %f %f %f %f %f %f\n",alpha[c0X   + c0Y*nx     + c0Z*nx*ny     + i*nx*ny*nz],alpha[c0X+1 + c0Y*nx     + c0Z*nx*ny     + i*nx*ny*nz],alpha[c0X   + (c0Y+1)*nx + c0Z*nx*ny     + i*nx*ny*nz],alpha[c0X   + c0Y*nx     + (c0Z+1)*nx*ny + i*nx*ny*nz],alpha[c0X+1 + (c0Y+1)*nx + c0Z*nx*ny     + i*nx*ny*nz],alpha[c0X   + (c0Y+1)*nx + (c0Z+1)*nx*ny + i*nx*ny*nz],alpha[c0X+1 + c0Y*nx     + (c0Z+1)*nx*ny + i*nx*ny*nz],alpha[c0X+1 + (c0Y+1)*nx + (c0Z+1)*nx*ny + i*nx*ny*nz]);
*/
/*          printf("%d %d %d %d %d %d %d %d\n",c0X+c0Y*nx+c0Z*nx*ny+i*nx*ny*nz,c0X+1+c0Y*nx+c0Z*nx*ny+i*nx*ny*nz,c0X+(c0Y+1)*nx+c0Z*nx*ny+i*nx*ny*nz,c0X+c0Y*nx+(c0Z+1)*nx*ny+i*nx*ny*nz,c0X+1+(c0Y+1)*nx+c0Z*nx*ny+i*nx*ny*nz,c0X+(c0Y+1)*nx+(c0Z+1)*nx*ny+i*nx*ny*nz,c0X+1+c0Y*nx+(c0Z+1)*nx*ny+i*nx*ny*nz,c0X+1+(c0Y+1)*nx+(c0Z+1)*nx*ny+i*nx*ny*nz);
*/        }

  }

  /* convolution */
  for (i=0;i<Ndim;i++){

     forward = fftw_plan_dft_r2c_3d(nz,ny,nx,alpha + i*N, aux, FFTW_ESTIMATE);
     fftw_execute(forward);

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
  fftw_free(aux);
  fftw_cleanup();
  free(outdim);
  free(alpha);
  return;

}
