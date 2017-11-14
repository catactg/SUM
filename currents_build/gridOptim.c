/* To be compile with : mex gridOptim.c -lfftw3 */
/* You should have installed fftw package first. See www.fftw.org */

#include <fftw3.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mex.h>
#include <matrix.h>

/* gamma = convol(alpha,fft3k) is identical to ifftn(fftn(alpha).*fft3k) in matlab except that the size of fft3k is not the same to take into account its symmetry (since it's supposed to be the fft of a real matrix)*/

/* BE CAREFUL : fft3k should be a real matrix and the fft of real matrix. This is the case if fft3k = fftn(A) where A is real and symmetric.
/* alpha should be real, so is gamma*/

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]){

  mwSize *outdim, *auxdim;
  int argu, Ndim, N, i, k, nx, ny, nz, nnx, Nin, Nout, c0X, c0Y, c0Z;
  double *fft3k_R, *taux, *cx, *cy, *tauy, *grille_origine, *grille_long, *alpha;
  double pas, pas3, rho000, rho100, rho010, rho001, rho110, rho101, rho011, rho111, deltaX, deltaY, deltaZ;
  fftw_complex *aux;
  fftw_plan forward, backward;


  /* CHECK FOR PROPER INPUT ARGUMENTS */
  if(nrhs != 7){
    mexErrMsgTxt("Seven inputs required.");
  }
  if(nlhs != 1){
    mexErrMsgTxt("Exactly one output is required.");
  }


  /* ------- the 0st input argument: cx ------- */
  argu = 0;
  if (mxGetNumberOfDimensions(prhs[argu]) != 2){
      mexErrMsgTxt("cx should be 2-dimensional array.");
  }
  /* cx should be of size 3*Npts*/
  if (mxGetM(prhs[argu]) != 3){
      mexErrMsgTxt("cx should have 3 rows");
  }
  Nin = mxGetN(prhs[argu]);
  cx = mxGetPr(prhs[argu]);

  /* ------- the first input argument: taux ------- */
  argu = 1;
  if (mxGetNumberOfDimensions(prhs[argu]) != 2){
      mexErrMsgTxt("taux should be 2-dimensional array.");
  }
  if (mxGetN(prhs[argu]) != Nin){
      mexErrMsgTxt("taux should have the same number of columns as cx");
  }
  Ndim = mxGetM(prhs[argu]);
  taux = mxGetPr(prhs[argu]);


  /* ------- the second input argument: cy ------- */
  argu = 2;
  if (mxGetNumberOfDimensions(prhs[argu]) != 2){
      mexErrMsgTxt("cy should be 2-dimensional array.");
  }
  /* cx should be of size 3*Npts*/
  if (mxGetM(prhs[argu]) != 3){
      mexErrMsgTxt("cy should have 3 rows");
  }
  Nout = mxGetN(prhs[argu]);
  cy = mxGetPr(prhs[argu]);



  /* ------- the third input argument: grille.long ------- */
  argu = 3;
  if (mxGetN(prhs[argu]) != 3 || mxGetM(prhs[argu]) != 1){
      mexErrMsgTxt("grille.long should be 1x3 array.");
  }
  grille_long = mxGetPr(prhs[argu]);
  nx = (int) grille_long[0];
  ny = (int) grille_long[1];
  nz = (int) grille_long[2];
  nnx = (int) (nx/2 + 1);
  N = nx*ny*nz;

  /* ------- the fourth input argument: grille.pas ------- */
  argu = 4;
  if( !mxIsDouble(prhs[argu]) || mxIsComplex(prhs[argu]) ||
      mxGetN(prhs[argu])*mxGetM(prhs[argu])!=1 ) {
    mexErrMsgTxt("Input grille.pas must be a scalar.");
  }
  /*  get the scalar input pas */
  pas = mxGetScalar(prhs[argu]);
  if (pas <= 0)
	  mexErrMsgTxt("Input grille.pas must be a positive number.");
  pas3 = pas*pas*pas;

  /* ------- the fifth input argument: grille.origine ------- */
  argu = 5;
  if (mxGetM(prhs[argu]) != 3 || mxGetN(prhs[argu]) != 1){
      mexErrMsgTxt("grille.origine should be 3x1 array.");
  }
  grille_origine = mxGetPr(prhs[argu]);


   /* the sixth input argument: fft(noyau) */
  argu = 6;
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
  outdim = mxMalloc(2*sizeof(*outdim));
  outdim[0] = Ndim;
  outdim[1] = Nout;
  auxdim = mxMalloc(sizeof(*auxdim));
  auxdim[0] = nx*ny*nz*Ndim;
  plhs[0] = mxCreateNumericArray(2, outdim, mxDOUBLE_CLASS, mxREAL);
  tauy = mxGetPr(plhs[0]);

  /********/
  /* MAIN */
  /********/
  /* allocation */
/*   alpha = (double *) malloc(sizeof(double) * nx*ny*nz*Ndim); */
/*   gamma = (double *) malloc(sizeof(double) * nx*ny*nz*Ndim);*/
  alpha = (double *) mxMalloc(nx*ny*nz*Ndim*sizeof(double));
  aux = (fftw_complex *) fftw_malloc(nnx*ny*nz*sizeof(fftw_complex));

  /* Tri-linear projection of the normals into the grid*/
  for (k=0;k<(nx*ny*nz*Ndim);k++){
      alpha[k] = 0;
  }

  for (k=0;k<Nin;k++){
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
          if (((c0X+1+(c0Y+1)*nx+(c0Z+1)*nx*ny+i*N)>=(nx*ny*nz*Ndim))||(c0X+c0Y*nx+c0Z*nx*ny+ i*N<0)){
             printf("ERREUR PROJECTION\n");
             printf("GRILLE: %f %d %d %d %f %f %f\n",pas,nx,ny,nz,grille_origine[0],grille_origine[1],grille_origine[2]);
             printf("c0: %d %d %d\n",c0X,c0Y,c0Z);
             printf("%d >= %d\n",c0X+1+(c0Y+1)*nx+(c0Z+1)*nx*ny+i*N,nx*ny*nz*Ndim);
             return;
		  }
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

  /* convolution with Fast Fourier Transform*/
  #pragma omp parallel for
  {
  for (i=0;i<Ndim;i++){

     forward = fftw_plan_dft_r2c_3d(nz,ny,nx, alpha + i*N, aux, FFTW_ESTIMATE);
     fftw_execute(forward);

     /*ATTENTION : fft3k should be a real matrix !*/
     for (k=0;k<(nnx*ny*nz);k++){
       aux[k][0] *= fft3k_R[k]/N;
       aux[k][1] *= fft3k_R[k]/N;
     }

     backward = fftw_plan_dft_c2r_3d(nz,ny,nx, aux, alpha + i*N, FFTW_ESTIMATE);
     fftw_execute(backward);

  }
  }


  /* interpolation on the points cy */
  for (k=0; k<Nout; k++){
		c0X = (int) ((cy[3*k]   - grille_origine[0])/pas);
        c0Y = (int) ((cy[1+3*k] - grille_origine[1])/pas);
        c0Z = (int) ((cy[2+3*k] - grille_origine[2])/pas);

        deltaX = cy[3*k]   - (grille_origine[0] + c0X*pas);
        deltaY = cy[1+3*k] - (grille_origine[1] + c0Y*pas);
        deltaZ = cy[2+3*k] - (grille_origine[2] + c0Z*pas);


        rho000 = (pas-deltaX) * (pas-deltaY)  * (pas-deltaZ);
        rho100 = deltaX       * (pas-deltaY)  * (pas-deltaZ);
        rho010 = (pas-deltaX) * deltaY        * (pas-deltaZ);
        rho001 = (pas-deltaX) * (pas-deltaY)  * deltaZ      ;
        rho110 = deltaX       * deltaY        * (pas-deltaZ);
        rho011 = (pas-deltaX) * deltaY        * deltaZ      ;
        rho101 = deltaX       * (pas-deltaY)  * deltaZ      ;
        rho111 = deltaX       * deltaY        * deltaZ      ;

        /*printf("%f %f %f %f %f %f %f %f\n",rho000,rho100,rho010,rho001,rho110,rho011,rho101,rho111);*/

        for (i=0;i<Ndim; i++){
           if (((c0X+1+(c0Y+1)*nx+(c0Z+1)*nx*ny+i*N)>=(nx*ny*nz*Ndim))||(c0X+c0Y*nx+c0Z*nx*ny+ i*N<0)){
             printf("ERREUR INTERPOLATION\n");
             return;
		  }
          tauy[i + k*Ndim] = (alpha[c0X   + c0Y*nx     + c0Z*nx*ny     + i*N]*rho000+
          					  alpha[c0X+1 + c0Y*nx     + c0Z*nx*ny     + i*N]*rho100+
          					  alpha[c0X   + (c0Y+1)*nx + c0Z*nx*ny     + i*N]*rho010+
          					  alpha[c0X   + c0Y*nx     + (c0Z+1)*nx*ny + i*N]*rho001+
          					  alpha[c0X+1 + (c0Y+1)*nx + c0Z*nx*ny     + i*N]*rho110+
          					  alpha[c0X   + (c0Y+1)*nx + (c0Z+1)*nx*ny + i*N]*rho011+
          					  alpha[c0X+1 + c0Y*nx     + (c0Z+1)*nx*ny + i*N]*rho101+
          					  alpha[c0X+1 + (c0Y+1)*nx + (c0Z+1)*nx*ny + i*N]*rho111)/pas3;
		}
 }

  fftw_destroy_plan(forward);
  fftw_destroy_plan(backward);
  fftw_free(aux);
  fftw_cleanup();
  mxFree(outdim);
  mxFree(auxdim);
  mxFree(alpha);
/*  mxFree(gamma);*/
  return;

}
