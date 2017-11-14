// 'exoShape' is a derived software by Stanley Durrleman, Copyright (C) INRIA (Asclepios team), All Rights Reserved, 2006-2009, version 1.0
// Based on MATCHINE v1.0 software.
// Copyright Université Paris Descartes
// Contributor: Joan Alexis GLAUNES (2006)
// alexis.glaunes@mi.parisdescartes.fr
//
// This software is a computer program whose purpose is to calculate an optimal 
// diffeomorphic transformation in 3D-space that allows to match two datasets 
// like points, curves or surfaces.
//
// This software is governed by the CeCILL-B license under French law and
// abiding by the rules of distribution of free software. You can use, 
// modify and/ or redistribute the software under the terms of the CeCILL-B 
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info". 
//
// As a counterpart to the access to the source code and rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author, the holder of the
// economic rights, and the successive licensors have only limited 
// liability. 
//
// In this respect, the user's attention is drawn to the risks associated
// with loading, using, modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean that it is complicated to manipulate, and that also
// therefore means that it is reserved for developers and experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or 
// data to be ensured and, more generally, to use and operate it in the 
// same conditions as regards security. 
//
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-B license and that you accept its terms.
// 
// Any use of this code should make reference to:
// - S. Durrleman, X. Pennec, A. Trouve and N. Ayache, Statistical Models of Sets of Curves and Surfaces based on Currents, Medical Image Analysis, (2009), DOI: 10.1016/j.media.2009.07.007
//  
// AND to at least one of the following references:
// - [for the use of 'landmark or measure matching']:
// J. Glaunes, L. Younes and A. Trouve, Diffeomorphic matching of distributions: A new approach for unlabelled point-sets and sub-manifolds matching, Proc. of the 2004 IEEE Computer Society Conference on Computer Vision and Pattern Recognition. (CVPR'04), Vol. 2, pp. 712--718, DOI: 10.1109/CVPR.2004.81 
// 
// - [for the use of 'curve matching']:
// J. Glaunes, A. Qiu, M. Miller, L. Younes, Large Deformation Diffeomorphic Metric Curve Mapping, International Journal of Computer vision, Springer, 2008, Vol. 80, No. 3, pp. 317--336, DOI: 10.1007/s11263-008-0141-9 
// 
// - [for the use of 'surface matching']:
// M. Vaillant and J. Glaunes, Surface Matching via Currents, Proc. of Information Processing in Medical Imaging (IPMI'05), Lecture Notes in Computer Science vol. 3565, Springer 2005, pp. 381--392

//-------------------------------------------------------------------
// Optimization routines based on Fast Gauss Transform use Figtree v0.9.3
// The code was written by Vlad Morariu, Vikas Raykar, and Changjiang Yang 
// and is copyrighted under the Lesser GPL: 
//
// Copyright (C) 2006-2010 Vlad Morariu and Vikas Raykar and Changjiang Yang 
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 or later.
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the GNU Lesser General Public License for more details. 
// You should have received a copy of the GNU Lesser General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, 
// MA 02111-1307, USA.  
//
// The author may be contacted via email at:
// morariu(at)umd.edu, vikas(at)umiacs(.)umd(.)edu, cyang(at)sarnoff(.)com
//-------------------------------------------------------------------


#include <fftw3.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <mex.h>

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

  /* OUTPUT ARGUMENT */
  plhs[0] = mxCreateNumericArray(Ndim, dim, mxDOUBLE_CLASS, mxREAL);

  /********/
  /* MAIN */
  /********/
  aux = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) *nnx*ny*nz);

  for (i=0;i<n_iter;i++){

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
  return;

}
