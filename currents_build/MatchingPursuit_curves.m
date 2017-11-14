% 'exoShape' is a derived software by Stanley Durrleman, Copyright (C) INRIA (Asclepios team), All Rights Reserved, 2006-2009, version 1.0
%--------------------------------------------------------------------------
% Based on MATCHINE v1.0 software.
% Copyright Université Paris Descartes
% Contributor: Joan Alexis GLAUNES (2006)
% alexis.glaunes@mi.parisdescartes.fr
% 
% This software is a computer program whose purpose is to calculate an optimal 
% diffeomorphic transformation in 3D-space that allows to match two datasets 
% like points, curves or surfaces.
% 
% This software is governed by the CeCILL-B license under French law and
% abiding by the rules of distribution of free software. You can use, 
% modify and/ or redistribute the software under the terms of the CeCILL-B 
% license as circulated by CEA, CNRS and INRIA at the following URL
% "http://www.cecill.info". 
% 
% As a counterpart to the access to the source code and rights to copy,
% modify and redistribute granted by the license, users are provided only
% with a limited warranty  and the software's author, the holder of the
% economic rights, and the successive licensors have only limited 
% liability. 
% 
% In this respect, the user's attention is drawn to the risks associated
% with loading, using, modifying and/or developing or reproducing the
% software by the user in light of its specific status of free software,
% that may mean that it is complicated to manipulate, and that also
% therefore means that it is reserved for developers and experienced
% professionals having in-depth computer knowledge. Users are therefore
% encouraged to load and test the software's suitability as regards their
% requirements in conditions enabling the security of their systems and/or 
% data to be ensured and, more generally, to use and operate it in the 
% same conditions as regards security. 
% 
% The fact that you are presently reading this means that you have had
% knowledge of the CeCILL-B license and that you accept its terms.
%--------------------------------------------------------------------------
 
% Any use of this code should make reference to:
% - S. Durrleman, X. Pennec, A. Trouve and N. Ayache, Statistical Models of Sets of Curves and Surfaces based on Currents, Medical Image Analysis, (2009), DOI: 10.1016/j.media.2009.07.007
% 
% AND to at least one of the following references:
% - [for the use of 'landmark or measure matching']:
% J. Glaunes, L. Younes and A. Trouve, Diffeomorphic matching of distributions: A new approach for unlabelled point-sets and sub-manifolds matching, Proc. of the 2004 IEEE Computer Society Conference on Computer Vision and Pattern Recognition. (CVPR'04), Vol. 2, pp. 712--718, DOI: 10.1109/CVPR.2004.81 
% 
% - [for the use of 'curve matching']:
% J. Glaunes, A. Qiu, M. Miller, L. Younes, Large Deformation Diffeomorphic Metric Curve Mapping, International Journal of Computer vision, Springer, 2008, Vol. 80, No. 3, pp. 317--336, DOI: 10.1007/s11263-008-0141-9 
% 
% - [for the use of 'surface matching']:
% M. Vaillant and J. Glaunes, Surface Matching via Currents, Proc. of Information Processing in Medical Imaging (IPMI'05), Lecture Notes in Computer Science vol. 3565, Springer 2005, pp. 381--392

function [Tx Tvx] = Mean_projectionWstar(gamma,grille,tau,etype)

  % this function takes a field gamma (4D image of vectors) and returns an approximation of its associated current alpha.
  % the approximation is computed thanks to an orthogonal matching pursuit on the dictionnary of the dirac currents located at the points of a fixed grid.
  % etype is the standard deviation that controls the quality of the approximation
  % tau is the percentage of standard deviation to be explained before exiting


if (~isequal(size(gamma(:,:,:,1)),grille.long))
 disp('Les donnees n ont pas des dimensions coherentes !!');
end

beta = zeros(size(gamma));

nx = grille.long(1);
ny = grille.long(2);
nz = grille.long(3);
nn = nx*ny*nz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% translation de la grille %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x y z] = ndgrid(0:(nx-1),0:(ny-1),0:(nz-1));
x = x*grille.pas + grille.origine(1);
y = y*grille.pas + grille.origine(2);
z = z*grille.pas + grille.origine(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul de l'approximation a tau% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
noy = ifftn(grille.fft3k);
gammaN = sum(gamma.^2,4);
maxG2 = max(gammaN(:));
indMax = 5000;
indX = zeros(1,indMax);
indY = zeros(1,indMax);
indZ = zeros(1,indMax);
ind = zeros(1,indMax);
K = zeros(indMax);
c = 0;
% disp(['maxG = ' num2str(maxG) ' ; ecart-type = ' num2str(etype)]);
while ((sqrt(maxG2) > tau*etype)&&(c<=indMax))
  c = c+1;
  %recherche du point realisant le maximum
  [indXaux indYaux indZaux] = ind2sub(size(gammaN),find(gammaN == maxG2));
  indX(c) = indXaux(1); indY(c) = indYaux(1); indZ(c) = indZaux(1);
  ind(c) = indX(c) + nx*(indY(c)-1) + nx*ny*(indZ(c)-1);
  %construction de la matrice d'interpolation
  Kaux = circshift(noy,[indX(c)-1,indY(c)-1,indZ(c)-1]);
  Kaux = Kaux(ind(1:c));
  K(c,1:c) = Kaux;
  K(1:c,c) = Kaux';
  
  %projection sur la sous-grille : resolution du systeme lineaire
  for i=1:3
      beta(ind(1:c) + (i-1)*nn) = K(1:c,1:c)\gamma(ind(1:c) + (i-1)*nn)';
  end
  residu = gamma -  convol(beta,grille.fft3k_d);
  %mise a jour des parametres d'arret
  gammaN = sum(residu.^2,4);
  maxG2 = max(gammaN(:));
end

if (c>indMax)
   Tx = -1;
   Tvx = -1;
   disp(['Attention : plus de ' num2str(indMax) ' points ! verifiez les conditions aux bords du domaine...']);
   return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% representation des donnees %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tx = zeros(3,2*c);
for j=1:c
  pos = [x(ind(j)) y(ind(j)) z(ind(j))]';
  tau = [beta(ind(j)) beta(ind(j)+nn) beta(ind(j)+2*nn)]';
  Tx(:,2*j-1) = pos - tau/2;
  Tx(:,2*j) = pos + tau/2;
end
Tvx = [1:2:(2*c-1);2:2:(2*c)];

end










