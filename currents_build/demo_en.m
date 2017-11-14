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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% demo / tutorial for matching code %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% In this example, we will match two sets composed of four pairs
% source/target of types surface, curve, unlabelled points and landmarks

% 1/ In a first variable "s", we input everything concerning the
% deformation part.

% s.x must be a 3*n array containing 3D coordinates of all source points.
% Here the surface is composed of two triangles (4 vertices), the source 
% curve is sampled with four points, the points set counts 2 points and there
% are 2 source landmarks, thus n=12
clear s
s.x = [0,1,1,0,1.5,1.6,1.8,2, 3,3.2,3.8, 4;
       0,0,1,1, 0 , 0 , 0 ,0, 0, 0 , 0 , 0 ;
       0,0,0,0, 0 , 0 , 0 ,0, 0, 0 , 0 , 0];
   
% s.sigmaV is the kernel size which defines the scale of deformation.
% In other words, the deformation kernel is a function of the euclidean
% distance between two points divided by sigmaV. One should fix a value
% corresponding to the range of coordinates of points. For example, 
% if data is inside a box of length D, sigmaV = 0.2*D or 0.5*D can be
% good choices.
s.sigmaV = .5;

% Other important parameters (see also match.m file)

s.rigidmatching = 0; % put = 1 to perform a rigid matching before elastic matching.

s.gammaR = 0;        % weight of deformation part in functional
% gammaR = 0 means we minimize only the matching term,
% but in this case the deformation will remain regular due to the space in
% which minimization is performed. To compute distances between
% objects, one should not set it to zero though.

s.numbminims = 1; % means we only perform one minimization. 
% Fixing for example numbminims=3 will tell the algorithm to do
% 3 consecutive minimizations, decreasing the size of the matching kernels
% sigmaW, sigmaI (see below) after each of them.

s.usefgt = 0; % setting 1 would mean we use Fast Gauss Transform for convolutions 
% (should be used when number of points involved is above 50 points approximately).




% 2/ In a second variable "target", we input everything concerning the
% targets and the matching methods used

clear target

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First target of type surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

target{1}.method = 'surfcurr';

% Put indices of triangles of SOURCE surface in target{1}.vx
% These indices refer to corresponding columns of s.x
target{1}.vx = [1,1;
                2,3;
                3,4];
% Put coordinates of vertices of TARGET surface in target{1}.y
target{1}.y = [.5,1,0
                0,1,.5
               .5,1,.5];
% Put indices of triangles of TARGET surface in target{1}.vy            
target{1}.vy = [1;
                2;
                3];
% scale of surface matching kernel. Small values give more precision, but
% the algorithm may get stuck in local minima. So a good choice is a value
% corresponding to the gross distance between the source and target
% objects. However, if you use several minimizations (variable "numbminims",
% see above), then it is possible to fix a small value since the algorithm
% will perform the first steps at larger scales.
target{1}.sigmaW = .5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second target of type curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

target{2}.method = 'curvecurr';

% Put indices of segments of SOURCE curve on target{2}.vx
% These indices refer to corresponding columns of s.x
target{2}.vx = [5,6,7;
                6,7,8];
% Put coordinates of points of TARGET curve in target{2}.y
target{2}.y = [1.5, 2 ,2.5;
                0 , 0 ,0.5;
               0.5,0.5, 1 ];
% Put indices of segments of TARGET curve in target{2}.vy            
target{2}.vy = [1,2;
                2,3];
% scale of curve matching kernel
target{2}.sigmaW = .5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Third target of type point set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

target{3}.method = 'measures';

% Put indices of SOURCE points in target{3}.vx
% These indices refer to corresponding columns of s.x
target{3}.vx = [9,10];
% Put coordiantes of TARGET points in target{3}.y
target{3}.y = [3 , 3.2,3.4;
              0.5, 0.7,0.5;
               0.5,0.5,0.5];
% scale of point matching kernel
target{3}.sigmaI = 2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fourth target of type landmarks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

target{4}.method = 'landmarks';

% Put indices of SOURCE points in target{4}.vx
% These indices refer to crresponding columns of s.x
target{4}.vx = [11,12];
% Put coordinates of TARGET points in target{3}.y
target{4}.y = [ 4 , 4;
               0.5,0.5;
               0.5,0.8];



           
% Running the program:
s = match(s,target);

% The output structure "s" contains the variable "X" 
% (trajectories of template points) et "mom" (momentum vectors)
% which parameterize the optimal deformation. Also "distIdPhi" gives the
% distance to Identity in the space of diffeomorphisms.
disp(['deformation cost D = D(id,phi) = ',num2str(s.distIdPhi)])
disp(' ')

% display results
s.show = {'phi','y','xrig','x'};
s.showtraj = 1;    % display trajectories of source points (in blue)
s.showmomtraj = 0; % display momentum vectors (green arrows)
s.showpoints = 1;
s.showlegend = 0;
clf
affiche(s);
view(-20,20)
zoom(1.7)

% The program flow.m allows to compute the image of any 3D points
% under the optimal deformation map

disp('The''image of points')
V = rand(3,4)
disp('under the deformation map is')
phiV = flow(s,V)


