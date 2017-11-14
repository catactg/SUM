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

function h = plotbyu(byufile,arg2,arg3)

% program to plot a byu surface on the current 3D graph
% usage : plotbyu(byufile)
%         plotbyu(byufile,rgbfile)
%         plotbyu(byufile,color)
%         plotbyu(byufile,'curv')
%         plotbyu(byufile,'curv',color)


fbyu = fopen(byufile,'r');

% read header
ncomponents = fscanf(fbyu,'%d',1);	% number of components
npoints = fscanf(fbyu,'%d',1);		% number of vertices
nfaces = fscanf(fbyu,'%d',1);		% number of faces
nedges = fscanf(fbyu,'%d',1);		% number of edges
fscanf(fbyu,'%d',2*ncomponents);	% components (ignored)

% read data
V = fscanf(fbyu,'%f',[3,npoints]);		% vertices
F = flipud(fscanf(fbyu,'%d',3*nfaces));		% faces

fclose(fbyu);

if nargin > 1
    if ischar(arg2)
        switch arg2
            case 'curv'
                C = meancurv(byufile);
                stdC = std(C);
                mnC = mean(C);
                lim = 1;
                minlim = mean(C)-std(C)*lim;
                maxlim = mean(C)+std(C)*lim;
                C(find(C<minlim)) = minlim;
                C(find(C>maxlim)) = maxlim;
                C = C - min(C);
                C = C / max(C);            
                C = shiftdim(reshape(C(abs(F)),3,prod(size(F))/3),-1);
                if nargin == 3
                    C = [arg3(1)*C;arg3(2)*C;arg3(3)*C];
                else
                    C = [C;C;C];
                end
                C = shiftdim(C,1);
            case 'nocurv'
                if nargin == 3
                    if size(arg3,1) == 3 && size(arg3,2) == 1
                        C = repmat(arg3,1,npoints);
                    elseif size(arg3,1) == 1 && size(arg3,2) == npoints
                        C = repmat(arg3,3,1);
                    else
                        C = arg3;
                    end 
                else
                    % default color is red
                    C = repmat([1;0;0],1,npoints);
                end 
            otherwise
                rgbfile = arg2;
                frgb = fopen(rgbfile,'r');
                C = fscanf(frgb,'%d',[3,npoints]);
                fclose(frgb);
        end
    elseif size(arg2,1) == 3 && size(arg2,2) == 1
        C = repmat(arg2,1,npoints);
    elseif size(arg2,1) == 1 && size(arg2,2) == npoints
        C = repmat(arg2,3,1);
    else
        C = arg2;
    end
else
    % default arg2 is red
    C = repmat([1;0;0],1,npoints);
end 

if min(min(min(C))) < max(max(max(C)))
    C = (C-min(min(min(C))))/(max(max(max(C)))-min(min(min(C))));
end

% plot polygons
hold on
ind = [find(F<0);nedges+1];
dind = diff(ind);

    F = abs(F);
    X = reshape(V(1,F),dind(1),nfaces);
    Y = reshape(V(2,F),dind(1),nfaces);
    Z = reshape(V(3,F),dind(1),nfaces);
    if ndims(C) ~= 3
        C = shiftdim(reshape(C(:,F),3,dind(1),nfaces),1);
    end
    h = patch(X,Y,Z,C);

% anti-bug visualisation
set(gcf,'Renderer','zbuffer')

axis equal
axis off
shading interp
camlight right
material dull

