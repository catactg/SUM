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

function s = shooting(s)

% geodesic shooting for point matching

T = s.T;
tau = 1/(T-1);
nx = size(s.X,2);
X = zeros(3,nx,T);
mom = zeros(3,nx,T);
X(:,:,1) = s.transmatrix * s.x + repmat(s.transvector,1,nx);
mom(:,:,1) = s.mom(:,:,1);
sigmaV = s.sigmaV;
sigmaV2 = sigmaV.^2;
stdV = s.stdV;
stdV2 = stdV.^2;
usefgt = 0;
usegrid = 0;
if isfield(s,'usefgt')&&s.usefgt
    usefgt = 1;
%     order = s.order;
%     K = s.K;
%     e = s.e;
      epsilon = s.epsilon;
elseif isfield(s,'usegrid')&&s.usegrid
    usegrid = 1;
    ratio = s.ratio;
end


if usegrid
    sourcegrid = cell(1,T);
    mini = min(X(:,:,1),[],2);
    maxi = max(X(:,:,1),[],2);
    sourcegrid{1} = setgrid(mini,maxi,1);
    for i=2:T
        sourcegrid{i} = sourcegrid{1};
    end
    disp(['source''s grid 1 has been set:  ' num2str([sourcegrid{1}.pas sourcegrid{1}.long sourcegrid{1}.origine'])]);
 end

dX = zeros(3,nx);
dmom = zeros(3,nx);
for t=1:T-1
    kernelsum(t);
    dkernelsum(t);
   
    X(:,:,t+1) = X(:,:,t) + tau * dX;
    mom(:,:,t+1) = mom(:,:,t) + (tau * 2 / sigmaV2(t)) * dmom;
    if usegrid
        updatesourcegrid(t+1)
    end
    
    kernelsum(t+1);
    dkernelsum(t+1);
    
    dX = .5 * dX;
    dmom = .5 * dmom;

    X(:,:,t+1) = X(:,:,t) + tau * dX;
    mom(:,:,t+1) = mom(:,:,t) + (tau * 2 / sigmaV2(t)) * dmom;
    if usegrid
        updatesourcegrid(t+1)
    end

    dX(:) = 0;
    dmom(:) = 0;
end

iter = 0;
method = 'shooting';

clear dX dmom
savestruct('s',who)

    function kernelsum(t)
        if ~(usefgt||usegrid)
            for m = 1:nx
                for l = 1:nx
                    argin = -( ...
                        (X(1,m,t)-X(1,l,t))^2 + ...
                        (X(2,m,t)-X(2,l,t))^2 + ...
                        (X(3,m,t)-X(3,l,t))^2)/sigmaV2(t);
argout = stdV2(t)*exp(argin);  %% BUILT IN KERNEL kerV, do not remove this comment
                    dX(1,m) = dX(1,m) + argout * mom(1,l,t);
                    dX(2,m) = dX(2,m) + argout * mom(2,l,t);
                    dX(3,m) = dX(3,m) + argout * mom(3,l,t);
                end
            end
        elseif usefgt
            dX = dX + stdV2(t)*figtree(3,X(:,:,t),mom(:,:,t),X(:,:,t),sigmaV(t),epsilon);
        else
            dX = dX + stdV2(t)*gridOptim(X(:,:,t),mom(:,:,t),X(:,:,t),sourcegrid{t}.long,sourcegrid{t}.pas,sourcegrid{t}.origine,sourcegrid{t}.fft3k_d);
        end
    end



    function dkernelsum(t)
        if ~(usefgt||usegrid)
            for m = 1:nx
                for l = 1:nx
                    argin = -( ...
                        (X(1,m,t)-X(1,l,t))^2 + ...
                        (X(2,m,t)-X(2,l,t))^2 + ...
                        (X(3,m,t)-X(3,l,t))^2)/sigmaV2(t);
argout = stdV2(t)*exp(argin);  %% BUILT IN KERNEL derV, do not remove this comment
                    argout = argout * (mom(1,l,t)*mom(1,m,t)+mom(2,l,t)*mom(2,m,t)+mom(3,l,t)*mom(3,m,t));
                    dmom(1,m) = dmom(1,m) + argout * (X(1,m,t)-X(1,l,t));
                    dmom(2,m) = dmom(2,m) + argout * (X(2,m,t)-X(2,l,t));
                    dmom(3,m) = dmom(3,m) + argout * (X(3,m,t)-X(3,l,t));
                end
            end
        else
            X1 = X(1,:,t);
            X2 = X(2,:,t);
            X3 = X(3,:,t);
            mom1 = mom(1,:,t);
            mom2 = mom(2,:,t);
            mom3 = mom(3,:,t);

            if usefgt
                res = stdV2(t)*figtree(3,X(:,:,t),[mom1;mom2;mom3;...
                    mom1.*X1;mom1.*X2;mom1.*X3;...
                    mom2.*X1;mom2.*X2;mom2.*X3;...
                    mom3.*X1;mom3.*X2;mom3.*X3]...
                    ,X(:,:,t),sigmaV(t),epsilon);
            else
                res = stdV2(t)*gridOptim(X(:,:,t),[mom1;mom2;mom3;...
                    mom1.*X1;mom1.*X2;mom1.*X3;...
                    mom2.*X1;mom2.*X2;mom2.*X3;...
                    mom3.*X1;mom3.*X2;mom3.*X3]...
                    ,X(:,:,t),sourcegrid{t}.long,sourcegrid{t}.pas,sourcegrid{t}.origine,sourcegrid{t}.fft3k_d);
            end

            temp = mom1.*res(1,:);
            dmom(1,:) = dmom(1,:) + temp.*X1;
            dmom(2,:) = dmom(2,:) + temp.*X2;
            dmom(3,:) = dmom(3,:) + temp.*X3;
            temp = mom2.*res(2,:);
            dmom(1,:) = dmom(1,:) + temp.*X1;
            dmom(2,:) = dmom(2,:) + temp.*X2;
            dmom(3,:) = dmom(3,:) + temp.*X3;
            temp = mom3.*res(3,:);
            dmom(1,:) = dmom(1,:) + temp.*X1;
            dmom(2,:) = dmom(2,:) + temp.*X2;
            dmom(3,:) = dmom(3,:) + temp.*X3;

            temp = mom1.*res(4,:);
            dmom(1,:) = dmom(1,:) - temp;
            temp = mom1.*res(5,:);
            dmom(2,:) = dmom(2,:) - temp;
            temp = mom1.*res(6,:);
            dmom(3,:) = dmom(3,:) - temp;

            temp = mom2.*res(7,:);
            dmom(1,:) = dmom(1,:) - temp;
            temp = mom2.*res(8,:);
            dmom(2,:) = dmom(2,:) - temp;
            temp = mom2.*res(9,:);
            dmom(3,:) = dmom(3,:) - temp;

            temp = mom3.*res(10,:);
            dmom(1,:) = dmom(1,:) - temp;
            temp = mom3.*res(11,:);
            dmom(2,:) = dmom(2,:) - temp;
            temp = mom3.*res(12,:);
            dmom(3,:) = dmom(3,:) - temp;

        end

    end



 function grille = setgrid(mini,maxi,tp,grille)
         pas = ratio * sigmaV(tp); %grid's step
         long = (maxi-mini)/pas + 3/ratio; %circonf??rence du tore
%           long = (long<=32)*32 + (long>32).*((long<=64)*64 + (long>64).*((long<=128)*128 + (long>128)*2.*ceil(long/2)));
         long = (long<=16)*16 + (long>16).*((long<=32)*32 + (long>32).*((long<=64)*64 + (long>64).*2.*ceil(long/2)));
         grille.pas = pas;
         grille.origine = mini - (long*pas-maxi+mini)/2;
         %mise a jour du noyau si necessaire
         if ~(isfield(grille,'long')&&(isequal(grille.long,long')))
            grille.long = long';
            grille.fft3k_d = noyau3D_PAIR(grille,sigmaV(tp));
         end
 end

 function p = changegrid(mini,maxi,t,grille)
         p = sum( (maxi > (grille.origine + grille.pas*grille.long' - sigmaV(t))) | (mini < (grille.origine + sigmaV(t))) ) ~= 0;
         if (~p)
           long = (maxi-mini)/grille.pas + 3/ratio; %circonference du tore
%             long = (long<=32)*32 + (long>32).*((long<=64)*64 + (long>64).*((long<=128)*128 + (long>128)*2.*ceil(long/2)));
           long = (long<=16)*16 + (long>16).*((long<=32)*32 + (long>32).*((long<=64)*64 + (long>64).*2.*ceil(long/2)));
           p = sum(long < grille.long') ~= 0;
         end
 end

    function updatesourcegrid(tp)
        mini = min(X(:,:,tp),[],2);
        maxi = max(X(:,:,tp),[],2);
        if (changegrid(mini,maxi,tp,sourcegrid{tp}))
            sourcegrid{tp} = setgrid(mini,maxi,tp,sourcegrid{tp});
            disp(['source''s grid ' num2str(tp) ' has been changed:  ' num2str([sourcegrid{tp}.pas sourcegrid{tp}.long sourcegrid{tp}.origine'])]);
        end
    end


end
