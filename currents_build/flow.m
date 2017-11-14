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


function [p traj] = flow(s,p,t1,t2)

% [p traj] = flow(s,p)
% flow points p with deformation map given in s
% p must be 3*N
% traj is the points trajectories along time
% s is the returned value of a call to match function
%
% p = flow(s,p,t1,t2) computes from time steps t1 to t2
% where t1 and t2 are integers between 0 and s.T
% t1=0 means we compute rigid transformation first
%
% p = flow(s,p,-1) is same as p = flow(s,p,s.T,0) (inverse flow)


kerV = 'exp(argin)';  %% BUILT IN KERNEL FLAG kerV, do not remove this comment
if ~strcmp(s.kerV,kerV)
    error('Matching was not performed with the current built-in kernel.')
end


X = 0;
mom = 0;
transmatrix = eye(3);
transvector = zeros(3,1);
T = 30;
sigmaV = 0;
sigmaV2 = 0;
normcoefV = ones(1,T);
stdV = ones(1,T);
stdV2 = stdV.^2;
tau = 0;


usefgt = 0;
% order = 0;
% K = 0;
% e = 0;
epsilon = 0;        % IFGT required precision
                    %replace order, K and e by single variable epsilon for figtree
usegrid = 0;
ratio = 0;
optim_verbosemode = 0;

loadstruct('s',who)

nx = size(X,2);

if (usefgt&&usegrid)
    error('fgt ou grille : il faut choisir !');
end



if ~strcmp(kerV,'exp(argin)') && usefgt
    warning('Current built-in kernel is not gaussian. Fast Gauss Transform disabled')
    usefgt = 0;
end

if nargin == 2
    t1 = 0;
    t2 = T;
end

if nargin == 3 && t1 == -1
        t1 = T;
        t2 = 0;
end

prefix = 'flow step t=';

N = size(p,2);

if t1 == 0
    if optim_verbosemode
        disp('performs affine transformation')
    end
    % direct affine transformation
    p = transmatrix * p + repmat(transvector,1,size(p,2));
    t1 = 1;
end

if t2 == 0
    doinvaffine = 1;
    t2 = 1;
else
    doinvaffine = 0;
end

if t1 > t2
    doinvflow = 1;
    X = X(:,:,T:-1:1);
    mom = - mom(:,:,T:-1:1);
    t1 = T-t1+1;
    t2 = T-t2+1;
else
    doinvflow = 0;
end

dp = zeros(size(p));

if optim_verbosemode
    if usefgt
        disp('use fgt for flow computation')
    elseif usegrid
        disp('use grids for flow computation');
    end
end

if usegrid
    sourcegrid = cell(1,T);
    for t=1:T
       mini = min([X(:,:,t) p],[],2);
       maxi = max([X(:,:,t) p],[],2);
       sourcegrid{t} = setgrid(mini,maxi,t);
    end
end

if (nargout == 2)
    traj = zeros([size(p), t2-t1+1]);
    traj(:,:,1) = p;
end

for t = t1:t2-1
    if optim_verbosemode
        if doinvflow
            disp(['flow step t=',num2str(T-t+1),' to ',num2str(T-t)])
        else
            disp(['flow step t=',num2str(t),' to ',num2str(t+1)])
        end
    end
    if usefgt
        dp = stdV2(t)*figtree(3,X(:,:,t),mom(:,:,t),p,sigmaV(t),epsilon);
    elseif usegrid
        dp = stdV2(t)*gridOptim(X(:,:,t),mom(:,:,t),p,sourcegrid{t}.long,sourcegrid{t}.pas,sourcegrid{t}.origine,sourcegrid{t}.fft3k_d);
    else
        dp(:) = 0;
        for l = 1:N
            for k = 1:nx
                argin = -( ...
                    (X(1,k,t)-p(1,l))^2 + ...
                    (X(2,k,t)-p(2,l))^2 + ...
                    (X(3,k,t)-p(3,l))^2)/sigmaV2(t);
argout = stdV2(t)*exp(argin);  %% BUILT IN KERNEL kerV, do not remove this comment
                dp(1,l) = dp(1,l) + argout*mom(1,k,t);
                dp(2,l) = dp(2,l) + argout*mom(2,k,t);
                dp(3,l) = dp(3,l) + argout*mom(3,k,t);
            end
        end
    end
    if usegrid
        updatesourcegrids(min([X(:,:,t+1),(p+tau*dp)],[],2),max([X(:,:,t+1), (p+tau*dp)],[],2),t+1);
    end

    %%%% BEGIN corrector scheme
    %%%% comment this block to use simple Euler scheme
    if usefgt
        dp = dp + stdV2(t+1)*figtree(3,X(:,:,t+1),mom(:,:,t+1),p+tau*dp,sigmaV(t+1),epsilon);
    elseif usegrid
        dp = dp + stdV2(t+1)*gridOptim(X(:,:,t+1),mom(:,:,t+1),p+tau*dp,sourcegrid{t+1}.long,sourcegrid{t+1}.pas,sourcegrid{t+1}.origine,sourcegrid{t+1}.fft3k_d);
    else
        ptemp = p + tau * dp;
        for l = 1:N
            for k = 1:nx
                argin = -( ...
                    (X(1,k,t+1)-ptemp(1,l))^2 + ...
                    (X(2,k,t+1)-ptemp(2,l))^2 + ...
                    (X(3,k,t+1)-ptemp(3,l))^2)/sigmaV2(t+1);
argout = stdV2(t+1)*exp(argin);  %% BUILT IN KERNEL kerV, do not remove this comment
                dp(1,l) = dp(1,l) + argout*mom(1,k,t+1);
                dp(2,l) = dp(2,l) + argout*mom(2,k,t+1);
                dp(3,l) = dp(3,l) + argout*mom(3,k,t+1);
            end
        end
    end
    dp = dp / 2;
    %%%% END corrector scheme
    p = p + tau * normcoefV(t) * dp;
    if (nargout == 2)
        traj(:,:,t+1) = p;
    end
    if usegrid
        updatesourcegrids(min([X(:,:,t+1) p],[],2),max([X(:,:,t+1) p],[],2),t+1);
    end
end

if doinvaffine
    if optim_verbosemode
        disp('performs inverse affine transformation')
    end
    % inverse affine transformation
    p = inv(transmatrix) * (p - repmat(transvector,1,size(p,2)));
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
         if optim_verbosemode
            disp(['source''s grid ' num2str(tp) ' has been set:  ' num2str([grille.pas grille.long grille.origine'])]);
         end
 end


    function updatesourcegrids(mini,maxi,tp)
        grille = sourcegrid{tp};
        paux = sum( (maxi > (grille.origine + grille.pas*grille.long' - sigmaV(t))) | (mini < (grille.origine + sigmaV(t))) ) ~= 0;
        if (~paux)
          long = (maxi-mini)/grille.pas + 3/ratio; %circonf??rence du tore
          long = (long<=16)*16 + (long>16).*((long<=32)*32 + (long>32).*((long<=64)*64 + (long>64).*2.*ceil(long/2)));
          paux = sum(long < grille.long') ~= 0;
        end
        if (paux)
          pas = ratio * sigmaV(tp); %grid's step
          long = (maxi-mini)/pas + 3/ratio; %circonf??rence du tore
          long = (long<=16)*16 + (long>16).*((long<=32)*32 + (long>32).*((long<=64)*64 + (long>64).*2.*ceil(long/2)));
          grille.pas = pas;
          grille.origine = mini - (long*pas-maxi+mini)/2;
          if ~(isfield(grille,'long')&&(isequal(grille.long,long')))
             grille.long = long';
             grille.fft3k_d = noyau3D_PAIR(grille,sigmaV(tp));
          end
          sourcegrid{tp} = grille;
          if optim_verbosemode
             disp(['source''s grid ' num2str(tp) ' has been changed:  ' num2str([sourcegrid{tp}.pas sourcegrid{tp}.long sourcegrid{tp}.origine'])]);
          end
        end
    end


end
