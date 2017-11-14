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
% - J. Glaunes, L. Younes and A. Trouve, Diffeomorphic matching of distributions: A new approach for unlabelled point-sets and sub-manifolds matching, Proc. of the 2004 IEEE Computer Society Conference on Computer Vision and Pattern Recognition. (CVPR'04), Vol. 2, pp. 712--718, DOI: 10.1109/CVPR.2004.81 


function s = measures(s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Measure matching - to be used with match.m     %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
method = 'measures';

timept = 1;
% y (target points) is required
if ~isfield(s,'y')
    error('input structure must contain field <y> (target points)');
end

% y must be 3D data
if size(s.y,1)~=3
    error('field <y> must be 3-rows matrix');
end

y = s.y;
ny = size(y,2);

if isfield(s,'x')~=1
    if isfield(s,'nx')~=1
        if isfield(s,'wx')~=1
            if isfield(s,'vx')~=1
                error('at least one of fields x, nx, wx, vx should be given')
            else
                nx = length(s.vx);
            end
        else
            nx = length(s.wx);
        end
    else
        nx = s.nx;
    end
else
    nx = size(s.x,2);
end

if isfield(s,'wx')~=1
    wx = ones(1,nx)/nx;
else
    wx = s.wx;
end

if isfield(s,'wy')~=1
    wy = ones(1,ny)/ny;
else
    wy = s.wy;
end

if isfield(s,'vx')~=1
    vx = 1:nx;
else
    vx = s.vx;
end

numbminims = 1;
sigmaI = 'auto';            % matching kernel size
decsigmaI = 2.5;            % decrease rate of sigmaI after each minimization

usefgt = 0;                 % if 1 use IFGT algorithm
testfgt = 0;
% order = 3;                  % IFGT order
% K = min([10,floor((nx+ny)/2)]);  % IFGT number of clusters
% e = 9;                      % IFGT cutoff ratio for targets
epsilon = 1e-3;             % IFGT required precision
                           %replace order, K and e by single variable epsilon for figtree

usegrid = 0;               % if 1 use optimisation with grid
testgrid = 0;
ratio = 0.2;               % ratio grid's step / sigmaW : the smaller, the smaller approximation's error

% update by input structure
loadstruct('s',who)
%check error
if (usefgt && usegrid)
  error('fgt ou grille : il faut choisir !');
end


% auto scale
sclsz = mean(max(y')-min(y'));
if strcmp(sigmaI,'auto')
    sigmaI = .5*sclsz;
end
sigmaI2 = sigmaI^2;

if usegrid
   premin = @fungridpremin;
else
premin = @funpremin;
end
setgrid = @funsetgrid;
changegrid = @funchangegrid;

matching = @funmatching;
gradmatching = @fungradmatching;
gradmatchingtime = @fungradmatchingtime;


if usefgt
    IdScalProd = @fgtIdScalProd;
    fungrad = @fgtgradmatching;
elseif usegrid
    IdScalProd = @gridIdScalProd;
    fungrad = @gridgradmatching;
else
    IdScalProd = @regIdScalProd;
    fungrad = @reggradmatching;
end

clear s
savestruct('s',who);

wz = [wx,-wy];

vxspec = 3 * vx;
vxspec = [vxspec-2;vxspec-1;vxspec];
vxspec= vxspec(:);

% multiscale settings
sigmaI = sigmaI*decsigmaI^numbminims;  %%% numbminims should be given by match.m
sigmaI2 = sigmaI^2;

% process this before each minimization
    function funpremin
        sigmaI = sigmaI/decsigmaI;
        sigmaI2 = sigmaI^2;
        if usefgt*testfgt
            essdir = regIdScalProd(y,wy,y,wy,sigmaI2);
            essfgt = fgtIdScalProd(y,wy,y,wy,sigmaI2);
            disp(['relative error fgt measures: ',num2str(abs((essdir-essfgt)/essdir))])
        end
    end

    function grille = fungridpremin(maxix,minix)
        sigmaI = sigmaI/decsigmaI;
        sigmaI2 = sigmaI^2;
        grille = setgrid(maxix,minix);
        if testgrid
            essdir = regIdScalProd(y,wy,y,wy,sigmaI2);
            essgrid = gridIdScalProd(y,wy,y,wy,grille);
            disp(['relative error grid surfcurr: ', num2str(abs(essdir-essgrid)/essdir)])
        end
    end


    function grille = funsetgrid(minix, maxix,grille)
         mini = min(min(y')',minix);
         maxi = max(max(y')',maxix);
         pas = ratio * sigmaI; %grid's step
         long = (maxi-mini)/pas + 10/ratio; %circonference du tore
         long = (long<=16)*16 + (long>16).*((long<=32)*32 + (long>32).*((long<=64)*64 + (long>64).*2.*ceil((long+1)/2)));
         grille.pas = pas;
         grille.origine = mini - (long*pas-maxi+mini)/2;
         %mise a jour du noyau si necessaire
         if ~(isfield(grille,'long')&&(isequal(grille.long,long')))
            grille.long = long';
            grille.fft3k_d = noyau3D_PAIR(grille,sigmaI);
            disp(['mise a jour du noyau : ' num2str(long(1)) ' ' num2str(long(2)) ' ' num2str(long(3))]);
         end
    end

    function p = funchangegrid(minix,maxix,grille)
         mini = min(min(y')',minix);
         maxi = max(max(y')',maxix);
         p = sum( (maxi > (grille.origine + grille.pas*grille.long' - sigmaI)) | (mini < (grille.origine + sigmaI)) ) ~= 0;
         if (~p)
           long = (maxi-mini)/grille.pas + 3/ratio; %circonference du tore
           long = (long<=16)*16 + (long>16).*((long<=32)*32 + (long>32).*((long<=64)*64 + (long>64).*2.*ceil(long/2)));
           p = sum(long < grille.long') ~= 0;
         end
    end


    function R = funmatching(phix,grille)
        phix = reshape(phix,3,prod(size(phix))/3);        
        z = [phix(:,vx),y];        
        if usegrid
            R = IdScalProd(z,wz,z,wz,grille);
        else
            R = IdScalProd(z,wz,z,wz,sigmaI2);
        end
    end


	function eta = fungradmatchingtime(phix,phixP,grille)
	 	z = [phix(:,vx),y];
		dz = [phixP(:,vx),phix(:,vx)];
		wd = [wx,-wx];
		if usegrid
			eta = IdScalProd(z,wz,dz,wd,grille);
		else
			eta = IdScalProd(z,wz,dz,wd,sigmaI2);
		end
	end


    function R = regIdScalProd(a,wa,b,wb,sigmaI2)
        R = 0;
        na = size(a,2);
        nb = size(b,2);
        for k = 1:na
            Rk = 0;
            for j = 1:nb
                argin = -((a(1,k)-b(1,j))^2+(a(2,k)-b(2,j))^2+(a(3,k)-b(3,j))^2)/sigmaI2;
argout = exp(argin);  %% BUILT IN KERNEL kerI, do not remove this comment
                Rk = Rk + wb(j) * argout;
            end
            R = R + wa(k) * Rk;
        end
    end

    function R = fgtIdScalProd(a,wa,b,wb,sigmaI2)
        resfgt = figtree(3,a,wa,b,sqrt(sigmaI2),epsilon);
        R = sum(wb.*resfgt);
    end


    function R = gridIdScalProd(a,wa,b,wb,grille)
        resgrid = gridOptim(a,wa,b,grille.long,grille.pas,grille.origine,grille.fft3k_d);
        R = sum(wb.*resgrid);
    end


    function eta = fungradmatching(phix,grille)
        if usegrid
            eta = fungrad(phix,grille);
        else
            eta = fungrad(phix,sigmaI2);
        end
    end


    function eta = reggradmatching(phix,sigmaI2)
        eta = zeros(size(phix));
        for ee = 1:nx
            loce = 3*(vx(ee)-1);
            eta1loc = 0;
            eta2loc = 0;
            eta3loc = 0;
            for f = 1:nx
                locf = 3*(vx(f)-1);
                argin = -(...
                    (phix(1+locf)-phix(1+loce))^2 + ...
                    (phix(2+locf)-phix(2+loce))^2 + ...
                    (phix(3+locf)-phix(3+loce))^2 )/sigmaI2;
argout = exp(argin);  %% BUILT IN KERNEL derI, do not remove this comment
                eta1loc = eta1loc + argout * wx(f) * (phix(1+locf)-phix(1+loce));
                eta2loc = eta2loc + argout * wx(f) * (phix(2+locf)-phix(2+loce));
                eta3loc = eta3loc + argout * wx(f) * (phix(3+locf)-phix(3+loce));
            end
            for f = 1:ny
                locf = 3*(f-1);
                argin = -(...
                    (y(1+locf)-phix(1+loce))^2 + ...
                    (y(2+locf)-phix(2+loce))^2 + ...
                    (y(3+locf)-phix(3+loce))^2 )/sigmaI2;
argout = exp(argin);  %% BUILT IN KERNEL derI, do not remove this comment
                eta1loc = eta1loc - argout * wy(f) * (y(1+locf)-phix(1+loce));
                eta2loc = eta2loc - argout * wy(f) * (y(2+locf)-phix(2+loce));
                eta3loc = eta3loc - argout * wy(f) * (y(3+locf)-phix(3+loce));
            end
            eta(1+loce) = 4 * wx(ee) * eta1loc / sigmaI2;
            eta(2+loce) = 4 * wx(ee) * eta2loc / sigmaI2;
            eta(3+loce) = 4 * wx(ee) * eta3loc / sigmaI2;
        end
    end

    function eta = fgtgradmatching(phix,sigmaI2)
        sz = size(phix);
        phix = phix(vxspec);
        phix = reshape(phix,3,length(vx));
        z = [phix,y];
        resfgt = figtree(3,z,[wz;specprod(z,wz,1)],phix,sigmaI,epsilon);
        eta = zeros(3,prod(sz)/3);
        eta(:,vx) = (4/sigmaI2) * (specprod(resfgt(2:4,:),wx,1) - specprod(phix,wx.*resfgt(1,:),1));
        eta = reshape(eta,sz);
    end


    function eta = gridgradmatching(phix,grille)
      sz = size(phix);
      phix = phix(vxspec);
      phix = reshape(phix,3,length(vx));
      z = [phix,y];
      resgrid = gridOptim(z,[wz;specprod(z,wz,1)],phix,grille.long,grille.pas,grille.origine,grille.fft3k_d);
      eta = zeros(3,prod(sz)/3);
      eta(:,vx) = (4/sigmaI2) * (specprod(resgrid(2:4,:),wx,1) - specprod(phix,wx.*resgrid(1,:),1));
      eta = reshape(eta,sz);
    end

end


 function alpha = projectionTriLineaire(grille,taux,cx)
 
    step = grille.pas;
    % projection trilineaire des moments aux points de la grille
    alpha = zeros(grille.long);
    nfx = size(cx,2);
    for k=1:nfx
      c0 = floor((cx(:,k) - grille.origine)/step) + 1;
      delta = cx(:,k) - (grille.origine + (c0-1)*step);
 
      rho000 = (step-delta(1)) * (step-delta(2)) * (step-delta(3));
      rho100 = delta(1)        * (step-delta(2)) * (step-delta(3));
      rho010 = (step-delta(1)) * delta(2)        * (step-delta(3));
      rho001 = (step-delta(1)) * (step-delta(2)) * delta(3);
      rho110 = delta(1)        * delta(2)        * (step-delta(3));
      rho011 = (step-delta(1)) * delta(2)        * delta(3);
      rho101 = delta(1)        * (step-delta(2)) * delta(3);
      rho111 = delta(1)        * delta(2)        * delta(3);
 
      alpha(c0(1),  c0(2),  c0(3)  ) = alpha(c0(1),  c0(2),  c0(3)  ) + rho000*taux(k);
      alpha(c0(1)+1,c0(2),  c0(3)  ) = alpha(c0(1)+1,c0(2),  c0(3)  ) + rho100*taux(k);
      alpha(c0(1),  c0(2)+1,c0(3)  ) = alpha(c0(1),  c0(2)+1,c0(3)  ) + rho010*taux(k);
      alpha(c0(1),  c0(2),  c0(3)+1) = alpha(c0(1),  c0(2),  c0(3)+1) + rho001*taux(k);
      alpha(c0(1)+1,c0(2)+1,c0(3)  ) = alpha(c0(1)+1,c0(2)+1,c0(3)  ) + rho110*taux(k);
      alpha(c0(1),  c0(2)+1,c0(3)+1) = alpha(c0(1),  c0(2)+1,c0(3)+1) + rho011*taux(k);
      alpha(c0(1)+1,c0(2),  c0(3)+1) = alpha(c0(1)+1,c0(2),  c0(3)+1) + rho101*taux(k);
      alpha(c0(1)+1,c0(2)+1,c0(3)+1) = alpha(c0(1)+1,c0(2)+1,c0(3)+1) + rho111*taux(k);
    end
 
 end



function C = specprod(A,B,m)

n = size(A,2);
a = size(A,1);
b = size(B,1);

A = reshape(A,a/m,m,n);
B = reshape(B,m,b/m,n);
C = zeros(a/m,b/m,n);

for p = 1:n
    C(:,:,p) = A(:,:,p)*B(:,:,p);
end

C = reshape(C,a*b/m/m,n);

end


