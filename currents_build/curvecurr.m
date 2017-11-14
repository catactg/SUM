% 'exoShape' is a derived software by Stanley Durrleman, Copyright (C) INRIA (Asclepios team), All Rights Reserved, 2006-2009, version 1.0
%--------------------------------------------------------------------------
% Based on MATCHINE v1.0 software.
% Copyright Universit� Paris Descartes
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
% - J. Glaunes, A. Qiu, M. Miller, L. Younes, Large Deformation Diffeomorphic Metric Curve Mapping, International Journal of Computer vision, Springer, 2008, Vol. 80, No. 3, pp. 317--336, DOI: 10.1007/s11263-008-0141-9 

function s = curvecurr(s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Curve matching via currents - to be used with match.m %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% contact : joan_glaunes@yahoo.fr

method = 'curvecurr';

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

% vx and vy (template and target faces) are required components
if isfield(s,'vx')~=1 | isfield(s,'vy')~=1
    error('input structure must contain fields vx and vy (template and target faces)');
end


vx = s.vx;
vy = s.vy;
nfx = size(vx,2);         % number of template faces
nfy = size(vy,2);         % number of target faces
wx = ones(1,nfx);                % weight of template faces
wy = ones(1,nfy);                 % weight of target faces


numbminims = 1;
sigmaW = 'auto';            % matching kernel size
decsigmaW = 2.5;            % decrease rate of sigmaW after each minimization

usefgt = 0;                 % if 1 use IFGT algorithm
testfgt = 0;
% order = 5;                 % IFGT order
% K = min([10,floor(ny/2)]);  % IFGT number of clusters
% e = 9;                      % IFGT cutoff ratio for targets
epsilon = 1e-3;             % IFGT required precision
                            %replace order, K and e by single variable epsilon for figtree

usegrid = 0;                % if 1 use grid optimisation
testgrid = 0;
ratio = 0.2;

% update by input structure
loadstruct('s',who)

%check error
if (usefgt && usegrid)
  error('fgt ou grille : il faut choisir !');
end


% auto scale
sclsz = mean(max(y')-min(y'));
if strcmp(sigmaW,'auto')
    sigmaW = .5*sclsz;
end
sigmaW2 = sigmaW^2;

% computation of centers and normals for each face of target
[cy,tauy] = compcurr(y,vy,wy);

if usegrid
   premin = @fungridpremin;
else
premin = @funpremin;
end
setgrid = @funsetgrid;
changegrid = @funchangegrid;

matching = @funmatching;
gradmatching = @fungradmatching;

clear s
savestruct('s',who);

% multiscale settings
sigmaW = sigmaW*decsigmaW^numbminims;   %%% numbminims should be given by match.m
sigmaW2 = sigmaW^2;

if usefgt
    WdScalProd = @fgtWdScalProd;
    gradWdNorm = @fgtgradWdNorm;
elseif usegrid
    WdScalProd = @gridWdScalProd;
    gradWdNorm = @gridgradWdNorm;
else
    WdScalProd = @regWdScalProd;
    gradWdNorm = @reggradWdNorm;
end

% process this before each minimization
    function funpremin
        sigmaW = sigmaW/decsigmaW;
        sigmaW2 = sigmaW^2;
        if usefgt*testfgt
            essdir = regWdScalProd(cy,tauy,cy,tauy,sigmaW2);
            essfgt = fgtWdScalProd(cy,tauy,cy,tauy,sigmaW2);
            disp(['relative error fgt curvecurr: ',num2str(abs((essdir-essfgt)/essdir))])
        end
    end

    function grille = fungridpremin(maxix,minix)
        sigmaW = sigmaW/decsigmaW;
        sigmaW2 = sigmaW^2;
        grille = setgrid(maxix,minix);
        if testgrid
            essdir = regWdScalProd(cy,tauy,cy,tauy,sigmaW2);
            essgrid = gridWdScalProd(cy,tauy,cy,tauy,grille);
            disp(['relative error grid curvecurr: ', num2str(abs(essdir-essgrid)/essdir)])
        end
    end

    function grille = funsetgrid(minix, maxix,grille)
         mini = min(min(y')',minix);
         maxi = max(max(y')',maxix);
         pas = ratio * sigmaW; %grid's step
         long = (maxi-mini)/pas + 3/ratio; %circonférence du tore
         long = (long<=32)*32 + (long>32).*((long<=64)*64 + (long>64).*((long<=128)*128 + (long>128)*2.*ceil(long/2)));
%           long = (long<=32)*32 + (long>32).*((long<=64)*64 + (long>64).*2.*ceil(long/2));
         grille.pas = pas;
         grille.origine = mini - (long*pas-maxi+mini)/2;

         %mise a jour du noyau si necessaire
         if ~(isfield(grille,'long')&&(isequal(grille.long,long')))
            grille.long = long';
            grille.fft3k_d = noyau3D_PAIR(grille,sigmaW);
         end
    end

    function p = funchangegrid(minix,maxix,grille)
         mini = min(min(y')',minix);
         maxi = max(max(y')',maxix);
         p = sum( (maxi > (grille.origine + grille.pas*grille.long' - sigmaW)) | (mini < (grille.origine + sigmaW)) ) ~= 0;
         if (~p)
           long = (maxi-mini)/grille.pas + 3/ratio; %circonférence du tore
           long = (long<=32)*32 + (long>32).*((long<=64)*64 + (long>64).*((long<=128)*128 + (long>128)*2.*ceil(long/2)));
%             long = (long<=32)*32 + (long>32).*((long<=64)*64 + (long>64).*2.*ceil(long/2));
           p = sum(long < grille.long') ~= 0;
         end
    end


    function R = funmatching(phix,grille)
        % Function:  compute energy
        %
        % Input:
        %        phix     coordinates of vertices on the first curve
        % Output:
        %        R        total current energy
        [cphix,tauphix] = compcurr(phix,vx,wx); % computes centers and normals for each face
        if usegrid
            R = WdScalProd([cphix,cy],[tauphix,-tauy],[cphix,cy],[tauphix,-tauy],grille);
        else
            R = WdScalProd([cphix,cy],[tauphix,-tauy],[cphix,cy],[tauphix,-tauy],sigmaW2);
        end
    end


    function g = fungradmatching(phix,grille)
        [cphix,tauphix] = compcurr(phix,vx,wx); % computation of centers and normals for each face
        if usegrid
           g = gradWdNorm(phix,cphix,tauphix,grille);
        else
           g = gradWdNorm(phix,cphix,tauphix,sigmaW2);
        end
    end


    function [ca,taua] = compcurr(a,ba,wa)
        %Function:
        %      computer segments' centers and tangents
        %
        % Input:
        %     a    coordinates of vertices, matrix  3*number of vertices
        %     ba   segments'connectivity, matrix 3 * number of segments
        % Output:
        %     ca   coordinates of center for each segment, array  3*number of
        %     segments
        %     taua   tangent directions for each segment, array 3*number of segments
        %
        % Variables:
        %      na   number of vertices
        %      nfa  number of segments

        nfa = size(ba,2);
        ca = zeros(3,nfa);
        taua = zeros(3,nfa);

        % computation of centers and normals for each face
        for f=1:nfa,
            b1 = 3*(ba(1,f)-1);
            b2 = 3*(ba(2,f)-1);
            ca(1,f) = (a(1+b1)+a(1+b2))/2;
            ca(2,f) = (a(2+b1)+a(2+b2))/2;
            ca(3,f) = (a(3+b1)+a(3+b2))/2;
            taua(1,f) = wa(f)*(a(1+b2)-a(1+b1));
            taua(2,f) = wa(f)*(a(2+b2)-a(2+b1));
            taua(3,f) = wa(f)*(a(3+b2)-a(3+b1));
        end

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Scalar product in W^* space %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function:
%             compute currents in W* space
% Input
%         ca      center of each triangle on the first surface
%         taua      normal of each triangle on the first surface
%         cb      center of each triangle on the second surface
%         taub      normal of each triangle on the second surface
% Output
%        R        currnets
%

    function R = regWdScalProd(ca,taua,cb,taub,sigmaW2)
        nfa = size(ca,2);
        nfb = size(cb,2);
        R = 0;
        for k = 1:nfa
            for j = 1:nfb
                argin = -((ca(1,k)-cb(1,j))^2+(ca(2,k)-cb(2,j))^2+(ca(3,k)-cb(3,j))^2)/sigmaW2;
argout = exp(argin);  %% BUILT IN KERNEL kerI, do not remove this comment
                R = R + argout * (taua(1,k)*taub(1,j)+taua(2,k)*taub(2,j)+taua(3,k)*taub(3,j));
            end
        end
    end


    function R = fgtWdScalProd(ca,taua,cb,taub,sigmaW2)
        R = sum(sum(taub.*figtree(3,ca,taua,cb,sqrt(sigmaW2),epsilon)));
    end


    function R = gridWdScalProd(ca,taua,cb,taub,grille)
        R = sum(sum(taub.*gridOptim(ca,taua,cb,grille.long,grille.pas,grille.origine,grille.fft3k_d)));
    end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% gradient of W^* norm %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function:  gradient of data attachment (curve current)
%
% Input
%        phi     coordinates of vertices on the first surface
%        cphi    coordiantes of the center for each triangle on the first
%        surface
%        tauphi    normals for each triangle on the first surface
%
% Output:
%        eta  gradient of current with respect to the coordinates
%        of vertices
%
% Variables
%        etaf_w       gradient on normals with respect to the coordinates
%        of vertices
%        etaf_dw      gradient on kernel with respect to the coordinates of
%        vertices
%

    function eta = reggradWdNorm(phi,cphi,tauphi,sigmaW2)
        eta = zeros(size(phi));
        etaloc = zeros(3,2);
        etaf_w= zeros(3,1);
        etaf_dw = zeros(3,1);

        for e = 1:nfx
            loce = 3*(e-1);
            etaf_w(1) = 0;
            etaf_w(2) = 0;
            etaf_w(3) = 0;
            etaf_dw(1) = 0;
            etaf_dw(2) = 0;
            etaf_dw(3) = 0;
            for f = 1:nfx
                locf = 3*(f-1);
                argin = -(...
                    (cphi(1+locf)-cphi(1+loce))^2 + ...
                    (cphi(2+locf)-cphi(2+loce))^2 + ...
                    (cphi(3+locf)-cphi(3+loce))^2 )/sigmaW2;
argout = exp(argin);  %% BUILT IN KERNEL kerI, do not remove this comment
                etaf_w(1) = etaf_w(1) + argout * tauphi(1+locf);
                etaf_w(2) = etaf_w(2) + argout * tauphi(2+locf);
                etaf_w(3) = etaf_w(3) + argout * tauphi(3+locf);
argout = exp(argin);  %% BUILT IN KERNEL derI, do not remove this comment
                dottau = tauphi(1+locf)*tauphi(1+loce) + tauphi(2+locf)*tauphi(2+loce) + tauphi(3+locf)*tauphi(3+loce);
                etaf_dw(1) = etaf_dw(1) + argout * dottau * (cphi(1+locf)-cphi(1+loce));
                etaf_dw(2) = etaf_dw(2) + argout * dottau * (cphi(2+locf)-cphi(2+loce));
                etaf_dw(3) = etaf_dw(3) + argout * dottau * (cphi(3+locf)-cphi(3+loce));
            end
            for f = 1:nfy
                locf = 3*(f-1);
                argin = -(...
                    (cy(1+locf)-cphi(1+loce))^2 + ...
                    (cy(2+locf)-cphi(2+loce))^2 + ...
                    (cy(3+locf)-cphi(3+loce))^2 )/sigmaW2;
argout = exp(argin);  %% BUILT IN KERNEL kerI, do not remove this comment
                etaf_w(1) = etaf_w(1) - argout * tauy(1+locf);
                etaf_w(2) = etaf_w(2) - argout * tauy(2+locf);
                etaf_w(3) = etaf_w(3) - argout * tauy(3+locf);
argout = exp(argin);  %% BUILT IN KERNEL derI, do not remove this comment
                dottau = tauy(1+locf)*tauphi(1+loce) + tauy(2+locf)*tauphi(2+loce) + tauy(3+locf)*tauphi(3+loce);
                etaf_dw(1) = etaf_dw(1) - argout * dottau * (cy(1+locf)-cphi(1+loce));
                etaf_dw(2) = etaf_dw(2) - argout * dottau * (cy(2+locf)-cphi(2+loce));
                etaf_dw(3) = etaf_dw(3) - argout * dottau * (cy(3+locf)-cphi(3+loce));
            end
            etaf_dw(1) = 2 * etaf_dw(1) / sigmaW2;
            etaf_dw(2) = 2 * etaf_dw(2) / sigmaW2;
            etaf_dw(3) = 2 * etaf_dw(3) / sigmaW2;

            etaloc(1) = -2*wx(e)*etaf_w(1) + etaf_dw(1);
            etaloc(2) = -2*wx(e)*etaf_w(2) + etaf_dw(2);
            etaloc(3) = -2*wx(e)*etaf_w(3) + etaf_dw(3);
            etaloc(4) = 2*wx(e)*etaf_w(1) + etaf_dw(1);
            etaloc(5) = 2*wx(e)*etaf_w(2) + etaf_dw(2);
            etaloc(6) = 2*wx(e)*etaf_w(3) + etaf_dw(3);

            loc1 = 3*(vx(1+2*(e-1))-1);
            loc2 = 3*(vx(2+2*(e-1))-1);
            eta(1+loc1) = eta(1+loc1) + etaloc(1);
            eta(2+loc1) = eta(2+loc1) + etaloc(2);
            eta(3+loc1) = eta(3+loc1) + etaloc(3);
            eta(1+loc2) = eta(1+loc2) + etaloc(4);
            eta(2+loc2) = eta(2+loc2) + etaloc(5);
            eta(3+loc2) = eta(3+loc2) + etaloc(6);
        end

    end

    function eta = fgtgradWdNorm(phi,cphi,tauphi,sigmaW2)
        cz = [cphi,cy];
        tauz = [tauphi,-tauy];

        resfgt = figtree(3,cz,[tauz;specprod(cz,tauz,1)],cphi,sigmaW,epsilon);

        eta_dw = (2/sigmaW2) * (specprod(resfgt(4:end,:),tauphi,3) - specprod(specprod(cphi,tauphi,1),resfgt(1:3,:),3));
        eta_w = 2*resfgt(1:3,:);

        eta = zeros(size(phi));
        for m = 1:nfx
            locm = 3*m-2:3*m;
            bm1 = 3*vx(1,m);
            locbm1 = bm1-2:bm1;
            eta(locbm1) = eta(locbm1)' + eta_dw(locm) - wx(m)*eta_w(locm);
            bm2 = 3*vx(2,m);
            locbm2 = bm2-2:bm2;
            eta(locbm2) = eta(locbm2)' + eta_dw(locm) + wx(m)*eta_w(locm);
        end

    end


    function eta = gridgradWdNorm(phi,cphi,tauphi,grille)

      cz = [cphi,cy];
      tauz = [tauphi,-tauy];

      resgrid = gridOptim(cz,[tauz; specprod(cz,tauz,1)],cphi,grille.long,grille.pas,grille.origine,grille.fft3k_d);

      eta_dw = (2/sigmaW2) * (specprod(resgrid(4:end,:),tauphi,3) - specprod(specprod(cphi,tauphi,1),resgrid(1:3,:),3));
      eta_w = 2*resgrid(1:3,:);

      eta = zeros(size(phi));
      for m = 1:nfx
            locm = 3*m-2:3*m;
            bm1 = 3*vx(1,m);
            locbm1 = bm1-2:bm1;
            eta(locbm1) = eta(locbm1)' + eta_dw(locm) - wx(m)*eta_w(locm);
            bm2 = 3*vx(2,m);
            locbm2 = bm2-2:bm2;
            eta(locbm2) = eta(locbm2)' + eta_dw(locm) + wx(m)*eta_w(locm);
      end
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
