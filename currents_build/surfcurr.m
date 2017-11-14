function s = surfcurr(s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Surface matching via currents - to be used with match.m %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

method = 'surfcurr';

timept = 1;
% y (target points) is required
if ~isfield(s,'y')
    error('input structure must contain field <y> (target points)');
end

y = s.y;

% y must be 3D data
if size(y,1)~=3
    error('field <y> must be 3-rows matrix');
end

ny = size(y,2);

% vx and vy (template and target faces) are required components
if isfield(s,'vx')~=1 | isfield(s,'vy')~=1
    error('input structure must contain fields vx and vy (template and target faces)');
end

vx = s.vx;
vy = s.vy;
nfx = size(vx,2);         % number of template faces
nfy = size(vy,2);         % number of target faces


numbminims = 1;
sigmaW = 'auto';            % matching kernel size
decsigmaW = 2.5;            % decrease rate of sigmaW after each minimization

usefgt = 0;                 % if 1 use IFGT algorithm
testfgt = 0;
order = 5;                 % IFGT order
K = min([10,floor(ny/2)]);  % IFGT number of clusters
e = 9;                      % IFGT cutoff ratio for targets

usegrid = 0;               % if 1 use optimisation with grid
testgrid = 0;
ratio = 0.2;               % ratio grid's step / sigmaW : the smaller, the smaller approximation's error

wx = ones(1,nfx);
wy = ones(1,nfy);

% update by input structure
loadstruct('s',who)

% auto scale
sclsz = mean(max(y')-min(y'));
if strcmp(sigmaW,'auto')
    sigmaW = .5*sclsz;
end
sigmaW2 = sigmaW^2;

% computation of centers and normals for each face of target
[cy,Ny] = compcurr(y,vy,wy);


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

clear s
savestruct('s',who);
%check error
if (usefgt && usegrid)
  error('fgt ou grille : il faut choisir !');
end

% multiscale settings
sigmaW = sigmaW*decsigmaW^numbminims;
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
            essdir = regWdScalProd(cy,Ny,cy,Ny,sigmaW2);
            essfgt = fgtWdScalProd(cy,Ny,cy,Ny,sigmaW2);
            disp(['relative error fgt surfcurr: ',num2str(abs((essdir-essfgt)/essdir))])
        end
    end

    function grille = fungridpremin(maxix,minix)
        sigmaW = sigmaW/decsigmaW;
        sigmaW2 = sigmaW^2;
        grille = setgrid(maxix,minix);
        if testgrid
            essdir = regWdScalProd(cy,Ny,cy,Ny,sigmaW2);
            essgrid = gridWdScalProd(cy,Ny,cy,Ny,grille);
            disp(['relative error grid surfcurr: ', num2str(abs(essdir-essgrid)/essdir)])
        end
    end


    function grille = funsetgrid(minix, maxix,grille)
         mini = min(min(y')',minix);
         maxi = max(max(y')',maxix);
         pas = ratio * sigmaW; %grid's step
         long = (maxi-mini)/pas + 3/ratio; %circonférence du tore
%           long = (long<=32)*32 + (long>32).*((long<=64)*64 + (long>64).*((long<=128)*128 + (long>128)*2.*ceil(long/2)));
         long = (long<=32)*32 + (long>32).*((long<=64)*64 + (long>64).*2.*ceil(long/2));

         grille.pas = pas;
         grille.origine = mini - (long*pas-maxi+mini)/2;

         %mise a jour du noyau si necessaire
         if ~(isfield(grille,'long')&&(isequal(grille.long,long')))
            grille.long = long';
            grille.fft3k_d = noyau3D_PAIR(grille,sigmaW);
         end
%           disp([grille.long grille.origine' mini' maxi' min(y') max(y')]);
    end

    function p = funchangegrid(minix,maxix,grille)
         mini = min(min(y')',minix);
         maxi = max(max(y')',maxix);
         p = sum( (maxi > (grille.origine + grille.pas*grille.long' - sigmaW)) | (mini < (grille.origine + sigmaW)) ) ~= 0;
         if (~p)
           long = (maxi-mini)/grille.pas + 3/ratio; %circonférence du tore
%             long = (long<=32)*32 + (long>32).*((long<=64)*64 + (long>64).*((long<=128)*128 + (long>128)*2.*ceil(long/2)));
           long = (long<=32)*32 + (long>32).*((long<=64)*64 + (long>64).*2.*ceil(long/2));
           p = sum(long < grille.long') ~= 0;
         end
    end



    function R = funmatching(phix,grille)
        [cphix,Nphix] = compcurr(phix,vx,wx); % computes centers and normals for each face
        if usegrid
           R = WdScalProd([cphix,cy],[Nphix,-Ny],[cphix,cy],[Nphix,-Ny],grille);
        else
           R = WdScalProd([cphix,cy],[Nphix,-Ny],[cphix,cy],[Nphix,-Ny],sigmaW2);
        end
    end


    function g = fungradmatching(phix,grille)
        [cphix,Nphix] = compcurr(phix,vx,wx); % computation of centers and normals for each face
        if usegrid
           g = gradWdNorm(phix,cphix,Nphix,grille);
        else
           g = gradWdNorm(phix,cphix,Nphix,sigmaW2);
        end
    end

	function eta = fungradmatchingtime(phix,phixP,grille)
		[cphix Nphix] = compcurr(phix,vx,wx);
		[cphixP NphixP] = compcurr(phixP,vx,wx);
	 	if usegrid
			eta = WdScalProd([cphix,cy],[Nphix,-Ny],[cphixP cphix],[NphixP,-Nphix],grille);
		else
			eta = WdScalProd([cphix,cy],[Nphix,-Ny],[cphixP cphix],[NphixP,-Nphix],sigmaW2);
		end
	end


    function [ca,Na] = compcurr(a,va,wa)

        nfa = size(va,2);
        ca = zeros(3,nfa);
        Na = zeros(3,nfa);

        % computation of centers and normals for each face
        v = zeros(1,9);
        for f = 1:nfa
            locf = 3*(f-1);
            for k = 1:3
                for j = 1:3
                    v(k+3*(j-1)) = a(k+3*(va(j+locf)-1));
                end
            end
            % c = (v1+v2+v3)/3;
            ca(1+locf) = (v(1)+v(4)+v(7))/3;
            ca(2+locf) = (v(2)+v(5)+v(8))/3;
            ca(3+locf) = (v(3)+v(6)+v(9))/3;

            % N = [(v2-v1)a(v3-v1)]/2;
            Na(1+locf) = wa(f) * ((v(5)-v(2))*(v(9)-v(3))-(v(6)-v(3))*(v(8)-v(2)))/2;
            Na(2+locf) = wa(f) * ((v(6)-v(3))*(v(7)-v(1))-(v(4)-v(1))*(v(9)-v(3)))/2;
            Na(3+locf) = wa(f) * ((v(4)-v(1))*(v(8)-v(2))-(v(5)-v(2))*(v(7)-v(1)))/2;
        end

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Scalar product in W^* space %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function R = regWdScalProd(ca,Na,cb,Nb,sigmaW2)
        nfa = size(ca,2);
        nfb = size(cb,2);
        R = 0;
        for k = 1:nfa
            for j = 1:nfb
                argin = -((ca(1,k)-cb(1,j))^2+(ca(2,k)-cb(2,j))^2+(ca(3,k)-cb(3,j))^2)/sigmaW2;
argout = exp(argin);  %% BUILT IN KERNEL kerI, do not remove this comment
                R = R + argout * (Na(1,k)*Nb(1,j)+Na(2,k)*Nb(2,j)+Na(3,k)*Nb(3,j));
            end
        end
    end

    function R = fgtWdScalProd(ca,Na,cb,Nb,sigmaW2)
%          tic
        R = sum(sum(Nb.*fgt(3,ca,Na,cb,sqrt(sigmaW2),order,K,e)));
%          toc
    end


   function R = gridWdScalProd(ca,Na,cb,Nb,grille)
      
%         gammaX = ifftn(fftn(alphaA(:,:,:,1)).*fft3k);
%         gammaY = ifftn(fftn(alphaA(:,:,:,2)).*fft3k);
%         gammaZ = ifftn(fftn(alphaA(:,:,:,3)).*fft3k);
%         aux = alphaB(:,:,:,1).*gammaX + alphaB(:,:,:,2).*gammaY + alphaB(:,:,:,3).*gammaZ;
%         R = sum(aux(:));

        R = sum(sum(Nb.*gridOptim(ca,Na,cb,grille.long,grille.pas,grille.origine,grille.fft3k_d)));
    end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% gradient of W^* norm %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function eta = reggradWdNorm(phi,cphi,Nphi,sigmaW2)

        eta = zeros(size(phi));
        etaloc = zeros(3);
        etaf_w= zeros(3,1);
        etaf_dw = zeros(3,1);
        v = zeros(1,9);

        for e = 1:nfx
            loce = 3*(e-1);
            etaf_w(1) = 0;
            etaf_w(2) = 0;
            etaf_w(3) = 0;
            etaf_dw(1) = 0;
            etaf_dw(2) = 0;
            etaf_dw(3) = 0;
            if wx(e)
            for f = 1:nfx
                locf = 3*(f-1);
                argin = -(...
                    (cphi(1+locf)-cphi(1+loce))^2 + ...
                    (cphi(2+locf)-cphi(2+loce))^2 + ...
                    (cphi(3+locf)-cphi(3+loce))^2 )/sigmaW2;
argout = exp(argin);  %% BUILT IN KERNEL kerI, do not remove this comment
                etaf_w(1) = etaf_w(1) + argout * Nphi(1+locf);
                etaf_w(2) = etaf_w(2) + argout * Nphi(2+locf);
                etaf_w(3) = etaf_w(3) + argout * Nphi(3+locf);
argout = exp(argin);  %% BUILT IN KERNEL derI, do not remove this comment
                dotN = Nphi(1+locf)*Nphi(1+loce) + Nphi(2+locf)*Nphi(2+loce) + Nphi(3+locf)*Nphi(3+loce);
                etaf_dw(1) = etaf_dw(1) + argout * dotN * (cphi(1+locf)-cphi(1+loce));
                etaf_dw(2) = etaf_dw(2) + argout * dotN * (cphi(2+locf)-cphi(2+loce));
                etaf_dw(3) = etaf_dw(3) + argout * dotN * (cphi(3+locf)-cphi(3+loce));
            end
            for f = 1:nfy
                locf = 3*(f-1);
                argin = -(...
                    (cy(1+locf)-cphi(1+loce))^2 + ...
                    (cy(2+locf)-cphi(2+loce))^2 + ...
                    (cy(3+locf)-cphi(3+loce))^2 )/sigmaW2;
argout = exp(argin);  %% BUILT IN KERNEL kerI, do not remove this comment
                etaf_w(1) = etaf_w(1) - argout * Ny(1+locf);
                etaf_w(2) = etaf_w(2) - argout * Ny(2+locf);
                etaf_w(3) = etaf_w(3) - argout * Ny(3+locf);
argout = exp(argin);  %% BUILT IN KERNEL derI, do not remove this comment
                dotN = Ny(1+locf)*Nphi(1+loce) + Ny(2+locf)*Nphi(2+loce) + Ny(3+locf)*Nphi(3+loce);
                etaf_dw(1) = etaf_dw(1) - argout * dotN * (cy(1+locf)-cphi(1+loce));
                etaf_dw(2) = etaf_dw(2) - argout * dotN * (cy(2+locf)-cphi(2+loce));
                etaf_dw(3) = etaf_dw(3) - argout * dotN * (cy(3+locf)-cphi(3+loce));
            end
            end
            etaf_w(1) = wx(e) * etaf_w(1);
            etaf_w(2) = wx(e) * etaf_w(2);
            etaf_w(3) = wx(e) * etaf_w(3);
            etaf_dw(1) = 4 * etaf_dw(1) / sigmaW2 / 3;
            etaf_dw(2) = 4 * etaf_dw(2) / sigmaW2 / 3;
            etaf_dw(3) = 4 * etaf_dw(3) / sigmaW2 / 3;
            for k = 1:3
                for j = 1:3
                    v(k+3*(j-1)) = phi(k+3*(vx(j+loce)-1));
                end
            end
            etaloc(1) = (v(5)-v(8))*etaf_w(3)-(v(6)-v(9))*etaf_w(2) + etaf_dw(1);
            etaloc(2) = (v(6)-v(9))*etaf_w(1)-(v(4)-v(7))*etaf_w(3) + etaf_dw(2);
            etaloc(3) = (v(4)-v(7))*etaf_w(2)-(v(5)-v(8))*etaf_w(1) + etaf_dw(3);
            etaloc(4) = (v(8)-v(2))*etaf_w(3)-(v(9)-v(3))*etaf_w(2) + etaf_dw(1);
            etaloc(5) = (v(9)-v(3))*etaf_w(1)-(v(7)-v(1))*etaf_w(3) + etaf_dw(2);
            etaloc(6) = (v(7)-v(1))*etaf_w(2)-(v(8)-v(2))*etaf_w(1) + etaf_dw(3);
            etaloc(7) = (v(2)-v(5))*etaf_w(3)-(v(3)-v(6))*etaf_w(2) + etaf_dw(1);
            etaloc(8) = (v(3)-v(6))*etaf_w(1)-(v(1)-v(4))*etaf_w(3) + etaf_dw(2);
            etaloc(9) = (v(1)-v(4))*etaf_w(2)-(v(2)-v(5))*etaf_w(1) + etaf_dw(3);
            loc1 = 3*(vx(1+loce)-1);
            loc2 = 3*(vx(2+loce)-1);
            loc3 = 3*(vx(3+loce)-1);
            eta(1+loc1) = eta(1+loc1) + etaloc(1);
            eta(2+loc1) = eta(2+loc1) + etaloc(2);
            eta(3+loc1) = eta(3+loc1) + etaloc(3);
            eta(1+loc2) = eta(1+loc2) + etaloc(4);
            eta(2+loc2) = eta(2+loc2) + etaloc(5);
            eta(3+loc2) = eta(3+loc2) + etaloc(6);
            eta(1+loc3) = eta(1+loc3) + etaloc(7);
            eta(2+loc3) = eta(2+loc3) + etaloc(8);
            eta(3+loc3) = eta(3+loc3) + etaloc(9);
        end
    end


    function eta = fgtgradWdNorm(phi,cphi,Nphi,sigmaW2)
        cz = [cphi,cy];
        Nz = [Nphi,-Ny];

        resfgt = fgt(3,cz,[Nz;specprod(cz,Nz,1)],cphi,sigmaW,order,K,e);

        eta_dw = (4/3/sigmaW2) * (specprod(resfgt(4:end,:),Nphi,3) - specprod(specprod(cphi,Nphi,1),resfgt(1:3,:),3));
        eta_w = resfgt(1:3,:);

        eta = zeros(prod(size(phi)),1);
        for m = 1:nfx
            locm = 3*m-2:3*m;
            vm1 = 3*vx(1,m); locvm1 = vm1-2:vm1;
            vm2 = 3*vx(2,m); locvm2 = vm2-2:vm2;
            vm3 = 3*vx(3,m); locvm3 = vm3-2:vm3;
            e1 = phi(locvm2)-phi(locvm3);
            e2 = phi(locvm3)-phi(locvm1);
            e3 = phi(locvm1)-phi(locvm2);
            eta(locvm1) = eta(locvm1)' + eta_dw(locm) + wx(m) * cross(e1,eta_w(locm));
            eta(locvm2) = eta(locvm2)' + eta_dw(locm) + wx(m) * cross(e2,eta_w(locm));
            eta(locvm3) = eta(locvm3)' + eta_dw(locm) + wx(m) * cross(e3,eta_w(locm));
        end
        eta = reshape(eta,size(phi));
    end


    function eta = gridgradWdNorm(phi,cphi,Nphi,grille)

       cz = [cphi,cy];
       Nz = [Nphi,-Ny];

%        tic
       resgrid = gridOptim(cz,[Nz; specprod(cz,Nz,1)],cphi,grille.long,grille.pas,grille.origine,grille.fft3k_d);
%        disp(['Gradient grille ' num2str(toc)]);
%        tic
%        resfgt = fgt(3,cz,[Nz;specprod(cz,Nz,1)],cphi,sigmaW,order,K,e);
%        disp(['Gradient fgt ' num2str(toc)]);
       eta_dw = (4/3/sigmaW2) * (specprod(resgrid(4:end,:),Nphi,3) - specprod(specprod(cphi,Nphi,1),resgrid(1:3,:),3));
       eta_w = resgrid(1:3,:);

       eta = zeros(prod(size(phi)),1);
       for m = 1:nfx
           locm = 3*m-2:3*m;
           vm1 = 3*vx(1,m); locvm1 = vm1-2:vm1;
           vm2 = 3*vx(2,m); locvm2 = vm2-2:vm2;
           vm3 = 3*vx(3,m); locvm3 = vm3-2:vm3;
           e1 = phi(locvm2)-phi(locvm3);
           e2 = phi(locvm3)-phi(locvm1);
           e3 = phi(locvm1)-phi(locvm2);
           eta(locvm1) = eta(locvm1)' + eta_dw(locm) + wx(m) * cross(e1,eta_w(locm));
           eta(locvm2) = eta(locvm2)' + eta_dw(locm) + wx(m) * cross(e2,eta_w(locm));
           eta(locvm3) = eta(locvm3)' + eta_dw(locm) + wx(m) * cross(e3,eta_w(locm));
       end
       eta = reshape(eta,size(phi));
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




%     function alpha = projectionTriLineaire(taux,cx,grille)
% 
%       n = size(taux,1);
%       step = grille.pas;
%       % projection trilineaire des moments aux points de la grille
%       alpha = zeros([grille.long n]);
%       nf = size(cx,2);
%       for k=1:nf
%         c0 = floor((cx(:,k) - grille.origine)/step) + 1;
%         delta = cx(:,k) - (grille.origine + (c0-1)*step);
% 
%         rho000 = (step-delta(1)) * (step-delta(2)) * (step-delta(3));
%         rho100 = delta(1)        * (step-delta(2)) * (step-delta(3));
%         rho010 = (step-delta(1)) * delta(2)        * (step-delta(3));
%         rho001 = (step-delta(1)) * (step-delta(2)) * delta(3);
%         rho110 = delta(1)        * delta(2)        * (step-delta(3));
%         rho011 = (step-delta(1)) * delta(2)        * delta(3);
%         rho101 = delta(1)        * (step-delta(2)) * delta(3);
%         rho111 = delta(1)        * delta(2)        * delta(3);
% 
%         tauxk = zeros(1,1,1,n);
%         tauxk(1,1,1,:) = taux(:,k);
% 
%         alpha(c0(1),  c0(2),  c0(3),  :) = alpha(c0(1),  c0(2),  c0(3),  :) + rho000*tauxk;
%         alpha(c0(1)+1,c0(2),  c0(3),  :) = alpha(c0(1)+1,c0(2),  c0(3),  :) + rho100*tauxk;
%         alpha(c0(1),  c0(2)+1,c0(3),  :) = alpha(c0(1),  c0(2)+1,c0(3),  :) + rho010*tauxk;
%         alpha(c0(1),  c0(2),  c0(3)+1,:) = alpha(c0(1),  c0(2),  c0(3)+1,:) + rho001*tauxk;
%         alpha(c0(1)+1,c0(2)+1,c0(3),  :) = alpha(c0(1)+1,c0(2)+1,c0(3),  :) + rho110*tauxk;
%         alpha(c0(1),  c0(2)+1,c0(3)+1,:) = alpha(c0(1),  c0(2)+1,c0(3)+1,:) + rho011*tauxk;
%         alpha(c0(1)+1,c0(2),  c0(3)+1,:) = alpha(c0(1)+1,c0(2),  c0(3)+1,:) + rho101*tauxk;
%         alpha(c0(1)+1,c0(2)+1,c0(3)+1,:) = alpha(c0(1)+1,c0(2)+1,c0(3)+1,:) + rho111*tauxk;
% 
% 
%       end
%       alpha = alpha / (step*step*step);
% 
% %        if (iter==0)
% %        [x y z] = ndgrid(0:(grille.long(1)-1),0:(grille.long(2)-1),0:(grille.long(3)-1));
% %  	    x = x*grille.pas + grille.origine(1);
% %  	    y = y*grille.pas + grille.origine(2);
% %  	    z = z*grille.pas + grille.origine(3);
% %  
% %          disp(size(max(alphaX,[],3)));
% %          figure;
% %          hold on
% %          quiver(x(:,:,1),y(:,:,1),alphaX(:,:,8),alphaY(:,:,8),0,'r');
% %  %          for k = 1:nf
% %  %  		  plot([cx(1,k)-taux(1,k)/2, cx(1,k)+taux(1,k)/2], [cx(2,k)-taux(2,k)/2, cx(2,k)+taux(2,k)/2],'-b');
% %  %            plot(cx(1,k),cx(2,k),'+r')
% %  %  	    end
% %  	    hold off
% %        end
%     end
% 
% 
%     function res = interpolationTriLineaire(gamma,cx,grille)
% 
%       step = grille.pas;
%       nf = size(cx,2);
%       % interpolation trilineaire aux points cx
%       res = zeros(size(gamma,4),nf);
% 
%       for k=1:nf
%         c0 = floor((cx(:,k) - grille.origine)/step) + 1;
%         delta = cx(:,k) - (grille.origine + (c0-1)*step);
% 
%         rho000 = (step-delta(1)) * (step-delta(2)) * (step-delta(3));
%         rho100 = delta(1)        * (step-delta(2)) * (step-delta(3));
%         rho010 = (step-delta(1)) * delta(2)        * (step-delta(3));
%         rho001 = (step-delta(1)) * (step-delta(2)) * delta(3);
%         rho110 = delta(1)        * delta(2)        * (step-delta(3));
%         rho011 = (step-delta(1)) * delta(2)        * delta(3);
%         rho101 = delta(1)        * (step-delta(2)) * delta(3);
%         rho111 = delta(1)        * delta(2)        * delta(3);
% 
%         res(:,k) = rho000*gamma(c0(1),  c0(2),  c0(3),  :) +...
%                    rho100*gamma(c0(1)+1,c0(2),  c0(3),  :) +...
%                    rho010*gamma(c0(1),  c0(2)+1,c0(3),  :) +...
%                    rho001*gamma(c0(1),  c0(2),  c0(3)+1,:) +...
%                    rho110*gamma(c0(1)+1,c0(2)+1,c0(3),  :) +...
%                    rho011*gamma(c0(1),  c0(2)+1,c0(3)+1,:) +...
%                    rho101*gamma(c0(1)+1,c0(2),  c0(3)+1,:) +...
%                    rho111*gamma(c0(1)+1,c0(2)+1,c0(3)+1,:);
%        end
%        res = res / (step*step*step);
%      end

