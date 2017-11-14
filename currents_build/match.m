function s = match(s,target)

% Generic matching algorithm
% s = match(s,target);
%
% contact : joan_glaunes@yahoo.fr

kerV = 'exp(argin)';  %% BUILT IN KERNEL FLAG kerV, do not remove this comment

% x (template points) is required
if ~isfield(s,'x')
    error('input structure must contain at least field <x> (template points)');
end

% x must be 3D data
if size(s.x,1)~=3
    error('field <x> must be 3-rows matrix');
end

x = s.x;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%        Set parameters         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nx = size(x,2);           % number of template points

transmatrix = eye(3);       % additional affine transformation
transvector = zeros(3,1);   % additional affine transformation

method = '';
if nargin==2 && isfield(s,'target')
    s = rmfield(s,'target');
end

greedymatching = 0;
%greedystep = 1;
greedystep = 0.5;

rigidmatching = 1;
affinematching = 0;
elasticmatching = 1;
initialize = 0;

T = 30;                     % time discretization
sigmaV = 'auto';            % deformation kernel size
powcoefV = 0;          % normalization coef for deformation kernel = sigmaV^powcoefV
stdV = 1;                 % multiplication factor kernel = stdV*exp(-|x-y|^2/sigmaV^2);
gammaR = 'auto';            % balance between deformation and matching
resetmode = 1;              % if 1 then reset trajectories and moments

elapsedtime = 0;

useoptim = 'adaptdesc';      % optimization method ( 'adaptdesc' or 'conjgrad')
optim_maxiter = 500;        % maximum number of iterations
optim_stepsize = 'auto';    % gradient descent step
optim_breakratio = 1e-4;    % ratio for termination criteria
optim_verbosemode = 0;      % if 1 display text information while processing
optim_dosaveiter = 0;       % if 1 save result in file savename after each iteration
savename = 'result';

usefgt = 0;                 % if 1 use IFGT algorithm
testfgt = 0;
order = 5;                  % IFGT order
K = min([10,floor(nx/2)]);  % IFGT number of clusters
e = 9;                      % IFGT cutoff ratio for targets

usegrid = 0;             % if 1 use grid optimisation
testgrid = 0;
ratio = 0.2;               % ratio grid's step / sigmaW : the smaller, the smaller approximation's error


numbminims = 1;
J = [];
energie = [];
attache = [];

X = 0;
mom = 0;
distIdPhi = 0;

XInit = 0;
momInit = 0;

% update by parameters from input structure
loadstruct('s',who);

% error check
if (usefgt && usegrid)
 disp('fgt ou grille : il faut choisir !');
 return;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%        Initialize variables for matching process      %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tau = 1/(T-1);

% auto scaling
sclsz = mean(max(x')-min(x'));
if strcmp(sigmaV,'auto')
    sigmaV = sclsz*.05.^((0:-1:1-T)/(1-T));
end
if strcmp(gammaR,'auto')
    gammaR = .01*sclsz^2;
end
if strcmp(optim_stepsize,'auto')
    optim_stepsize = .5/sclsz;
end
clear sclsz

% kernel size
if length(sigmaV)==1
    sigmaV = sigmaV*ones(1,T);
end
if length(stdV)==1
    stdV = stdV*ones(1,T);
end
sigmaV2 = sigmaV.^2;
stdV2 = stdV.^2;
normcoefV = sigmaV.^powcoefV;


% initialize targets and get their function handles
if nargin==2
    ntargets = length(target);
    targetsgrid = cell(1,ntargets);
    for k = 1:ntargets
        if isfield(target{k},'weight')
            targetkweight = target{k}.weight;
        else
            targetkweight = 1;
        end
        target{k}.numbminims = numbminims;
        target{k} = feval(target{k}.method,target{k});
        target{k}.weight = targetkweight;
    end
    clear targetkweight

elseif exist('method')
    target{1} = feval(method,s);
    ntargets = 1;
    target{1}.weight = 1;
else
    ntargets = 0;
end


% functions for grid optimization; if needed
if usegrid
    mini = min(x,[],2);
    maxi = max(x,[],2);
    sourcegrid = cell(1,T);
    for t=1:T
       sourcegrid{t} = setgrid(mini,maxi,t);
%         disp(['source''s grid ' num2str(t) ' has been set:  ' num2str([sourcegrid{t}.pas sourcegrid{t}.long sourcegrid{t}.origine'])]);
    end
end
clear mini maxi


% set matching functionals and gradients

    function R = matching(phix)
        R = 0;
        phix = reshape(phix,3,nx);
        for k = 1:ntargets
            if isfield(target{k},'usegrid')&&(target{k}.usegrid)
                R = R + target{k}.weight * target{k}.matching(phix(:),targetsgrid{k});
            else
                R = R + target{k}.weight * target{k}.matching(phix(:));
            end
        end
        phix = phix(:);
    end

    function g = gradmatching(phix)
        g = zeros(size(phix));
        for k = 1:ntargets
            if isfield(target{k},'usegrid')&&target{k}.usegrid
               g = g + target{k}.weight * target{k}.gradmatching(phix,targetsgrid{k});
            else
               g = g + target{k}.weight * target{k}.gradmatching(phix);
            end
        end
    end

    function premin
        for k = 1:ntargets
            if (isfield(target{k},'usegrid')&&(target{k}.usegrid))
                pts = x(:,target{k}.vx(:));
                mini = min(pts')';
                maxi = max(pts')';
                targetsgrid{k} = target{k}.premin(mini,maxi);
%                 disp(['target ' num2str(k) ': grid has been set:  ' num2str([targetsgrid{k}.pas targetsgrid{k}.long targetsgrid{k}.origine'])]);
            elseif isfield(target{k},'premin')
                target{k}.premin();
            end
        end
    end

% build optimization options structure
optvars = who('optim_*');
optim.method = useoptim;
for v = 1:length(optvars)
    var = optvars{v};
    eval([regexprep(var,'optim_','optim.'),' = ',var,';']);
end
optim.saveiter = @saveiter;
clear optvars var v

% prior affine transformation
xrig = transmatrix * x + repmat(transvector,1,nx);


phix = zeros(3*nx,1);

dX = zeros(3*nx,1);
eta = 0;
deta = 0;

if ~strcmp(kerV,'exp(argin)') && usefgt
    warning('Current built-in kernel is not gaussian. Fast Gauss Transform disabled')
    usefgt = 0;
end

if usefgt
    kernelsum = @fgtkernelsum;
    compdeta = @fgtcompdeta;
elseif usegrid
    kernelsum = @gridkernelsum;
    compdeta = @gridcompdeta;
else
    kernelsum = @regkernelsum;
    compdeta = @regcompdeta;
end

if (usefgt&&testfgt)||(usegrid&&testgrid)
        a = rand(1,nx);
        essdir = zeros(1,nx);
        for m = 1:nx
            locm = 3*(m-1);
            for l = 1:nx
                locl = 3*(l-1);
                argin = -( ...
                    (x(1+locm)-x(1+locl))^2 + ...
                    (x(2+locm)-x(2+locl))^2 + ...
                    (x(3+locm)-x(3+locl))^2)/sigmaV2(1);
                argout = exp(argin);  %% BUILT IN KERNEL kerV, do not remove this comment
                essdir(m) = essdir(m) + argout * a(l);
            end
        end
        if usefgt
          essfgt = fgt(3,x,a,x,sigmaV(1),order,K,e);
          disp(['erreur relative fgt diffeos: ',num2str(mean(abs((essdir-essfgt)./essdir)))])
        else
          essgrid = gridOptim(x,a,x,sourcegrid{1}.long,sourcegrid{1}.pas,sourcegrid{1}.origine,sourcegrid{1}.fft3k_d);
          disp(['erreur relative grid diffeos: ',num2str(mean(abs((essdir-essgrid)./essdir)))])
        end
end

gradformula = compdeta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%              main                   %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
if rigidmatching
    rigidmatch;
end
if greedymatching
    greedymatch;
end
if affinematching
    affinematch;
end
if elasticmatching
    diffeomatch;
end
elapsedtime = elapsedtime + toc;

clear gradformula Jcur energieCur attacheCur dX deta eta optim meth
clear functions
savestruct('s',who)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% scale/rotation/translation matching %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function rigidmatch
        % recenter
        meanxrig = mean(xrig,2);
        xrig = xrig - repmat(meanxrig,1,nx);
        transvector = transvector - meanxrig ;

        q = [1;0;0;0;meanxrig]; % first 4 digits = quaternion, last 3 = translation vector / scale

        [q,Jcur] = feval(useoptim,q,@rigidfunct,@rigidgrad,optim);
        J = [J,Jcur(:)'];

        [xrig,Rq] = comprigid(q);
        transmatrix = Rq * transmatrix;
        transvector = Rq * transvector + [q(5);q(6);q(7)];
        phix = xrig;

        function [J energie attache] = rigidfunct(q)
            % functional for rigid matching algorithm
            phi = comprigid(q);
            % compute matching term
            J = matching(phi(:));
	    energie = 0;
	    attache = 0;
        end

        function G = rigidgrad(q)
            a = q(1);
            ncross = [0,-q(4),q(3);q(4),0,-q(2);-q(3),q(2),0];
            ndot = [q(2),q(3),q(4);q(2),q(3),q(4);q(2),q(3),q(4)];
            [phi,Rq] = comprigid(q);
            % computes gradient wrt points %%%
            Gp = gradmatching(phi(:));
            Gp = reshape(Gp,3,nx);
            % computes gradient wrt rigid motion params %%%
            G = zeros(7,1);
            ncrossxrig = ncross*xrig;
            G(1) = 2 * sum(dot(a*xrig + ncrossxrig , Gp));
            G(2:4) = 2 * sum(cross(ncrossxrig,Gp) + a*cross(xrig,Gp) + (ndot*xrig).*Gp,2);
            G(5:7) = sum(Gp,2);
        end

        function [phi,Rq] = comprigid(q)
            Rq = [q(1)^2+q(2)^2-q(3)^2-q(4)^2 , 2*(q(2)*q(3)-q(1)*q(4)) , 2*(q(2)*q(4)+q(1)*q(3));...
                2*(q(2)*q(3)+q(1)*q(4)) , q(1)^2-q(2)^2+q(3)^2-q(4)^2 , 2*(q(3)*q(4)-q(1)*q(2));...
                2*(q(2)*q(4)-q(1)*q(3)) , 2*(q(3)*q(4)+q(1)*q(2)) , q(1)^2-q(2)^2-q(3)^2+q(4)^2];
            phi = Rq * xrig + repmat([q(5);q(6);q(7)],1,nx);
        end

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% affine matching %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function affinematch
        % recenter
        meanxrig = mean(xrig,2);
        xrig = xrig - repmat(meanxrig,1,nx);
        transvector = transvector - meanxrig ;
        scalex = max(max(xrig')-min(xrig'));

        A = [1;0;0;0;1;0;0;0;1;meanxrig/scalex]; % first 9 digits = linear transfo, last 3 = translation vector / scale

        [A,Jcur] = feval(useoptim,A,@affinefunct,@affinegrad,optim);
        J = [J,Jcur(:)'];

        [xrig,a,b] = compaffine(A);
        transmatrix = a * transmatrix;
        transvector = a * transvector + b;
        phix = xrig;

        function J = affinefunct(A)
            % functional for affine matching algorithm
            phi = compaffine(A);
            % compute matching term
            J = matching(phi(:));
        end

        function G = affinegrad(A)
            phi = compaffine(A);
            % computes gradient wrt points %%%
            Gp = gradmatching(phi(:));
            Gp = reshape(Gp,3,nx);
            % computes gradient wrt affine motion params %%%
            G = zeros(12,1);
            g1 = Gp*xrig';
            G(1:9) = g1(:);
            G(10:12) = sum(Gp,2) * scalex;
        end

        function [phi,a,b] = compaffine(A)
            a = reshape(A(1:9),3,3);
            b = scalex * A(10:12);
            phi = a * xrig +  repmat(b,1,nx);
        end

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% greedy descent %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function greedymatch
        X = zeros(3*nx,T);
        X(:,1) = xrig(:);
        mom = zeros(3*nx,T);
        tmax = 1;
        momX0 = [zeros(1,3*nx),xrig(:)'];

        if 0
            [Tech,momX] = ode45(@greedyfct,tmax*tau*(0:T-1),momX0);
            for t = 1:length(Tech)
                dmomX = greedyfct(Tech(t),momX(t,:)');
                mom(:,t) = dmomX(1:3*nx);
            end
        else
            momX = zeros(T,length(momX0));
            momX(1,:) = momX0;
            J = zeros(1,T);
            J(1) = matching(momX(1,3*nx+1:end));
            for m = 1:numbminims
                premin;
                for t = 1+(m-1)*(T-1)/numbminims:m*(T-1)/numbminims
                    if optim.verbosemode
                        disp(['t=',num2str(t)])
                    end
                    dmomX = greedyfct(t,momX(t,:))';
                    mom(:,t) = dmomX(1:3*nx);
                    momX(t+1,:) = momX(t,:) + tau * dmomX;
                    %max(abs(dmomX))
                    momX(t+1,:) = momX(t,:) + (.5*tau) * (dmomX+greedyfct(t+1,momX(t+1,:))');
                    J(t+1) = matching(momX(t+1,3*nx+1:end));
                end
            end
        end

        mom = reshape(mom,3,nx,T)*tmax;
        X = reshape(momX(:,3*nx+1:end)',3,nx,T);
        phix = X(:,:,T);

        function dmomX = greedyfct(t,momXt)
            dmomX = zeros(6*nx,1);
            dmomX(1:3*nx) = - greedystep * gradmatching(momXt(3*nx+1:end));
            dX(:) = 0;
            X(:,1) = momXt(3*nx+1:end);
            feval(kernelsum,1,dmomX(1:3*nx));
            dmomX(3*nx+1:end) = dX(:);
        end
    end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%    diffeo matching         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function diffeomatch

        % reset trajectories and momentum
        if initialize
            disp('initialize');
            X = reshape(XInit,3*nx,T);
            mom = reshape(momInit,3*nx,T);
        elseif resetmode
            elapsedtime = 0;
            X = zeros(3*nx,T);
            X(:,1) = xrig(:);
            mom = zeros(3*nx,T);
            J = [];
	        energie = [];
	        attache = [];
        else
            X = reshape(X,3*nx,T);
            mom = reshape(mom,3*nx,T);
        end

        if strcmp(useoptim,'bfgs')
            mom = reshape(mom,3*nx*T,1);
        end

        % temp variables
        eta = zeros(3*nx,T);
        deta = zeros(3*nx,1);

        for minimstep = 1:numbminims
            premin;
            % call to optimization routine
            [mom,Jcur,energieCur,attacheCur] = feval(useoptim,mom,@diffeofunct,@diffeograd,optim);
            if strcmp(useoptim,'bfgs')
                mom = reshape(mom,3*nx,T);
            end
            J = [J,Jcur(:)'];
	        energie = [energieCur,energie(:)'];
	        attache = [attacheCur,attache(:)'];
        end

        % compute trajectories
        comptraj(mom);
        distIdPhi = sqrt(energy(mom));
        mom = reshape(mom,3,nx,T);
        X = reshape(X,3,nx,T);
        phix = X(:,:,T);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% functional for diffeomatch %%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [J,energie,attache] = diffeofunct(mom)
            comptraj(mom);
            phix = X(:,T);

%              disp('bounding box 1, diffeofunct');
%              phixAux = reshape(phix,[3,4032]);
%              size(phixAux)
%              disp(min(phixAux(:,target{1}.vx(:))')');
%              disp(max(phixAux(:,target{1}.vx(:))')');
%              disp(' ');

            energie = energy(mom);
            attache = matching(phix);
            J = gammaR * energie + attache;
%            fprintf('energy = %e,  attache= %e, J = %e\n',energy(mom),matching(phix),J);
        end

        function E = energy(mom)
            % compute energy term
            E = 0;
            for k = 1:nx
                lock = 3*(k-1);
                for t = 1:T-1
                    E = E + (mom(1+lock,t)+mom(1+lock,t+1))*(X(1+lock,t+1)-X(1+lock,t))...
                        + (mom(2+lock,t)+mom(2+lock,t+1))*(X(2+lock,t+1)-X(2+lock,t))...
                        + (mom(3+lock,t)+mom(3+lock,t+1))*(X(3+lock,t+1)-X(3+lock,t));
                end
            end
            E = .5*E;  %%% because we use mom(t)+mom(t+1) above
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% gradient for diffeomatch  %%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function G = diffeograd(mom)

             if strcmp(useoptim,'bfgs')
                mom = reshape(mom,3*nx,T);
             end
            comptraj(mom);
            phix = X(:,T);

%              disp('bounding box 1, diffeograd');
%              phixAux = reshape(phix,[3,4032]);
%              size(phixAux)
%              disp(min(phixAux(:,target{1}.vx(:))')');
%              disp(max(phixAux(:,target{1}.vx(:))')');
%              disp(' ');


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% computes eta at t=1
            %%% ie gradient of W^* norm wrt points
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            eta(:,T) = gradmatching(phix);

            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% computes eta(t)
            %%%%%%%%%%%%%%%%%%%%%%%%%%

            % computes eta(t) with backward integration centered scheme
            for t=T:-1:2
                feval(gradformula,t,mom);
                eta(:,t-1)=eta(:,t) + (tau*2/sigmaV2(t)) * normcoefV(t) * deta;
                feval(gradformula,t-1,mom);
                eta(:,t-1)=eta(:,t) + (tau/sigmaV2(t)) * normcoefV(t) * deta; % factor 2 cancels
                deta(:) = 0;
            end
%              disp('fin schema d''integration');
            G = 2*gammaR*mom + eta;
            if strcmp(useoptim,'bfgs')
                G = reshape(G,3*nx*T,1);
            end
        end

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function saveiter(mom,Jcur)
        savestruct('s',who)
        s.J = [s.J,Jcur(:)'];
        if prod(size(mom))==3*nx*T
        s.mom = reshape(mom,3,nx,T);
        s.X = reshape(X,3,nx,T);
        s.phix = s.X(:,:,T);
        end
        save(savename,'s');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function regcompdeta(t,mom)

        for m = 1:nx
            locm = 3*(m-1);
            for ll = 1:nx
                locl = 3*(ll-1);
                argin = -((X(locm+1,t)-X(locl+1,t))^2 + (X(locm+2,t)-X(locl+2,t))^2 + (X(locm+3,t)-X(locl+3,t))^2)/sigmaV2(t);
argout = stdV2(t)*exp(argin);  %% BUILT IN KERNEL derV, do not remove this comment
                argout = - argout * ( ...
                    (eta(locl+1,t)*mom(locm+1,t) + eta(locl+2,t)*mom(locm+2,t) + eta(locl+3,t)*mom(locm+3,t)) + ...
                    (eta(1+locm,t)*mom(1+locl,t) + eta(2+locm,t)*mom(2+locl,t) + eta(3+locm,t)*mom(3+locl,t)) + ...
                    2*gammaR * (mom(1+locl,t)*mom(1+locm,t) + mom(2+locl,t)*mom(2+locm,t) + mom(3+locl,t)*mom(3+locm,t)));

                deta(1+locm) = deta(1+locm) + argout * (X(1+locm,t)-X(1+locl,t));
                deta(2+locm) = deta(2+locm) + argout * (X(2+locm,t)-X(2+locl,t));
                deta(3+locm) = deta(3+locm) + argout * (X(3+locm,t)-X(3+locl,t));
            end
        end
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function fgtcompdeta(t,mom)
        Xt = reshape(X(:,t),3,nx);
        X1t = Xt(1,:);
        X2t = Xt(2,:);
        X3t = Xt(3,:);
        etat = reshape(eta(:,t),3,nx);
        eta1t = etat(1,:);
        eta2t = etat(2,:);
        eta3t = etat(3,:);
        momt = reshape(mom(:,t),3,nx);
        mom1t = momt(1,:);
        mom2t = momt(2,:);
        mom3t = momt(3,:);
        deta = reshape(deta,3,nx);

%         tic
        resfgt = fgt(3,Xt,[eta1t;eta2t;eta3t;mom1t;mom2t;mom3t;...
            eta1t.*X1t;eta1t.*X2t;eta1t.*X3t;...
            eta2t.*X1t;eta2t.*X2t;eta2t.*X3t;...
            eta3t.*X1t;eta3t.*X2t;eta3t.*X3t;...
            mom1t.*X1t;mom1t.*X2t;mom1t.*X3t;...
            mom2t.*X1t;mom2t.*X2t;mom2t.*X3t;...
            mom3t.*X1t;mom3t.*X2t;mom3t.*X3t],Xt,sigmaV(t),order,K,e);
%             toc

        resfgt = stdV2(t)*resfgt;

        temp = mom1t.*resfgt(1,:);
        deta(1,:) = deta(1,:) - temp.*X1t;
        deta(2,:) = deta(2,:) - temp.*X2t;
        deta(3,:) = deta(3,:) - temp.*X3t;
        temp = mom2t.*resfgt(2,:);
        deta(1,:) = deta(1,:) - temp.*X1t;
        deta(2,:) = deta(2,:) - temp.*X2t;
        deta(3,:) = deta(3,:) - temp.*X3t;
        temp = mom3t.*resfgt(3,:);
        deta(1,:) = deta(1,:) - temp.*X1t;
        deta(2,:) = deta(2,:) - temp.*X2t;
        deta(3,:) = deta(3,:) - temp.*X3t;

        temp = resfgt(4,:);
        dum = temp.*eta1t;
        deta(1,:) = deta(1,:) - dum.*X1t;
        deta(2,:) = deta(2,:) - dum.*X2t;
        deta(3,:) = deta(3,:) - dum.*X3t;
        dum = 2*gammaR*temp.*mom1t;
        deta(1,:) = deta(1,:) - dum.*X1t;
        deta(2,:) = deta(2,:) - dum.*X2t;
        deta(3,:) = deta(3,:) - dum.*X3t;

        temp = resfgt(5,:);
        dum = temp.*eta2t;
        deta(1,:) = deta(1,:) - dum.*X1t;
        deta(2,:) = deta(2,:) - dum.*X2t;
        deta(3,:) = deta(3,:) - dum.*X3t;
        dum = 2*gammaR*temp.*mom2t;
        deta(1,:) = deta(1,:) - dum.*X1t;
        deta(2,:) = deta(2,:) - dum.*X2t;
        deta(3,:) = deta(3,:) - dum.*X3t;

        temp = resfgt(6,:);
        dum = temp.*eta3t;
        deta(1,:) = deta(1,:) - dum.*X1t;
        deta(2,:) = deta(2,:) - dum.*X2t;
        deta(3,:) = deta(3,:) - dum.*X3t;
        dum = 2*gammaR*temp.*mom3t;
        deta(1,:) = deta(1,:) - dum.*X1t;
        deta(2,:) = deta(2,:) - dum.*X2t;
        deta(3,:) = deta(3,:) - dum.*X3t;


        temp = mom1t.*resfgt(7,:);
        deta(1,:) = deta(1,:) + temp;
        temp = mom1t.*resfgt(8,:);
        deta(2,:) = deta(2,:) + temp;
        temp = mom1t.*resfgt(9,:);
        deta(3,:) = deta(3,:) + temp;

        temp = mom2t.*resfgt(10,:);
        deta(1,:) = deta(1,:) + temp;
        temp = mom2t.*resfgt(11,:);
        deta(2,:) = deta(2,:) + temp;
        temp = mom2t.*resfgt(12,:);
        deta(3,:) = deta(3,:) + temp;

        temp = mom3t.*resfgt(13,:);
        deta(1,:) = deta(1,:) + temp;
        temp = mom3t.*resfgt(14,:);
        deta(2,:) = deta(2,:) + temp;
        temp = mom3t.*resfgt(15,:);
        deta(3,:) = deta(3,:) + temp;




        temp = resfgt(16,:);
        dum = temp.*eta1t;
        deta(1,:) = deta(1,:) + dum;
        dum = 2*gammaR*temp.*mom1t;
        deta(1,:) = deta(1,:) + dum;

        temp = resfgt(17,:);
        dum = temp.*eta1t;
        deta(2,:) = deta(2,:) + dum;
        dum = 2*gammaR*temp.*mom1t;
        deta(2,:) = deta(2,:) + dum;

        temp = resfgt(18,:);
        dum = temp.*eta1t;
        deta(3,:) = deta(3,:) + dum;
        dum = 2*gammaR*temp.*mom1t;
        deta(3,:) = deta(3,:) + dum;



        temp = resfgt(19,:);
        dum = temp.*eta2t;
        deta(1,:) = deta(1,:) + dum;
        dum = 2*gammaR*temp.*mom2t;
        deta(1,:) = deta(1,:) + dum;

        temp = resfgt(20,:);
        dum = temp.*eta2t;
        deta(2,:) = deta(2,:) + dum;
        dum = 2*gammaR*temp.*mom2t;
        deta(2,:) = deta(2,:) + dum;

        temp = resfgt(21,:);
        dum = temp.*eta2t;
        deta(3,:) = deta(3,:) + dum;
        dum = 2*gammaR*temp.*mom2t;
        deta(3,:) = deta(3,:) + dum;



        temp = resfgt(22,:);
        dum = temp.*eta3t;
        deta(1,:) = deta(1,:) + dum;
        dum = 2*gammaR*temp.*mom3t;
        deta(1,:) = deta(1,:) + dum;

        temp = resfgt(23,:);
        dum = temp.*eta3t;
        deta(2,:) = deta(2,:) + dum;
        dum = 2*gammaR*temp.*mom3t;
        deta(2,:) = deta(2,:) + dum;

        temp = resfgt(24,:);
        dum = temp.*eta3t;
        deta(3,:) = deta(3,:) + dum;
        dum = 2*gammaR*temp.*mom3t;
        deta(3,:) = deta(3,:) + dum;

        deta = deta(:);
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function gridcompdeta(t,mom)
        Xt = reshape(X(:,t),3,nx);
        X1t = Xt(1,:);
        X2t = Xt(2,:);
        X3t = Xt(3,:);
        etat = reshape(eta(:,t),3,nx);
        eta1t = etat(1,:);
        eta2t = etat(2,:);
        eta3t = etat(3,:);
        momt = reshape(mom(:,t),3,nx);
        mom1t = momt(1,:);
        mom2t = momt(2,:);
        mom3t = momt(3,:);
        deta = reshape(deta,3,nx);

%         tic
        resgrid = gridOptim(Xt,[eta1t;eta2t;eta3t;mom1t;mom2t;mom3t;...
            eta1t.*X1t;eta1t.*X2t;eta1t.*X3t;...
            eta2t.*X1t;eta2t.*X2t;eta2t.*X3t;...
            eta3t.*X1t;eta3t.*X2t;eta3t.*X3t;...
            mom1t.*X1t;mom1t.*X2t;mom1t.*X3t;...
            mom2t.*X1t;mom2t.*X2t;mom2t.*X3t;...
            mom3t.*X1t;mom3t.*X2t;mom3t.*X3t],Xt,sourcegrid{t}.long,sourcegrid{t}.pas,sourcegrid{t}.origine,sourcegrid{t}.fft3k_d);
%         toc


        resgrid = stdV2(t)*resgrid;

        temp = mom1t.*resgrid(1,:);
        deta(1,:) = deta(1,:) - temp.*X1t;
        deta(2,:) = deta(2,:) - temp.*X2t;
        deta(3,:) = deta(3,:) - temp.*X3t;
        temp = mom2t.*resgrid(2,:);
        deta(1,:) = deta(1,:) - temp.*X1t;
        deta(2,:) = deta(2,:) - temp.*X2t;
        deta(3,:) = deta(3,:) - temp.*X3t;
        temp = mom3t.*resgrid(3,:);
        deta(1,:) = deta(1,:) - temp.*X1t;
        deta(2,:) = deta(2,:) - temp.*X2t;
        deta(3,:) = deta(3,:) - temp.*X3t;

        temp = resgrid(4,:);
        dum = temp.*eta1t;
        deta(1,:) = deta(1,:) - dum.*X1t;
        deta(2,:) = deta(2,:) - dum.*X2t;
        deta(3,:) = deta(3,:) - dum.*X3t;
        dum = 2*gammaR*temp.*mom1t;
        deta(1,:) = deta(1,:) - dum.*X1t;
        deta(2,:) = deta(2,:) - dum.*X2t;
        deta(3,:) = deta(3,:) - dum.*X3t;

        temp = resgrid(5,:);
        dum = temp.*eta2t;
        deta(1,:) = deta(1,:) - dum.*X1t;
        deta(2,:) = deta(2,:) - dum.*X2t;
        deta(3,:) = deta(3,:) - dum.*X3t;
        dum = 2*gammaR*temp.*mom2t;
        deta(1,:) = deta(1,:) - dum.*X1t;
        deta(2,:) = deta(2,:) - dum.*X2t;
        deta(3,:) = deta(3,:) - dum.*X3t;

        temp = resgrid(6,:);
        dum = temp.*eta3t;
        deta(1,:) = deta(1,:) - dum.*X1t;
        deta(2,:) = deta(2,:) - dum.*X2t;
        deta(3,:) = deta(3,:) - dum.*X3t;
        dum = 2*gammaR*temp.*mom3t;
        deta(1,:) = deta(1,:) - dum.*X1t;
        deta(2,:) = deta(2,:) - dum.*X2t;
        deta(3,:) = deta(3,:) - dum.*X3t;


        temp = mom1t.*resgrid(7,:);
        deta(1,:) = deta(1,:) + temp;
        temp = mom1t.*resgrid(8,:);
        deta(2,:) = deta(2,:) + temp;
        temp = mom1t.*resgrid(9,:);
        deta(3,:) = deta(3,:) + temp;

        temp = mom2t.*resgrid(10,:);
        deta(1,:) = deta(1,:) + temp;
        temp = mom2t.*resgrid(11,:);
        deta(2,:) = deta(2,:) + temp;
        temp = mom2t.*resgrid(12,:);
        deta(3,:) = deta(3,:) + temp;

        temp = mom3t.*resgrid(13,:);
        deta(1,:) = deta(1,:) + temp;
        temp = mom3t.*resgrid(14,:);
        deta(2,:) = deta(2,:) + temp;
        temp = mom3t.*resgrid(15,:);
        deta(3,:) = deta(3,:) + temp;


        temp = resgrid(16,:);
        dum = temp.*eta1t;
        deta(1,:) = deta(1,:) + dum;
        dum = 2*gammaR*temp.*mom1t;
        deta(1,:) = deta(1,:) + dum;

        temp = resgrid(17,:);
        dum = temp.*eta1t;
        deta(2,:) = deta(2,:) + dum;
        dum = 2*gammaR*temp.*mom1t;
        deta(2,:) = deta(2,:) + dum;

        temp = resgrid(18,:);
        dum = temp.*eta1t;
        deta(3,:) = deta(3,:) + dum;
        dum = 2*gammaR*temp.*mom1t;
        deta(3,:) = deta(3,:) + dum;


        temp = resgrid(19,:);
        dum = temp.*eta2t;
        deta(1,:) = deta(1,:) + dum;
        dum = 2*gammaR*temp.*mom2t;
        deta(1,:) = deta(1,:) + dum;

        temp = resgrid(20,:);
        dum = temp.*eta2t;
        deta(2,:) = deta(2,:) + dum;
        dum = 2*gammaR*temp.*mom2t;
        deta(2,:) = deta(2,:) + dum;

        temp = resgrid(21,:);
        dum = temp.*eta2t;
        deta(3,:) = deta(3,:) + dum;
        dum = 2*gammaR*temp.*mom2t;
        deta(3,:) = deta(3,:) + dum;



        temp = resgrid(22,:);
        dum = temp.*eta3t;
        deta(1,:) = deta(1,:) + dum;
        dum = 2*gammaR*temp.*mom3t;
        deta(1,:) = deta(1,:) + dum;

        temp = resgrid(23,:);
        dum = temp.*eta3t;
        deta(2,:) = deta(2,:) + dum;
        dum = 2*gammaR*temp.*mom3t;
        deta(2,:) = deta(2,:) + dum;

        temp = resgrid(24,:);
        dum = temp.*eta3t;
        deta(3,:) = deta(3,:) + dum;
        dum = 2*gammaR*temp.*mom3t;
        deta(3,:) = deta(3,:) + dum;

        deta = deta(:);
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% updates trajectories %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function comptraj(mom)
        for t = 1:T-1
            feval(kernelsum,t,mom);
            X(:,t+1) = X(:,t) + tau * normcoefV(t) * dX;
            if usegrid
                updatesourcegrids;
            end
            feval(kernelsum,t+1,mom);
            dX = .5*dX;
            X(:,t+1) = X(:,t) + tau * normcoefV(t) * dX;
            dX(:) = 0;
            if usegrid
                updatesourcegrids;
            end
            updatetargetgrids;
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function regkernelsum(t,mom)
        for mm = 1:nx
            locm = 3*(mm-1);
            for ll = 1:nx
                locl = 3*(ll-1);
                argin = -( ...
                    (X(1+locm,t)-X(1+locl,t))^2 + ...
                    (X(2+locm,t)-X(2+locl,t))^2 + ...
                    (X(3+locm,t)-X(3+locl,t))^2)/sigmaV2(t);
argout = stdV2(t)*exp(argin);  %% BUILT IN KERNEL kerV, do not remove this comment
                dX(1+locm) = dX(1+locm) + argout * mom(1+locl,t);
                dX(2+locm) = dX(2+locm) + argout * mom(2+locl,t);
                dX(3+locm) = dX(3+locm) + argout * mom(3+locl,t);
            end
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function fgtkernelsum(t,mom)
	[aux rad] = fgt(3,X(:,t),mom(:,t),X(:,t),sigmaV(t),order,K,e);
        dX = dX + stdV2(t)*aux;
%	fprintf('%e ',rad);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function gridkernelsum(t,mom)
        Xt = reshape(X(:,t),3,nx);
        momt = reshape(mom(:,t),3,nx);
        aux = gridOptim(Xt,momt,Xt,sourcegrid{t}.long,sourcegrid{t}.pas,sourcegrid{t}.origine,sourcegrid{t}.fft3k_d);
        dX = dX + stdV2(t)*aux(:);
    end


%%%%%%%%%%%%%%%%%%%%%%%
%%% functions for grid optimization %%%%%
%%%%%%%%%%%%%%%%%%%%%%%
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
           long = (maxi-mini)/grille.pas + 3/ratio; %circonf??rence du tore
%             long = (long<=32)*32 + (long>32).*((long<=64)*64 + (long>64).*((long<=128)*128 + (long>128)*2.*ceil(long/2)));
           long = (long<=16)*16 + (long>16).*((long<=32)*32 + (long>32).*((long<=64)*64 + (long>64).*2.*ceil(long/2)));
           p = sum(long < grille.long') ~= 0;
         end
 end

    function updatesourcegrids
        X = reshape(X,3,nx,T);
        for tp=1:T
            mini = min(X(:,:,tp),[],2);
            maxi = max(X(:,:,tp),[],2);
            if (changegrid(mini,maxi,tp,sourcegrid{tp}))
                sourcegrid{tp} = setgrid(mini,maxi,tp,sourcegrid{tp});
%                  disp(['source''s grid ' num2str(tp) ' has been changed:  ' num2str([sourcegrid{tp}.pas sourcegrid{tp}.long sourcegrid{tp}.origine'])]);
            end
        end
        X = reshape(X,3*nx,T);
    end

    function updatetargetgrids
        X = reshape(X,3,nx,T);
        for k=1:ntargets
            if isfield(target{k},'usegrid')&&(target{k}.usegrid)
                mini = min(X(:,target{k}.vx(:),T),[],2);
                maxi = max(X(:,target{k}.vx(:),T),[],2);
                if (target{k}.changegrid(mini,maxi,targetsgrid{k}))
                    targetsgrid{k} = target{k}.setgrid(mini,maxi,targetsgrid{k});
%                      disp(['target ' num2str(k) ': grid has been changed.  ' num2str([targetsgrid{k}.pas targetsgrid{k}.long targetsgrid{k}.origine'])]);
                end
            end
        end
        X = reshape(X,3*nx,T);
    end

end



