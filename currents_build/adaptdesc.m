function [x,J,energie,attache] = adaptdesc(x,f,g,options,varargin)

% Steepest descent with adaptive stepsize
%
% f functional
% g gradient of f
% x variable
% in structure options:
%   maxiter = maximum number of iterations; default 500;
%   stepsize = initial stepsize of descent; default 1/norm(initial gradient)
%   stepmult = increase ratio of stepsize after each minimization success;
%     default 1.2
%   stepdiv = decrease ratio of stepsize after each minimization failure;
%     default 2
%   breakratio = ratio for termination criteria; default 1e-4
%   loopbreak = maximum number of loops in each iteration
%   verbosemode if 1 then display evaluation of f at each step
%
% varargin additional arguments passed to f and g
%
% returns x variable at the end of minimization
% and J record of function evaluations
%
% contact : Joan Glaun??? - joan_glaunes@yahoo.fr

maxiter = 500;
[J(1) energie(1) attache(1)] = feval(f,x,varargin{:});
G = feval(g,x,varargin{:});
stepsize = J(1)/sum(G(:).^2);
stepmult = 1.2;
stepdiv = 2;
breakratio = 1e-4;
loopbreak = 40;
verbosemode = 0;

loadstruct('options');

stepsize = J(1)/sum(G(:).^2);  % stepsize was overloaded by options structure, we reset it.

for iter = 1:maxiter
    G = feval(g,x,varargin{:});
    stepsize = stepsize * stepdiv;
    minimtest = 0;
    loop = 0;
    while(~minimtest & (loop<loopbreak)) 
        stepsize = stepsize / stepdiv;
        xnew = x - stepsize * G; 
        [J(iter+1) energie(iter+1) attache(iter+1)] = feval(f,xnew,varargin{:});
        minimtest = (J(iter+1) < J(iter));
        loop = loop + 1;
    end
    x = xnew;
    stepsize = stepsize * stepmult;
    if verbosemode
        disp(['iteration ',num2str(iter),...
            '     functional = ',num2str(J(iter)),...
            '     regularity = ',num2str(energie(iter)),...
            '     attache = ', num2str(attache(iter)),...
            '     stepsize = ', num2str(stepsize)]);
    end
    if dosaveiter
        feval(saveiter,x,J)
    end

%     fprintf('loop= %d, loopbreak= %d\n',loop,loopbreak);
% 	fprintf('ratio = %e, breakratio = %e\n',(J(iter)-J(iter+1))/(J(1)-J(iter+1)), breakratio);
    
    % termination criteria
    if (J(iter)-J(iter+1)) < breakratio*(J(1)-J(iter+1)) || loop==loopbreak
% 	fprintf('loop= %d, loopbreak= %d\n',loop,loopbreak);
% 	fprintf('ratio = %e, breakratio = %e\n',(J(iter)-J(iter+1))/(J(1)-J(iter+1)), breakratio);
        return
    end
end

if verbosemode
    disp('maximum number of iterations exceeded');
end
