% 'exoShape" is a derived software by Stanley Durrleman, Copyright (C) INRIA (Asclepios team), All Rights Reserved, 2006-2009, version 1.0
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

function makewrl(wrlfile,s)

% build vrml file from structure s

fwrlout = fopen(wrlfile,'w');
% header
fprintf(fwrlout,'#VRML V2.0 utf8\n');
fprintf(fwrlout,'Background { skyColor 1 1 1 }\n');
%fprintf(fwrlout,'Group {\n');
%fprintf(fwrlout,'children [\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% display parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

templatebyu = 0;                % use byu file as template data
show = {'x','phi','y'};         % elements to be plotted
xcolor = [0;0;1];               % color for template data
xrigcolor = [1;1;1];            % color for xrig
targetbyu = 0;                  % use byu file as target data
ycolor = [.1;.5;0];             % color for target data
phicolor = [1;0;0];             % color for deformation of template
showtraj = 0;                   % plot trajectories
showgrid = 0;
targetcolors = [1,0,0,1,1,0,1,0;...
    0,1,0,1,1,1,0,0;...
    0,0,1,1,0,1,1,0];
gridsize = 10;
transparent = 1;                % transparent display
showanim = 1;                   % show animation of deformation


% update by variables of structure s
loadstruct('s')
loadstruct('target{1}')

if showgrid
    fprintf(fwrlout,'Shape {\n');
    fprintf(fwrlout,'appearance Appearance {\n');
    fprintf(fwrlout,'material Material {\n');
    fprintf(fwrlout,'diffuseColor 1.0 1.0 1.0\n');
    fprintf(fwrlout,'}\n');
    fprintf(fwrlout,'}\n');
    fprintf(fwrlout,'geometry IndexedLineSet {\n');
    fprintf(fwrlout,'coord ');
    if showanim
        fprintf(fwrlout,'DEF CD_3 ');
    end
    fprintf(fwrlout,'Coordinate {\n');
    fprintf(fwrlout,'point [\n');
    gridstep = max(max(x) - min(x))/gridsize;
    g1 = min(x(1,:)):gridstep:max(x(1,:));
    s1 = length(g1);
    g2 = min(x(2,:)):gridstep:max(x(2,:));
    s2 = length(g2);
    g3 = min(x(3,:)):gridstep:max(x(3,:));
    s3 = length(g3);
    g1 = repmat(g1',1,s2*s3);
    g2 = kron(g2,ones(s3,s1));
    g3 = kron(ones(s1,s2),g3);
    g1 = g1(:)';
    g2 = g2(:)';
    g3 = g3(:)';
    g = [g1;g2;g3];
    gev = zeros([size(g),T]);
    gev(:,:,1) = flow(s,g,0,1);
    for t = 1:T-1
        gev(:,:,t+1) = flow(s,gev(:,:,t),t,t+1);
    end
    fprintf(fwrlout,'%f %f %f,\n',gev(:,:,1));
    fprintf(fwrlout,']\n}\n');
    fprintf(fwrlout,'coordIndex [\n');
    indices = 0:s1*s2*s3-1;
    eval(['fprintf(fwrlout,''',strcat(repmat('%d, ',1,s1)),' -1,\n'',indices);'])
    indices = reshape(indices,[s1,s2,s3]);indices = shiftdim(indices,1);indices = indices(:)';
    eval(['fprintf(fwrlout,''',strcat(repmat('%d, ',1,s3)),' -1,\n'',indices);'])
    indices = reshape(indices,[s3,s1,s2]);indices = shiftdim(indices,1);indices = indices(:)';
    eval(['fprintf(fwrlout,''',strcat(repmat('%d, ',1,s2)),' -1,\n'',indices);'])
    fprintf(fwrlout,']\n');
    fprintf(fwrlout,'}\n');
    fprintf(fwrlout,'},\n');
end


if showtraj
    for k = 1:ntargets
        fprintf(fwrlout,'Shape {\n');
        fprintf(fwrlout,'appearance Appearance {\n');
        fprintf(fwrlout,'material Material {\n');
        fprintf(fwrlout,'emissiveColor %f %f %f\n',targetcolors(:,k));
        fprintf(fwrlout,'}\n');
        fprintf(fwrlout,'}\n');
        fprintf(fwrlout,'geometry IndexedLineSet {\n');
        fprintf(fwrlout,'coord ');
        if showanim
            fprintf(fwrlout,'DEF CD_%d ',k+1);
        end
        fprintf(fwrlout,'Coordinate {\n');
        fprintf(fwrlout,'point [\n');
        ech = [ones(1,T),2:T];
        if showanim
            ech1 = ones(1,T);
        else
            ech1 = 1:T;
        end
        nvx = length(s.target{k}.vx);
        fprintf(fwrlout,'%f %f %f,\n',[...
            reshape(squeeze(X(1,s.target{k}.vx,ech1))',1,nvx*T);...
            reshape(squeeze(X(2,s.target{k}.vx,ech1))',1,nvx*T);...
            reshape(squeeze(X(3,s.target{k}.vx,ech1))',1,nvx*T)]);...
            fprintf(fwrlout,']\n}\n');
        fprintf(fwrlout,'coordIndex [\n');
        %eval(['fprintf(fwrlout,''',strcat(repmat('%d, ',1,T)),' -1,\n'',reshape(0:floor(nvx*.1)*T-1,T,floor(nvx*.1)));'])
        eval(['fprintf(fwrlout,''',strcat(repmat('%d, ',1,T)),' -1,\n'',reshape(0:nvx*T-1,T,nvx));'])
        fprintf(fwrlout,']\n');
        fprintf(fwrlout,'}\n');
        fprintf(fwrlout,'},\n');
    end
end

for i = 1:length(show)
    it = show{i};
    fprintf(fwrlout,'Transform {\n');
    fprintf(fwrlout,'children [\n');
    fprintf(fwrlout,'Shape {\n');
    fprintf(fwrlout,'geometry IndexedFaceSet {\n');
    fprintf(fwrlout,'solid FALSE\n');

    % write vertices data

    fprintf(fwrlout,'coord ');

    switch(it)
        case 'y'
            if targetbyu
                [V,F] = readbyu(targetbyu);
            else
                V = y;
                F = vy;
            end
        case 'x'
            if templatebyu
                [V,F] = readbyu(templatebyu);
            else
                V = x;
                F = vx;
            end
        case 'xrig'
            if templatebyu
                [V,F] = readbyu(templatebyu);
            else
                V = x;
                F = vx;
            end
            V = transmatrix * V + repmat(transvector,1,size(V,2));
        case 'phi'
            if templatebyu
                [V,F] = readbyu(templatebyu);
            else
                V = x;
                F = vx;
            end
            V = transmatrix * V + repmat(transvector,1,size(V,2));
            Vev = zeros([size(V),T]);
            Vev(:,:,1) = V;
            for t = 1:T-1
                Vev(:,:,t+1) = flow(s,Vev(:,:,t),t,t+1);
            end
            V = Vev(:,:,T);
            if showanim
                fprintf(fwrlout,'DEF CD_1 ');
            end
    end

    fprintf(fwrlout,'Coordinate {\n');
    fprintf(fwrlout,'point [');
    fprintf(fwrlout,'\n%f %f %f,',V(:,1:end-1));
    fprintf(fwrlout,'\n%f %f %f',V(:,end));

    % write faces data
    fprintf(fwrlout,'\n]\n}\ncoordIndex [');
    fprintf(fwrlout,'\n%d %d %d -1,',F(:,1:end-1)-1); % "-1" : indices start at 0
    fprintf(fwrlout,'\n%d %d %d -1',F(:,end)-1);
    fprintf(fwrlout,'\n]\n');

    % write color data
    clr = eval([it,'color']);
    fprintf(fwrlout,'}\n');
    fprintf(fwrlout,'appearance Appearance {\nmaterial Material {\n');
    fprintf(fwrlout,'diffuseColor %f %f %f\n',clr);
    if transparent
        fprintf(fwrlout,'transparency 0.5\n');
    end
    fprintf(fwrlout,'}\n}\n}\n]\n}\n');
end

if sum(strcmp(show,'phi')) & showanim
    fprintf(fwrlout,'DEF COD_INT1 CoordinateInterpolator {\n');
    fprintf(fwrlout,'key [ \n');
    fprintf(fwrlout,'%f ',.75*(0:T-1)/(T-1));
    fprintf(fwrlout,']\n');
    fprintf(fwrlout,'keyValue [\n\n');
    fprintf(fwrlout,'%f %f %f,\n',Vev);
    fprintf(fwrlout,']\n}\n');
    fprintf(fwrlout,'DEF TIMER1 TimeSensor {\n');
    fprintf(fwrlout,'cycleInterval 10\n');
    fprintf(fwrlout,'loop TRUE\n}\n');
    fprintf(fwrlout,'ROUTE TIMER1.fraction_changed TO COD_INT1.set_fraction\n');
    fprintf(fwrlout,'ROUTE COD_INT1.value_changed TO CD_1.set_point\n');
end

if showtraj & showanim
    for k = 1:ntargets
        fprintf(fwrlout,'DEF COD_INT%d CoordinateInterpolator {\n',k+1);
        fprintf(fwrlout,'key [ \n');
        fprintf(fwrlout,'%f ',.75*(0:T-1)/(T-1));
        fprintf(fwrlout,']\n');
        fprintf(fwrlout,'keyValue [\n\n');
        nvx = length(s.target{k}.vx);
        for t=1:T
            fprintf(fwrlout,'%f %f %f,\n',[...
                reshape(squeeze(X(1,s.target{k}.vx,ech(t:T+t-1)))',1,nvx*T);...
                reshape(squeeze(X(2,s.target{k}.vx,ech(t:T+t-1)))',1,nvx*T);...
                reshape(squeeze(X(3,s.target{k}.vx,ech(t:T+t-1)))',1,nvx*T)]);
        end
        fprintf(fwrlout,']\n}\n');
        fprintf(fwrlout,'DEF TIMER%d TimeSensor {\n',k+1);
        fprintf(fwrlout,'cycleInterval 10\n');
        fprintf(fwrlout,'loop TRUE\n}\n');
        fprintf(fwrlout,'ROUTE TIMER%d.fraction_changed TO COD_INT%d.set_fraction\n',[k+1 k+1]);
        fprintf(fwrlout,'ROUTE COD_INT%d.value_changed TO CD_%d.set_point\n',[k+1 k+1]);
    end
end

if showgrid & showanim
    fprintf(fwrlout,'DEF COD_INT3 CoordinateInterpolator {\n');
    fprintf(fwrlout,'key [ \n');
    fprintf(fwrlout,'%f ',.75*(0:T-1)/(T-1));
    fprintf(fwrlout,']\n');
    fprintf(fwrlout,'keyValue [\n\n');
    fprintf(fwrlout,'%f %f %f,\n',gev);
    fprintf(fwrlout,']\n}\n');
    fprintf(fwrlout,'DEF TIMER3 TimeSensor {\n');
    fprintf(fwrlout,'cycleInterval 10\n');
    fprintf(fwrlout,'loop TRUE\n}\n');
    fprintf(fwrlout,'ROUTE TIMER3.fraction_changed TO COD_INT3.set_fraction\n');
    fprintf(fwrlout,'ROUTE COD_INT3.value_changed TO CD_3.set_point\n');
end

%fprintf(fwrlout,']\n}');
fclose(fwrlout);
