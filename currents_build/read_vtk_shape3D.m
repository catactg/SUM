function [v,f,array] = read_vtk_shape3D(filename, arrayname)

array = [];

%only ASCII
data = textread(filename,'%s');

c=1;

while ~strcmp(data{c}, 'POINTS') && c<numel(data)
    c=c+1;
end;

c=c+1;
npts = str2num(data{c});
v = zeros(npts,3);
c=c+2; %skip datatype
for i=1:npts
    v(i,1) = str2num(data{c});
    c=c+1;
    v(i,2) = str2num(data{c});
    c=c+1;
    v(i,3) = str2num(data{c});
    c=c+1;
end;    

while ~strcmp(data{c}, 'POLYGONS') && c<numel(data)
    c=c+1;
end;
c=c+1;

nfaces = str2num(data{c});
c=c+2;

f = zeros(nfaces,3);
for i=1:nfaces
    c=c+1;
    f(i,1) = str2num(data{c})+1;
    c=c+1;
    f(i,2) = str2num(data{c})+1;
    c=c+1;
    f(i,3) = str2num(data{c})+1;
    c=c+1;
end;    

if exist('arrayname','var')
    while ~strcmp(data{c}, arrayname) && c<numel(data)
        c=c+1;
    end;
    c=c+2;
    nscalars = str2num(data{c});
    array = zeros(nscalars,1);
    c=c+1;
    for i=1:nscalars
        c=c+1;
        array(i) = str2num(data{c});
    end;    
end;

