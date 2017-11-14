function save_vtk_shape3D(filename,v,f, cellscalars)
%scalars = [] - no scalars
%scalars.name = name of the array
%scalars.data = array itself

fid = fopen(filename,'w','b');
fprintf(fid,'# vtk DataFile Version 3.0\x0A');
fprintf(fid,'vtk output\x0A');
fprintf(fid,'ASCII\x0A');
fprintf(fid,'DATASET POLYDATA\x0A');

str = sprintf('POINTS %d float\x0A',size(v,1));
fprintf(fid,str);
%fwrite(fid,reshape(v',numel(v),1),'float');

for i=1:size(v,1)
    fprintf(fid,'%f %f %f \x0A',v(i,1),v(i,2),v(i,3));
end;

str = sprintf('POLYGONS %d %d\x0A',size(f,1),numel(f)+size(f,1));
fprintf(fid,str);

%fwrite(fid,reshape(f1',numel(f1),1),'int');
for i=1:size(f,1)
   fprintf(fid,'3 %d %d %d \x0a',f(i,1)-1,f(i,2)-1,f(i,3)-1);
end;

%scalars
if ~isempty(cellscalars)
    str = sprintf('CELL_DATA %d\x0A',size(f,1));
    fprintf(fid,str);
    str = sprintf('FIELD FieldData 1\x0A',size(f,1));
    fprintf(fid,str);
    str = sprintf('%s 1 %d short\x0A',cellscalars.name, size(f,1));
    fprintf(fid,str);
    
    for i=1:numel(cellscalars.data)
        fprintf(fid,'%d ',cellscalars.data(i));
    end;
end;


fclose(fid);


