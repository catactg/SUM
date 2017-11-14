function savestruct(sname,field)

% save workspace variables in structure 
% savestruct('s') only updates fields which are already defined in s 
%
% savestruct('s',field) specify the list of fields

if nargin==1
   field = evalin('caller',['fieldnames(',sname,')']);
end

for n = 1:length(field)
   fld = field{n};
   if ~strcmp(sname,fld)  % prevent from saving structure inside itself..
       evalin('caller',[sname,'.',fld,'=eval(''',fld,''',''eval(''''',sname,'.',fld,''''',''''0;'''')'');']);
   end
end

