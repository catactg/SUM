function loadstruct(sname,vars)

% each field of structure is loaded as a workspace variable

if evalin('caller',['isstruct(',sname,')'])
    if nargin == 1
        vars = evalin('caller',['fieldnames(',sname,')']);
    end
    for n = 1:length(vars)
        fld = vars{n};
        if evalin('caller',['isfield(',sname,',''',fld,''')']) && ~evalin('caller',['isempty(',sname,'.',fld,')'])
            evalin('caller',[fld,' = getfield(',sname,',''',fld,''');']);
        end
    end
end
