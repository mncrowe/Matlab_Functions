function name = append_number(name, ext)
% appends number to a filepath to avoid duplicates
%
% - name: filename excluding extension
% - ext: file extension

arguments
    name char
    ext char
end

save_num = '';
save_index = 0;

while exist([name save_num '.' ext], 'file') == 2
    save_index = save_index + 1;
    save_num = ['_' num2str(save_index)];
end

name=[name save_num '.' ext];

end