function gpu = GPUdisp()
% Displays properties about available GPU

gpu = gpuDevice;

disp(' ')
disp(['Name: ' gpu.Name])
disp(['Total memory: ' num2str(round(10 * gpu.TotalMemory / 1024^3) / 10) 'Gb'])
disp(['Available memory: ' num2str(round(1000 * gpu.AvailableMemory / gpu.TotalMemory) / 10) '%'])
disp(' ')

if nargout == 0
    gpu = [];
end

end