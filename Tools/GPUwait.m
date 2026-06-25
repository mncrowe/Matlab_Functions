function GPUwait(options)
% Waits until the specified 'available_fraction' of the GPU memory is free.
% memory is checked every 'check_interval' seconds and waiting stops
% automatically after 'timeout' seconds.

arguments
    options.timeout            (1,1) double  = inf   % measured in seconds
    options.available_fraction (1,1) double  = 0.8   % between 0 and 1
    options.check_interval     (1,1) double  = 60    % measured in seconds
    options.print_output       (1,1) logical = true  % display output to screen
end

time = tic;

fraction = get_fraction();

while toc(time) < options.timeout && fraction < options.available_fraction

    pause(options.check_interval)

    fraction = get_fraction();

    if options.print_output
        disp(['Wait time: ' num2str(round(10 * toc(time)) / 10) ' seconds,' ... 
            ' Available memory: ' num2str(round(1000 * fraction) / 10) '%'])
    end

end

end

function fraction = get_fraction()

    gpu = gpuDevice;
    fraction = gpu.AvailableMemory / gpu.TotalMemory;

end