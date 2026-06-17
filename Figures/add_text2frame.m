function add_text2frame(frame_name, text, T, position, save_name, sig_figs)
% Adds text to frames in the specified location
%
% - frame_name: location and name of frames, do not include trailing numbers
% - text: string with text to include
% - T: vector of numbers to add to end of string, optional (default: [])
% - position: position of text, from top left corner (default: [0 0])
% - save_name: location and name of output frames (default: frame_name)
% - sig_figs: number of significant figures for displaying T (default: 3)

arguments
    frame_name char
    text char
    T (:,:) double        = []
    position (2,1) double = [0, 0]
    save_name char        = frame_name
    sig_figs (1,1) int64  = 3
end

N = length(dir([frame_name '*']));
d = 1 + floor(log(N) / log(10));

if nargin > 2
    if length(T) ~= 1 & length(T) ~= N
        error('Length of T inconsistent with number of frames')
    end
end

for i = 1:N

    file_name = [frame_name '_' pad_zeros(i,d) '.png'];
    savename = [save_name '_' pad_zeros(i,d) '.png'];

    if length(T) == N
        text_string = [text sig_fig_str(T(i), sig_figs)];
    else
        text_string = [text sig_fig_str(T, sig_figs)];
    end

    [I, map] = imread(file_name);
    I = reshape(map(I + 1, :), [size(I) 3]);
    
    RGB = insertText(I, position', text_string, 'BoxOpacity', 0, ...
        'TextColor', 'black', 'FontSize', 18);
    
    imwrite(RGB, savename)

end

end