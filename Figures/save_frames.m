function save_frames(f, savename, options)
% saves 3D field as N frames f(:,:,i) for i = {1,..,N}
%
% - f: f(x,y,i), with x as horizontal axis (column)
% - savename: string, name of saved images, suffixed by '_00i', excluding file extension (default: 'frame')
% - ext: string, file extension (default: 'png')
% - colormap: colormap, matrix of size [N 3] describing colour scale (default: cmap())
% - scale: vector [v1 v2], fix the ends of the colorbar to these values (default: [-M M], M = max(|f|))
% - NaN_col: color of NaN values, [R G B] vector (default: [0.5 0.5 0.5], grey)
% - add_contour: adds contours to frames (default: false)
% - num_contour: number of contours (default: 5)
% - col_contour: contour colour, [R G B] vector (default: [1 1 1], black)
% - text: character array of text to add to figure (default: '')
% - text_col: text colour (default: 'black')
% - text_pos: position of text from top left corner (default: [0 0])
% - text_size: text size (default: 18)
% - text_num: array of numbers to add to text string (default: [])
% - text_sigfig: number of significant figures for text_num (default: 3)

arguments
    f (:,:,:)                 double
    savename                  char        = 'frame'
    options.ext               char        = 'png'
    options.colormap    (:,3) double      = cmap()
    options.scale       (1,:) double      = max(max(max(abs(f)))) * [-1 1];
    options.NaN_col     (1,3) double      = [0.5 0.5 0.5]
    options.add_contour (1,1) logical     = false
    options.num_contour (1,1) double      = 5
    options.col_contour (1,3) double      = [0 0 0]
    options.text              char        = ''
    options.text_col                      = 'black'
    options.text_pos    (1,2) double      = [0 0]
    options.text_size   (1,1) double      = 18
    options.text_num    (1,:) double      = []
    options.text_sigfig (1,1) int64       = 3 
end

s = size(f); N = s(3);
d = 1 + floor(log(N) / log(10));

for i = 1:N
    fi = squeeze(f(:,:,i));

    if numel(options.text_num) == N
        num = sig_fig_str(options.text_num(i), options.text_sigfig);
    else
        num = '';
    end

    text = [options.text num];

    save_image(fi, [savename '_' pad_zeros(i,d)], options.ext, ...
        colormap=options.colormap, scale=options.scale, NaN_col=options.NaN_col, ...
        add_contour=options.add_contour, num_contour=options.num_contour, ...
        col_contour=options.col_contour, text=text, text_col=options.text_col, ...
        text_pos=options.text_pos, text_size=options.text_size)
end

end

