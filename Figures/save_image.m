function save_image(field, savename, ext, options)
% Saves the given field as an image 'savename.ext' using specified colormap
%
% - field: field(x,y), with x as horizontal axis (column)
% - savename: string, name of saved image, excluding file extension
% - ext: string, file extension (default: 'png')
% - colormap: colormap, matrix of size [N 3] describing colour scale (default: cmap())
% - scale: vector [v1 v2], fix the ends of the colorbar to these values
% - NaN_col: color of NaN values, [R G B] vector (default: [0.5 0.5 0.5], grey)
% - add_contour: adds contours to image (default: false)
% - num_contour: number of contours (default: 5)
% - col_contour: contour colour, [R G B] vector (default: [1 1 1], black)
% - text: character array of text to add to figure (default: '')
% - text_col: text colour (default: 'black')
% - text_pos: position of text from top left corner (default: [0 0])
% - text_size: text size (default: 18)

% ----------------------------------------------------------------------------
% Note: Enter field' for traditional matrix approach with first index as row.
% ----------------------------------------------------------------------------

arguments
    field               (:,:) double      = gallery("circul", 128)
    savename                  char        = 'img'
    ext                       char        = 'png'
    options.colormap    (:,3) double      = cmap()
    options.scale       (1,:) double      = [min(min(field)) max(max(field))]
    options.NaN_col     (1,3) double      = [0.5 0.5 0.5]
    options.add_contour (1,1) logical     = false
    options.num_contour (1,1) double      = 5
    options.col_contour (1,3) double      = [0 0 0]
    options.text              char        = ''
    options.text_col          {CheckCol}  = 'black'
    options.text_pos    (1,2) double      = [0 0]
    options.text_size   (1,1) double      = 18
end

% Set values outside scale limts to min/max value:
field = min(max(field, options.scale(1)), options.scale(2));

% Normalise values between 0 and 1:
field = squeeze((field' - options.scale(1)) / (options.scale(2) - options.scale(1)));

% Flip field:
field = field(end:-1:1, :);

% Create RGB array by interpolating colormap onto normalised field:
field_rgb = interp1(linspace(0, 1, length(options.colormap)), options.colormap, field, 'nearest');

% Replace NaN entries with NaN colour:
field_rgb(isnan(field_rgb)) = kron(options.NaN_col, ones(1, sum(isnan(field), "all")));

% Add contours to image:
if options.add_contour

    % Define coutour heights:
    n = options.num_contour;
    contours = 0.5/n:1/n:1-0.5/n;

    % Get array of x and y positions along contours:
    C = max(1, round(contourc(field', contours)));
    C(1, :) = min(size(field, 1), C(1, :));
    C(2, :) = min(size(field, 2), C(2, :));

    % Add contours to image:
    I = sub2ind([size(field) 3], [C(1,:) C(1,:) C(1, :)], [C(2,:) C(2,:) C(2, :)], kron([1 2 3], ones(1, length(C(1,:)))));
    field_rgb(I) = kron(options.col_contour, ones(1, length(C(1,:))));
end

% Add text to image:
if isa(options.text_col, "double"); options.text_col = reshape(options.text_col, [1 3]); end
if ~isempty(options.text)
    field_rgb = insertText(field_rgb, options.text_pos, options.text, 'BoxOpacity', 0, ...
        'TextColor', options.text_col, 'FontSize', options.text_size);
end

% Save figure:
imwrite(field_rgb, [savename '.' ext]);

end

% Validation for text colour input:
function CheckCol(col)
    if ~isa(col, "char") && ~(isa(col, "double") && numel(col) == 3)
        eid = 'Colour:DataType';
        msg = 'text_col should be either a character array (e.g. ''black'') or 1 x 3 vector.';
        error(eid, msg)
    end
end


