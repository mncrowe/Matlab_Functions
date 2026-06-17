function save_figure(filename,ext,figure_no,renderer,resolution)
% Saves open matlab figure windows as files, may save multiple figures with different parameters
%
% - filename: name of file(s), text input (default: 'fig')
% - ext: extension, e.g. 'png', 'eps', text input (default: 'png')
% - figure_no: figure number(s) to save, integer (default: 1)
% - renderer: renderer(s) used to create file:
%           - 0: set based on file extension
%           - 1: opengl (default for non-eps and non-pdf)
%           - 2: vector (default for eps and pdf)
% - resolution: resolution(s) in DPI, opengl only, integer (default: screen resolution)
%
% -------------------------------------------------------------------------
%
% Note: text inputs should be character arrays (e.g. 'figure_1'), strings
% (e.g. "figure_1") or arrays of strings (e.g. ["figure_1", "figure_2"].
% Integer inputs may be arrays of integers which match the length of
% filename or of length 1. Parameters of length 1 will apply to all figures
% generated.

arguments
    filename   (1,:) string                                    = 'fig'
    ext        (1,:) string {MatchLength(filename,ext)}        = 'png'
    figure_no  (1,:) int64  {MatchLength(filename,figure_no)}  = 1
    renderer   (1,:) int64  {MatchLength(filename,renderer)}   = 0
    resolution (1,:) int64  {MatchLength(filename,resolution)} = 0
end

N = length(filename);

if length(ext) == 1; ext = repmat(ext, [1 N]); end
if length(figure_no) == 1; figure_no = repmat(figure_no, [1 N]); end
if length(renderer) == 1; renderer = repmat(renderer, [1 N]); end
if length(resolution) == 1; resolution = repmat(resolution, [1 N]); end

disp(' ')

for n = 1:N

    if renderer(n) == 0
        if ext(n) == "pdf" || ext(n) == "eps"
            render = '-vector';
        else
            render = '-opengl';
        end
    end

    if renderer(n) == 1; render = '-opengl'; end
    if renderer(n) == 2; render = '-vector'; end

    res = ['-r' num2str(resolution(n))];

    if ext(n) == "jpg"
        format_type = '-djpeg';
    else
        if ext(n) == "eps"
            format_type = '-depsc';
        else
            format_type = ['-d' char(ext(n))];
        end
    end

    fname = [char(filename(n)) '.' char(ext(n))];
    fnum = ['-f' num2str(figure_no(n))];

    disp('Creating figure:')
    disp([' - Figure: ' fnum])
    disp([' - Name: ' fname])
    disp([' - Type: ' format_type])
    disp([' - Renderer: ' render])
    disp([' - Resolution: ' res])
    disp(' ')
    print(fnum, fname, format_type, render, res)

end

end

function MatchLength(V1, V2)
    if (length(V2) ~= length(V1)) & (length(V2) ~= 1)
        eid = 'Size:LengthMismatch';
        msg = 'Length of input must be 1 or must match the length of filename.';
        error(eid, msg)
    end
end