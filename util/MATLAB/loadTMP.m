% Load an ImageStack TMP file, an up to 4-dimensional array of integer or
% floating point data.
%
% output permutes the dimensions so that the output is:
% [height width channels frames]
%
% Format specification is here:
% https://github.com/mit-gfx/ImageStack/blob/master/src/FileTMP.cpp
function output = loadTMP(filename)

type_strings = {'*float32', '*float64', '*uint8', '*int8', ...
    '*uint16', '*int16', '*uint32', '*int32', '*uint64', '*int64'};

fid = fopen(filename, 'r');
header = fread(fid, 5, '*int32');

width = header(1);
height = header(2);
frames = header(3);
channels = header(4);
type_code = header(5);

if width < 1 || height < 1 || frames < 1 || channels < 1 || ...
        type_code < 0 || type_code > 9
    fclose(fid);
    error('Invalid header: %d %d %d %d %d', ...
        width, height, frames, channels, type_code);
end

num_elements = width * height * frames * channels;
type_string = type_strings{type_code + 1}; % MATLAB is one-based.

output = fread(fid, num_elements, type_string);
fclose(fid);

% Input is stored row major.
output = reshape(output, [width, height, frames, channels]);

% Permute dimensions so that each image is column major, and frames comes
% last.
output = permute(output, [2 1 4 3]);