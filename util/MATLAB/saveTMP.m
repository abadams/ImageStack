% Write an up to 4-dimensional matrix 'a' to 'filename' in the TMP format.
%
% a should be in the MATLAB convention: [height width channels frames]
% a will be permuted to the TMP file convention:
%   [width height frames channels]
%
% Format specification is here:
% https://github.com/abadams/ImageStack/blob/master/src/FileTMP.cpp
function [] = saveTMP(a, filename)

if ndims(a) > 4
   error('TMP only supports up to 4 dimensions.');
end

switch class(a)
    case 'single'
        type_code = int32(0);
    case 'double'
        type_code = int32(1);
    case 'uint8'
        type_code = int32(2);
    case 'int8'
        type_code = int32(3);
    case 'uint16'
        type_code = int32(4);
    case 'int16'
        type_code = int32(5);
    case 'uint32'
        type_code = int32(6);
    case 'int32'
        type_code = int32(7);
    case 'uint64'
        type_code = int32(8);
    case 'int64'
        type_code = int32(9);
    otherwise
        error('TMP only supports single, double, and {u}int{8,16,32,64} numeric matrices.');
end

% Pick a type string for MATLAB's fwrite() based on the type code.
type_strings = {'*float32', '*float64', '*uint8', '*int8', ...
    '*uint16', '*int16', '*uint32', '*int32', '*uint64', '*int64'};
type_string = type_strings{type_code + 1};

% Undo permutation from loadTMP:
% Turn MATLAB [height, width, channels, frames] back into TMP format's
% [width, height, frames, channels].
a = ipermute(a, [2 1 4 3]);

% Now grab its dimensions.
width = int32(size(a, 1));
height = int32(size(a, 2));
frames = int32(size(a, 3));
channels = int32(size(a, 4));

% Turn it into a single giant vector.
a = a(:);

% Write the file.
fid = fopen(filename, 'w');
fwrite(fid, [width height frames channels type_code], '*int32');
fwrite(fid, a, type_string);
fclose(fid);
