% Test TMP file reading and writing by saving to the given filename,
% reading it back, and checking for equality.
function ok = testTMP(temp_filename)

% Some decent size that's not too big.
% Pick odd sizes on purpose.
height = 767;
width = 1021;
channels = 3;
frames = 7;

sz = [height, width, channels, frames];

fprintf('Testing single...');
a = single(rand(sz));
saveTMP(a, temp_filename);
b = loadTMP(temp_filename);

if ~isequal(a, b)
    fprintf('FAILED\n');
    ok = false;
    return;
else
    fprintf('passed.\n');
end

fprintf('Testing double...');
a = double(rand(sz));
saveTMP(a, temp_filename);
b = loadTMP(temp_filename);

if ~isequal(a, b)
    fprintf('FAILED\n');
    ok = false;
    return;
else
    fprintf('passed.\n');
end

fprintf('Testing uint8...');
a = uint8(rand(sz));
saveTMP(a, temp_filename);
b = loadTMP(temp_filename);

if ~isequal(a, b)
    fprintf('FAILED\n');
    ok = false;
    return;
else
    fprintf('passed.\n');
end

fprintf('Testing int8...');
a = int8(rand(sz));
saveTMP(a, temp_filename);
b = loadTMP(temp_filename);

if ~isequal(a, b)
    fprintf('FAILED\n');
    ok = false;
    return;
else
    fprintf('passed.\n');
end

fprintf('Testing uint16...');
a = uint16(rand(sz));
saveTMP(a, temp_filename);
b = loadTMP(temp_filename);

if ~isequal(a, b)
    fprintf('FAILED\n');
    ok = false;
    return;
else
    fprintf('passed.\n');
end

fprintf('Testing int16...');
a = int16(rand(sz));
saveTMP(a, temp_filename);
b = loadTMP(temp_filename);

if ~isequal(a, b)
    fprintf('FAILED\n');
    ok = false;
    return;
else
    fprintf('passed.\n');
end

fprintf('Testing uint32...');
a = uint32(rand(sz));
saveTMP(a, temp_filename);
b = loadTMP(temp_filename);

if ~isequal(a, b)
    fprintf('FAILED\n');
    ok = false;
    return;
else
    fprintf('passed.\n');
end

fprintf('Testing int32...');
a = int32(rand(sz));
saveTMP(a, temp_filename);
b = loadTMP(temp_filename);

if ~isequal(a, b)
    fprintf('FAILED\n');
    ok = false;
    return;
else
    fprintf('passed.\n');
end

fprintf('Testing uint64...');
a = uint64(rand(sz));
saveTMP(a, temp_filename);
b = loadTMP(temp_filename);

if ~isequal(a, b)
    fprintf('FAILED\n');
    ok = false;
    return;
else
    fprintf('passed.\n');
end

fprintf('Testing int64...');
a = int64(rand(sz));
saveTMP(a, temp_filename);
b = loadTMP(temp_filename);

if ~isequal(a, b)
    fprintf('FAILED\n');
    ok = false;
    return;
else
    fprintf('passed.\n');
end

fprintf('All tests PASSED.');
ok = true;