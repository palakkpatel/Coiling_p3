function [faces, vertices] = stlread(filename)
% Read the .stl file and build the double connected linked list.
% Input: filename is a string. It should be 'name.extention'.
% Output:
%   vertices: v_cnt * 3 matrix, each row is a vertex
%   faces: f_cnt * 3 matrix. Each row is the vertex_id of the triangluar face.

% Open the file, assumes STL ASCII format.
fid = fopen(filename, 'r');
if fid == -1 
    error('File could not be opened, check name or path.')
end

% Total number of vertices, edges and faces.
v_cnt = 0;
f_cnt = 0;

% Read the data. Get the count of vertices, edges and faces.
while feof(fid) == 0
    % Read one line of the file.
    tline = fgetl(fid);
    fword = sscanf(tline, '%s');
    % A new face
    if strncmpi(fword, 'f', 1) == 1
        f_cnt = f_cnt + 1;
    else
        if strncmpi(fword, 'v', 1) == 1
            v_cnt = v_cnt + 1;
        end
    end
end

% Allocate memory for outputs.
vertices = zeros(v_cnt, 3);
faces = zeros(f_cnt, 3);

% Go back to the beginning of the file.
status = fseek(fid, 0, 'bof');
if status < 0
    error('fseek() fails!\n')
end

% Temporary index for vertices, edges and faces.
v_idx = 0;
f_idx = 0;

% Temporary vertex.
tempv = zeros(1,3);

% Whether we are inside a face.
% in_face = 0: outside any face
% in_face = 1: first vertex of a face
% in_face = 2: second vertex
% in_face = 3: third vertex
in_face = 0;

% Index of current vertex, could be less than v_idx.
cv = 0;

% Read the file again and fill in the matrices.
while feof(fid) == 0
    % Read one line of the file.
    tline = fgetl(fid);
    fword = sscanf(tline, '%s');
    if strncmpi(fword, 'f', 1) == 1
        f_idx = f_idx + 1;
        in_face = 1;
    end
    if strncmpi(fword, 'v', 1) == 1
        tempv = sscanf(tline, '%*s %f %f %f');
        cv = -1;
        for ii = v_idx:-1:1
            % Check whether this vertex has appeared before.
            if cmpv(tempv, vertices(ii,:)) == 1
                cv = ii;
                break;
            end
        end
        if cv < 0
            v_idx = v_idx + 1;
            vertices(v_idx, :) = tempv;
            cv = v_idx;
        end
        faces(f_idx,in_face) = cv;
        in_face = in_face + 1;
    end % if strncmpi(fword, 'v', 1) == 1
end

vertices = vertices(1:v_idx, :);
% Close the file.
fclose(fid);
