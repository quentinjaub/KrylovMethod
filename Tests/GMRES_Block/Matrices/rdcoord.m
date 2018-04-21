function A = rdcoord(infilename)
%
% A = rdcoord(filename)     reads a sparse matrix from a coordinate text file
%                           with the following format:
%
%                            M   N   nz
%                            i_1   j_1   [val_1   [val_img_1]]
%                            i_2   j_2   [val_2   [val_img_2]]
%                            .      .      .          .
%                            .      .      .          .
%                            .      .      .          .
%                            i_nz  j_nz  [val_nz  [val_img_nz]]
%
%
% with one nonzero entry per line.  If only "i" and "j" values appear on each 
% line, then the matrix is assumed to be pattern-only.  One extra floating 
% point value represents a real matrix, and two floating point values represent
% a complex matrix.
%
%
% Author: R. Pozo 
%         (http://math.nist.gov/pozo)

fid = fopen(infilename, 'rt');
if (fid == -1) 
    error('File not found.');
end;

line = fgets(fid);
v = sscanf(line,'%d');
m = v(1);  n = v(2); nz=v(3);

% bogus return

if (nz < 1) 
    return;
end;

% get first line to see if pattern, real, or complex
%
% sscanf() returns a column vector...
%
line = fgets(fid);
A = sscanf(line, '%f');
num_items = size(A,1);
shape = [1 num_items];

T = zeros(nz,num_items);
T(1,:)= A';

% get rest of lines
for i=2:nz,
    line = fgets(fid);
    t = sscanf(line,'%f', num_items);
    T(i,:) =  t';
end

fclose(fid);

% Now construct the sparse matrix...


% pattern only
if (num_items == 2)     
    A = sparse(T(:,1), T(:,2), zeros(nz,1) , m , n);

% real values
elseif (num_items == 3)
	if (max(abs(T(:,3))) == 0.0 )		% really a pattern matrix?
										% if so, convert to nonzero values,
										% otherwise MATLAB will not record
										% sparsity structure information
		A = sparse(T(:,1), T(:,2), ones(nz,1), m , n);
	else
    	A = sparse(T(:,1), T(:,2), T(:,3), m , n);
	end;

% complex values
elseif (num_items == 4)
    A = sparse(T(:,1), T(:,2), T(:,3) + T(:,4)*sqrt(-1), m , n);

end

