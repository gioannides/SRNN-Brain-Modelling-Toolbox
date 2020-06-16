function  [R,D] = breadthdist(CIJ)

% input:  
%           CIJ = connection/adjacency matrix
% output: 
%           R   = reachability matrix
%           D   = distance matrix

% This function is potentially less memory-hungry than 'reachdist.m',
% particularly if the characteristic path length is rather large.
%
% Olaf Sporns, Indiana University, 2002/2007

N = size(CIJ,1);

D = zeros(N);
for i=1:N
   D(i,:) = breadth(CIJ,i);
end;

% replace zeros with 'Inf's
D(D==0) = Inf;

% construct R
R = double(D~=Inf);

