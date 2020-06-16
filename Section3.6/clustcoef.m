function [gamma,gammaG] = clustcoef(CIJ)

% input:  
%           CIJ    = binary connection/adjacency matrix
% output: 
%           gamma  = clustering coefficient for each vertex
%           gammaG = mean clustering coefficient for entire graph
%
% Compute clustering coefficient.
% CIJ must be binary and must not have any self-connections.
% Note: Watts/Strogatz use un-directed graphs, here we use directed graphs.
% Find the immediate neighbors of a given vertex (both in and out), 
% then determine how many connections exist between them out 
% of all possible.  Note: immediate neighbors are all those vertices that
% either send or receive an edge from the target vertex.
%
% Olaf Sporns, Indiana University, 2002/2007

% ensure CIJ is binary...
CIJ = double(CIJ~=0);

% initialize
N = size(CIJ,1);
gamma = zeros(1,N);

% loop over all vertices
for v=1:N
   [nb] = find(CIJ(v,:) + CIJ(:,v)');
   lnb = length(nb);
   if (lnb>1)
      gamma(v) = sum(sum(CIJ(nb,nb)))./(lnb^2-lnb);
   end;
   if (lnb==1)
      gamma(v) = 0;
   end;
end;

% get average over all nodes
% Note: for some applications the median may be more appropriate
gammaG = mean(gamma);
