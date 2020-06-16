function SWI = SmallWorldIndex(CIJ)
% Computes the small-world index of the graph with connection matrix CIJ
% Self-connections are ignored, as are cyclic paths

N = length(CIJ); % number of nodes
K = sum(sum(CIJ))/length(CIJ); % average degree

[G,CC] = clustcoef(CIJ); % clustering coefficient
[RR,DD] = breadthdist(CIJ); % distance matrix
DD = DD+diag(inf(1,length(DD))); % ignore self-connections and cycles
PL = charpath(DD); % characteristic path length

% Small world index
CCs = CC/(K/N);
PLs = PL/(log(N)/log(K));
SWI = CCs/PLs;