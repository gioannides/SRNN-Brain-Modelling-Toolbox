function [ newCIJ ] = ExtendNetwork(CIJ, n)
% Expands a network by adding the new n to the bottom right
    if size(CIJ) == 0
        newCIJ = n;
    else
        [i,j] = size(n);
        [I,J] = size(CIJ);
        newCIJ = [[CIJ,zeros(I,j)];[zeros(i,J),n]];
    end
end

