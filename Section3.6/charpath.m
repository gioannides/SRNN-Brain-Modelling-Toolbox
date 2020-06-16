function  [lambda,ecc,radius,diameter] = charpath(D)

% input:  
%           D          distance matrix
% output: 
%           lambda     characteristic path length
%           ecc        eccentricity (for each vertex)
%           radius     radius of graph
%           diameter   diameter of graph
%
% Characteristic path length is calculated as the global mean of the
% distance matrix D, not taking into account any 'Infs' but including the
% distances on the main diagonal.
%
% Olaf Sporns, Indiana University, 2002/2007
% =========================================================================

% Mean of finite entries of D(G)
lambda = sum(sum(D(D~=Inf)))/length(nonzeros(D~=Inf));

% Eccentricity for each vertex (note: ignore 'Inf') 
ecc = max(D.*(D~=Inf),[],2);

% Radius of graph
radius = min(ecc);  % but what about zeros?

% Diameter of graph
diameter = max(ecc);
