function [ MI ] = Mutual_information( X,S )

  % Find mutual information MI
  
  N = size(S, 1);
  
  MI = find_entropy(S(X, :)) + find_entropy(S(setdiff(1:N, X), :)) - find_entropy(S);

end

