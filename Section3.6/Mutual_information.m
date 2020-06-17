function [ MI ] = Mutual_information( X,S )

  % Find mutual information MI
  
  N = size(S, 1);
  
  MI = get_entropy(S(X, :)) + get_entropy(S(setdiff(1:N, X), :)) - get_entropy(S);

end

