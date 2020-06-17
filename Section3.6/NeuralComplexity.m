function [ DynComplexity ] = NeuralComplexity(S)
  
  % Find integration
  H_X = 0;
  S = aks_diff(S);
  S = aks_diff(S);
  N = size(S, 1);
  
  for i = 1:N    
    H_X = H_X + get_entropy(S(i, :));
  end

  Integration = H_X - get_entropy(S);
  
  % Find Interaction Complexity
  
  %disp(Integration);
  
  IC = 0;
  
  for i = 1:N
    IC = IC + Mutual_information(i,S);
  end

  %disp(IC);
  
  IC = IC - Integration;

  DynComplexity = IC;
  
  
end

