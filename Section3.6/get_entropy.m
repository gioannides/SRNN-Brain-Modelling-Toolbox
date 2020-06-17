function entropy = get_entropy(S)

  N = size(S, 1);
  entropy = 0.5*log( ((2*pi*exp(1))^N) .* det(cov(S')) + eps );

end
