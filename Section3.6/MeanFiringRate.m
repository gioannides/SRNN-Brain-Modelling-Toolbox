function [MFR] = MeanFiringRate(firings,ws,ds)
   %ws = 50; % window size - must be greater than or equal to Dmax
   %ds = 20; % slide window by ds
   Tmax = max(firings(:,1));
   MFR = zeros(1,ceil(Tmax/ds));   
   % MF(i,j) is the firing rate for the ith largest module for the
   % ws data points up to but not including j
  for m=1:8
    for j=1:ds:Tmax 
        MFR(m,ceil(j/ds)) = sum(firings(:,1) >= j-ws & firings(:,1) < j...
            & firings(:,2) <= m*100 & firings(:,2) >= 100*(m-1)+1 ) / ws ; 
    end
   end
end

