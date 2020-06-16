function [layer,v_] = IzNeuronUpdate(layer,i,t,Dmax,v_)
% Updates membrane potential v and reset rate u for neurons in layer i
% using Izhikevich's neuron model. Dmax is the maximum
% conduction delay

%%

%%

dt = 0.1526531;
% Calculate current from incoming spikes
for j=1:length(layer)
   S = layer{i}.S{j};
   if ~isempty(S)
      firings = layer{j}.firings;
      if ~isempty(firings)
         % Find incoming spikes (taking account of propagation delays)
         delay = layer{i}.delay{j};
         F = layer{i}.factor{j};
         % Sum current from incoming spikes
         k = size(firings,1);
         while (k>0 && firings(k,1)>t-(Dmax+1))
            spikes = (delay(:,firings(k,2))==t-firings(k,1));
            if ~isempty(layer{i}.I(spikes))
               layer{i}.I(spikes) = layer{i}.I(spikes)+S(spikes,firings(k,2))*F;
            end
            k = k-1;
         end;
         % Don't let I go below zero (shunting inhibition)
         % layer{i}.I = layer{i}.I.*(layer{i}.I > 0);
      end
   end
end
% Update v and u using Izhikevich's model in increments of dt

for k=1:1/dt
   v = layer{i}.v;
   u = layer{i}.u;
   layer{i}.v = v+(dt*(0.04*v.^2+5*v+140-u+layer{i}.I));
   layer{i}.u = u+(dt*(layer{i}.a.*(layer{i}.b.*v-u)));
   % Reset neurons that have spiked
   fired = find(layer{i}.v>=30); % indices of spikes
   v_ = layer{i}.v;
   if ~isempty(fired)
      layer{i}.firings = [layer{i}.firings ; t+0*fired, fired];
      layer{i}.v(fired) = layer{i}.c(fired);
      layer{i}.u(fired) = layer{i}.u(fired)+layer{i}.d(fired);
   end  
    
end
end