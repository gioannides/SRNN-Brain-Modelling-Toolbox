% Created by Eugene M. Izhikevich, February 25, 2003
% STDP by Gaia Silvestri 2013
% The randomly connected spiking neural network undergoes plastic changes

% Excitatory neurons    Inhibitory neurons
Ne=800;                 Ni=200;
n=Ne+Ni;
re=rand(Ne,1);          ri=rand(Ni,1);
a=[0.02*ones(Ne,1);     0.1*ones(Ni,1)];
b=[0.2*ones(Ne,1);      0.2*ones(Ni,1)];
c=[-65*ones(Ne,1);      -65*ones(Ni,1)];
d=[8*ones(Ne,1);        2*ones(Ni,1)];
S=[0.5*rand(Ne+Ni,Ne),  -rand(Ne+Ni,Ni)];

p=0.7;                  % 30% of zero-weight synapses when p=0.7
totzeros=0;             % Count of zero-weight synapses
plasticity=20;          % Plastic changes happen after the first 20 ms of runtime
runtime=1000;           % Runtime in ms
v=-65*ones(Ne+Ni,1);    % Initial values of v
u=b.*v;                 % Initial values of u
firings=[];             % spike timings
LFP=zeros(runtime,3);
ADJ=zeros(runtime,n);
count=0;
count2=0;


%MAKE THE NETWORK SPARSELY CONNECTED 
M=rand(n,n);
for i=1:n
    for jjj=1:n
        if M(i,jjj)>p
             S(i,jjj)=0;
             totzeros=totzeros+1;
        end
        if i==jjj
             S(i,jjj)=0; %no self connections
             totzeros=totzeros+1;
        end
    end
end
SPARSE=S;

for t=1:runtime                    % simulation of runtime ms
  I=[5*randn(Ne,1);2*randn(Ni,1)]; % thalamic input
  fired=find(v>=30);               % indices of spikes
  for jj=1:numel(fired)
      ADJ(t,fired(jj))=1;          % RECORD WHICH NEURONS FIRE AT EACH MS
  end
  if t>plasticity        %t>20 for STDP on, t>runtime for STDP off
      for j=1:n
        count2=count2+1;
        if ADJ(t,j)==1   %neuron j fired at time t
          count=count+1;
          for x=1:20     %range of ms over which synaptic strengthening and weakening occur                           
              for k=1:n                  
                if ADJ(t-x,k)==1  % find out which neurons fired in the previous 20 ms              
                  % synaptic weight of pair undergoes larger increment if dT is smaller and negative
                  S(j,k)=S(j,k).*(1+(0.9*exp(x/20.0))); 
                  % synaptic weight of pair undergoes larger decrement is dT is smaller and positive
                  S(k,j)=S(k,j).*(1+(-0.925*exp(-x/20.0))); 
                  % set a maximum value for synaptic weights
                  if S(j,k)>2.0    
                      S(j,k)=2.0;
                  end
                  if S(j,k)<-2.0
                      S(j,k)=-2.0;
                  end
                  if S(k,j)>2.0
                      S(k,j)=2.0;
                  end
                  if S(k,j)<-2.0
                      S(k,j)=-2.0;
                  end
                end
              end
          end
        end
      end
  end
  
  firings=[firings; t+0*fired,fired]; % left column is time, right is index of neuron firing out of 1000. 
  v(fired)=c(fired);
  u(fired)=u(fired)+d(fired);
  I=I+sum(S(:,fired),2);
  v=v+0.5*(0.04*v.^2+5*v+140-u+I); % step 0.5 ms
  v=v+0.5*(0.04*v.^2+5*v+140-u+I); % for numerical stability
  u=u+a.*(b.*v-u);                 % update u
  LFP(t,1)=sum(v(1:Ne,1));         % sum of voltages of excitatory per ms 
  LFP(t,2)=sum(v((Ne+1):n,1));     % sum of voltages of inhibitory per ms
  LFP(t,3) = sum(v(:));            % sum of all voltages for each ms
end

%PLOT FIRINGS
figure(1);
plot(firings(:,1),firings(:,2),'.'); 
title('Spike raster - RS,FS STDP in a random network');
%PLOT LFP
figure(2);
plot(LFP);
legend('all','inhibitory','excitatory');
title('LFP - RS,FS STDP in a random network');