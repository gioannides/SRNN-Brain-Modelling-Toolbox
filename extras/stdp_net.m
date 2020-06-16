rand('seed',1);
%Spiking network with axonal conduction delays and STDP
% Created by Eugene M.Izhikevich.                February 3, 2004
% Modified to allow arbitrary delay distributions.  April 16,2008
M=100;                 % number of synapses per neuron
D=20;                  % maximal conduction delay 
% excitatory neurons   % inhibitory neurons      % total number 
Ne=800;                Ni=200;                   N=Ne+Ni;
a=[0.02*ones(Ne,1);    0.1*ones(Ni,1)];
d=[   8*ones(Ne,1);    2*ones(Ni,1)];
sm=10;                 % maximal synaptic strength

x0  = load('CzChannel80Epochs.mat');
x = x0.dataset{1,1}(:,:,1);
data = double(x);
Fs = 128;

data = interp(data(1:end),50); % increase data samples by 10 times

dt = 0.15;
%Storage variables for synapse integration  
IPSC = zeros(N,1); %post synaptic current 
h = zeros(N,1); 
r = zeros(N,1);
hr = zeros(N,1);
JD = zeros(N,1);

% post=ceil([N*rand(Ne,M);Ne*rand(Ni,M)]); 
% Take special care not to have multiple connections between neurons
delays = cell(N,D);
for i=1:Ne
    p=randperm(N);
    post(i,:)=p(1:M);
    for j=1:M
        delays{i, ceil(D*rand)}(end+1) = j;  % Assign random exc delays
    end;
end;
for i=Ne+1:N
    p=randperm(Ne);
    post(i,:)=p(1:M);
    delays{i,1}=1:M;                    % all inh delays are 1 ms.
end;

s=[6*ones(Ne,M);-5*ones(Ni,M)];         % synaptic weights
sd=zeros(N,M);                          % their derivatives
step = 2;

% Make links at postsynaptic targets to the presynaptic weights
pre = cell(N,1);
aux = cell(N,1);
for i=1:Ne
    for j=1:D
        for k=1:length(delays{i,j})
            pre{post(i, delays{i, j}(k))}(end+1) = N*(delays{i, j}(k)-1)+i;
            aux{post(i, delays{i, j}(k))}(end+1) = N*(D-1-j)+i; % takes into account delay
        end;
    end;
end;
  
tr = 2;
td = 20;

STDP = zeros(N,1001+D);
v = -65*ones(N,1);                      % initial values
u = 0.2.*v;                             % initial values
firings=[-D 0];                         % spike timings

%% Target signal  COMMENT OUT TEACHER YOU DONT WANT, COMMENT IN TEACHER YOU WANT. 
target_signal = data;

%%
k = min(size(target_signal)); %used to get the dimensionality of the approximant correctly.  Typically will be 1 unless you specify a k-dimensional target function.  
z = zeros(k,1);  %initial approximant
BPhi = zeros(N,k); %initial decoder.  Best to keep it at 0.  
tspike = zeros(5*length(target_signal),2);  %If you want to store spike times, 
ns = 0; %count toal number of spikes

%E = (2*rand(N,k)-1)*Q;  %Weight matrix is OMEGA0 + E*BPhi'; 
Pinv = eye(N)*10; %initial correlation matrix, coefficient is the regularization constant as well 
current = zeros(length(target_signal),k);  %store the approximant 
RECB = zeros(length(target_signal),5); %store the decoders 
REC = zeros(length(target_signal),10); %Store voltage and adaptation variables for plotting 
increase = 0;

%for sec=1:1                      % simulation of 1 day
  for i=1:length(target_signal)                          % simulation of 1 sec
    I=zeros(N,1);        
    I(ceil(N*rand))=20+z;                 % random thalamic input 
    fired = find(v>=30);                % indices of fired neurons
    v_ = v;
    v(fired)=-65;  
    u(fired)=u(fired)+d(fired);
    STDP(fired,i+D)=0.1;
    for k=1:length(fired)
      sd(pre{fired(k)})=sd(pre{fired(k)})+STDP(N*i+aux{fired(k)});
    end;
    firings=[firings;i*ones(length(fired),1),fired];
    k=size(firings,1);
    temp = I-I;
    while firings(k,1)>i-D
      del=delays{firings(k,2),i-firings(k,1)+1};
      ind = post(firings(k,2),del);
      I(ind)=I(ind)+s(firings(k,2), del)';
      temp(ind) = temp(ind) + s(firings(k,2), del)';
      sd(firings(k,2),del)=sd(firings(k,2),del)-1.2*STDP(ind,i+D)';
      k=k-1;
    end;
    %synapse for double exponential
    IPSC = IPSC*exp(-dt/tr) + h*dt;
    h = h*exp(-dt/td) + temp*(length(fired)>0)/(tr*td);  %Integrate the current

    r = r*exp(-dt/tr) + hr*dt; 
    hr = hr*exp(-dt/td) + (v_>=30)/(tr*td);
    
    z = BPhi'*r; %approximant 
    err = z - target_signal(:,i); %error 
    %% RLS 
     if mod(i,step)==1 
           cd = Pinv*r;
           BPhi = BPhi - (cd*err');
           Pinv = Pinv -((cd)*(cd'))/( 1 + (r')*(cd));
     end 
    
    v=v+0.5*((0.04*v+5).*v+140-u+I);    % for numerical 
    v=v+0.5*((0.04*v+5).*v+140-u+I);    % stability time 
    u=u+a.*(0.2*v-u);                   % step is 0.5 ms
    STDP(:,i+D+1)=0.95*STDP(:,i+D);     % tau = 20 ms
    REC(i,:) = [v(1:5)',u(1:5)'];  
    current(i,:) = z';
    RECB(i,:)=BPhi(1:5);
    
    if mod(i,9600)==1 
        drawnow
        figure(1)
        plot(target_signal(1:i),'k','LineWidth',2), hold on
        plot(current(1:i),'r--','LineWidth',2), hold off
        xlabel('Timestep')
        ylabel('Amplitude, $\mu$V','Interpreter','latex')
        legend('Target EEG Signal', 'Modelled EEG');
        set(gca,'FontSize',15);
        grid on;
        figure(3)
        plot(firings(:,1),firings(:,2),'k.', 'MarkerSize',0.05);
        xlabel('Timestep');
        ylabel('Neuron index');
        ylim([0,N]);
    end   
  end;
  plot(firings(:,1),firings(:,2),'.');
  axis([0 1000 0 N]); drawnow;
  STDP(:,1:D+1)=STDP(:,1001:1001+D);
  ind = find(firings(:,1) > 1001-D);
  firings=[-D 0;firings(ind,1)-1000,firings(ind,2)];
  s(1:Ne,:)=max(0,min(sm,0.01+s(1:Ne,:)+sd(1:Ne,:)));
  sd=0.9*sd;
