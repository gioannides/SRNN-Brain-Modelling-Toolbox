function [net] = createNet_version2(N,p,target_signal,flag,BPhi,E,OMEGA,dt)

%s.t = 15; %Total time in ms
%s.dt = s.t/length(target_signal); %Integration time step in ms 
s.dt = dt;
s.n_timesteps = length(target_signal); %Time steps
s.N = N;  %Number of neurons 
s.p=p;

%%
rng(2);
N_ = s.N; % total number of nodes in all modules
K = 6; % K nearest neighbors on either side
beta = s.p; %rewiring probability beta.
modules = {};
M = 4; % number of modules
n_over_m = N_/M;
no_edges = n_over_m*K;
A = zeros(N_*K,3);

for i=0:M-1
    modules{i+1} = WattsStrogatz(n_over_m,K,beta);
    col_1_2 = table2array(modules{i+1}.Edges(:,1)) + n_over_m*(i);
    col_3 = table2array(modules{i+1}.Edges(:,2));
    A(i*no_edges+1:(i+1)*no_edges, 1:2) = col_1_2;
    A(i*no_edges+1:(i+1)*no_edges, 3) = col_3;
end

s.S = zeros(N_,1);

%disp(size(s.S));

for i=1:(N_*K)    
    idx_1 = A(i,1);
    idx_2 = A(i,2);
    s.S(idx_1,idx_2) = A(i,3);
    s.S(idx_2,idx_1) = A(i,3);
    if idx_1 == idx_2
        s.S(idx_1,idx_2) = 0;
    end
end

p = 0.1;
s.S = RewireNetwork(s.S,M,p);

%%

%K = 5; % K nearest neighbors on either side
%beta = s.p; %rewiring probability beta.

%watts = WattsStrogatz(N,K,beta);
%plot(watts);
%s.S = zeros(N);

s.Pinv = eye(s.N)*25; %initial correlation matrix, coefficient is the regularization constant as well 

%A = table2array(watts.Edges); % the weights of each neuron

%for i=1:(N*K)    
%    idx_1 = A(i,1);
%    idx_2 = A(i,2);
%    s.S(idx_1,idx_2) = A(i,3);
%end

%%


%s.C = 250;  %capacitance
%s.C = 100;  %RS
s.C = zeros(N,1);
s.C(1:3*(N/4)) = 100; % IB
s.C(3*(N/4)+1:end) = 100; % RS
%s.C(3*(N/4)+1:end) = 150; % RS

%s.b = -2;  %constant (1/R), RS
s.b = zeros(N,1);
s.b(1:3*(N/4)) = -2; % IB
s.b(3*(N/4)+1:end) = -2; % RS
%s.b(3*(N/4)+1:end) = 5; % RS

%s.gain_izhikevich = 2.5;  %k parameter for Izhikevich, gain on v
s.gain_izhikevich = zeros(N,1);
s.gain_izhikevich(1:3*(N/4)) = 0.7;  %k parameter for Izhikevich, gain on v, IB, RS
s.gain_izhikevich(3*(N/4)+1:end) = 0.7;
%s.gain_izhikevich(3*(N/4)+1:end) = 1.2;

%s.vpeak = 30;  % spike cutoff value
%s.vpeak = 35;  % RS spike cutoff value
s.vpeak = zeros(N,1);
s.vpeak(1:3*(N/4)) = 35; %  IB
s.vpeak(3*(N/4)+1:end) = 35; % RS
%s.vpeak(3*(N/4)+1:end) = 50; % RS

%s.vreset = -65; % reset voltage value
%s.vreset = -50; % RS reset voltage value
s.vreset = zeros(N,1);
s.vreset(1:3*(N/4)) = -50; %  IB
s.vreset(3*(N/4)+1:end) = -50; % RS
%s.vreset(3*(N/4)+1:end) = -56; % RS

%s.vr = -60;   %resting membrane potential RS
s.vr = zeros(N,1);
s.vr(1:3*(N/4)) = -60; % IB
s.vr(3*(N/4)+1:end) = -60; % RS
%s.vr(3*(N/4)+1:end) = -75; % RS

%s.vt = s.vr + 40 -(s.b/s.gain_izhikevich); %instanteneous threshold potential
%s.vt = -40; %RS
s.vt = zeros(N,1);
s.vt(1:3*(N/4)) = -40; % IB 
s.vt(3*(N/4)+1:end) = -40; % RS
%s.vt(3*(N/4)+1:end) = -45; % RS

s.u = zeros(N,1);  %recovery variable
s.a = zeros(N,1);
s.a(1:3*(N/4)) = 0.03; % IB
s.a(3*(N/4)+1:end) = 0.03; %recovery time constant RS
%s.a(3*(N/4)+1:end) = 0.03; %recovery time constant RS
%s.a = 0.03; %RS recovery time constant

%s.d = 200; %outward minus inward currents activated during the spike and affecting the after-spike behavior
%s.d = 100; %RS outward minus inward currents activated during the spike and affecting the after-spike behavior
s.d = zeros(N,1);
s.d(1:3*(N/4)) = 100; % IB
s.d(3*(N/4)+1:end) = 100; % RS
%s.d(3*(N/4)+1:end) = 130; % RS

%s.tr = 2;  %synaptic rise time 
s.tr = 8;  %synaptic rise time RS, IB

s.td = 20; %decay time 
s.G =15*10^3; %Gain on the static matrix with 1/sqrt(N) scaling weights.  Note that the units of this have to be in pA. 
           % The greater the G the greater the chaos 

s.CHAOS_DURATION = 5000; %In milliseconds
%s.chaos_period = round(s.CHAOS_DURATION/s.dt); %time before starting RLS, gets the network to chaotic attractor 
s.chaos_period = s.CHAOS_DURATION; %time before starting RLS, gets the network to chaotic attractor 
s.Q =15*10^3; %Gain on the rank-k perturbation modified by RLS.  Note that the units of this have to be in pA   

s.Irh = 0.25*s.gain_izhikevich.*(s.vt-s.vr).^2;
%Storage variables for synapse integration  
s.I_postsynaptic = zeros(N,1); %post synaptic current 
s.h = zeros(N,1); 
s.r = zeros(N,1);
s.hr = zeros(N,1);
s.JD = zeros(N,1);
%-----Initialization---------------------------------------------
s.v = s.vr+(s.vpeak-s.vr).*rand(s.N,1); %initial distribution 
s.v_ = s.v; %These are just used for Euler integration, previous time step storage

if flag==0
    %s.OMEGA = s.G*(randn(s.N,s.N)).*(rand(s.N,s.N)<s.p)/(s.p*sqrt(s.N)); %Static weight matrix.
    s.OMEGA = s.G*s.S; %Static weight matrix
    %s.OMEGA = s.OMEGA>0;
    %s.OMEGA = s.OMEGA*s.G;
else
    s.OMEGA = OMEGA;
end
%set the row average weight to be zero, explicitly.

s.k = min(size(target_signal)); %used to get the dimensionality of the approximant correctly.  Typically will be 1 unless you specify a k-dimensional target function.    

if flag==0
    for i = 1:1:s.N
        QS = find(abs(s.OMEGA(i,:))>0);
        s.OMEGA(i,QS) = s.OMEGA(i,QS) - sum(s.OMEGA(i,QS))/length(QS);
    end
end

if flag == 0
    s.BPhi = zeros(s.N,s.k); %initial decoder.  Best to keep it at 0.
else
    s.BPhi = BPhi;
end

s.ns = 0; %counter of total number of spikes
s.BIAS = zeros(N,1);
%s.BIAS(1:3*(N/4))=350;
s.BIAS(3*(N/4)+1:end) = 50;
s.BIAS(3*(N/4)+1:end) = 50;
%s.BIAS = 1000; %Bias current, note that the Rheobase is around 950 or something.  I forget the exact formula for this but you can test it out by shutting weights and feeding constant currents to neurons 
if flag==0
    s.E = (2*rand(s.N,s.k)-1)*s.Q;  %Weight matrix is OMEGA0 + E*BPhi';
else
    s.E = E;
end
%%
s.track_WEIGHTS_std = zeros(length(target_signal)/2,1);
s.track_WEIGHTS_avg = zeros(length(target_signal)/2,1);
%%

s.z = zeros(s.k,1);  %initial approximated target signal
s.tspike = zeros(N*s.n_timesteps,2);  %storage of spike times 
s.step = 2; %optimize with RLS only every 10 steps 



s.current = zeros(s.n_timesteps,s.k);  %store the approximated target signal
s.RECB = zeros(s.n_timesteps,5); %store the decoders 
s.REC = zeros(s.n_timesteps,10); %Store voltage and adaptation variables for plotting 

s.simulation_end = round(length(target_signal)/s.dt); %end simulation at this time step and keep the RSNN keep spiking what it learned

net = s;

end