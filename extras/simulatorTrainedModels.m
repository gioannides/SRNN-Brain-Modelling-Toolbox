%% Force Method with Izhikevich Network 
clear all
close all
clc 
%m = load('../IZplots/CzChannelEpoch1_FITTED_RS_neuronsmat.mat');
%m = load('../IZplots/Czfullforce.mat');
%m = load('SWN/fullforceIB.mat');
m = load('matlab.mat');
m = m.task;
%x = m.x;
%save('force.mat','x');
rng(1);
%t = m.t; %Total time in ms
dt = m.dt; %Integration time step in ms 
n_timesteps = 500; %Time steps
N = m.N;  %Number of neurons 

%% Izhikevich Parameters
C = zeros(N,1);
C(1:(N/2)) = 100; % IB
C((N/2)+1:end) = 100; % RS
%C = 100;  %capacitance

b = m.b;  %constant (1/R)

gain_izhikevich = m.gain_izhikevich;  %k parameter for Izhikevich, gain on v

vpeak = m.vpeak;  % spike cutoff value

vreset = m.vreset; % reset voltage value

vr = m.vr;   %resting membrane potential

vt = m.vt;

u = m.u;  %recovery variable

a = m.a; %recovery time constant

d = m.d; %outward minus inward currents activated during the spike and affecting the after-spike behavior

tr = m.tr;  %synaptic rise time 
td = m.td; %decay time 
%p = m.p; %small-worldness
G =m.G; %Gain on the static matrix with 1/sqrt(N) scaling weights.  Note that the units of this have to be in pA. 
           % The greater the G the greater the chaos 
           


Q =m.Q; %Gain on the rank-k perturbation modified by RLS.  Note that the units of this have to be in pA 
           
Irh = m.Irh; 

%Storage variables for synapse integration  
I_postsynaptic = zeros(N,1); %post synaptic current 
h =zeros(N,1);  
r = zeros(N,1); 
hr = zeros(N,1); 
JD = zeros(N,1); 

%-----Initialization---------------------------------------------
v = m.vr+(m.vpeak-m.vr).*rand(m.N,1); %initial distribution ; %initial distribution 
v_ = v; %These are just used for Euler integration, previous time step storage;



k = m.k; %used to get the dimensionality of the approximant correctly.  Typically will be 1 unless you specify a k-dimensional target function.  

OMEGA = m.OMEGA; %Static weight matrix.  

z = m.z;  %initial approximated target signal

BPhi = m.BPhi; %initial decoder.  Best to keep it at 0.  
tspike = zeros(5*n_timesteps,2);  %storage of spike times 

ns = m.ns; %counter of total number of spikes
BIAS = m.BIAS; %Bias current, note that the Rheobase is around 950 or something.  I forget the exact formula for this but you can test it out by shutting weights and feeding constant currents to neurons 
E = m.E;  %Weight matrix is OMEGA0 + E*BPhi'; 
%%

Pinv = m.Pinv; %initial correlation matrix, coefficient is the regularization constant as well 
step = m.step; %optimize with RLS only every 10 steps

current = zeros(n_timesteps,k);  %store the approximated target signal
RECB =m.RECB; %store the decoders 
REC = m.REC; %Store voltage and adaptation variables for plotting 
i=1;

%% SIMULATION
CHAOS_DURATION = n_timesteps; %In timesteps
%tspikes = zeros(CHAOS_DURATION,N);
waiting = waitbar(0,'Please wait...');
s = clock;
ns = 0;
for j = 1:1:CHAOS_DURATION
    %% EULER INTEGRATE
    I = I_postsynaptic + E*z  + BIAS;  %postsynaptic current 
    for pp=1:1:(round(1/m.dt))
    v = v + dt*((gain_izhikevich.*(v-vr).*(v-vt) - u + I))./C ; % v(t) = v(t-1)+dt*v'(t-1)
    end
    u = u + dt*(a.*(b.*(v_-vr)-u)); %same with u, the v_ term makes it so that the integration of u uses v(t-1), instead of the updated v(t)
    %% 
    index = find(v>=vpeak);
    %tspikes(j,index) = 1;
    if ~isempty(index)
        JD = sum(OMEGA(:,index),2); %compute the increase in current due to spiking  
        ns = ns + length(index);
    end

    %synapse for single exponential 
    if tr == 0 
        I_postsynaptic = I_postsynaptic*exp(-dt/td)+  JD*(~isempty(index))/(td);
        r = r *exp(-dt/td) + (v>=vpeak)/td;
    else

        %synapse for double exponential
        I_postsynaptic = I_postsynaptic*exp(-dt/tr) + h*dt;
        h = h*exp(-dt/td) + JD*(~isempty(index))/(tr*td);  %Integrate the current

        r = r*exp(-dt/tr) + hr*dt; 
        hr = hr*exp(-dt/td) + (v>=vpeak)/(tr*td);
    end


    z = BPhi'*r; %approximated 

    %% Store, and plot.  
    u = u + d.*(v>=vpeak);  %implements set u to u+d if v>vpeak, component by component. 
    v = v+(vreset-v).*(v>=vpeak); %implements v = c if v>vpeak add 0 if false, add c-v if true, v+c-v = c
    v_ = v;  % sets v(t-1) = v for the next itteration of loop
    REC(i,:) = [v(1:5)',u(1:5)'];  
    current(j,:) = z';
    RECB(i,:)=BPhi(1:5);
    
    %%%begin estimate remaining time
     if j==1
      is = etime(clock,s);
      esttime = is * CHAOS_DURATION;
     end
     waiting = waitbar(j/CHAOS_DURATION,waiting,['Network simulating. Please wait...']);
     %%%end estimate remaining time
     
      if mod(round(j/(dt*10)),5)==0
        drawnow
        figure(1)
        plot(current(1:j),'b','LineWidth',2), hold off
        xlabel('Timestep')
        ylabel('Amplitude, $\mu$V','Interpreter','latex')
        legend('Modelled EEG - IB');
        %ylim([-6 12]);
        %figure(2)
        %plot(m.tspike(1:ns,2),m.tspike(1:ns,1),'k.')
        %ylim([0,N])
    end   
    
end

delete(waiting);

data = decimate(current,1/m.dt); % increase data samples by 10 times
figure();
plot(data, 'r');
xlabel('Time (s)');
ylabel('Amplitude');

Fs = 128;
NFFT = length(current);
Y = fft(current,NFFT);
F = ((0:1/NFFT:1-1/NFFT)*Fs).';
magnitudeY = abs(Y);        % Magnitude of the FFT

figure()
dB_mag=mag2db(magnitudeY);
plot(F(1:NFFT/2),dB_mag(1:NFFT/2),'r');
title('Magnitude response of Approximated EEG signal');
ylabel('Magnitude(dB)');
xlabel('Frequency in Hz')

%figure()
%Fs = 128;
%NFFT = length(m.target_signal(1:round(1*n_timesteps)));
%C = fft(m.target_signal(1:round(1*m.n_timesteps)),NFFT);
%F = ((0:1/NFFT:1-1/NFFT)*Fs).';
%magnitudeC = abs(C);        % Magnitude of the FFT

%dB_mag=mag2db(abs(magnitudeC));
%plot(F(1:NFFT/2),dB_mag(1:NFFT/2),'b');
%title('Magnitude response of Actual EEG signal');
%ylabel('Magnitude(dB)');
%xlabel('Frequency in Hz')

figure();
Fs = 128;
N = length(current);
q = fft(current,N);
ff = 0:Fs/N:Fs-Fs/N;
ff = ff(1:floor(N/2)+1);
q = q(1:floor(N/2)+1);
plot(ff,log(abs(q)/N)); hold off;
Fs = 128;
%N = length(m.target_signal(1:round(1*m.n_timesteps)));
%q = fft(m.target_signal(1:round(1*m.n_timesteps)),N);
%ff = 0:Fs/N:Fs-Fs/N;
%ff = ff(1:floor(N/2)+1);
%q = q(1:floor(N/2)+1);
%plot(ff,log(abs(q)/N)); hold off;

xlabel('Frequency (Hz)');
ylabel('Log Power Spectrum Density');
legend({'Approximated EEG','Real EEG'});

%total_spikes_perstep = zeros(CHAOS_DURATION,1);
%for j = 1:1:CHAOS_DURATION
%    total_spikes_perstep(j) = sum(tspikes(j,m.N-80));
%end
%disp(sum(total_spikes_perstep));

%[R,P] = corrcoef(current,m.target_signal(1:length(current)));
%disp(R(2));