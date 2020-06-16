%% Force Method with Izhikevich Network 
clear all
close all
clc 
warning('off');

x0  = load('CzChannel80Epochs.mat');

x = x0.dataset{1,1}(:,:,1);
data = double(x);
Fs = 128;

% Uncomment to bandlimit EEG to certain frequency range
%data = bandpass(double(data1),[30 100],Fs);

data = interp(data(1:end),50); % increase data samples by 50 times

dt = 0.1526531; %Integration time step in ms 
n_timesteps = length(data); %Time steps

N =2000; %Number of neurons
p = 0.1; %rewiring probability
rng('default');
rng(2);

%% Target signal  COMMENT OUT TEACHER YOU DONT WANT, COMMENT IN TEACHER YOU WANT. 
target_signal = data;

flag = 0;

task = createNet(N,p,target_signal,flag,0,0,0,dt);
target = createNet(N,p,target_signal,flag,0,0,0,dt);
scaler = 1;
prev_err = 0;

epochs = 2;

%% SIMULATION

%task = induceChaos(task);
%target = induceChaos(target);

%waiting = waitbar(0,'Please wait...');
%s = clock;
step=2;
%i=1;
tic
%ilast = i ;

lambda = 1; %forget factor

%error_plot = [];
%options = odeset('Events',@events);
errors_train = zeros(ceil(n_timesteps),1);
simulation_end = length(target_signal(1:round(3*(end/4))));
stimulus = pinknoise(length(target_signal))';%10*(1:1:length(target_signal));
errors_validation = zeros(ceil(n_timesteps),1);

ADJ=zeros(n_timesteps,N);
count=0;
count2=0;

plasticity=20;   % Plastic changes happen after the first 20 ms of runtime
% St maximum and minimum synaptic weights
max_weight = 100;
min_weight = -100;

record_spikes = zeros(N,task.n_timesteps);

waiting = waitbar(0,'Please wait...');
s = clock;

%Record number of synapse updates
weight_updates_STDP = 0;
tracker = zeros(task.chaos_period,1);

% Carry out STDP
for i = 1:1:task.chaos_period
     %EULER integrate TASK net
   task.I = task.I_postsynaptic + task.E*task.z  + task.BIAS + target_signal(:,i) + stimulus(:,i);  %postsynaptic current
   
    task.v = task.v_ + task.dt*((task.gain_izhikevich.*(task.v-task.vr).*(task.v-task.vt) - task.u + task.I))./task.C ; % v(t) = v(t-1)+dt*v'(t-1)
    
    task.u = task.u + task.dt*(task.a.*(task.b.*(task.v_-task.vr)-task.u));
    
    task.index = find(task.v>=task.vpeak);
    
    for jj=1:numel(task.index)
          ADJ(i,task.index(jj))=1;          % RECORD WHICH NEURONS FIRE AT EACH MS
    end

        if i>plasticity        %t>plasticity period for STDP on, t>runtime for STDP off
          for j=1:N
            count2=count2+1;
            if ADJ(i,j)==1   %neuron j fired at time t
              count=count+1;
              for xx=1:20     %range of ms over which synaptic strengthening and weakening occur                           
                  for k=1:N                  
                    if ADJ(i-xx,k)==1  % find out which neurons fired in the previous plasticity ms              
                      % synaptic weight of pair undergoes larger increment if dT is smaller and negative
                      task.OMEGA(j,k)=task.OMEGA(j,k).*(1+(0.9*exp(xx/20))); 
                      % synaptic weight of pair undergoes larger decrement is dT is smaller and positive
                      task.OMEGA(k,j)=task.OMEGA(k,j).*(1+(-0.925*exp(-xx/20))); 
                      weight_updates_STDP = weight_updates_STDP + 1;
                      % set a maximum value for synaptic weights
                      if task.OMEGA(j,k)>max_weight    
                          task.OMEGA(j,k)=max_weight;
                          weight_updates_STDP = weight_updates_STDP -1;
                      end
                      if task.OMEGA(j,k)<min_weight
                          task.OMEGA(j,k)=min_weight;
                          weight_updates_STDP = weight_updates_STDP - 1;
                      end
                      if task.OMEGA(k,j)>max_weight
                          task.OMEGA(k,j)=max_weight;
                          weight_updates_STDP = weight_updates_STDP - 1;
                      end
                      if task.OMEGA(k,j)<min_weight
                          task.OMEGA(k,j)=min_weight;
                          weight_updates_STDP = weight_updates_STDP -1;
                      end
                    end
                  end
              end
            end
          end
        end
    tracker(i,:) = weight_updates_STDP;
    weight_updates_STDP = 0;
    if mod(i,5000)==0
        plot(tracker(1:i,:));
        xlabel('Timestep');
        ylabel('Number of synpase updates');
    end
        
    if ~isempty(task.index)
        task.JD = sum(task.OMEGA(:,task.index),2); %compute the increase in current due to spiking  
        task.tspike(task.ns+1:task.ns+length(task.index),:) = [task.index,0*task.index+i];
        
        task.ns = task.ns + length(task.index); 
    end

    %synapse for single exponential 
    if task.tr == 0 
        task.I_postsynaptic = task.I_postsynaptic*exp(-task.dt/task.td)+  task.JD*(~isempty(task.index))/(task.td);
        task.r = task.r *exp(-task.dt/task.td) + (task.v>=task.vpeak)/task.td;
    else
        %synapse for double exponential
        task.I_postsynaptic = task.I_postsynaptic*exp(-task.dt/task.tr) + task.h*task.dt;
        task.h = task.h*exp(-task.dt/task.td) + task.JD*(~isempty(task.index))/(task.tr*task.td);  %Integrate the current

        task.r = task.r*exp(-task.dt/task.tr) + task.hr*task.dt; 
        task.hr = task.hr*exp(-task.dt/task.td) + (task.v>=task.vpeak)/(task.tr*task.td);
    end
    
    task.u = task.u + task.d.*(task.v>=task.vpeak);  %implements set u to u+d if v>vpeak, component by component. 
    task.v = task.v+(task.vreset-task.v).*(task.v>=task.vpeak); %implements v = c if v>vpeak add 0 if false, add c-v if true, v+c-v = c
    task.v_ = task.v;  % sets v(t-1) = v for the next itteration of loop
    
    
    waiting = waitbar(i/task.chaos_period,waiting,['Please wait...']);
end

delete(waiting);

waiting = waitbar(0,'Please wait...');
s = clock;

for i = 1:1:target.chaos_period
     %EULER integrate TASK net
   target.I = target.I_postsynaptic + target.E*target.z  + target.BIAS + target_signal(:,i) + stimulus(:,i);  %postsynaptic current
   
    target.v = target.v_ + target.dt*((target.gain_izhikevich.*(target.v-target.vr).*(target.v-target.vt) - target.u + target.I))./target.C ; % v(t) = v(t-1)+dt*v'(t-1)
    
    target.u = target.u + target.dt*(target.a.*(target.b.*(target.v_-target.vr)-target.u));
    
    target.index = find(target.v>=target.vpeak);
    
    for jj=1:numel(target.index)
          ADJ(i,target.index(jj))=1;          % RECORD WHICH NEURONS FIRE AT EACH MS
    end

        if i>plasticity        %t>plasticity period for STDP on, t>runtime for STDP off
          for j=1:N
            count2=count2+1;
            if ADJ(i,j)==1   %neuron j fired at time t
              count=count+1;
              for xx=1:20     %range of ms over which synaptic strengthening and weakening occur                           
                  for k=1:N                  
                    if ADJ(i-xx,k)==1  % find out which neurons fired in the previous plasticity ms              
                      % synaptic weight of pair undergoes larger increment if dT is smaller and negative
                      target.OMEGA(j,k)=target.OMEGA(j,k).*(1+(0.9*exp(xx/20))); 
                      % synaptic weight of pair undergoes larger decrement is dT is smaller and positive
                      target.OMEGA(k,j)=target.OMEGA(k,j).*(1+(-0.925*exp(-xx/20))); 
                      weight_updates_STDP = weight_updates_STDP + 1;
                      % set a maximum value for synaptic weights
                      if target.OMEGA(j,k)>max_weight    
                          target.OMEGA(j,k)=max_weight;
                          weight_updates_STDP = weight_updates_STDP -1;
                      end
                      if target.OMEGA(j,k)<min_weight
                          target.OMEGA(j,k)=min_weight;
                          weight_updates_STDP = weight_updates_STDP - 1;
                      end
                      if target.OMEGA(k,j)>max_weight
                          target.OMEGA(k,j)=max_weight;
                          weight_updates_STDP = weight_updates_STDP - 1;
                      end
                      if target.OMEGA(k,j)<min_weight
                          target.OMEGA(k,j)=min_weight;
                          weight_updates_STDP = weight_updates_STDP -1;
                      end
                    end
                  end
              end
            end
          end
        end
    %tracker(i,:) = weight_updates_STDP;
    %weight_updates_STDP = 0;
    %if mod(i,1)==0
    %     figure(10);
    %    plot(tracker(1:i,:));
    %    xlabel('Timestep');
    %    ylabel('Number of synpase updates');
    %end
        
    if ~isempty(target.index)
        target.JD = sum(target.OMEGA(:,target.index),2); %compute the increase in current due to spiking  
        target.tspike(target.ns+1:target.ns+length(target.index),:) = [target.index,0*target.index+i];
        
        target.ns = target.ns + length(target.index); 
    end

    %synapse for single exponential 
    if target.tr == 0 
        target.I_postsynaptic = target.I_postsynaptic*exp(-target.dt/target.td)+  target.JD*(~isempty(target.index))/(target.td);
        target.r = target.r *exp(-target.dt/target.td) + (target.v>=target.vpeak)/target.td;
    else
        %synapse for double exponential
        target.I_postsynaptic = target.I_postsynaptic*exp(-target.dt/target.tr) + target.h*target.dt;
        target.h = target.h*exp(-target.dt/target.td) + target.JD*(~isempty(target.index))/(target.tr*target.td);  %Integrate the current

        target.r = target.r*exp(-target.dt/target.tr) + target.hr*target.dt; 
        target.hr = target.hr*exp(-target.dt/target.td) + (target.v>=target.vpeak)/(target.tr*target.td);
    end
    
    target.u = target.u + target.d.*(target.v>=target.vpeak);  %implements set u to u+d if v>vpeak, component by component. 
    target.v = target.v+(target.vreset-target.v).*(target.v>=target.vpeak); %implements v = c if v>vpeak add 0 if false, add c-v if true, v+c-v = c
    target.v_ = target.v;  % sets v(t-1) = v for the next itteration of loop
    
    
    waiting = waitbar(i/target.chaos_period,waiting,['Please wait...']);
end

delete(waiting);
flag=1;


store_mats = zeros(N, ceil(length(target_signal)/2));
LFP=zeros(task.n_timesteps,3);

%community_1 = zeros(N*task.n_timesteps,2);
%community_2 = zeros(N*task.n_timesteps,2);
%community_3 = zeros(N*task.n_timesteps,2);
%community_4 = zeros(N*task.n_timesteps,2);
r = randi([0 100], N, 1);
tic
for num=1:1:epochs
    fprintf('EPOCH: %.1f \n', num);
    %if num==epochs
        waiting = waitbar(0,'Please wait...');
        s = clock;
        close all;
        task = createNet(N,p,target_signal,flag,task.BPhi, task.E, task.OMEGA,dt);
        %target = createNet(N,p,target_signal,flag,target.BPhi, target.E, target.OMEGA,dt);
        %task = induceChaos(task);
        %target = induceChaos(target);
        %break;
    %end
    flag = 1;
    error_ratio = [];
    
    for i = 1:1:(task.n_timesteps) 
        
    
    if mod(i, step)==1 && num==epochs
        task.BPhi = store_mats(:,floor(i/step)+1);
        %task.BPhi = mean(store_mats,2);
    end
    
    if mod(i, step)==0 && num<epochs
        store_mats(:,i/step) = task.BPhi;
        %task = createNet(N,p,target_signal,flag,task.BPhi, task.E, task.OMEGA,dt);
        %target = createNet(N,p,target_signal,flag,target.BPhi, target.E, target.OMEGA,dt);
    end
    %% EULER INTEGRATE target
    target.I = target.I_postsynaptic + target.E*task.z  + target.BIAS + target_signal(:,i) + stimulus(:,i);  %postsynaptic current
    %[t1,v] = ode23s(@(t1,v) odevfcn(t1,v,u,I,C,vr,vt,gain_izhikevich),[i,i+1],v);
    %v = v(end,:).';
    %[t2,u] = ode23s(@(t2,u) odeufcn(t2,u,a,b,v_,vr),[i,i+1],u);
    %u = u(end,:).';
    %if num ~= epochs
        %target.v_ = target.v;
        %for pp=1:1:(1/target.dt)
        target.v = target.v_ + target.dt*((target.gain_izhikevich.*(target.v-target.vr).*(target.v-target.vt) - target.u + target.I))./target.C ; % v(t) = v(t-1)+dt*v'(t-1)
        %target.v2 = target.v + 0.5*target.dt*target.v;
        %target.v3 = target.v + 0.5*target.dt*target.v2;
        %target.v4 = target.v + target.dt*target.v3;
        %end

        %target.v = target.v_ + (target.dt*(target.v+2*target.v2+2*target.v3+target.v4))/6;
        %for pp=1:1:(1/target.dt)
        %target.u_ = target.u;    
        target.u = target.u + target.dt*(target.a.*(target.b.*(target.v_-target.vr)-target.u)); %same with u, the v_ term makes it so that the integration of u uses v(t-1), instead of the updated v(t)
        %target.u2 = target.u + 0.5*target.dt*target.u;
        %target.u3 = target.u + 0.5*target.dt*target.u2;
        %target.u4 = target.u + target.dt*target.u3;
        %end
        %target.u = target.u_ + (target.dt*(target.u2+2*target.u2+2*target.u3+target.u4))/6;
        target.index = find(target.v>=target.vpeak);
        %FR(i) = sum(target.index>0);
        if ~isempty(target.index)
            target.JD = sum(target.OMEGA(:,target.index),2); %compute the increase in current due to spiking  
            target.tspike(target.ns+1:target.ns+length(target.index),:) = [target.index,0*target.index+target.dt*i];

            target.ns = target.ns + length(target.index); 
        end

        %synapse for single exponential 
        if target.tr == 0 
            target.I_postsynaptic = target.I_postsynaptic*exp(-target.dt/target.td)+  target.JD*(~isempty(target.index))/(target.td);
            target.r = target.r *exp(-target.dt/target.td) + (target.v>=target.vpeak)/target.td;
        else
            %synapse for double exponential
            target.I_postsynaptic = target.I_postsynaptic*exp(-target.dt/target.tr) + target.h*target.dt;
            target.h = target.h*exp(-target.dt/target.td) + target.JD*(~isempty(target.index))/(target.tr*target.td);  %Integrate the current

            target.r = target.r*exp(-target.dt/target.tr) + target.hr*target.dt; 
            target.hr = target.hr*exp(-target.dt/target.td) + (target.v>=target.vpeak)/(target.tr*target.td);
        end
         target.z = target.BPhi'*target.r; %approximated
    %end
     
     
    %%
    %EULER integrate TASK net
    task.I = task.I_postsynaptic + task.E*task.z  + task.BIAS + stimulus(:,i);  %postsynaptic current
    %[t1,v] = ode45(@(t1,v) odevfcn(t1,v,u,I,C,vr,vt,gain_izhikevich),[i:1/30000:i+1],v);
    %v = v(end,:).';
    %[t2,u] = ode45(@(t2,u) odeufcn(t2,u,a,b,v_,vr),[i:1/30000:i+1],u);
    %u = u(end,:).';
    
    %for pp=1:1:(1/task.dt)
    %task.v_ = task.v;
    task.v = task.v_ + task.dt*((task.gain_izhikevich.*(task.v-task.vr).*(task.v-task.vt) - task.u + task.I))./task.C ; % v(t) = v(t-1)+dt*v'(t-1)
    %task.v2 = task.v + 0.5*task.dt*task.v;
    %task.v3 = task.v + 0.5*task.dt*task.v2;
    %task.v4 = task.v + task.dt*task.v3;
    
    %task.v = task.v_ + (task.dt*(task.v+2*task.v2+2*task.v3+task.v4))/6;
    %end
    
    %for pp=1:1:(1/target.dt)
    %task.u_ = task.u;    
    task.u = task.u + task.dt*(task.a.*(task.b.*(task.v_-task.vr)-task.u)); %same with u, the v_ term makes it so that the integration of u uses v(t-1), instead of the updated v(t)
    %task.u2 = task.u + 0.5*task.dt*task.u;
    %task.u3 = task.u + 0.5*task.dt*task.u2;
    %task.u4 = task.u + task.dt*task.u3;
    
    %task.u = task.u_ + (task.dt*(task.u2+2*task.u2+2*task.u3+task.u4))/6;
    %end
    
     
    task.index = find(task.v>=task.vpeak);
    %if num==epochs
    %    if i == 1
    %        qq = 1;
    %    else
    %        qq = (i-1)*N+(i-1);
    %    end
    %    for b=1:length(task.index)
    %        if task.index(b)<=750
    %            community_1(qq,1) = task.index(b);
    %            community_1(qq,2) = i;
    %        end
    %        if task.index(b)>750 && task.index(b)<=1500
    %            community_2(qq,1) = task.index(b);
    %            community_2(qq,2) = i;
    %        end
    %        if task.index(b)>1500 && task.index(b)<=2250
    %            community_3(qq,1) = task.index(b);
    %            community_3(qq,2) = i;
    %        end
    %        if task.index(b)>2250
    %            community_4(qq,1) = task.index(b);
    %            community_4(qq,2) = i;
    %        end
    %        qq = qq+1;
    %    end
    %end
    if ~isempty(task.index)
        task.JD = sum(task.OMEGA(:,task.index),2); %compute the increase in current due to spiking  
        if num==epochs
            task.tspike(task.ns+1:task.ns+length(task.index),:) = [task.index,0*task.index+i];

            task.ns = task.ns + length(task.index); 
        end
    end

    %synapse for single exponential 
    if task.tr == 0 
        task.I_postsynaptic = task.I_postsynaptic*exp(-task.dt/task.td)+  task.JD*(~isempty(task.index))/(task.td);
        task.r = task.r *exp(-task.dt/task.td) + (task.v>=task.vpeak)/task.td;
    else
        %synapse for double exponential
        task.I_postsynaptic = task.I_postsynaptic*exp(-task.dt/task.tr) + task.h*task.dt;
        task.h = task.h*exp(-task.dt/task.td) + task.JD*(~isempty(task.index))/(task.tr*task.td);  %Integrate the current

        task.r = task.r*exp(-task.dt/task.tr) + task.hr*task.dt; 
        task.hr = task.hr*exp(-task.dt/task.td) + (task.v>=task.vpeak)/(task.tr*task.td);
    end
     task.z = task.BPhi'*task.r; %approximant
    %%   
     if num ~= epochs
        task.err = task.z - rand(min(size(target.r)))*target.z - target_signal(:,i) - stimulus(:,i); %error
        errors_train(i) = task.err;
     end
      
     if num==epochs
        task.err = task.z - rand(min(size(target.r)))*target.z - target_signal(:,i) - stimulus(:,i); %error
        errors_validation(i) = task.err;
     end
    
     %% RLS 
     if mod(i,step)==0 && num<epochs
        %if i > imin 
             %if i < simulation_end 
               task.cd = task.Pinv*task.r; %correlation matrix, Pinv, times the synapse for double exponential, r
               task.BPhi = task.BPhi - (task.cd*task.err'); %BPhi is the weight matrix, cd is the 
               task.Pinv = (1/lambda)*task.Pinv -(1/lambda)*((task.cd)*(task.cd'))/( 1 + (task.r')*(task.cd));
               %disp(task.err/prev_err);
               task.track_WEIGHTS_std(i/step) = std(task.BPhi);
               task.track_WEIGHTS_avg(i/step) = mean(task.BPhi);
               error_ratio = [error_ratio abs(task.err/prev_err)];
             %else
             %    break;
             %end
       prev_err = task.err;
       % end 
     end
     %task.z = task.BPhi'*tanh(task.r); %approximant
     %tas.err_plus = task.z - rand(min(size(target.r)))*target.z - target_signal(:,i); %error 

    if num == epochs
        record_spikes(task.index,i) = 1;
    end


    %% Store, and plot TARGET.  
    target.u = target.u + target.d.*(target.v>=target.vpeak);  %implements set u to u+d if v>vpeak, component by component. 
    target.v = target.v+(target.vreset-target.v).*(target.v>=target.vpeak); %implements v = c if v>vpeak add 0 if false, add c-v if true, v+c-v = c
    target.v_ = target.v;  % sets v(t-1) = v for the next itteration of loop
    target.REC(i,:) = [target.v(1:5)',target.u(1:5)'];  
    target.current(i,:) = target.z';
    target.RECB(i,:)=target.BPhi(1:5);
    
    
    %% Store, and plot TASK.  
    task.u = task.u + task.d.*(task.v>=task.vpeak);  %implements set u to u+d if v>vpeak, component by component. 
    task.v = task.v+(task.vreset-task.v).*(task.v>=task.vpeak); %implements v = c if v>vpeak add 0 if false, add c-v if true, v+c-v = c
    task.v_ = task.v;  % sets v(t-1) = v for the next itteration of loop
    LFP(i,1)=sum(task.v(1:end));         % sum of voltages of excitatory per ms 
    task.REC(i,:) = [task.v(1:5)',task.u(1:5)'];  
    task.current(i,:) = task.z';
    task.RECB(i,:)=task.BPhi(1:5);


   if mod(i,9600)==0 
        drawnow
        figure(1)
        plot(target_signal(1:i),'k','LineWidth',2), hold on
        plot(task.current(1:i),'g--','LineWidth',2), hold off
        xlabel('Timestep')
        ylabel('Amplitude, $\mu$V','Interpreter','latex')
        legend('Target EEG Signal', 'Modelled EEG');
        set(gca,'FontSize',15);
        grid on;
        if num==epochs
            figure(2)
            plot(errors_validation(1:i), 'LineWidth',2);
            xlabel('Timestep');
            ylabel('Error');
            grid on;
        end
        figure(7);
        %plot((1:1:i)*dt/1000,task.RECB(1:1:i,:))
        %if num==epochs
        %    figure(3)
        %    for b=1:task.ns
        %        if task.tspike(b,1)<=750
        %            plot(task.tspike(b,2),task.tspike(b,1), 'r.', 'MarkerSize',0.01);hold on;
        %        end
        %        if task.tspike(b,1)<=1500 && task.tspike(b,1)>750
        %            plot(task.tspike(b,2),task.tspike(b,1), 'g.', 'MarkerSize',0.01);hold on;
        %        end
        %        if task.tspike(b,1)<=2250 && task.tspike(b,1)>1500
        %            plot(task.tspike(b,2),task.tspike(b,1), 'b.', 'MarkerSize',0.01);hold on;
        %        end
        %        if task.tspike(b,1)<=3000 && task.tspike(b,1)>2250
        %            plot(task.tspike(b,2),task.tspike(b,1), 'k.', 'MarkerSize',0.01);hold on;
        %        end
        %    end
        plot(task.tspike(1:task.ns,2),task.tspike(1:task.ns,1),'k.', 'MarkerSize',0.01);
        %if mod(i,task.n_timesteps)==0 && num==epochs
        %community_1 = community_1(any(community_1,2),:);
        %community_2 = community_2(any(community_2,2),:);
        %community_3 = community_3(any(community_3,2),:);
        %community_4 = community_4(any(community_4,2),:);
        %figure(2);
        %xlim([0 task.n_timesteps]);
        %ylim([0,task.N]);
        %plot(community_1(1:end,2),community_1(1:end,1),'r.', 'LineWidth',0.01); hold on;
        %plot(community_2(1:end,2),community_2(1:end,1),'g.', 'LineWidth',0.01)
        %plot(community_3(1:end,2),community_3(1:end,1),'b.', 'LineWidth',0.01)
        %plot(community_4(1:end,2),community_4(1:end,1),'k.', 'LineWidth',0.01); hold off;
        xlabel('Timestep');
        ylabel('Neuron index');
        %end
        %end
        set(gca,'FontSize',15);hold off;
        figure(4)
        plot(task.track_WEIGHTS_std(1:end),'LineWidth',2);
        xlabel('Timestep');
        ylabel('$\Delta$W (s.d. in Weight updates)','interpreter','latex');
        figure(5)
        plot(task.track_WEIGHTS_avg(1:end),'LineWidth',2);
        xlabel('Timestep');
        ylabel('$\Delta$W (Mean Weight update)','interpreter','latex');
        %PLOT LFP
        figure(6);
        plot(LFP(1:i), 'Color', [0.8500, 0.3250, 0.0980], 'Linewidth', 2);
        legend('Collective Voltage of Neurons');
        ylabel('mV');
        xlabel('Timestep');
        set(gca,'FontSize',15);
        grid on;
    end   
    
     %%%begin estimate remaining time
     if i ==1
      is = etime(clock,s);
      esttime = is * task.n_timesteps;
     end
     
     waiting = waitbar(i/task.n_timesteps,waiting,['Please wait...']);
     %%%end estimate remaining time
     %task.D = abs(target_signal(:,gg:1:i)'-task.current(gg:1:i,:)).^2;
     %task.error_metric = sqrt(sum(task.D(:))/numel(target_signal(:,gg:1:i)'));
     %disp(task.error_metric);
     %task.error_plot(i) = task.error_metric; 
     
    end
    delete(waiting);
end
toc
delete(waiting);
summation = 0;
for i=1:1:length(errors_validation)
    summation = summation + errors_validation(i)^2;
end
rmse = sqrt(mean(summation));
fprintf('RMSE: %.1f \n', rmse);

figure();
histfit(errors_train);
pd = fitdist(errors_train,'Normal');

Fs = 128;
approx = decimate(task.current(1:end),50);
NFFT = length(approx);
Y = fft(approx,NFFT);
F = ((0:1/NFFT:1-1/NFFT)*Fs).';
magnitudeY = abs(Y);        % Magnitude of the FFT

figure()
dB_mag=mag2db(magnitudeY);
plot(F(1:NFFT/2),dB_mag(1:NFFT/2),'r');
title('Magnitude response of Approximated EEG signal');
ylabel('Magnitude(dB)');
xlabel('Frequency in Hz')
set(gca,'FontSize',15);

figure()
Fs = 128;
NFFT = length(x);
C = fft(x,NFFT);
F = ((0:1/NFFT:1-1/NFFT)*Fs).';
magnitudeC = abs(C);        % Magnitude of the FFT

dB_mag=mag2db(abs(magnitudeC));
plot(F(1:NFFT/2),dB_mag(1:NFFT/2),'b');
title('Magnitude response of Actual EEG signal');
ylabel('Magnitude(dB)');
xlabel('Frequency in Hz')
set(gca,'FontSize',15);

figure();
Fs = 128;
N = length(approx);
q = fft(approx,N);
ff = 0:Fs/N:Fs-Fs/N;
ff = ff(1:floor(N/2)+1);
q = q(1:floor(N/2)+1);
plot(ff,log(abs(q)/N), 'Linewidth', 2); hold on;
Fs = 128;
N = length(x);
q = fft(x,N);
ff = 0:Fs/N:Fs-Fs/N;
ff = ff(1:floor(N/2)+1);
q = q(1:floor(N/2)+1);
plot(ff,log(abs(q)/N), 'Linewidth', 2); hold off;

xlabel('Frequency (Hz)');
ylabel('Log Power Spectrum Density');
legend({'Modelled EEG','Target EEG'})
set(gca,'FontSize',15);
grid on;

tspike = task.tspike(task.tspike(:,2)~=0,:); 
M = tspike(tspike(:,2)>2992); 
AverageFiringRate = length(M)/(task.N*(2992));
%figure()
%Z = eig(task.OMEGA+task.E*task.BPhi'); %eigenvalues after learning (THE CONNECTION MATRIX)
%Z2 = eig(task.OMEGA); %eigenvalues before learning 
%
%plot(Z2,'r.'), hold on 
%plot(Z,'k.') 
%legend('Pre-Learning','Post-Learning')
%xlabel('Real \lambda')
%ylabel('Imaginary \lambda')
%beep on;