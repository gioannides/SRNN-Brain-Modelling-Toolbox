function [w] = induceChaos(w)

waiting = waitbar(0,'Please wait...');
s = clock;

for j = 1:1:w.chaos_period
    %% EULER INTEGRATE
    w.I = w.I_postsynaptic + w.E*w.z  + w.BIAS;  %postsynaptic current
    for pp=1:(1/w.dt)
        w.v = w.v + w.dt*(( w.gain_izhikevich.*(w.v-w.vr).*(w.v-w.vt) - w.u + w.I))./w.C ; % v(t) = v(t-1)+dt*v'(t-1)
    end
    w.u = w.u + w.dt*(w.a.*(w.b.*(w.v_-w.vr)-w.u)); %same with u, the v_ term makes it so that the integration of u uses v(t-1), instead of the updated v(t)
    %% 
    w.index = find(w.v>=w.vpeak);
    if ~isempty(w.index)
        w.JD = sum(w.OMEGA(:,w.index),2); %compute the increase in current due to spiking  
        %tspike(ns+1:ns+length(index),:) = [index,0*index+dt*i];
        
        %ns = ns + length(index); 
    end

    %synapse for single exponential 
    if w.tr == 0 
        w.I_postsynaptic = w.I_postsynaptic*exp(-w.dt/w.td)+  w.JD*(~isempty(w.index))/(w.td);
        w.r = w.r *exp(-w.dt/w.td) + (w.v>=w.vpeak)/w.td;
    else

        %synapse for double exponential
        w.I_postsynaptic = w.I_postsynaptic*exp(-w.dt/w.tr) + w.h*w.dt;
        w.h = w.h*exp(-w.dt/w.td) + w.JD*(~isempty(w.index))/(w.tr*w.td);  %Integrate the current

        w.r = w.r*exp(-w.dt/w.tr) + w.hr*w.dt; 
        w.hr = w.hr*exp(-w.dt/w.td) + (w.v>=w.vpeak)/(w.tr*w.td);
    end


    w.z = w.BPhi'*w.r; %approximated 

    %% Store, and plot.  
    w.u = w.u + w.d.*(w.v>=w.vpeak);  %implements set u to u+d if v>vpeak, component by component. 
    w.v = w.v+(w.vreset-w.v).*(w.v>=w.vpeak); %implements v = c if v>vpeak add 0 if false, add c-v if true, v+c-v = c
    w.v_ = w.v;  % sets v(t-1) = v for the next itteration of loop
    %REC(i,:) = [v(1:5)',u(1:5)'];  
    %current(j,:) = z';
    %RECB(i,:)=BPhi(1:5);
    
    %%%begin estimate remaining time
     if j ==1
      is = etime(clock,s);
      esttime = is * w.chaos_period;
     end
     waiting = waitbar(j/w.chaos_period,waiting,['Network approaching chaotic attractor. Please wait...']);
     %%%end estimate remaining time
    
end

delete(waiting);

end