close all;
clc;
clear all;

%Uncomment this to try out different rewiring probabilities
%p_ = 1:1:10;
p_ = 1;
% Define number of timesteps in upsampled EEG signal
time = 19200;
% Define number of excitatory neurons
Ne_ = 1600;
% Define number of inhibitory neurons
Ni_ = 400;
% Set the seed for reproducibility
rng(8);
% Store the errors of each EEG channel in this array to plot a boxplot
boxes = [];
window_size = 50; % window size (it must be greater than or equal to Dmax)
slide_window = 20; % slide window
mkdir('results');
mkdir('./results/Connectivity');
mkdir('./results/MFR');
for p = p_
    % Build the topology of modular SWN
    beta = p/10;
    [CIJ,Ne_per_module] = BuildTopology(beta, Ne_, Ni_);
    % store the connectivity matrix
    clf;
    spy(CIJ,'k.');
    xlabel('');
    saveas( gcf,['./results/Connectivity/CM_' num2str(p/10) '.png'], 'png' );
    % Create the TASK network
    layer = ConnectNetwork(CIJ,Ne_,Ni_);   
    % Create the TARGET network
    layer2 = ConnectNetwork(CIJ,Ne_,Ni_);  
    % Carry out modified-Full-Force
    [firings, firings2, boxes] = Simulate(layer,layer2,time,Ne_,Ni_,Ne_per_module);
    % Estimate MFR of excitatory neurons of TASK network
    MFR = MeanFiringRate(firings,window_size,slide_window, Ne_per_module);    
    % Estimate Dynamical Complexity of TASK network
    complexity = NeuralComplexity(MFR);
    fprintf('Dynamical Complexity: %.5f \n', complexity);
    
    %Raster plots of excitatory firings of TASK network  
    figure();
    subplot(2,1,1)
    plot(firings(:,1),firings(:,2),'k.', 'Markersize', 0.1);
    title(strcat('\beta_{0}=',num2str(p/10)));
    xlabel(['Timestep'])
    xlim([0 time])
    ylabel('Neuron Index')
    ylim([0 Ne_+1])
    
    %Raster plots of inhibitory firings of TASK network (Not needed but can do)
    
    %subplot(3,1,2)
    %plot(firings2(:,1),firings2(:,2),'k.', 'Markersize', 0.1);
    %title(strcat('p=',num2str(p/10)));
    %xlabel(['Time (ms)'])
    %xlim([0 time])
    %ylabel('Neuron Index')
    %ylim([0 200+1])
    
    % Plot the Mean Firing Rates of each excitatory module of TASK network
    subplot(2,1,2);
    plot(1:slide_window:time,MFR);
    xlabel(['Timestep']);
    ylabel('Mean Firing Rate');
    saveas(gcf,strcat('./results/MFR/_MFR_',num2str(p/10),'.png'),'png');
end

   
