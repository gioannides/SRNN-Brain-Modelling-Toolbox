function [layer] =  ConnectNetwork(CIJ,Ne_,Ni_)
% Constructs two layers of Izhikevich neurons and connects them together
% excitatory (layer 1)
NE = Ne_;
ME = 1;
% inhibitory (layer 2)
NI = Ni_;
MI = 1;
% layer construction
layer{1} = ConstructLayer('excitatory',NE,ME); % excitatory
layer{2} = ConstructLayer('inhibitory',NI,MI); % inhibitory
% layer conversion parameters
e = 1;
i = 2;
ex = 1:NE;
in = NE+1:NE+NI;

%---Setup Layers---
%excitatory-to-excitatory
layer = SetupLayer(layer,e,e,CIJ(ex,ex),17,randi([0,20],NE));
%excitatory-to-inhibitory
layer = SetupLayer(layer,e,i,CIJ(ex,in).*rand(NE,NI),50,ones(NE,NI));
%inhibitory-to-excitatory
layer = SetupLayer(layer,i,e,CIJ(in,ex).*-rand(NI,NE),2,ones(NI,NE));
%inhibitory-to-inhibitory
layer = SetupLayer(layer,i,i,CIJ(in,in).*-rand(NI,NI),1,ones(NI,NI));
end

function [layer2] = SetupLayer(layer,i,j,S,factor,delay)
    layer2 = layer;
    layer2{j}.S{i} = S';
    layer2{j}.factor{i} = factor;
    layer2{j}.delay{i} = delay';
end

function [ layer ] = ConstructLayer(type,N,M )
    r = rand(N,M);
    if strcmp(type,'excitatory')
        a = 0.02*ones(N,M);
        b = 0.2*ones(N,M);
        c = -65+15*r.^2;
        d = 8-6*r.^2; 
    elseif strcmp(type,'inhibitory')
        a = 0.02+0.08*r;
        b = 0.25-0.05*r;
        c = -65*ones(N,M);
        d = 2*ones(N,M);
    end
    layer.rows = N;
    layer.columns = M;
    layer.a = a;
    layer.b = b;
    layer.c = c;
    layer.d = d;
end

