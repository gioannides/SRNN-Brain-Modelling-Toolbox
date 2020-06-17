function [ CIJ, Ne_per_module ] = BuildTopology(p, Ne_, Ni_)
%BuildTopology creates the topology of the ex with rewiring probability p
rng(8);
CIJ = [];
% create 8 modules of Ne_/8 excitatory neurons 
% (with randConnections randomly assigned 1-way connections)
randConnections = 2000;
Modules = 8;
for i=1:Modules
    n = RandomNetwork(Ne_/Modules,randConnections);
    CIJ = ExtendNetwork(CIJ,n);
end
% rewire
CIJ = RewireNetwork(CIJ,Modules,p);
% create 1 module of Ni_ inhibitory neurons which is disconnected
CIJ = ExtendNetwork(CIJ,zeros(Ni_,Ni_));

% excitatory-to-inhibitory
% each inhibitory neuron has connections from four excitatory neurons
% (all from the same module)
% For the thesis, 100 was used for Ne_/8 exists
for j= Ne_+1:Ne_+Ni_
    mdule = randi([0,Modules-1]);
    rs = randi([1,Ne_/Modules],1,4);
    while (length(unique(rs)) ~= length(rs))
        rs = randi([1,Ne_/Modules],1,4);
    end
    i = mdule*(Ne_/Modules)*ones(1,4) + rs;
    CIJ(i,j) = 1;
end
assert(sum(sum(CIJ(1:Ne_,Ne_+1:Ne_+Ni_))) == Ne_);

% -- DIFFUSE -- 
% every inhibitory neuron projects to every neuron in the whole network
CIJ(Ne_+1:Ne_+Ni_,1:Ne_+Ni_) = 1;
Ne_per_module = Ne_/Modules;
%For thesis results set 'Ne_per_module = 100';
end

