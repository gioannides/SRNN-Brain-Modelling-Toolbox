function [ CIJ ] = BuildTopology(p, Ne_, Ni_)
%BuildTopology creates the topology of the ex with rewiring probability p
rng(8);
CIJ = [];
Conns = 2000;
% create 8 modules of excitatory neurons 
% (with Conns randomly assigned 1-way connections)
for i=1:8
    n = NetworkDirectedRandom(Ne_/8,Conns);
    CIJ = ExpandNetwork(CIJ,n);
end
% rewire
CIJ = RewireNetwork(CIJ,8,p);
% create 1 module of inhibitory neurons (disconnected)
CIJ = ExpandNetwork(CIJ,zeros(Ni_,Ni_));

% excitatory-to-inhibitory
% each inhibitory neuron has connections from four excitatory neurons
% (all from the same module)
for j= Ne_+1:Ne_+Ni_
    md = randi([0,7]);
    rs = randi([1,100],1,4);
    while (length(unique(rs)) ~= length(rs))
        rs = randi([1,100],1,4);
    end
    i = md*100*ones(1,4) + rs;
    CIJ(i,j) = 1;
end
assert(sum(sum(CIJ(1:Ne_,Ne_+1:Ne_+Ni_))) == Ne_);

% -- DIFFUSE -- 
% every inhibitory neuron projects to every neuron in the whole network
CIJ(Ne_+1:Ne_+Ni_,1:Ne_+Ni_) = 1;
end

