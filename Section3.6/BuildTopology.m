function [ CIJ, Ne_per_module ] = BuildTopology(p, Ne_, Ni_)
%BuildTopology creates the topology of the ex with rewiring probability p
rng(8);
CIJ = [];
% create 8 modules of Ne_/8 excitatory neurons 
% (with 1000 randomly assigned 1-way connections)

for i=1:8
    n = NetworkDirectedRandom(Ne_/8,2000);
    CIJ = ExpandNetwork(CIJ,n);
end
% rewire
CIJ = RewireNetwork(CIJ,8,p);
% create 1 module of 200 inhibitory neurons (disconnected)
CIJ = ExpandNetwork(CIJ,zeros(Ni_,Ni_));

% excitatory-to-inhibitory
% each inhibitory neuron has connections from four excitatory neurons
% (all from the same module)
% For the thesis, 100 was used for Ne_/8 exists
for j= Ne_+1:Ne_+Ni_
    md = randi([0,7]);
    rs = randi([1,Ne_/8],1,4);
    while (length(unique(rs)) ~= length(rs))
        rs = randi([1,Ne_/8],1,4);
    end
    i = md*(Ne_/8)*ones(1,4) + rs;
    CIJ(i,j) = 1;
end
%disp(sum(sum(CIJ(1:2250,2251:3000))));
assert(sum(sum(CIJ(1:Ne_,Ne_+1:Ne_+Ni_))) == Ne_);

% -- DIFFUSE -- 
% every inhibitory neuron projects to every neuron in the whole network
CIJ(Ne_+1:Ne_+Ni_,1:Ne_+Ni_) = 1;
Ne_per_module = Ne_/8;
%For thesis results set 'Ne_per_module = 100';
end

