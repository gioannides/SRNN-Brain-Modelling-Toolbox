function CIJ = RandomNetwork(N,randConns)
% Creates a random network with N nodes and C randomly assigned one-way connections
CIJ = zeros(N,N);
for i = 1:randConns
    f = 0;
    t = 0;
    while f == t || CIJ(f,t) == 1
        f = randi(N);
        t = randi(N);
    end
    CIJ(f,t) = 1;
    %Uncomment this line for mixed Undirected network (remove if statement for fully undirected topology)
    %if rand(1) > 0.1
        %CIJ(t,f) = 1;
    %end
end