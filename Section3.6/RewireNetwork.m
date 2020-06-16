function [ CIJ ] = RewireNetwork(CIJ,m,p)
%RewireNetwork rewires a network CIJ of m modules with probability p
if p == 0
    return
end
[N,M] = size(CIJ);
assert(N==M);
msize = N/m;
for i = 1:N
    for j = 1:M
        % if intra-community connection
        if CIJ(i,j) == 1 && floor((i-1)/msize) == floor((j-1)/msize)
            % we rewire with probability p
            if rand < p
                % we lose the intra-connection
                CIJ(i,j) = 0;
                % current module
                cm = floor((i-1)/msize);
                % possible target modules
                tms = [0:(cm-1),(cm+1):(m-1)];
                % target module
                tm = tms(randi(m-1));
                % target node
                t = tm*msize + randi(msize);                
                while(CIJ(i,t) == 1)
                    t = tm*msize + randi(msize);
                end
                CIJ(i,t) = 1;
            end
        end
    end
end



end

