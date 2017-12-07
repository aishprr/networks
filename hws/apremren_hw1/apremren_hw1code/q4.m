
function q4(N)
s = [];
t = [];
% N = number of nodes
for i = 0:N-1
    for e = -2:2
        if (e ~= 0)
            s(end + 1) = i + 1;
            t(end + 1) = mod(i + e, N) + 1;
        end
    end
end

rng(0, 'twister')
figure
G = digraph(s,t);
adj = adjacency(G);
subplot(3,1,1)
plot(G, 'Layout', 'force');
xlabel('Original Graph')

oave = avePathLength(adj);
ocoeff = clustCoeff(adj);

x = logspace(-4, 0, 30);
yave = [];
ycoeff = [];
for elem = x
    totalave = 0;
    totalcoeff = 0;
    for i = 1:4
        [newS, newT] = rewireEdge(s, t, N, elem);
        G2 = digraph(newS,newT);
        adj2 = adjacency(G2);
        totalave = totalave + avePathLength(adj2);
        totalcoeff = totalcoeff + clustCoeff(adj2);
    end
    yave(end + 1) = (totalave / 4) / oave;
    ycoeff(end + 1) = (totalcoeff / 4) / ocoeff;
end


subplot(3,1,2)
loglog(x, yave, '-s');
xlabel('Prob (p)');
ylabel('Ratio(Char Path Length)')

subplot(3,1,3)
loglog(x, ycoeff, '-s');
xlabel('Prob (p)');
ylabel('Ratio(Clustering Coefficient)')
end

function [newS, newT] = rewireEdge(s, t, N, p)

randStream = [];
for i = 1:N
    if (rand() < p)
        randStream(end + 1) = 1;
    else
        randStream(end + 1) = 0;
    end
end

origG = digraph(s,t);
origAdj = adjacency(origG);
newS = [];
newT = [];

for i = 1:length(s)
    newG = digraph(newS, newT);
    S = s(i);
    T = t(i);
    
    if ((randStream(S) == 1) && (addOffset(S, 1, N) == T ||  addOffset(S, 2, N) == T))
        % I want to try to rewire it
        retarget = randi(N-1);
        if (retarget >= S)
            retarget = addOffset(retarget, 1, N);
        end
        if (findedge(origG, S, retarget) == 0)
            % if it is not in the old graph, then try to add it to this one
            if (findedge(newG, S, retarget) ~=0)
                % if the other forward edge from this already caught this
                % target, then set it to the original one -- that one could
                % not have caught this one also, since it had to have 
                % checked against the original graph
                newS(end + 1) = S;
                newT(end + 1) = T;
            else
                newS(end + 1) = S;
                newT(end + 1) = retarget;
            end
        else
            newS(end + 1) = S;
            newT(end + 1) = T;
        end
    else
        newS(end + 1) = S;
        newT(end + 1) = T;
    end
end
end

function [res] = addOffset(a, add, N)
    if (a + add > N)
        res = a + add - N;
    elseif (a + add <= 0)
        res = N + a + add;
    else
        res = a + add;
    end
end
