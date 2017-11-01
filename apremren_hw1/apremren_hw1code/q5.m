function [a] = q5(m, n)

done = 0;
initial = randi([0 1], m,m);
initialAdj = triu(initial) + triu(initial, 1)';

origMatrix = zeros(n);
origMatrix(1:m, 1:m) = initialAdj(1:m, 1:m);
G = graph(origMatrix);

for i = 1:m
    if (findedge(G, i, i) ~= 0)
        G = rmedge(G, i, i);
    end
end

adj = adjacency(G);

while ~done
    % go through and find an unconnected node
    un = -1;
    for i = 1:n
        if (degree(G, i) == 0)
            un = i;
        end
    end
    % found my unconnected node
    if (un == -1)
        % you are done, no more unconnected
        done = 1;
        break;
    end
    
    for i = 1:n
        % the degree of this one should be 0 here
        if (degree(G, un) >= m) 
            break;
            
        end
        if (i == un || (findedge(G, i, un) ~= 0))
            continue;
        end
        prob = degree(G, i) / sum(sum(adj));
        
        if (rand() <= prob)
            % connect the two
            G = addedge(G, un, i, 1);
        end
        adj = adjacency(G);
    end
end

edges = numedges(G);
newEdges = 0;
newGraph = addnode(graph(), n);
newAdj = adjacency(newGraph);
p = 0.05;

while (newEdges < edges)
   
    s = randi([1 n],1,1);
    t = randi([1 n],1,1);
    if (s == t || (rand() > p))
        continue
    end
    if (findedge(newGraph, s, t) == 0)
        newGraph = addedge(newGraph, s, t, 1);
        newEdges = newEdges + 1;
    end
end
figure
plot(G)
degPlot(G);
figure
plot(newGraph)
degPlot(newGraph);

clustCoeff1 = clustCoeff(adj)
clustCoeff1 = clustCoeff(adjacency(newGraph))

end

function degPlot(G)
    n = numnodes(G);
    degList= degree(G);
    histogram(degList,'BinMethod', 'integers')
    xx = unique(degList);
    M = max(xx)
    m = min(xx)
    for ii=1:numel(xx)
      xpdf(ii) = xx(ii);
      ypdf(ii) = length(find(xx(ii)==degList)) / length(degList);
    end
    figure
    plot(xpdf,ypdf)
end