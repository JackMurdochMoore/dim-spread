% A correlated version of the BA network model.
% 
% Inclusivity model with inclusivity parameter r, introduced in
% Moore et al., "Inclusvity enhances robustness and efficiency of social
% networks", 2021.
% 
% N  -   final network size
% m  -   number of edges added per incoming node
% m0 -  growth begins from a fully connected network with m0 nodes (require
%       m0 >= m) 
% r  -  inclusivity parameter; each neighbour of an incoming node is chosen
%       from among those within network distance r of an established
%       neighbour (after r steps of a random walk starting at a
%       randomly chosen established neighbour)
% 
% Jack Murdoch Moore, 2021-09
% 
function A = inclusivity(N, m, m0, r)

assert((N >= m0) && (m0 >= m) && (r >= 1), 'Require (N >= m0) && (m0 >= m) && (r >= 1).');

A = logical(spalloc(N, N, 2*m0*(m0 - 1) + 2*(N - m0)*m));

A(1:m0, 1:m0) = ones(m0, m0) - eye(m0);

for ii = (m0 + 1):N
    nodes = 1:(ii-1);
    allDegrees = full(sum(A(nodes, nodes)));%Compute row vector containing degrees. I.e., Compute the sum of each column of A, i.e., the degree of each node.
    if (sum(allDegrees) == 0)%If no established nodes has any links then randomly choose an established node.
        allDegrees = ones(size(allDegrees));
        disp('No nonzero degrees.');
    end
    node = datasample(nodes, 1, 'Weights', allDegrees);
    
    targetNodes = node; numNSF = 1;%Keep track of the neighbours added so far.
    
    for numNSF = 2:m
        node = datasample(targetNodes, 1);
        while any(targetNodes == node)
            node = datasample(targetNodes, 1);
            for kk = 1:r
                node = datasample(find(A(node, nodes)), 1);
            end
        end
        targetNodes(numNSF) = node;%Keep track of the chosen node.
    end
    
    for node = targetNodes
        A(ii, node) = 1; A(node, ii) = 1;%Connect new node to chosen old node.
    end
    
end
end