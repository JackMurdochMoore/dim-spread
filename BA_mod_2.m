% Generate BA scale-free network
% 
% Moore et al. (2024), "Network spreading from network dimension"
%
% Jack Moore, 2024
%
function A = BA_mod_2(N, m, m0)

A = spalloc(N, N, m0*(m0 - 1) + 2*(N - m0)*m);
A(1:m0, 1:m0) = 1;
A(1:m0, 1:m0) = A(1:m0, 1:m0) - eye(m0);

for ii = (m0 + 1):N
    nodes = 1:(ii-1);
    degSeq = full(sum(A(nodes, nodes)));
    if (sum(degSeq) == 0); degSeq = ones(size(degSeq)); end
    nodesToAdd = datasample(nodes, m, 'Weights', degSeq, 'Replace', false);
    for node = nodesToAdd
        A(ii, node) = 1; A(node, ii) = 1;%Connect new node to chosen old node.
    end
end

%A=full(A);