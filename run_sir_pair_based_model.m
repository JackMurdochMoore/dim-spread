% Adapted from
% pair_based_model.m by K. J. Sharkey (2010),
% which is available in the supplementary material for
% “Deterministic epidemic models on contact networks: correlations and
% unbiological terms”, Theoretical Population Biology, Volume 79, Issue 4,
% June 2011, Pages 115-129.
% 
% 
%
% Adapted by Jack Moore, 2024
% for
% Moore et al. (2024), "Network spreading from network dimension"  
%
function [tt, nnS, nnI] = run_sir_pair_based_model(A, lam, gam, tt)
[~, Y] = sir_pair_based_model(A, lam, gam, tt);
N = size(A, 1);
S = Y(:, 1:N);
I = Y(:, (N + 1):(2*N));
nnS = sum(S, 2);
nnI = sum(I, 2);
end

function [tt, Y] = sir_pair_based_model(A, lam, gam, tt)
%Solves the pair-based model for an arbitrary transmission network.
%A is a sparse matrix where A_{ij} is the strength of contact from site j
%to site i.
%lam is the infection rate.
%gam is the removal rate.
%tt is the list of time points at which to evaluate the solution.
%Returned is a matrix of length (tt) rows where each column specifies the
%probability time course for the susceptible/infectious states of each
%individual and the states of all IS SS and II pairs on the
%network.

A=sparse(A);%Ensure that transmission matrix is sparse

G=spones(A);%Contact network.

%Symmetrise contact network.
G_A=G.*(1-G');
H=G+G_A';

%Obtain information on the presence of closed and open triples.
Q=min(H^2,1);Q=Q-diag(diag(Q));%Nodes joined by a path-length of 2.
H_c=Q.*H;%Nodes joined by a path-length of 1 and 2.
H_o=Q-H_c;%Nodes joined by a path length of 2 but not 1.

%Set initial conditions for pairs.
F_AB=H;
F_AA=tril(H);

%Determine the locations of each pair.
N=length(A(:,1));
W_AB=find(reshape(F_AB,N^2,1));
W_AA=find(reshape(F_AA,N^2,1));
d_AB=length(W_AB);
d_AA=length(W_AA);

%Set initial state of infectious and susceptible populations.
I_vec=ones(N,1)/N;
S_vec=ones(N,1) - I_vec;

%Generate sparse diagonal matrices for multiplication
I_mat=spdiags(I_vec,0,N,N);
S_mat=spdiags(S_vec,0,N,N);

IS=I_mat*F_AB*S_mat;
SS=S_mat*F_AA*S_mat;%Symmetry means that we only need half the matrix.
II=I_mat*F_AA*I_mat;

Y0=[S_vec;I_vec;IS(W_AB);SS(W_AA);II(W_AA)];

options=odeset('abstol',1e-4); %Improves numerical stability for some networks.

[tt,Y]=ode45(@model_function,tt,Y0,options,A,H_c,H_o,lam,gam,N,W_AA,W_AB,d_AA,d_AB);

end

%______________________________________________________
function [dY] = model_function(t,Y0,T,H_c,H_o,lam,gam,N,W_AA,W_AB,d_AA,d_AB)

Y0=max(Y0,0);%Ensures that numerical errors do not cause negative values.

IS=spalloc(N,N,d_AB);
SS=spalloc(N,N,d_AB);
II=spalloc(N,N,d_AB);

%Generate vectors of individuals and matrices of pairs.
S=Y0(1:N);
I=Y0(N+1:2*N);
IS(W_AB)=Y0(2*N+1:2*N+d_AB);
SS(W_AA)=Y0(2*N+d_AB+1:2*N+d_AB+d_AA);
II(W_AA)=Y0(2*N+d_AB+d_AA+1:2*N+d_AB+2*d_AA);
SS=SS+SS'; %Regenerate symmetric matrices.
II=II+II';

%Diagonal matrices of reciprocal values.
inv_S=spdiags(spfun(@inve,S),0,N,N);
inv_I=spdiags(spfun(@inve,I),0,N,N);

%_______________________________________________________
%Kirkwood-type closure for triples.
%Here r denotes "right" and l denoted "left" indicating the
%relevant direction of network contact.

R=T.*IS;

%IrSrS triple.
Q=inv_I*(IS.*H_c);
IrSrS=(R'*Q).*(inv_S*SS*inv_S);%Closed part.
IrSrS=IrSrS+R'*H_o.*(inv_S*SS);%Addition of open part.

%IrSlI triple.
Q=(II.*H_c)*inv_I;
IrSlI=(inv_I*IS*inv_S).*(Q*R);%Closed part.
IrSlI=IrSlI+IS*inv_S.*(H_o*R);%Addition of open part.

SrSlI=IrSrS';%Due to symmetrising of networks.
IrSrI=IrSlI';

%_________________________________________________________
%Pair-based model equations
dT=sum(R)';
dS=-lam*dT;
dI=lam*dT-gam*I;
dSS=-lam*IrSrS-lam*SrSlI;
dIS=lam*IrSrS-lam*IrSlI-lam*R-gam*IS;
dII=lam*IrSlI+lam*IrSrI+lam*R+lam*R'-2*gam*II;

dY=[dS;dI;dIS(W_AB);dSS(W_AA);dII(W_AA)];

end

%_________________________________________________________
function M_out=inve(M_in)
%Determines the reciprocal value.
M_out=M_in.^(-1);
end