delimiterIn = ',';
headerlinesIn = 1;
G_route = "Output/Debug/G.csv";
G_data = importdata(G_route,delimiterIn, headerlinesIn);
G_i = G_data.data(:,1);G_j = G_data.data(:,2);G_ij = G_data.data(:,3);
clear G_data;
clear G_route;
G_i = G_i + ones(size(G_i));G_j = G_j + ones(size(G_j));
G = sparse(G_i,G_j,G_ij);
G_sparsity =  size(G_ij);
G_sparsity = G_sparsity(1)/(max(G_i)^2);
clear G_i; clear G_j; clear G_ij;
disp("Not preconditioned magnitudes\n")
disp(sprintf('Sparsity is %.2e \n',G_sparsity))
disp(sprintf('Norm 2 condition number is %.4e\n', cond(full(G))))
disp(sprintf('Norm infinity condition number is %.4e\n', cond(full(G),inf)))
G_ndiag = G-speye(size(G));
delta = max_outlier(0.0,G_ndiag);
%disp(sprintf('Max outlier is %.4e\n',1.0*delta))
disp("Computing delta\n")
fun = @(x)(max_outlier(x,G_ndiag));
delta = fzero(fun,abs(delta));
delta = 2*delta;
E = G_ndiag;
E(G_ndiag>=-delta) = -1.0*G_ndiag(G_ndiag>=-delta);
E(G_ndiag<-delta) = delta;
%plot(sum(G+E,2))
%disp(sprintf('Norm 2 condition number is %.4e\n', cond(full(G+E))))
%disp(sprintf('Norm infinity condition number is %.4e\n', cond(full(G+E),inf)))
disp("Computing preconditioner\n")
C = -G_ndiag -E;
%Inv_approx = speye(size(G)) + C;
%for i = 2:10
%    Inv_approx = speye(size(G)) + C*Inv_approx;
%end
Inv_approx = speye(size(C))+C;
for l = 1:10
    Inv_approx = Inv_approx*(speye(size(C))+C^(2^l));
end
disp(sprintf('Norm 2 condition number E matrix is %.4e\n', cond(full(Inv_approx*G))))
disp(sprintf('Norm infinity condition number E matrix is %.4e\n', cond(full(Inv_approx*G),inf)))
%eps = max(abs(sum(G,2)))*1.1;
%Inv_approx = speye(size(G)) + C/eps;
%for i = 2:10
%    Inv_approx = speye(size(G)) + (C/eps)*Inv_approx;
%end
%Inv_approx = Inv_approx/eps;
%disp(sprintf('Norm 2 condition number Increasing Diagonal is %.4e\n', cond(full(Inv_approx*G))))
%disp(sprintf('Norm infinity condition number Increasing Diagonal is %.4e\n', cond(full(Inv_approx*G),inf)))