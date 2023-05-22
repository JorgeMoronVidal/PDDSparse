n_subd = 10;
n_k_pi = 64;
[NS,NKPI] = meshgrid(n_subd,n_k_pi);
NKT = 2*NS.*(NS-1).*NKPI-3*(NS-1).^2;
figure;
surf(NS.^2,NKPI,NKT);
set(gca,"zscale","log")
%set(gca,"yscale","log")
%set(gca,"xscale","log")
set(gca,"ColorScale","log")
xlabel("Number of Subdomains")
ylabel("Number of Knots per Interface")
zlabel("Total number of knots")



