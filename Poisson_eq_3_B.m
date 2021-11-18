#! octave-qf
%p36.m (Modified by Jorge Moron) - Poisson eq. on [-1,1]x[-1,1] with nonzero BC's
%Set up grid and 2D Laplacian, boundary points included:
args = argv();
id = args{1};
file = sprintf("Input/Buffer/parameters_%s.csv", id);
table = csvread(file);
parameters = table(:,1);
clear table;
file = sprintf("Input/Buffer/positions_%s.csv", id);
table = csvread(file);
x_knot = table(:,1);
y_knot = table(:,2);
clear table;
N = 40;
x_stencil = parameters(1) + (parameters(3)-parameters(1))*(0:39)/(N-1);
y_stencil = parameters(2) + (parameters(4)-parameters(2))*(0:39)/(N-1);
[Dx,Dy,x,y] = cheb(N,x_stencil,y_stencil);
[xx,yy] = meshgrid(x,y); xx = xx(:); yy = yy(:);
D2x = Dx^2; D2y = Dy^2; I = eye(N+1); L = kron(D2x,I) + kron(I,D2y);
%Impose boundary conditions and -f function by replacing appropriate rows of L:
b = find(xx==parameters(1) | xx == parameters(3) | yy==parameters(2) | yy == parameters(4));
% boundary pts
L(b,:) = zeros(4*N,(N+1)^2); 
L(b,b) = eye(4*N);
omegax = 0.23; omegay = 0.49; omegapx = 0.331; omegapy = 0.667;
rhs = -pi*pi*(omegax*omegax + omegay*omegay)*sin(omegax*pi*xx + omegay*pi*yy);
rhs += -pi*pi*(omegapx*omegapx + omegapy*omegapy)*cos(omegapx*pi*xx + omegapy*pi*yy);
%rhs(b) = sin(omegax*pi*xx(b) + omegay*pi*yy(b)) + cos(omegapx*pi*xx(b) + omegapy*pi*yy(b));
rhs(b) = 0*xx(b);
% Solve Poisson equation, reshape to 2D, and plot:
u = L\rhs; uu = reshape(u,N+1,N+1);
xx = reshape(xx,N+1,N+1);
yy = reshape(yy,N+1,N+1);
file = sprintf("Input/Buffer/B_%s.csv", id);
B = interp2(xx,yy,uu,x_knot,y_knot,'cubic');
save("-ascii",file,"B");
