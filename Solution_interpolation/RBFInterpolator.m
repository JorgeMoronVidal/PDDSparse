function yev = RBFInterpolator(xi, RHS,xev)
%[max_er,cond] = testRBFInterpolator([0, rand(1,40),1], @(x)1+ sin(3*x).*exp(x))
%Construimos la matriz de inteprolacion
xi = xi(~isnan(xi));
RHS = RHS(~isnan(RHS));
N = length(xi);
%Nodos de Chebychev
%xj = NaN(1,N-2);
%for k = 1:N-2
%    xj(k) =0.5* (1 + cos((2*k -1)*pi/(2 * N)));
%end
%xj = [1, xj, 0]
xj = xi;
rbf = 'mq';


c2j = ones(N,1)*(0.5*0.5);
MAT = NaN(N);

for k = 1:N
    MAT(k,:) = feval(rbf,xi(k),xj,c2j);
end

condi = cond(MAT)
%Calcular alphas

alpha = MAT\RHS;
%Evaluar en un grid fino
Nev = length(xev);
yev = NaN(Nev,1);
for i = 1:Nev
    esto = feval(rbf,xev(i),xj,c2j);
    yev(i) = sum(alpha.*feval(rbf,xev(i),xj,c2j));
end
%Error

%err = yev - yex;
%errel = err./yev;
%merr = max(err);

%figure
%subplot(1,2,1), plot(xev,yev,'r-',xev,yex,'b-'), %legend('interpolated value','exact value')
%hold on, for k = 1:N, plot(xi(k), sum(alpha'.*feval(rbf,xi(k),xj,c2j)),'ko'),end
%subplot(1,2,2), plot(xev,errel,'k-')
%hold on, for k = 1:N, plot(xi(k),0.0,'ro', xj(k),0.0, 'bo'),end

return
%multricuadrica
function res = mq(x,xj,c2j)
    res = sqrt((x-xj).*(x-xj) + c2j);
return
%multicuadrica inversa
function res = imq(x,xj,c2j)
    res = 1./sqrt((x-xj).*(x-xj) + c2j);
return