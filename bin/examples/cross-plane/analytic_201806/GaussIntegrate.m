% original file: lgwt.m
% Written by Greg von Winckel - 02/25/2004
%
% This script is for computing definite integrals using Legendre-Gauss 
% Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval
% [a,b] with truncation order N
%
% Suppose you have a continuous function f(x) which is defined on [a,b]
% which you can evaluate at any x in [a,b]. Simply evaluate it at all of
% the values contained in the x vector to obtain a vector f. Then compute
% the definite integral using sum(f.*w);
%

function [point,weight] = GaussIntegrationPoints(a,b,N)
% input:
%       a,b constitute the interval [a,b]
%       nGaussIntegralPoint means how many Legendre-Gauss points are required
% output:
%       point is a nGaussIntegralPoint*1 vector which is made up of the coordinates of Legendre-Gauss points
%       weight is the correspoind weight with respect to the point

% a = -1;b = 1;nGaussIntegralPoint = 15;      %% just for testing the function
N=N-1;
N1=N+1; N2=N+2;
% xu=linspace(-1,1,N1)';

% Initial guess
%y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);
%y=cos((2*(N:-1:0)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2); %modified by Yongxing Shen for an ascending output
y = cos((2*(N:-1:0)'+1)*pi/(2*N+2));

% Legendre-Gauss Vandermonde Matrix
L=zeros(N1,N2);
L(:,1)=1;
% Derivative of LGVM
Lp=zeros(N1,1); 

% Compute the zeros of the N+1 Legendre Polynomial
% using the recursion relation and the Newton-Raphson method
y0=2;

% Iterate until new points are uniformly within epsilon of old points
while max(abs(y-y0))>eps       
   % Lp(:,1)=0;
    
    L(:,2)=y;
   % Lp(:,2)=1;
    
    for k=2:N1
        L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
    end
    Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);   
    
    y0=y;
    y=y0-L(:,N2)./Lp;
end
% Linear map from[-1,1] to [a,b]
point=(a*(1-y)+b*(1+y))/2;      

% Compute the weights
weight=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;
end
%% the way by using defining symbol variable
% syms x;
% P = 1/(2^N*factorial(N))*diff((x^2-1)^N,N);  % Rodrigues' Formula
% dP = diff(P);
% point = double(solve(P));
% point = sort(point);
% weight = zeros(size(point));
% % weight = 2./((1-point.^2).*subs(dP,x,point).^2);
% for i = 1:length(point)
%     weight(i) = 2/((1-point(i)^2)*subs(dP,x,point(i))^2);
% end
