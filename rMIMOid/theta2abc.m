% Computes the state space model matrices
% Autor: Rodrigo Alvite Romano
%
function [a,b,c] = theta2abc(theta,l,m,A_D,C_D)

p = length(l);
n = sum(l);

D = zeros(n,p);
B = zeros(n,m);
G = zeros(p);

for jp = 1:p
    ct = 0;
    for i = 1:p
        D(ct+1:ct+l(i),jp) = theta{i}(1 + l(i)*(jp-1):l(i)*jp);
        ct = ct + l(i);
    end
    
    if(jp > 1), G(jp,1:jp-1) = theta{jp}((m+p)*l(jp)+1:end); end
end

for jm = 1:m
    ct = 0;
    for i = 1:p
        B(ct+1:ct+l(i),jm) = theta{i}(1 + l(i)*p + l(i)*(jm-1):l(i)*(p+jm));
        ct = ct + l(i);
    end
end

c = (eye(p) - G)\C_D;
b = B;
a = A_D + D*c;