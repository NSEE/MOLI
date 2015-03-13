

clear all;
close all;
clc

%% Plant
A = [.8,-.4,.2;...
	0,.3,-.5;...
	0,0,.5];

B = [0,0;...
	0,-.6;...
	.5,0];

C = [.5,.5,0;...
	0,0,1];

W = [.055;.04;.045];

V = [.025;.03];

%% I/O data

N = 103;
x = zeros(size(A,1),N);
y = zeros(size(C,1),N);

rng(1);
u = rand(size(B,2),N);

rng(2);
w = 0*randn(size(W,2),N);
v = 0*randn(size(V,2),N);


for k = 1:N-1
	x(:,k+1) = A*x(:,k) + B*u(:,k) + W*w(:,k);
	y(:,k) = C*x(:,k) + V*v(:,k);
end

y = y';
u = u';

%% Parameter estimation

m0 = ss(A,B,C,zeros(2),1);

[a,b,c] = moli(y, u, [2,1], poly([0 0]));
ssx = ss(a,b,c,zeros(2), 1);

[a,b,c,P] = rmoli(y, u, [2,1], poly([0 0]));
ssr = ss(a,b,c,zeros(2), 1);

dt = iddata(y,u,1);
step(m0,ssx,pem(dt,3),ssr)
