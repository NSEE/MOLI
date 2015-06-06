%
% test3rdsys_single_realiz.m
%
% Compara vers?o recursiva e batelada do algoritmo baseados na parametriza??o MOLI
% baseado no resultado de uma ?nica realizacao

clear all;
%close all;
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

m0 = ss(A,B,C,zeros(2),1);

%% I/O data

N = 103;
x = zeros(size(A,1),N);
y = zeros(size(C,1),N);

rng(15);
u = randn(size(B,2),N);
w = randn(size(W,2),N);
v = randn(size(V,2),N);


for k = 1:N-1
	x(:,k+1) = A*x(:,k) + B*u(:,k) + W*w(:,k);
	y(:,k) = C*x(:,k) + V*v(:,k);
end

y = y';
u = u';

%% Parameter estimation

l = [2,1];	%list of observability indices

tic
[a,b,c] = moli(y, u, l, poly([0 0]));
ssx = ss(a,b,c,zeros(2), 1);
toc

tic
[a,b,c,P,theta_cell,A_D,C_D] = rmoli(y, u, l, poly([0.15 0.15]));
ssr = ss(a,b,c,zeros(2), 1);
toc

tic
[theta_cell_zoft,A_D_cell,C_D_cell,Pzoft] = rmoli_zoft(y, u, l, 1);
toc

rmoli_eig = zeros(length(A),N);
rmolizoft_eig = zeros(length(A),N);
theta = cell(size(C,1));
detP = zeros(N,1);
detPzoft = zeros(N,1);


for k=2:N
	theta = theta_cell(:,k);
	rmoli_eig(:,k) = sort(eig(theta2abc(theta,l,size(B,2),A_D,C_D)));
	rmolizoft_eig(:,k) = sort(eig(theta2abc(theta_cell_zoft(:,k),l,size(B,2),...
		A_D_cell{k},C_D_cell{k})));
	detP(k) = det(blkdiag(P{1,k},P{2,k}));
	detPzoft(k) = det(blkdiag(Pzoft{end,1,k},Pzoft{end,2,k}));
end

% Recursive Subspace Identification parameters
n = 3;          % number of states
p = n + 2;      % past window size
f = n;          % past window size
lambda = 1;		% forgetting factor
ireg = 1e-6;    % initial regularisation
rlsopts = struct('ireg',[ireg ireg ireg],'lambda',[lambda lambda lambda],'reg',0);

% Start Recursive Subspace Identification
idopts = struct('method','varx','weight',1,'ltv',1,'noD',1,'past',0,'Kalm',0);
% 'Kalm' = 1 demands more computation effort
%default: struct('method','varx','weight',0,'ltv',0,'noD',0,'past',0,'Kalm',0);
tic
[Ak,Bk,Ck,Dk,Kk,err1,eigA1,dampA1] = rpbsid(u,y,f,p,n,[],idopts,rlsopts);
toc
sspbsid = ss(Ak,Bk,Ck,Dk, 1);

%% Report results
sort([eig(A) eig(ssx.A) eig(ssr.A) eig(Ak)])

real_poles = ones(N,1)*sort(eig(A))';
moli_poles = ones(N,1)*sort(eig(ssx.A))';


figure(1);
for i=1:length(A)
	subplot(3,1,i)
	plot((1:N)', real_poles(:,i),'-k','Linewidth',1.1); hold on;
	plot((1:N)', real(rmolizoft_eig(i,:))','r--','Linewidth',1.6);%,'color',[.3 .3 .3],'Linewidth',3);
	plot((1:N)', real(eigA1(i,:))','-.','color',[.6 .6 .6],'Linewidth',1.6);
	plot((1:N)', real(rmoli_eig(i,:))','x:','color',[0 0 1],'Linewidth',1.6);
	set(gca,'Fontsize',12);
	xlim([0,N])
	grid on;
end

legend('actual','rMOLI-ZOFT','rPBSID','rMOLI','Orientation','Horizontal',4);
xlabel('Samples');

% ---------
figure(2)
semilogy(detP,':','color',[.3 .3 .3],'Linewidth',1.6);
hold on;
semilogy(detPzoft,'r--','Linewidth',1.6);
set(gca,'Fontsize',12);
ylabel('det({\it P_k})');
xlabel('Samples');
xlim([0,N])
grid on;

%%
%step(m0,ssx,ssr,sspbsid)
