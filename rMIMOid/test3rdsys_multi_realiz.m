%
% test3rdsys_multi_realiz.m
%
% Compara vers?o recursiva e batelada do algoritmo baseados na parametriza??o MOLI
% baseado no resultado de uma ?nica realizacao

clear all;
close all;
clc

%% Primer configurations
Nruns = 200;
N = 103;

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
real_poles = ones(N,1)*sort(eig(A))';

%%
% PBSID - Recursive Subspace Identification settings
n = 3;          % number of states
p = n + 2;      % past window size
f = n;          % past window size
lambda = 1;		% forgetting factor
ireg = 1e-6;    % initial regularisation

rlsopts = struct('ireg',[ireg ireg ireg],'lambda',[lambda lambda lambda],'reg',0);
idopts = struct('method','varx','weight',1,'ltv',1,'noD',1,'past',0,'Kalm',0);
%default: struct('method','varx','weight',0,'ltv',0,'noD',0,'past',0,'Kalm',0);

%% Data generation & model estimation

moli_eig = zeros(Nruns,N,length(A));
rmoli_eig = zeros(Nruns,N,length(A));
rpbsid_eig = zeros(Nruns,N,length(A));
detP = zeros(N,Nruns);

for r = 1:Nruns
	
	rng(r);	% set seed
	
	% I/O data
	x = zeros(size(A,1),N);
	y = zeros(size(C,1),N);
	
	u = randn(size(B,2),N);	
	w = randn(size(W,2),N);
	v = randn(size(V,2),N);
	
	
	for k = 1:N-1
		x(:,k+1) = A*x(:,k) + B*u(:,k) + W*w(:,k);
		y(:,k) = C*x(:,k) + V*v(:,k);
	end
	
	y = y';
	u = u';
	
	% Parameter estimation
	
	l = [2,1];
	% MOLI
	[a,b,c] = moli(y, u, l, poly([0 0]));
	moli_eig(r,:,:) = ones(N,1)*sort(eig(a))';
	%ssx = ss(a,b,c,zeros(2), 1);
	
	% rMOLI
	[ar,br,cr,P,theta_cell,A_D,C_D] = rmoli(y, u, l, poly([0 0]));
	%ssr = ss(ar,br,cr,zeros(2), 1);
	
	theta = cell(size(C,1));
	for k=1:N
		for i = 1:size(C,1), theta{i} = theta_cell{i,k}; end
		rmoli_eig(r,k,:) = sort(eig(theta2abc(theta,l,size(B,2),A_D,C_D)));
		detP(k,r) = det(blkdiag(P{1,k},P{2,k}));
	end
	
% 	figure(3);
% 	for i=1:length(A)
% 		subplot(3,1,i)
% 		plot((1:N)', real_poles(:,i),'-k','Linewidth',1.1); hold on;
% 		plot((1:N)', real(rmoli_eig(r,:,i))','r--','Linewidth',0.5);
% 		set(gca,'Fontsize',12);
% 		xlim([30,N])
% 		grid on;
% 	end
	
	% rPBSID
	[Ak,Bk,Ck,Dk,Kk,err1,eigA1,dampA1] = rpbsid(u,y,f,p,n,[],idopts,rlsopts);
	rpbsid_eig(r,:,:) = eigA1';
	%sspbsid = ss(Ak,Bk,Ck,Dk,1);

end

%% Report results

figure(1);

moli_pole_mean = zeros(N,length(A));
rmoli_pole_mean = zeros(N,length(A));
rpbsid_pole_mean = zeros(N,length(A));

moli_pole_std = zeros(N,length(A));
rmoli_pole_std = zeros(N,length(A));
rpbsid_pole_std = zeros(N,length(A));

for i=1:length(A)
	
	for k=1:N
		moli_pole_mean(k,i) = mean(moli_eig(:,k,i));
		rmoli_pole_mean(k,i) = mean(rmoli_eig(:,k,i));
		rpbsid_pole_mean(k,i) = mean(rpbsid_eig(:,k,i));
		
		moli_pole_std(k,i) = std(real(moli_eig(:,k,i)),1);
		rmoli_pole_std(k,i) = std(real(rmoli_eig(:,k,i)),1);
		rpbsid_pole_std(k,i) = std(real(rpbsid_eig(:,k,i)),1);
	end
	
	subplot(3,1,i)
	plot((1:N)', real_poles(:,i),'-r','Linewidth',1.1); hold on;
	%plot((1:N)', real(moli_pole_mean(:,i)),':','color',[.3 .3 .3],'Linewidth',3);
	
	plot((1:N)', real(rpbsid_pole_mean(:,i)),'--','color',[.6 .6 .6],'Linewidth',1.7);
	plot((1:N)', real(rpbsid_pole_mean(:,i)) + 3*rpbsid_pole_std(:,i),...
		'--','color',[.6 .6 .6],'Linewidth',1);
	plot((1:N)', real(rpbsid_pole_mean(:,i)) - 3*rpbsid_pole_std(:,i),...
		'--','color',[.6 .6 .6],'Linewidth',1);
	
	plot((1:N)', real(rmoli_pole_mean(:,i)),'-.','color',[.3 .3 .3],'Linewidth',1.7);
	plot((1:N)', real(rmoli_pole_mean(:,i)) + 3*rmoli_pole_std(:,i),...
		'-.','color',[.3 .3 .3],'Linewidth',1);
	plot((1:N)', real(rmoli_pole_mean(:,i)) - 3*rmoli_pole_std(:,i),...
		'-.','color',[.3 .3 .3],'Linewidth',1);

	set(gca,'Fontsize',12);
	xlim([1,25])
	grid on;
end

%legend('actual','batch MOLI','rPBSID','rMOLI','Orientation','Horizontal',4);
xlabel('Samples');

% ---------
% figure(2)
% 
% detPmean = mean(detP,2)';
% detPstd = std(detP,1,2)';
% 
% semilogy(detPmean + 3*detPstd,'--','color',[.3 .3 .3],'Linewidth',1); hold on;
% semilogy(detPmean,'-','color',[0 0 0],'Linewidth',1.4);
% set(gca,'Fontsize',12);
% ylabel('det({\it P_k})');
% xlabel('Samples');
% xlim([0,N])
% grid on;

