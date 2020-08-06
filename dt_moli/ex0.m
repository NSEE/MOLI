%
% CT-MOLI example
%

clc;

%% Data-generating model

% Use one to update the model in each run
if(0 || ~exist('m0','var'))
    nu = 2; ny = 2;
    m0 = drss(4,nu,ny);
    %[A0,B0,C0,~] = ssdata(m0);
    % Avoid improper models
    m0.D = zeros(ny,nu);
end

%% Generate data

N = 500;
try
    u_ = prbs(12,20);
    u = [u_(1:N),u_(end-N+1:end)];
catch
    error("N and PRBS parameters not compatible");
end
% Sampling period
Ts = 1;
m0.Ts = Ts;
% Time vector
t = Ts*((0:N-1)');
% Output data
ye0 = lsim(m0,u,t);
ye = ye0 + randn(size(u))*[0.05,0;0,0.01];


%% MOLI parameterization & filter settings
% List of observability indices
l = [2;2];



%%


wc = linspace(.1,1,11)*pi/Ts;
zeta = [0.3,.5, .7, .9];
sfactor = 2;
opteta = zoft(ye,u,l,Ts,@evalBFR,wc,zeta,sfactor);
alpha = genalphapoly(opteta,l,Ts);
[a,b,c,~] = moli(ye,u,l,alpha);
m = ss(a,b,c,zeros(ny,nu),Ts);


figure(3)
step(m0,m)
sort([eig(m0.A), eig(a)])

figure(4)
ysim = lsim(m,u);
subplot(2,1,1)
plot([ye(:,1),ysim(:,1)])
subplot(2,1,2)
plot([ye(:,2),ysim(:,2)])