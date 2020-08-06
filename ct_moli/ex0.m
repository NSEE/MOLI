%
% CT-MOLI example
%

%% Data-generating model

% Use one to update the model in each run
if(0 || ~exist('m0','var'))
    nu = 2; ny = 2;
    m0 = rss(4,nu,ny);
    %[A0,B0,C0,~] = ssdata(m0);
    % Avoid improper models
    m0.D = zeros(ny,nu);
end

%% Generate data

N = 500;
try
    u_ = prbs(12,40);
    u = [u_(1:N),u_(end-N+1:end)];
catch
    error("N and PRBS parameters not compatible");
end
% Sampling period
Ts = 0.1;
% Time vector
t = Ts*((0:N-1)');
% Output data
ye = lsim(m0,u,t);


%% MOLI parameterization & filter settings
% List of observability indices
l = [2;2];



%%

% wc = 0.1*2*pi/Ts;
% zeta = 0.5;
% alphact = [1 2*zeta*wc wc^2];
% [A,B,C,epsilon] = ctmoli(ye, u, Ts, l, alphact);
% [eig(m0.A), eig(A)]

% Create data structure;
dte.y = ye; dte.u = u; dte.Ts = Ts;

wc = logspace(log10(1/40*Ts),log10(1/4*Ts),15);

[Ab,Bb,Cb,alphab] = butterzoft(ye,u,Ts,l,wc,ye);
mb = ss(Ab,Bb,Cb,zeros(ny,nu));
zeta = [.6,.7,1];
[Az,Bz,Cz,alphaz] = ctmolizoft(dte,l,wc,zeta,@evalBFRc,dte);
mz = ss(Az,Bz,Cz,zeros(ny,nu));


figure(3)
step(m0,mb,mz,c2d(m0,Ts))
sort([eig(m0.A), eig(Ab)])

figure(4)
yb = lsim(mb,u,t);
yz = lsim(mz,u,t);
subplot(2,1,1)
plot([ye(:,1),yb(:,1),yz(:,1)])
subplot(2,1,2)
plot([ye(:,2),yb(:,2),yz(:,2)])