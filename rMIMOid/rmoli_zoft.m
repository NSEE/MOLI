function [theta_cell,A_D_cell,C_D_cell,P] = rmoli_zoft(y, u, l, T)%, Omega0)
% ---------------------------------------------------------------------------------------
% Function to recursively estimate state-space multivariable model parameters. In an
% inner loop the plant parameter vector (\theta) is recursively updated. The filter
% tuning is addressed in an outer loop, which is slower than the inner one.
%
% Author: Rodrigo A. Romano - Mai/2015
%
% ---------------------------------------------------------------------------------------
% function [theta_cell,A_D_cell,C_D_cell,P] = rmoli_zoft(y, u, l, T)
%
% ---------------------------------------------------------------------------------------
% Inputs:
%
% - y: output data matrix (N x p, where "p" is the number of outputs)
% - u: input data matrix (N x m, where "m" is the number of inputs)
% - l: list of observability indices
% - T: sampling period
%
% Outputs:
%
% - The parameter vector relative to the ith output computed in iteration k is
% provided through theta_cell{i,k}.
% - Design model matrices: A_D and C_D
% - P is a cell array with de covariance matrices relative to the parameter
% estimation. P{i,k} stores the covariance relative to the ith output computed in
% iteration k
% ---------------------------------------------------------------------------------------

%% Verification of the input arguments

if(nargin < 4), error('too few input arguments...'); end

[N, p] = size(y);
[~, m] = size(u);
n_ = max(l);

if(length(l) ~= p)
	error('The length of l must be equal to the number of outputs');
end

%if(length(alpha) ~= n_+1), error('alpha order mismatch!'); end


%% ZOFT design variables
Talpha = 30;	% sampling instant multiple to update the filter
wc = (1/T)*abs(log([.3 .25 .1 .15]));
zeta = 1;

% Initialize set of test points (Omega) using wc and zeta vectors
[alphacand, Omega] = genalphapoly(wc, zeta, 2, l, T);
n_eta = size(Omega,1);

% MOLI matrices
diagM_T_set = cell(n_eta,1);
C_D_set = cell(n_eta,1);
A_D_set = cell(n_eta,1);
a = cell(n_eta,1);
b = cell(n_eta,1);
c = cell(n_eta,1);
filtercoeff = zeros(n_eta,n_);

for r = 1:n_eta
	[diagM_T_set{r}, C_D_set{r}, A_D_set{r}, filtercoeff(r,:)] =...
		comp_MOLI_arrays(alphacand(r,:),m,l);
end

%% Initial conditions for recursive prameter estimation

P = cell(n_eta,p,N);	% Covariance matrices
R1 = cell(p,1);	% Parameter covariance matrices
R2 = cell(p,1);	% Output variance
theta_cell = cell(p,N);
A_D_cell = cell(N,1);
C_D_cell = cell(N,1);
theta = cell(n_eta,p);	% Parameter cell array

ini = n_;
for i = 1:p
	theta_cell{i,ini-1} = zeros((m+p)*l(i) + i - 1,1);
	R1{i} = zeros((m+p)*l(i) + i - 1); %1e-3*eye((m+p)*l(i) + i - 1);%
	R2{i} = 1e-1;
	
	for r = 1:n_eta
		P{r,i,ini-1} = eye((m+p)*l(i) + i - 1)*10^6;
		theta{r,i} = theta_cell{i,ini-1};
	end
end

%% Recursive parameter estimation

z = [y,u];
zf = zeros(N,m+p,n_eta);
xi_k = zeros((m+p)*n_,n_eta);

for k = ini:N
	
	for r = 1:n_eta % -- suitable for parallel computing
	
		% Identifier state recursive computation
		zf(k+1,:,r) = filtercoeff(r,:)*zf(k-n_+1:k,:,r) + z(k,:);
		for ii=1:(p+m)
			xi_k(1+(ii-1)*n_:ii*n_,r) = zf(k-n_+1:k,ii,r);
		end
		
		for i = 1:p
			if(l(i) == n_)
				phi_i_k = [xi_k(:,r);y(k,1:i-1)];
			else
				phi_i_k = [diagM_T_set{r}{i}*xi_k(:,r);y(k,1:i-1)];
			end
			
			% Recursive LS - Kalman filter algorithm
			aux = phi_i_k'*P{r,i,k-1}*phi_i_k + R2{i};
			K_k = (P{r,i,k-1}*phi_i_k)/aux;
			theta{r,i} = theta{r,i} + K_k*(y(k,i)-phi_i_k'*theta{r,i});
			P{r,i,k} = P{r,i,k-1} + R1{i} -...
				(P{r,i,k-1}*(phi_i_k*phi_i_k')*P{r,i,k-1})/aux;
		end
		
		if(~rem(k,Talpha))
			 [a{r},b{r},c{r}] = theta2abc(theta(r,:),l,m,A_D_set{r},C_D_set{r});
		end
		
	end
	
	theta_cell(:,k) = theta(n_eta,:);
	A_D_cell{k} = A_D_set{n_eta};
	C_D_cell{k} = C_D_set{n_eta};

	if(~rem(k,Talpha))
		Omega(n_eta,:) = barycenter(y(k-Talpha+1:k,:),u(k-Talpha+1:k,:),l,...
			Omega, alphacand, a, b, c);
		optalpha = genalphapoly(Omega(n_eta,:),l,T);
		
		if((mod(max(l),2) == 0) && any(mod(l,2)))
			if(~isreal(roots(optalpha)))
				optalpha = push_poly_roots(optalpha);
			end
		end
		
		alphacand(n_eta,:) = optalpha;
		[diagM_T_set{n_eta}, C_D_set{n_eta}, A_D_set{n_eta}, filtercoeff(n_eta,:)] =...
			comp_MOLI_arrays(optalpha,m,l);
	end
	
end

%% --------------------------------------------------------------------------------------
function opteta = barycenter(y, u, l, Omega, alphacand, A, B, C)

% Zero-order oracle (derivative-free) filter tuning
% opteta = barycenter(y, u, l, Omega, alphacand, A, B, C)
%
%%

p = length(l);
m = size(u,2);
nii = size(alphacand,1);

J = zeros(nii,1);
aux = zeros(size(Omega));  




for ii = 1:nii	% For each guess, compute the performance index J
				% parallel computing might be explored
	
	ys = dltisim(A{ii},B{ii},C{ii},zeros(p,m),u);
	vafys = zeros(p,1);
	% Employ function handle!
	for j=1:p, vafys(j) = 100 - vaf(y(:,j),ys(:,j)); end

	J(ii) = sum(vafys)/p;
end

mu = 40/min(J);

for ii = 1:nii, aux(ii,:) = Omega(ii,:)*exp(-mu*J(ii)); end

% Barycenter
opteta = sum(aux,1)./sum(exp(-mu*J));	


%% --------------------------------------------------------------------------------------
function pushedpoly = push_poly_roots(pol)
    Q = cplxpair(roots(pol));
    Q(1) = abs(Q(1))^2;
    Q(2) = real(Q(2));
    pushedpoly = poly(Q);

%% --------------------------------------------------------------------------------------
function [diagM_T, C_D, A_D, filtercoeff] = comp_MOLI_arrays(alpha,m,l)

n_ = max(l);
p = length(l);
n = sum(l);
m_plus_p_times_n_ = (m+p)*n_;

filtercoeff = fliplr(-alpha(2:end));
A_ = [ [zeros(1,n_-1);eye(n_-1)], filtercoeff'];

A_D = zeros(n);
C_D = zeros(p,n);
M = cell(p,1);
diagM_T = cell(p,1);
ct = 0;

for i = 1:p
    
    if(l(i) ~= n_)
        M{i} = zeros(n_,l(i));
       
        if(~exist('alpha_roots','var')), alpha_roots = cplxpair(roots(alpha)); end;
		
		R = rem(l(i),2);	%R = 1 if the jth observability index is odd
		
		if(R && (~isreal(alpha_roots(end))))
			disp('Unable to build a design model using the specified alpha polynomial.');
			error('For such choice of l, alpha must have at least one real root.');
		end       
        
		alphai_roots = alpha_roots(1:l(i));
		if(R), alphai_roots(end) = alpha_roots(end); end
		
        alphai = poly(alphai_roots);
        
        alpha_div_alphai = deconv(alpha,alphai);
		for ii = 1:l(i)
            M{i}(ii:(ii+length(alpha_div_alphai)-1),ii) = flipud(alpha_div_alphai');
		end
		
		% Parameter estimation adjust matrix
		diagM_T{i} = zeros((m+p)*l(i),m_plus_p_times_n_);
		for r = 1:(m+p), diagM_T{i}(1+(r-1)*l(i):r*l(i),1+(r-1)*n_:r*n_) = M{i}'; end;
		
        A_D(ct+1:ct+l(i),ct+1:ct+l(i)) = ...
            [[zeros(1,l(i)-1); eye(l(i)-1)], flipud(-alphai(2:end)')];
        
	else	% ni = n_
        A_D(ct+1:ct+l(i),ct+1:ct+l(i)) = A_;
    end 
    
    ct = ct + l(i);
    C_D(i,ct) = 1;
end