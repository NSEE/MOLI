function [a,b,c,P] = rmoli(y, u, l, alpha)
% ---------------------------------------------------------------------------------------
% Function to recursively estimate state-space multivariable model parameters. 
%
% Author: Rodrigo A. Romano - Mar/2015
%
% ---------------------------------------------------------------------------------------
% function theta = rmoli(y, u, l, alpha)
%
% ---------------------------------------------------------------------------------------
% Inputs:
%
% - y: output data matrix (N x p, where "p" is the number of outputs)
% - u: input data matrix (N x m, where "m" is the number of inputs)
% - l: list of observability indices
% - alpha: filter characteristic polynomial
%
% Outputs:
%
% - State space model matrices a, b and c

% ---------------------------------------------------------------------------------------
% Rev0 - Mar/2015: .

[N, p] = size(y);
[~, m] = size(u);
n_ = max(l);
n = sum(l);
m_plus_p_times_n_ = (m+p)*n_;

%% Verification of the input arguments

if(nargin < 4), error('too few input arguments...'); end

if(length(l) ~= p)
	error('The length of l must be equal to the number of outputs');
end

if(length(alpha) ~= n_+1), error('alpha order mismatch!'); end

%% Some necessary matrices

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

%% Recursive parameter estimation

n_ = length(alpha)-1;
ini = n_;
z = [y,u];
zf = zeros(size(z));
xi_k = zeros(m_plus_p_times_n_,1);

% Initial conditions

P = cell(p,N);	% Covariance matrices
R1 = cell(p,1);	% Parameter covariance matrices
R2 = cell(p,1);	% Output variance
theta = cell(p,N);	% Parameter cell array

for i = 1:p
	P{i,ini-1} = eye((m+p)*l(i) + i - 1)*10^6;
	theta{i,ini-1} = zeros((m+p)*l(i) + i - 1,1);
	R1{i} = zeros((m+p)*l(i) + i - 1); %eye((m+p)*l(i) + i - 1);
	R2{i} = 1e-9;
end
	
for k = ini:N
	for ii=1:(p+m)	% Identifier recursive state computation
		xi_k(1+(ii-1)*n_:ii*n_) = zf(k-n_+1:k,ii);
		zf(k+1,ii) = filtercoeff*zf(k-n_+1:k,ii) + z(k,ii);
	end
	
	for i = 1:p
		if(l(i) == n_)
			phi_i_k = [xi_k;y(k,1:i-1)];
		else
			phi_i_k = [diagM_T{i}*xi_k;y(k,1:i-1)];
		end
		
		aux = phi_i_k'*P{i,k-1}*phi_i_k + R2{i};
		K_k = (P{i,k-1}*phi_i_k)/aux;
		theta{i,k} = theta{i,k-1} + K_k*(y(k,i)-phi_i_k'*theta{i,k-1});
		P{i,k} = P{i,k-1} - (P{i,k-1}*(phi_i_k*phi_i_k')*P{i,k-1})/aux + R1{i};
	end
end

%% Computation of the state space model matrices

D = zeros(n,p);
B = zeros(n,m);
G = zeros(p);

for jp = 1:p
    ct = 0;
    for i = 1:p
        D(ct+1:ct+l(i),jp) = theta{i,k}(1 + l(i)*(jp-1):l(i)*jp);
        ct = ct + l(i);
    end
    
    if(jp > 1), G(jp,1:jp-1) = theta{jp,k}((m+p)*l(jp)+1:end); end
end

for jm = 1:m
    ct = 0;
    for i = 1:p
        B(ct+1:ct+l(i),jm) = theta{i,k}(1 + l(i)*p + l(i)*(jm-1):l(i)*(p+jm));
        ct = ct + l(i);
    end
end

c = (eye(p) - G)\C_D;%((eye(p) - G)^-1)*C_D;
b = B;
a = A_D + D*c;



