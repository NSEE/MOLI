function [a,b,c,epsilon] = moli(y, u, l, alpha)
% ---------------------------------------------------------------------------------------
% Function to estimate the state-space multivariable model parameters. 
%
%
% ---------------------------------------------------------------------------------------
% function [a,b,c] = moli(y, u, l, alpha)
%
% ---------------------------------------------------------------------------------------
% Inputs:
%
% - y: output data matrix (N x p, where "p" is the number of outputs)
% - u: input data matrix (N x m, where "m" is the number of inputs)
% - l: list of observability indices
% - alpha: identifier characteristic polynomial
%
% Outputs:
%
% - State space model matrices a, b and c
% - epsilon: model residuals

% ---------------------------------------------------------------------------------------
% Author: Rodrigo A. Romano
%
% Mai.2013 - Correction in the computation of alpha_i polynomial roots(alphai_roots).
% Mar.2014 - renamed to moli.
% Out.2014 - Tikhonov regularization & Direct filter.
% Out.2015 - Estimation based on the range [n_,N] and other minor changes.

%%
[N, p] = size(y);
[~, m] = size(u);
n_ = max(l);
n = sum(l);
m_plus_p_times_n_ = (m+p)*n_;

%% Verification of the input arguments

if(nargin < 4), error('too few input arguments...'); end

if(length(l) ~= p), error('The length of l must be equal to the number of outputs'); end

if(length(alpha) ~= n_+1), error('alpha order mismatch!'); end

%% Some necessary matrices

A_ = [ [zeros(1,n_-1);eye(n_-1)], flipud(-alpha(2:end)')];

A_I = zeros(m_plus_p_times_n_);
B_I = zeros(m_plus_p_times_n_,m);
D_I = zeros(m_plus_p_times_n_,p);

for j = 1:(m+p)
    A_I(1+(j-1)*n_:j*n_,1+(j-1)*n_:j*n_) = A_';
    if (j <= p)
        D_I(j*n_,j) = 1;
    else
        B_I(j*n_,j-p) = 1;
    end
end

A_D = zeros(n);
C_D = zeros(p,n);
M = cell(p,1);
ct = 0;

for j = 1:p
    
    if(l(j) ~= n_)
        M{j} = zeros(n_,l(j));
       
        if(~exist('alpha_roots','var')), alpha_roots = cplxpair(roots(alpha)); end;
		
		R = rem(l(j),2);	%R = 1 if the jth observability index is odd
		
		if(R && (~isreal(alpha_roots(end))))
			disp('Unable to build a design model using the specified alpha polynomial.');
			error('For such choice of l, alpha must have at least one real root.');
		end       
        
		alphai_roots = alpha_roots(1:l(j));
		if(R), alphai_roots(end) = alpha_roots(end); end
		
        alphai = poly(alphai_roots);
        
        alpha_div_alphai = deconv(alpha,alphai);
        for ii = 1:l(j)
            M{j}(ii:(ii+length(alpha_div_alphai)-1),ii) = flipud(alpha_div_alphai');
        end
        
        A_D(ct+1:ct+l(j),ct+1:ct+l(j)) = ...
            [[zeros(1,l(j)-1); eye(l(j)-1)], flipud(-alphai(2:end)')];
        
	else	% ni = n_
        A_D(ct+1:ct+l(j),ct+1:ct+l(j)) = A_;
    end 
    
    ct = ct + l(j);
    C_D(j,ct) = 1;
end

%% Identifier state computation

% disp('Direct Filter')

% Phi = zeros(N-n_,m_plus_p_times_n_);
% Zf = filter([0, 1],alpha,[y,u]);
% jPhi = n_:n_:(m+p)*n_;
% for j = 0:n_-1
% 	Phi(1+j:end,jPhi-j) = Zf(1+n_:end-j,:);
% end
% 
% y = y(1+n_:end,:);

Phi = zeros(N-n_,m_plus_p_times_n_);
Zf = filter([0, 1],alpha,[y,u]);
jPhi = n_:n_:(m+p)*n_;
for j = 0:n_-1
	Phi(1:end,jPhi-j) = Zf(1-j+n_:end-j,:);
end

y = y(1+n_:end,:);

%% Parameter estimation
% sqrt_lambda = 1e-3;% .1 (TWR) %1 (IEEETrans)
theta = cell(p,1);  % Parameter matrix

if(nargout >= 4), epsilon = zeros(size(y)); end

for j = 1:p
	if(j == 1)
		if(l(j) == n_)
			theta{j} = Phi \ y(:,j);
			if(nargout >= 4), epsilon(:,j) = y(:,j) - Phi*theta{j}; end
		else
			diagM = [];     % Parameter estimation adjust matrix
			for r = 1:(m+p), diagM = blkdiag(diagM,M{j}); end;
			
			theta{j} = (Phi*diagM) \ y(:,j);
			if(nargout >= 4), epsilon(:,j) = y(:,j) - (Phi*diagM)*theta{j}; end
		end
	else
		if(l(j) == n_)
			
 			theta{j} = [Phi, y(:,1:j-1)] \ y(:,j);
			
% 			theta{j} = [[Phi, y(:,1:j-1)]; sqrt_lambda*eye(l(j)*(m+p)+j-1)] \...
% 				[y(:,j); zeros(l(j)*(m+p)+j-1,1)];
			
			if(nargout >= 4)
				epsilon(:,j) = y(:,j) - [Phi, y(:,1:j-1)]*theta{j};
			end
			
		else
			diagM = [];     % Parameter estimation adjust matrix
			for r = 1:(m+p), diagM = blkdiag(diagM,M{j}); end;
			
 			theta{j} = [Phi*diagM, y(:,1:j-1)] \ y(:,j);
			
% 			theta{j} = [[Phi*diagM, y(:,1:j-1)]; sqrt_lambda*eye(l(j)*(m+p)+j-1)] \...
% 				[y(:,j); zeros(l(j)*(m+p)+j-1,1)];

			if(nargout >= 4)
				epsilon(:,j) = y(:,j) - [Phi*diagM, y(:,1:j-1)]*theta{j};
			end
		end
	end
		
end


%% Computation of the state space model matrices

D = zeros(n,p);
B = zeros(n,m);
G = zeros(p);

for jp = 1:p
    ct = 0;
    for j = 1:p
        D(ct+1:ct+l(j),jp) = theta{j}(1 + l(j)*(jp-1):l(j)*jp);
        ct = ct + l(j);
    end
    
    if(jp > 1), G(jp,1:jp-1) = theta{jp}((m+p)*l(jp)+1:end); end
end

for jm = 1:m
    ct = 0;
    for j = 1:p
        B(ct+1:ct+l(j),jm) = theta{j}(1 + l(j)*p + l(j)*(jm-1):l(j)*(p+jm));
        ct = ct + l(j);
    end
end

c = (eye(p) - G)\C_D;
b = B;
a = A_D + D*c;




