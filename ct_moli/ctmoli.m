function [A,B,C,epsilon] = ctmoli(y, u, Ts, l, alpha, zeroBel)
% -----------------------------------------------------------------------------------
% Function to estimate the prameters of continuous-time state-space multivariable 
% model parameters. The user defined matrix A is assumed to be composed of companion
% matrices and C is block diagonal constructed from vectors of the form [0 ... 0 1]. 
%
% -----------------------------------------------------------------------------------
% function [A,B,C,epsilon] = ctmoli(y, u, Ts, l, alpha, zeroBel)
%
% -----------------------------------------------------------------------------------
% Inputs:
%
% - y: output data matrix (N x ny, where "ny" is the number of outputs)
% - u: input data matrix (N x nu, where "nu" is the number of inputs)
% - Ts: sampling period
% - l: list of observability indices
% - alpha: Morse observer characteristic polynomial
% - zeroBel: defines which coefficients of B are zero (non-zero entries set to "0")
% 
%
% Outputs:
%
% - Model matrices A, B and C of the state-space continuous time model
%
%	dx/dt = Ax(t) + Bu(t)
%	y(t) = Cx(t)
%
% - epsilon: model residuals

% -----------------------------------------------------------------------------------
%
% Author:	Rodrigo A. Romano
%

%% Verification of the input arguments

if(nargin < 5), error('Too few input arguments...'); end

[N, ny] = size(y);
if(length(l) ~= ny)
	error('The length of l must be equal to the number of outputs');
end

[~, nu] = size(u);
n_ = max(l);
if(length(alpha) ~= n_+1), error('alpha order mismatch!'); end

% Test if alpha is monic
if(alpha(1) ~= 1), alpha = alpha./alpha(1); end

nx = sum(l);
n_I = (nu+ny)*n_;

restrictB = 0;
if(nargin > 5)
	if(size(zeroBel,1) == nx && size(zeroBel,2) == nu), restrictB = 1;
    else, disp(['Matrix zeroBcoeffs must be nx x nu',...
		'No restrictions will be imposed on the coefficients of matrix B']);	
	end
end

%% Some necessary matrices

A_ = [ [zeros(1,n_-1);eye(n_-1)], flipud(-alpha(2:end)')];


% ctrans_times_c = zeros(n_);
% ctrans_times_c(end,end) = 1;
% 
% Wo = dlyap(A_', ctrans_times_c);
% try
% 	T = inv(chol(Wo,'lower'))'; % REVER T
% 	
% 	C_ = [zeros(1,n_-1), 1]*T;
% 	A_ = T\A_*T;
% 	
% 	C_t = C_';
% 	A_t = A_';
% 
% catch
% 	T = eye(n_);
% 	C_ = [zeros(1,n_-1), 1];
% 	
% 	C_t = C_';
% 	A_t = A_';
% 	disp('xxxx')
% end


% Observer state equation matrices

A_I = zeros(n_I);
B_I = zeros(n_I,nu);
D_I = zeros(n_I,ny);

for j = 1:(nu+ny)
	
% 	A_I(1+(j-1)*n_:j*n_,1+(j-1)*n_:j*n_) = A_t;
% 	
%     if (j <= ny), D_I(1+(j-1)*n_:j*n_,j) = C_t;
% 	else B_I(1+(j-1)*n_:j*n_,j-ny) = C_t;
%     end
	
    A_I(1+(j-1)*n_:j*n_,1+(j-1)*n_:j*n_) = A_';
    if (j <= ny)
        D_I(j*n_,j) = 1;%alpha(end);
    else
        B_I(j*n_,j-ny) = 1;%alpha(end);
    end
end

%% Identifier state computation

Xi = lsim(ss(A_I,[D_I, B_I],eye(n_I),zeros(n_I,ny+nu)),...
	[y,u],0:Ts:(N-1)*Ts,[],'zoh');

% Remove data dominated by initial conditions
y = y(1+n_:end,:);
Xi = Xi(1+n_:end,:);

%% Computation of Mi, A_D, C_D and selCol (if any entry of B shall be set to zero)

A_D = zeros(nx);
C_D = zeros(ny,nx);
M = cell(ny,1);
if(restrictB), selcolXi = cell(ny,1); end
ct = 0;

for j = 1:ny
    
	if(l(j) ~= n_)
        M{j} = zeros(n_,l(j));
       
        if(~exist('alpha_roots','var')), alpha_roots = cplxpair(roots(alpha)); end
		
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
    
	
	if(restrictB)
		selcolXi{j} = [(1:l(j)*ny)';...
			l(j)*ny + find(reshape(zeroBel(ct+1:ct+l(j),:), l(j)*nu, 1) == 0)];
	end
	
	
    ct = ct + l(j);
    C_D(j,ct) = 1;%alpha(end)
end


%% Parameter estimation
theta = cell(ny,1);  % Parameter matrix
if(nargout >= 4), epsilon = zeros(size(y)); end

llssolver = 'default';%'leastn';%

switch llssolver
	case 'Tik'
		%disp('Tikhonov reg');
		sqrt_lambda = 1e-1;
        for j = 1:ny
            if(j == 1)
                if(l(j) == n_), theta{j} = [Xi, y(:,1:j-1)] \ y(:,j);
                else
                    diagM = [];     % Parameter estimation adjust matrix
                    for r = 1:(nu+ny), diagM = blkdiag(diagM,M{j}); end
                    
                    theta{j} = [Xi*diagM, y(:,1:j-1)] \ y(:,j);
                end
            else
                if(l(j) == n_)
                    theta{j} = [[Xi, y(:,1:j-1)];...
                        sqrt_lambda*eye(l(j)*(nu+ny)+j-1)] \...
                        [y(:,j); zeros(l(j)*(nu+ny)+j-1,1)];
                    
                else
                    diagM = [];     % Parameter estimation adjust matrix
                    for r = 1:(nu+ny), diagM = blkdiag(diagM,M{j}); end
                    
                    theta{j} = [[Xi*diagM, y(:,1:j-1)];...
                        sqrt_lambda*eye(l(j)*(nu+ny)+j-1)] \...
                        [y(:,j); zeros(l(j)*(nu+ny)+j-1,1)];
                end
            end
            
        end
		
	case 'leastn'
		%disp('Least norm solution');
        for j = 1:ny
            if(j == 1)
                if(l(j) == n_)
                    theta{j} = [Xi, y(:,1:j-1)] \ y(:,j);
                else
                    diagM = [];     % Parameter estimation adjust matrix
                    for r = 1:(nu+ny), diagM = blkdiag(diagM,M{j}); end
                    
                    theta{j} = [Xi*diagM, y(:,1:j-1)] \ y(:,j);
                end
            else
                if(l(j) == n_)
                    [Q,R] = qr([Xi, y(:,1:j-1)]);
                    theta{j} = R\(Q'*y(:,j));
                else
                    diagM = [];     % Parameter estimation adjust matrix
                    for r = 1:(nu+ny), diagM = blkdiag(diagM,M{j}); end
                    
                    [Q,R] = qr([Xi*diagM, y(:,1:j-1)]);
                    theta{j} = R\(Q'\y(:,j));
                end
            end
            
        end
		
	case 'default'
		% default
		if(restrictB) % Test --- if entries of B shall be set to zero
			for j = 1:ny
				theta{j} = zeros(l(j)*(ny+nu) + j-1,1);
				if(j == 1)
					if(l(j) == n_)
						theta{j}([selcolXi{j}]) = Xi(:,selcolXi{j})\y(:,j);
						if(nargout >= 4), epsilon(:,j) = y(:,j) - Xi*theta{j}; end
					else
						diagM = [];     % Parameter estimation adjust matrix
						for r = 1:(nu+ny), diagM = blkdiag(diagM,M{j}); end
						XidiagMi = Xi*diagM;
						theta{j}([selcolXi{j}]) = XidiagMi(:,selcolXi{j})\y(:,j);
						if(nargout >= 4)
							epsilon(:,j) = y(:,j) - (Xi*diagM)*theta{j};
						end
					end
				else
					if(l(j) == n_)				
						theta{j}([selcolXi{j};l(j)*(ny+nu) + (1:j-1)']) = ...
							[Xi(:,selcolXi{j}), y(:,1:j-1)]\y(:,j);
						if(nargout >= 4)
							epsilon(:,j) = y(:,j) - [Xi, y(:,1:j-1)]*theta{j};
						end
					else
						diagM = [];     % Parameter estimation adjust matrix
						for r = 1:(nu+ny), diagM = blkdiag(diagM,M{j}); end
						XidiagMi = Xi*diagM;
						theta{j}([selcolXi{j};l(j)*(ny+nu) + (1:j-1)']) = ...
							[XidiagMi(:,selcolXi{j}), y(:,1:j-1)]\y(:,j);
						if(nargout >= 4)
							epsilon(:,j) = y(:,j) - [Xi*diagM, y(:,1:j-1)]*theta{j};
						end
					end
				end
				
			end
		else
			for j = 1:ny
				if(j == 1)
					if(l(j) == n_)
						theta{j} = Xi \ y(:,j);
						if(nargout >= 4), epsilon(:,j) = y(:,j) - Xi*theta{j}; end
					else
						diagM = [];     % Parameter estimation adjust matrix
						for r = 1:(nu+ny), diagM = blkdiag(diagM,M{j}); end
						
						theta{j} = Xi*diagM \ y(:,j);
						if(nargout >= 4)
							epsilon(:,j) = y(:,j) - (Xi*diagM)*theta{j};
						end
					end
				else
					if(l(j) == n_)	
						theta{j} = [Xi, y(:,1:j-1)]\y(:,j);
						if(nargout >= 4)
							epsilon(:,j) = y(:,j) - [Xi, y(:,1:j-1)]*theta{j};
						end
						
					else
						diagM = [];     % Parameter estimation adjust matrix
						for r = 1:(nu+ny), diagM = blkdiag(diagM,M{j}); end
						
						theta{j} = [Xi*diagM, y(:,1:j-1)]\y(:,j);
						
						if(nargout >= 4)
							epsilon(:,j) = y(:,j) - [Xi*diagM, y(:,1:j-1)]*theta{j};
						end
						
					end
				end
				
			end
		end
end



%% Computation of the state space model matrices A, B and C

D = zeros(nx,ny);
B = zeros(nx,nu);
G = zeros(ny);

for jp = 1:ny
    ct = 0;
    for j = 1:ny
        D(ct+1:ct+l(j),jp) = theta{j}(1 + l(j)*(jp-1):l(j)*jp);
        ct = ct + l(j);
    end
    
    if(jp > 1), G(jp,1:jp-1) = theta{jp}((nu+ny)*l(jp)+1:end); end
end

for jm = 1:nu
    ct = 0;
    for j = 1:ny
        B(ct+1:ct+l(j),jm) =...
            theta{j}(1 + l(j)*ny + l(j)*(jm-1):l(j)*(ny+jm));
        ct = ct + l(j);
    end
	
end

C = (eye(ny) - G)\C_D;
A = A_D + D*C;


