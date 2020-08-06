%% --------------------------------------------------------------------------------------
function opteta = zoft(varargin)

% Zero-order oracle (derivative-free) filter tuning
% Usage:
% opteta = zoft(y,u,l,T,J_handle,wczetacand)
% or
% opteta = zoft(y,u,l,T,J_handle,wc,zeta,sfactor)
%
%
%%

% Output data
y = varargin{1};
% Input data
u = varargin{2};
% List of observability indices
l = varargin{3};
% Sampling period
T = varargin{4};
% Objective function handle
J_handle = varargin{5};

[~, nu] = size(u);
ny = length(l);

if(nargin == 6), [alphacand, wczetacand] = genalphapoly(varargin{6},l,T);
else
	[alphacand, wczetacand] = genalphapoly(varargin{6},varargin{7},varargin{8},l,T);
end

nii = size(alphacand,1);

J = zeros(nii,1);

% Employ function handle!
aux = zeros(size(wczetacand));  

for ii = 1:nii	% For each guess, compute the performance index J
	
	if((mod(max(l),2) == 0) && any(mod(l,2)))
		if(~isreal(roots(alphacand(ii,:))))
			alphacand(ii,:) = push_poly_roots(alphacand(ii,:));
		end
	end
	
	[a,b,c] = moli(y, u, l, alphacand(ii,:));
% 	ys = lsim(ss(a,b,c,zeros(p,m),T),u);
% 	
% 	vafys = zeros(p,1);
% 	for j=1:p, vafys(j) = 100 - vaf(y(:,j),ys(:,j)); end

    dtv.y = y; dtv.u = u;
	J(ii) = J_handle(dtv,ss(a,b,c,zeros(ny,nu),T)); %sum(vafys)/p;
end

mu = 40/min(J);

for ii = 1:nii, aux(ii,:) = wczetacand(ii,:)*exp(-mu*J(ii)); end

% Barycenter
opteta = sum(aux,1)./sum(exp(-mu*J));	

end



%% --------------------------------------------------------------------------------------
function pushedpoly = push_poly_roots(pol)
    Q = cplxpair(roots(pol));
    Q(1) = abs(Q(1))^2;
    Q(2) = real(Q(2));
    pushedpoly = poly(Q);
end