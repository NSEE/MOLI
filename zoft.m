%% --------------------------------------------------------------------------------------
function opteta = zoft(varargin)

% Zero-order oracle (derivative-free) filter tuning
% opteta = zoft(y,u,l,T,wczetacand)
% or
% opteta = zoft(y,u,l,T,wc,zeta,sfactor)
%
%
%%

y = varargin{1};
u = varargin{2};
l = varargin{3};
T = varargin{4};

[~, m] = size(u);
p = length(l);

if(nargin == 5), [alphacand, wczetacand] = genalphapoly(varargin{5},l,T);
else
	[alphacand, wczetacand] = genalphapoly(varargin{5},varargin{6},varargin{7},l,T);
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
	ys = lsim(ss(a,b,c,zeros(p,m),T),u);
	
	vafys = zeros(p,1);
	for j=1:p, vafys(j) = 100 - vaf(y(:,j),ys(:,j)); end

	J(ii) = sum(vafys)/p;
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