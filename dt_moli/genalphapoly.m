%% --------------------------------------------------------------------------------------
function [alphacand, wczetacand] = genalphapoly(varargin)

%%
if(nargin >= 5)	% Generate alphacand and wczetacand given wc,zeta,l and T
	wc = varargin{1};
	zeta = varargin{2};
	sfactor = varargin{3};
	l = varargin{4};
	T = varargin{5};
	
	ct = 0;
	nwc = length(wc);
	nzeta = length(zeta);
	nfactor = length(sfactor);
	alphacand = zeros(nzeta*nwc*nfactor,1 + max(l));
	wczetacand = zeros(nzeta*nwc*nfactor,3);
	
	for ii = 1:nwc
		for jj = 1:nzeta
			for kk = 1:nfactor
				ct = ct + 1;
				
				if(zeta(jj) < 1)
					z1 = exp(wc(ii)*(-zeta(jj)+1i*sqrt(1-zeta(jj)^2))*T);
					alphacand(ct,end-2:end) = poly([z1,conj(z1)]);
				else
					if(sum(l) == length(l))
						alphacand(ct,end-1:end) = [1 -exp(-wc(ii)*T)];
					else
						z1 = exp(wc(ii)*(-zeta(jj)+sqrt(zeta(jj)^2-1))*T);
						z2 = exp(wc(ii)*(-zeta(jj)-sqrt(zeta(jj)^2-1))*T);
						alphacand(ct,end-2:end) = poly([z1,z2]);
					end
				end
				
				jalpha = 1;
				while (jalpha <= max(l)-2)
					
					if(max(l)-2-jalpha >= 1)
						z = z1.^sfactor(kk);
						alphacand(ct,end-3-jalpha:end) = ...
							conv(alphacand(ct,end-jalpha-1:end), poly([z,conj(z)]));
						z1 = z;
						jalpha = jalpha + 2;
					else
						alphacand(ct,end-2-jalpha:end) = ...
							conv(alphacand(ct,end-jalpha-1:end),...
							[1 -abs(z1)^sfactor(kk)]);
						jalpha = jalpha + 1;
					end
				end
				
				wczetacand(ct,:) = [wc(ii),zeta(jj),sfactor(kk)];
			end
		end
	end
	
	%% Plot the roots of the alpha candidates
% 	if((nzeta ~= 1) && (sum(l) >= 4))	%<=== B2B Porto plot
% 		for jj=1:nzeta
% 			ialphacand = jj:nzeta:nzeta*nwc';
% 			for iaux = 1:nwc, alphar(iaux,:) = roots(alphacand(ialphacand(iaux),:))'; end;
% 			
% 			plot(real(alphar(:,1)),imag(alphar(:,1)),'s-.','color',[.5 .5 .5]); hold on;
% 			plot(real(alphar(:,2)),imag(alphar(:,2)),'s-.','color',[.5 .5 .5]);
% 		end
% 	end
	
elseif(nargin == 3)	% Generate alphacand given wczetacand,l and T
	etacand = varargin{1};
	l = varargin{2};
	T = varargin{3};
	
	neta = length(etacand(:,1));
	alphacand = zeros(neta,1 + max(l));
	
	for ii = 1:neta
		
		if(etacand(ii,2) < 1)
			z1 = exp(etacand(ii,1)*(-etacand(ii,2)+1i*sqrt(1-etacand(ii,2)^2))*T);
			alphacand(ii,end-2:end) = poly([z1,conj(z1)]);
		else
			if(sum(l) == length(l))
				alphacand(ii,end-1:end) = [1 -exp(-etacand(ii,1)*T)];
			else
				z1 = exp(etacand(ii,1)*(-etacand(ii,2)+sqrt(etacand(ii,2)^2-1))*T);
				z2 = exp(etacand(ii,1)*(-etacand(ii,2)-sqrt(etacand(ii,2)^2-1))*T);
				alphacand(ii,end-2:end) = poly([z1,z2]);
			end
		end
		
		jalpha = 1;
		while (jalpha <= max(l)-2)
			
			if(max(l)-2-jalpha >= 1)
				z = z1.^etacand(ii,3);
				alphacand(ii,end-3-jalpha:end) = ...
					conv(alphacand(ii,end-jalpha-1:end), poly([z,conj(z)]));
				z1 = z;
				jalpha = jalpha + 2;
			else
				alphacand(ii,end-2-jalpha:end) = ...
					conv(alphacand(ii,end-jalpha-1:end),...
					[1 -abs(z1)^etacand(ii,3)]);
				jalpha = jalpha + 1;
			end
		end
		
	end
	wczetacand = etacand;
else
	error('Wrong usage of function genalphapoly - contact the authors...');
end