%% --------------------------------------------------------------------------------------
function [alphacand, wczetacand] = genctalpha(varargin)

%%
if(nargin == 4)	% Generate alphacand and wczetacand given wc,zeta and l
	wc = varargin{1};
	zeta = varargin{2};
	sfactor = varargin{3};
	l = varargin{4};
	
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
				
				% Dominating poles
				if(zeta(jj) < 1)	% Complex poles
					s1 = wc(ii)*(-zeta(jj)+1i*sqrt(1-zeta(jj)^2));
					alphacand(ct,end-2:end) = poly([s1,conj(s1)]);
				else
					if(sum(l) == length(l))	% Single real pole
						alphacand(ct,end-1:end) = [1 wc(ii)];
					else					% Overdamped poles - exception
						Delta = sqrt(zeta(jj)^2-1);
						s1 = wc(ii)*(-zeta(jj) + Delta);
						s2 = wc(ii)*(-zeta(jj) - Delta);
						alphacand(ct,end-2:end) = poly([s1,s2]);
					end
				end
				
				% Other poles to complete required filter order
				jalpha = 1;
				wcaux = wc(ii);
				while (jalpha <= max(l)-2)
					wcaux = wcaux*sfactor(kk);
					
					if(max(l)-2-jalpha >= 1)
						
						if(zeta(jj) < 1)	% Complex poles
							s = wcaux*(-zeta(jj)+1i*sqrt(1-zeta(jj)^2));
							alphacand(ct,end-3-jalpha:end) = ...
								conv(alphacand(ct,end-jalpha-1:end), poly([s,conj(s)]));
						else
							% Overdamped poles - exception
							Delta = sqrt(zeta(jj)^2-1);
							s1 = wcaux*(-zeta(jj) + Delta);
							s2 = wcaux*(-zeta(jj) - Delta);
							alphacand(ct,end-3-jalpha:end) = ...
								conv(alphacand(ct,end-jalpha-1:end), poly([s1,s2]));
						end
						jalpha = jalpha + 2;
					
					else	% alpha is odd
						alphacand(ct,end-2-jalpha:end) = ...
							conv(alphacand(ct,end-jalpha-1:end),...
							[1 wcaux]);
						jalpha = jalpha + 1;
					end
				end
				
				wczetacand(ct,:) = [wc(ii),zeta(jj),sfactor(kk)];
			end
		end
	end
	
	%% Plot the roots of the alpha candidates
% 	if((nzeta ~= 1) && (sum(l) >= 4))	%<=== Debug plot
% 		for jj=1:nzeta
% 			ialphacand = jj:nzeta:nzeta*nwc';
% 			for iaux = 1:nwc, alphar(iaux,:) = roots(alphacand(ialphacand(iaux),:))'; end;
% 			
% 			plot(real(alphar(:,1)),imag(alphar(:,1)),'s-.','color',[.5 .5 .5]); hold on;
% 			plot(real(alphar(:,2)),imag(alphar(:,2)),'s-.','color',[.5 .5 .5]);
% 		end
% 	end

elseif(nargin == 2)	% Generate alphacand given wczetacand and l
	etacand = varargin{1};
	l = varargin{2};
	
	neta = length(etacand(:,1));
	alphacand = zeros(neta,1 + max(l));
	
	for ii = 1:neta
		% Dominating poles
		if(etacand(ii,2) < 1)	% Complex poles
			s1 = etacand(ii,1)*(-etacand(ii,2)+1i*sqrt(1-etacand(ii,2)^2));
			alphacand(ii,end-2:end) = poly([s1,conj(s1)]);
		else
			if(sum(l) == length(l))	% Single real pole
				alphacand(ii,end-1:end) = [1 etacand(ii,1)];
			else		% Overdamped poles - exception
				Delta = sqrt(etacand(ii,2)^2-1);
				s1 = etacand(ii,1)*(-etacand(ii,2) + Delta);
				s2 = etacand(ii,1)*(-etacand(ii,2) - Delta);
				alphacand(ii,end-2:end) = poly([s1,s2]);
			end
		end
		
		% Other poles to complete required filter order
		jalpha = 1;
		wcaux = etacand(ii,1);
		while (jalpha <= max(l)-2)
			wcaux = wcaux*etacand(ii,3);
			
			if(max(l)-2-jalpha >= 1)
				
				if(etacand(ii,2) < 1)	% Complex poles
					s = wcaux*(-etacand(ii,2)+1i*sqrt(1-etacand(ii,2)^2));
					alphacand(ii,end-3-jalpha:end) = ...
						conv(alphacand(ii,end-jalpha-1:end), poly([s,conj(s)]));
				else
					% Overdamped poles - exception
					Delta = sqrt(etacand(ii,2)^2-1);
					s1 = wcaux*(-etacand(ii,2) + Delta);
					s2 = wcaux*(-etacand(ii,2) - Delta);
					alphacand(ii,end-3-jalpha:end) = ...
						conv(alphacand(ii,end-jalpha-1:end), poly([s1,s2]));
				end	
				
				jalpha = jalpha + 2;
			else		% alpha is odd
				alphacand(ii,end-2-jalpha:end) = ...
					conv(alphacand(ii,end-jalpha-1:end),[1 wcaux]);
				jalpha = jalpha + 1;
			end
		end
		
	end
	wczetacand = etacand;

else
	
	error('Wrong usage of function genalphapoly - contact the authors...');
end