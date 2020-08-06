% function [A,B,C,alpha] = butterzoft(y,u,Ts,l,wc,yv)
% Derivative-free filter tuning - Butterworth filter prototype
%
%
%% --------------------------------------------------------------------------------------
function [A,B,C,alpha] = butterzoft(y,u,Ts,l,wc,yv)


[~, nu] = size(u);
[N, ny] = size(y);

if(nargin < 5), Ts = 1; end

nii = length(wc);

J = zeros(nii,1);
aux = size(J); 
bfrys = zeros(1,ny);
n_ = max(l);

norm_contfiltpoles = exp(1i*(pi*(1:2:2*n_-1)/(2*n_)+pi/2)).';
alphacand = zeros(nii,n_+1);

%%
for ii = 1:nii

	alphacand(ii,:) = real(poly(norm_contfiltpoles.*wc(ii)));
	
	if((mod(max(l),2) == 0) && any(mod(l,2)))
		if(~isreal(roots(alphacand(ii,:))))
			alphacand(ii,:) = fixrealp(alphacand(ii,:));
		end
	end
	
	[a,b,c] = ctmoli(y, u, Ts, l, alphacand(ii,:));
	ys = lsim(ss(a,b,c,zeros(ny,nu)),u,0:Ts:(N-1)*Ts);


% 	J(ii) = min(norm(yv-ys,2)/...
%         norm(yv-kron(mean(yv),ones(length(yv),1)),2),1);

	
	for iny=1:ny
		bfrys(iny) = min(norm(yv(:,iny)-ys(:,iny))/...
         norm(yv(:,iny)-mean(yv(:,iny))),1);
	end
	J(ii) = sum(bfrys)/ny;
end
	
	mu = min(8/min(J));

for ii = 1:nii
	aux(ii) = wc(ii)*exp(-mu*J(ii));
end

barywc = sum(aux)./sum(exp(-mu*J));

alpha = real(poly(norm_contfiltpoles*barywc));


%%

if((mod(max(l),2) == 0) && any(mod(l,2)))
    if(~isreal(roots(alpha)))
        alpha = fixrealp(alpha);
    end
end
[A,B,C] = ctmoli(y, u, Ts, l, alpha);



%% Plot barycenter and curiosity points

plotbary = 0;

if(plotbary)
	figure(3);
	plot(wc,J,'s--','color',[.4 .4 .4],'Linewidth',1.5,...
		'Markersize',12); grid; hold on;
	set(gcf,'Position',[1031 276 540 400]);
	set(gca,'Position',[0.11 0.1 .87 .87],'Fontsize',13);

	ys = lsim(ss(A,B,C,zeros(ny,nu)),u,0:Ts:(N-1)*Ts);
	
% 	baryJ = min(norm(yv-ys,2)/...
% 		norm(yv-kron(mean(yv),ones(length(yv),1)),2),1);
	
	for iny=1:ny
		bfrys(iny) = min(norm(yv(:,iny)-ys(:,iny))/...
         norm(yv(:,iny)-mean(yv(:,iny))),1);
	end
	
	baryJ = sum(bfrys)/ny;
	
	semilogx(barywc,baryJ,'rx','Linewidth',3,'Markersize',16);
	legend('Evaluated frequencies','Barycenter');
	
	J(ii+1) = baryJ;
	
	ylabel('Index of merit \itJ'); xlabel('\omega_c (rad/s)')
	axis([0.1 1.05*wc(end) 0.99*min(J) 1.005*max(J)])
	hold off;

end

	
end

%% --------------------------------------------------------------------------------------
function fixedp = fixrealp(pol)
    Q = cplxpair(roots(pol));
    Q(1) = abs(Q(1))^2;
    Q(2) = real(Q(2));
    fixedp = poly(Q);
end