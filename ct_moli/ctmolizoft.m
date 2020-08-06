
%
% --------------------------------------------------------------------------------------
function [A,B,C,alpha] = ctmolizoft(dte,l,wcvec,zetavec,J_handle,varargin)

% function [A,B,C,alpha] = ctmolizoft(dte,l,wcvec,zetavec,dtv,J_handle,varargin)
% Zero-order Oracle Filter Tuning - Underdamped filter prototype
%

%   - J_handle: filter tuning criterion function handle
%   + varargin: function handle input arguments

%%

%dte == y,u,Ts

nu = size(dte.u,2);
ny = size(dte.y,2);

sfactor = 2;
[alphacand, etacand] = genctalpha(wcvec,zetavec,sfactor,l);


nii = size(alphacand,1);

J = zeros(nii,1);
aux = zeros(size(etacand));


%%
for ii = 1:nii
    
    if((mod(max(l),2) == 0) && any(mod(l,2)))
        if(~isreal(roots(alphacand(ii,:))))
            alphacand(ii,:) = fixrealp(alphacand(ii,:));
        end
    end
    
    [a,b,c] = ctmoli(dte.y, dte.u, dte.Ts, l, alphacand(ii,:));
    
    J(ii) = J_handle(varargin{:},ss(a,b,c,zeros(ny,nu))); %% sum(bfrys)/ny;
    
end

mu = min(10/std(J),250);

for ii = 1:nii, aux(ii,:) = etacand(ii,:)*exp(-mu*J(ii)); end

opteta = sum(aux,1)./sum(exp(-mu*J));
alpha = genctalpha(opteta,l);

%%

if((mod(max(l),2) == 0) && any(mod(l,2)))
    if(~isreal(roots(alpha)))
        alpha = fixrealp(alpha);
    end
end
[A,B,C] = ctmoli(dte.y, dte.u, dte.Ts, l, alpha);


%% Plot barycenter and curiosity points

plotbary = 0;

if(plotbary)
    figure(7)
    subplot(121)
    [Cx,h] = contour(wcvec,zetavec,reshape(J,length(zetavec),[]));
    text_handle = clabel(Cx,h);
    set(text_handle,'BackgroundColor',[1 1 .6],...
        'Edgecolor',[.7 .7 .7])
    hold on;
    plot(etacand(:,1),etacand(:,2),'k+');
    plot(opteta(1),opteta(2),'b*');
    xlabel('Frequency \omega_c');
    ylabel('Damping ratio \zeta');
    
    [~, imin] = min(J);
    bestalpha = genctalpha(etacand(imin,:),l);
    plot(etacand(imin,1),etacand(imin,2),'ro');
    plot(real(roots(bestalpha)),imag(roots(bestalpha)),'hc');
    % 	optalpha = bestalpha;
    ylim([.9*min(zetavec),1.1*max(zetavec)]);
    hold off;
    
    subplot(122);
    plot(J,'+k'); hold on;
    plot(imin,J_handle(varargin{:},ss(A,B,C,zeros(ny,nu))),'b*');
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