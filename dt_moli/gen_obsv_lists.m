function [list_set, checkoblist] = gen_obsv_lists(p, nmax)
% Generate all possible lists of observability indexes of linear 
% dynamic models up to order nmax
%
% Sintaxe:
% [list_set, checkoblist] = gen_obsv_lists(p, nmax)
%
% list_set: cell with the set of lists. The lists of different order
% are stored in separate cell. Each list is stored in a differentrow.
%
% checkoblist: a matrix with the amount of lists (second row) of each
% order.


% Author: Rodrigo A. Romano - Apr, 2014
%

vecn = p:nmax;
list_set = cell(length(vecn),1);
ind_list = zeros(length(vecn),1);

% amountofmodels = 0;
checkoblist = [vecn; zeros(1,length(vecn))];

for r = 1:nmax-p+1
	for s = 1:nmax-p+1
		for v = 1:nmax-p+1
			
			l = [r,s,v];
			n = sum(l);
			if(n > nmax), continue; end;
			
			ind_list(n-p+1) = ind_list(n-p+1) + 1;		
			list_set{n-p+1,1}(ind_list(n-p+1),1:p) = l;
			
			checkoblist(2,n-p+1) = checkoblist(2,n-p+1)+1;
			if(exist('amountofmodels','var'))
				amountofmodels = amountofmodels + 1;
				disp(int2str(amountofmodels));
			end
		end
	end
end

end

