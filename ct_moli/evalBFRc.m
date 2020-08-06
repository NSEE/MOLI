%

function BFRc = evalBFRc(dtv,linear_model)
 
% function BFRc = evalBFRc(dtv,model)
%
%    Compute model performance using a criterion based on the
% Best Fit Rate (BFR), for p = 2 it is equivalent to the RMSE
% (Root Mean Squared Error).
%

N = size(dtv.u,1);
yhat = lsim(linear_model,dtv.u,0:dtv.Ts:(N-1)*dtv.Ts);

p = 2;
BFRc = min(norm(dtv.y - yhat,p)/...
	norm(dtv.y - kron(ones(N,1),mean(dtv.y)),p),1);
