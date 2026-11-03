function R = grangerLinearF(src, tgt, lag)
% grangerLinearF  Classical linear Granger test src to tgt (F test on least-squares fits).
%
% CALLS     None (uses fcdf from the Statistics and Machine Learning Toolbox)
% CALLED BY runCase1a
%
% INPUT
%   src, tgt   numeric vectors of equal length
%   lag        number of lags, scalar
% OUTPUT
%   R   struct with F, p, lag
src = src(:); tgt = tgt(:); N = numel(tgt);
y  = tgt(lag+1:N);
X0 = ones(numel(y),1);
for j = 1:lag, X0 = [X0, tgt(lag+1-j:N-j)]; end
X1 = X0;
for j = 1:lag, X1 = [X1, src(lag+1-j:N-j)]; end
r0 = y - X0*(X0\y); r1 = y - X1*(X1\y);
df2 = numel(y) - size(X1,2);
R.F = ((r0'*r0 - r1'*r1)/lag) / ((r1'*r1)/df2);
R.p = 1 - fcdf(R.F, lag, df2);
R.lag = lag;
end
