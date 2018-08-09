% This function prewhitens a time series by fitting autoregressive
% coefficients describing an AR(p) process such that the process has a
% spectrum that closely matches the spectrum of the time series. As
% described in Percival and Walden (Ch 11), these coefficients are then
% used to compute a whitened time series: 
% X(t) = X(t-1)*a(1) + ... + X(t-p)*a(p) + W(t), t = 1...n 
% where X(t) is the time series we observe and W(t) is the whitenened time 
% series we can recover by applying the coefficients "in reverse", allowing
% us to recover n-p samples of W(t).
% The current implementation uses ARfit() to estimate the lag coefficients
% of the given process from the frequency domain.
% 
% IN:
% ts: input time series. must be a vector in the current implementation
% 'p': (default 4) order of the AR(p) process to fit
%
% OUT:
% ws: whitened time series
%
% TO DO:
% - could be generalized to prewhitening by any sort of autoregressive
%   process; i.e. we could force prewhitening by a pink noise process.
%   Might consider making particular functions or need to think of a
%   general way to implement prewhitening by restricting the set of
%   coefficients a(1)...a(p)
%
% Adrian Tasistro-Hart, adrianraph-at-gmail.com, 09.08.2018

function ws = prewhiten(ts,varargin)

% parse inputs
parser = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addRequired(parser,'ts',@isnumeric)
addParameter(parser,'p',4,validScalarPosNum);

parse(parser,ts,varargin{:})

ts = parser.Results.ts;
p = parser.Results.p;

% number of samples
n = length(ts);
% make sure have enough samples
assert(n>p,'need enough samples for given order')
% make column
ts = ts(:);

% fit AR(p) process to time series
% ab = arburg(ts,p);
% a = -ab(2:end);
% % make column
% a = a(:);

% fit autoregressive coefficients from spectrum
[pxx,w] = pmtm(detrend(ts),2,[],1);
a = ARfit(p,w,pxx,1/2);

% with AR(p) parameters, subtract weighted observations from ts to generate
% ws
ws = zeros(n,1);
for ii = n:-1:p+1
    ws(ii) = ts(ii) - sum(flipud(a).*ts(ii-p:ii-1));
end
% truncate ws
ws = ws(p+1:end);


end