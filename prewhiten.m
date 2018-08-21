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
% 'process': (default 'ar') 'ar' or 'pink'; type of spectrum to fit to
%   background spectrum of time series
% 'ncoeff': (default 4) number of coefficients to fit either to
%   AR(p) process or to pink noise process
% 'nw': (default 2) time half bandwidth product for computing spectrum from
%   which the lag coefficients are estimated
% 'a': user-specified lag-coefficients
%
% OUT:
% ws: whitened time series
% a: lag coefficients used for prewhitening
%
% TO DO:
% - could be generalized to prewhitening by any sort of autoregressive
%   process; i.e. we could force prewhitening by a pink noise process.
%   Might consider making particular functions or need to think of a
%   general way to implement prewhitening by restricting the set of
%   coefficients a(1)...a(p)
%
% Adrian Tasistro-Hart, adrianraph-at-gmail.com, 21.08.2018

function [ws,a] = prewhiten(ts,varargin)

% parse inputs
parser = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addRequired(parser,'ts',@isnumeric)
addParameter(parser,'process','ar',@ischar)
addParameter(parser,'ncoeff',4,validScalarPosNum);
addParameter(parser,'nw',2,validScalarPosNum);
addParameter(parser,'a',[],@isnumeric);

parse(parser,ts,varargin{:})

ts = parser.Results.ts;
process = parser.Results.process;
ncoeff = parser.Results.ncoeff;
nw = parser.Results.nw;
a = parser.Results.a;

% validate process
process = validatestring(process,{'ar','pink'});

% number of samples
n = length(ts);
% make sure have enough samples
assert(n>ncoeff,'need enough samples for given order')
% make column
ts = ts(:);

% fit AR(p) process to time series
% ab = arburg(ts,p);
% a = -ab(2:end);
% % make column
% a = a(:);

% only estimate coefficients if user has not specified them
if isempty(a)
    % estimate spectrum from time series (not going to detrend, assume user
    % has done so already)
    [pxx,w] = pmtm(ts,nw,[],1);

    switch process
        case 'ar'
            a = ARfit(ncoeff,w,pxx,1/2);
        case 'pink'
            A = pinkfit(w,pxx);
            a = pinkcoeff(A,ncoeff,'filter','ar');
            % convert from matlab convention to Percival Walden convention
            a = -a(2:end);
    end

else
    ncoeff = length(a);
end

% with AR(p) parameters, subtract weighted observations from ts to generate
% ws
ws = zeros(n,1);
for ii = n:-1:ncoeff+1
    ws(ii) = ts(ii) - sum(flipud(a).*ts(ii-ncoeff:ii-1));
end
% truncate ws
ws = ws(ncoeff+1:end);


end