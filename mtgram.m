% This function computes the multi-taper spectrogram of a given signal.
% IN:
% y: signal
% nw: time-halfbandwidth product
% window: size of window in samples for sliding window
% 'noverlap': size of overlap in samples
% 'f': frequencies
% 'fs': sampling frequency
% 'ax': handle to axis in which to plot
%
% OUT:
% pxx: power spectral density estimate
% w: frequencies
% t: time axis with windows centered at each point in time
%
% TO DO:
% - finish commenting and cleaning up
%
% Adrian Tasistro-Hart, adrianraph-at-gmail.com, 04.08.2018

function [pxx,w,t] = mtgram(y,nw,window,varargin)

%% parse
parser = inputParser;
addRequired(parser,'y',@isnumeric);
addRequired(parser,'nw',@isscalar);
addRequired(parser,'window',@isscalar);
addOptional(parser,'noverlap',[],@isscalar);
addOptional(parser,'f',[],@isnumeric);
addOptional(parser,'fs',1,@isscalar);
addOptional(parser,'axis',[],@ishandle);

parse(parser,y,nw,window,varargin{:});

y = parser.Results.y;
nw = parser.Results.nw;
window = parser.Results.window;
noverlap = parser.Results.noverlap;
f = parser.Results.f;
fs = parser.Results.fs;
ax = parser.Results.axis;
    
%% set dynamic defaults and validate

% if only one column of data, make sure its a column
if size(y,1) == 1
    y = y';
end

% number of data points
n = size(y,1);
% make sure that window is smaller than the length of the data
assert(window < n,'window is larger than length of data')
% make window an integer
window = double(window);

% if no value for overlp, make 50%
if isempty(noverlap)
    noverlap = window/2;
% otherwise make sure that overlap is less than the window
else
    assert(noverlap < window, 'overlap must be less than window length')
end

% set default number of frequencies
if isempty(f)
    % Frequency axis increment is the sampling frequency divided by number of
    % samples
    fi = fs/window;
    % Frequency axis goes from 0 to 1/2 the sampling frequency
    f = 0:fi:fs/2;
end

% recalculate overlap given discrete number of slides and given window
% width
nslides = floor((n-noverlap)/(window-noverlap));  % number of slides
noverlap = ceil((nslides*window-n)/(nslides-1));

%% window and compute psds
Y = zeros(window,nslides);
% do all but last windows
for j = 1:nslides-1
    lidx = round( (j-1)*(window-noverlap)+1 );
    ridx = round( j*window-(j-1)*(noverlap) );
    Y(:,j) = detrend(y(lidx:ridx));
end
% do last window
Y(:,end) = detrend(y(end-window+1:end));
% compute multitaper estimates
[pxx,w] = pmtm(Y,nw,f,fs);
% get time axis
dt = 1/fs;
tw = (window-noverlap) * dt; % spacing of time axis
t1 = (window-1)/2 * dt;
t2 = ( (n-1) - (window-1)/2 ) * dt;
t = t1:tw:t2;
t(end) = t2;

%% plot
if isempty(ax)
    figure
else
    axes(ax)
end
surf(w,t,10*log10(abs(pxx')+eps),'edgecolor','none')
view(2)
xlim([min(f) max(f)])
ylim([min(t) max(t)])
set(gca,'xscale','log')

colormap(jet)

end