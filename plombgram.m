% This function computes the Lomb-Scargle spectrogram of a given signal.
% 
% IN:
% x: time series coordinates
% y: time series values
% window: size of window in samples for sliding window
% 'noverlap': size of overlap in samples
% 'f': frequencies at which to sample
% 'ax': handle to axis in which to plot
% 'plotit': (true) whether to plot the spectrogram
%
% OUT:
% pxx: power spectral density estimate
% w: frequencies
% t: time axis with windows centered at each point in time
%
% TO DO:
% - 
%
% Adrian Tasistro-Hart, adrianraph-at-gmail.com, 02.09.2018

function [pxx,w,t] = plombgram(x,y,window,varargin)

%% parse
parser = inputParser;
addRequired(parser,'x',@isnumeric);
addRequired(parser,'y',@isnumeric);
addRequired(parser,'window',@isscalar);
addParameter(parser,'noverlap',[],@isscalar);
addParameter(parser,'f',[],@isnumeric);
addParameter(parser,'axis',[],@ishandle);
addParameter(parser,'plotit',true,@islogical);

parse(parser,x,y,window,varargin{:});

x = parser.Results.x;
y = parser.Results.y;
window = parser.Results.window;
noverlap = parser.Results.noverlap;
f = parser.Results.f;
ax = parser.Results.axis;
plotit = parser.Results.plotit;
    

%% set dynamic defaults and validate

n = length(y);  % number of samples
assert(length(x) == n, 'x and y must be same length')

% make column
y = y(:);

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

% recalculate overlap given discrete number of slides and given window
% width
nslides = floor((n-noverlap)/(window-noverlap));  % number of slides
noverlap = ceil((nslides*window-n)/(nslides-1));

% set default number of frequencies
if isempty(f)
    nfft = 2*window;
else
    nfft = length(f);
end


%% window and compute psds

% array of windowed signals
Y = zeros(window,nslides);
% array of coordinates for points in each window
X = zeros(window,nslides);
% do all windows except last
for ii = 1:nslides-1
    lidx = round( (ii-1)*(window-noverlap)+1 );
    ridx = round( ii*window-(ii-1)*(noverlap) );
    Y(:,ii) = detrend(x(lidx:ridx),y(lidx:ridx));
    X(:,ii) = x(lidx:ridx);
end
% do last window
Y(:,end) = detrend(x(end-window+1:end),y(end-window+1:end));
X(:,end) = x(end-window+1:end);

% compute Lomb-Scargle estimates
pxx = zeros(nfft,nslides);
for ii = 1:nslides
    [pxx(:,ii),w] = plomb(Y(:,ii),X(:,ii),f);
end

% get time axis; need center of each group of coordinates
t = mean(X);


%% plot

if plotit
    
    if isempty(ax)
        figure
    else
        axes(ax)
    end

    surf(w,t,10*log10(abs(pxx')+eps),'edgecolor','none')
    view(2)
    xlim([min(w) max(w)])
    ylim([min(t) max(t)])
    set(gca,'xscale','log')

    colormap(jet)

end


end