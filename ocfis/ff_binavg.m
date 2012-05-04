function [ybin, xbin, n, s,  Y, X] = ff_binavg(x,y,db)
% ff_binavg.m - bin average y(x) onto ybin(xbin)
%
% use:   [ybin,xbin,n,s] = ff_binavg(x,y,db)
% input:
%    y(x) - data (vel/prof)
%    x    - data axis (time/depth)
%    db   - bin interval in x units
%
% output:
%    ybin - binned data
%    xbin - new binned axis
%    n    - number of points in each bin
%    s    - std of data in each bin
%
% example:
%    [ybin, n, s] = ff_binavg(time, uvel, 1)
%
% other m-files required: matlab's histc.m
%
% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     28-May-2009
% modified: 04-May-2012
%
% obs: modified from bindata.m
%

% Make them all column vectors.
x = x(:);
y = y(:);

% Cut the corners.
x_min = ceil(min(x));
x_max = floor(max(x));

f = find(x >= x_min & x <= x_max);
x = x(f);
y = y(f);

xbin = x_min:db:x_max;

% Using find.
len_x = length(xbin);
for kk=1:len_x-1
        f = find(x >= xbin(kk) & x <= xbin(kk+1));
        Y(kk) = mean(y(f));
        X(kk) = (1/db) * (xbin(kk) + xbin(kk+1));
end

% NaNs check.
idx  = find(~isnan(y));
if (isempty(idx))
    b = nan * xbin;
    n = nan * xbin;
    s = nan * xbin;
end

x = x(idx);
y = y(idx);

binwidth = diff(xbin);
xx = [xbin(1)-binwidth(1)/2 xbin(1:end-1) + binwidth/2 xbin(end)+binwidth(end)/2];

% shift bins so the interval is "( ]" instead of "[ )".
bins = xx + max(eps, eps * abs(xx));

[n, bin] = histc(x, bins, 1);
n = n(1:end-1);
n(n==0) = NaN;

idx  = find(bin > 0);
sum  = full(sparse(bin(idx), idx * 0 + 1, y(idx)));
sum  = [sum; zeros(length(xbin) - length(sum), 1) * nan];
ybin = sum./n;

sum  = full(sparse(bin(idx), idx * 0 + 1, y(idx).^2));
sum  = [sum; zeros(length(xbin) - length(sum), 1) * nan];
s    = sqrt(sum./(n-1) - ybin.^2.*n./(n-1) );

xbin = xbin + 1/db;  % FIXME