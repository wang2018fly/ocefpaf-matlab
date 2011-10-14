function map = cbathy(m);
% cbathy.m - colormap for bathymetry
%
% use:  map = cbathy(m);
% input:
%    input - colormap length (64)
%
% output:
%    output - colormap
%
% example:
%    cm1 = cbathy(255);
%
% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     22-May-2009
% modified: 22-May-2009
%
% obs: original matfile was only 64 length instead of 255
%

map = [...
     8   241   255
     8   237   255
     8   234   255
     8   231   255
     8   228   255
     8   225   255
     8   222   255
     8   218   255
     8   215   255
     8   212   255
     8   209   255
     8   204   255
     8   199   255
     9   194   255
     9   190   255
     9   185   255
     9   180   255
     9   175   255
     9   171   255
     9   166   255
     9   161   255
     9   160   255
     9   158   255
     9   157   255
     9   156   255
     9   154   255
     9   153   255
     9   152   255
     9   150   255
     9   149   255
     9   148   255
     9   146   255
     9   145   255
     9   144   255
     9   142   255
     9   141   255
    10   140   255
    10   134   255
    10   129   255
    10   124   255
    10   119   255
    10   114   255
    10   109   255
    10   103   255
    10    98   255
    10    93   255
    10    88   255
    11    83   255
    11    77   255
    10    72   247
     9    67   240
     8    62   232
     8    57   225
     7    51   217
     6    46   210
     5    41   202
     5    36   195
     4    31   188
     3    25   180
     2    20   173
     2    15   165
     1    10   158
     0     5   150
     0     0   143];

map = map./255;

if nargin > 0
  if m ~=64
    map = interp1(1:64,map,linspace(1,64,m));
  end
end
