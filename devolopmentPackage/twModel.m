function M = twModel(fs)  

%td
%  Add noise option
%  add variable noise
%  add constant offset


% x = sin(linspace(1,10,20));
% y = sin(linspace(2,11,20));
% z = cos(linspace(2,11,20));

% q&d data gen
q = sin(linspace(2,11,fs)+.5);
r = sin(linspace(2,11,fs)+1);
s = sin(linspace(2,11,fs)+1.5);
t = sin(linspace(2,11,fs)+2);

% noisX = x + .5 * rand(1, length(x));


M = [q; r; s; t];


% corr(x,y)
% corr(x,x)
% corr(x, noisX)
% 
% corrcoef(x,noisX)
% corrcoef(x,y)
% corrcoef(x,2*y)
% corrcoef(x,1+y)
% corrcoef(x,cos(y))

end