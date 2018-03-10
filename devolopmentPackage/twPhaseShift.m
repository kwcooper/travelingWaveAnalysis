%find phase shift 

x = sin(linspace(1,10,100));
noisX = x + .5 * rand(1, length(x));
plot(noisX);

w1 = x;
w2 = noisX;
% calculate phase offset in rad then convert to deg
offsetRad = acos(dot(w1,w2)/(norm(w1)*norm(w2)));
offsetDeg = offsetRad*360/(2*pi);



% fit sin wave
y = noisX'

yu = max(y);
yl = min(y);
yr = (yu-yl);                               % Range of ‘y’
yz = y-yu+(yr/2);
zx = x(yz .* circshift(yz,[0 1]) <= 0);     % Find zero-crossings
per = 2*mean(diff(zx));                     % Estimate period
ym = mean(y);                               % Estimate offset
fit = @(b,x)  b(1).*(sin(2*pi*x./b(2) + 2*pi/b(3))) + b(4);    % Function to fit
fcn = @(b) sum((fit(b,x) - y).^2);                              % Least-Squares cost function
s = fminsearch(fcn, [yr;  per;  -1;  ym])                       % Minimise Least-Squares
xp = linspace(min(x),max(x));
figure(1)
plot(x,y,'b',  xp,fit(s,xp), 'r')