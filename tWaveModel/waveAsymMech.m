%x = 1:100;
%plot(sin(x/100))
%%
Fs = 10000;                  % samples per second
dt = 1/Fs;                   % seconds per sample
duration = 1;                % seconds
t = (0:dt:duration-dt)';     % seconds
%% Define wave:
Fc = 6;                      % Hz
theta = cos(2*pi*Fc*t);
figure; plot(t,theta);
xlabel('time (in seconds)');
title('Theta');

%%
% make activation pulse
%  to do: turn into function, negitave pulse, 
z = zeros(1,length(t));
startPulse = 250;
numPulse = 6;
pulseWidth = 500; 
pulseHeight = .5;
for i = startPulse:floor(length(t)/numPulse):length(t)
    % define pulse
    st = linspace(0,1,pulseWidth);
    pulse = sin(2*pi*.5*st) * pulseHeight; % Generate Sine pulse
    w = floor(pulseWidth/2);
    if  (i - w) <= 0 % check if first pulse
       z(i:i+w) = pulse(length(pulse)/2:end);
    elseif (i + w) > length(z)
        z(i-w:i-1) = pulse(1:length(pulse)/2);
    else
       z(i-w:i+w-1) = pulse;
    end
end
figure; plot(t,z);

figure; plot(t,z'+theta);

%figure; plot(st,x);



