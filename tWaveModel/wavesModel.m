




thetaEC = nan(numCells, 1);


P_EC = pi/2;
n = 1;
for t = 1:100
th(t,:) = .5 * sin((t + (P_EC * n))) + .5;
end
figure; plot(th);


num_neurons = 5;
P_EC = pi/2;

EC_A = [1; 1; 1; 1; 1];

%iterate through time
for t = 1:10
  % iterate over neurons
  for i = 1:size(EC_A,1)
    EC_AT(i,t) = EC_A(i) * (.5 * sin((t + (P_EC - (pi/size(EC_A,1) * i)))) + .5);
  end
end

figure; 
for i = 1:size(EC_AT, 1)
  subplot(size(EC_AT,1), 1, i); plot(EC_AT(i,:))
end


nn = 5;
for i = 1:5
(pi / nn) * i
end

180 / num_neurons