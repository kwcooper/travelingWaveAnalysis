% tw grab asymetry
function [asmScores] = twComputeAsym(root)

asmScores = [];
for i = 1:size(root.b_lfp,2)
  root.active_lfp = i;
  [~,~,~,asymScores] = thetaAsymmetryAnalysis(root);
  asmScores(1:size(asymScores,1), end+1) = asymScores;
end

%mean(asmScores)

end

