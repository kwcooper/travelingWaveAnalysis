function [asmScores] = twComputeAsym(root)
% Wrapper for thetaAsymmetryAnalysis 
% Returns list of scores 
% 180710 kwc


asmScores = [];
for i = 1:size(root.b_lfp,2)
  root.active_lfp = i;
  [~,~,~,asymScores] = thetaAsymmetryAnalysis(root);
  asmScores(1:size(asymScores,1), end+1) = asymScores;
end

%mean(asmScores)

end

