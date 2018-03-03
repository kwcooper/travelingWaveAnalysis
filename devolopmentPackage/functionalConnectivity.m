%functional connectivity

M = twModel();

M = origData;


for i = 1:size(M,1)
  for j = 1:size(M,1)
    cc = corrcoef(M(i,:), M(j,:));
    corrMat(i,j) = cc(1,2);
  end
end

figure; imagesc(corrMat)
title('DFC: Rio')

for i = 1:4
  chan = linspace(1,10,20)

end




  





























