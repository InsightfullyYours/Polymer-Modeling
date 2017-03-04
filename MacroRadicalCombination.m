function resultDist=MacroRadicalCombination(xaxis,signal1,signal2)
%this is a subfunction of Hamielec that calculates the macroradical
%recombination of any 2 input distributions.  It assumes the distributions
%are correct (number fraction)).
%we again employ nested for loops; the outer to cycle through r and the
%inner to cycle through s.
%input should be small MW to large MW for all.

%determine the largest chain to go to.
maxMonomerNum=max(xaxis)-min(xaxis)+1;
minMonomerNum=min(xaxis);
%create variables
resultDist=zeros(size(signal1,1),1);

%cycle through n(r)
for r=2:1:maxMonomerNum
    %cycle through all smaller molecules
    thesum=0;
    for s=1:1:r-1
        thesum=thesum + s.*(r-s).*signal1(s,:).*signal2(r-s,:);
    end
    resultDist(r,1)=thesum;
end

%make sure it is a number fraction
resultDist=resultDist./repmat(sum(resultDist,1),size(resultDist),1);
end
