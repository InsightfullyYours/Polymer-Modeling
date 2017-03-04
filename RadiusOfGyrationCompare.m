function [FNaxisRg estimates]=RadiusOfGyrationCompare(linRg2b2,BranchRgs,FNaxisUse,FNcompare)

start_point = 5.*rand(1, 1);
options=optimset('Display','iter','MaxIter',10000,'MaxFunEvals',10000,'TolX',1E-40,'TolFun',1E-40);
estimates = fminsearch(@RadiusCompare, start_point,options);
fitted=BranchRgs+estimates;

for i=1:1:size(BranchRgs,2)
    %we will "round" to nearest Rg on the linear Rg axis.
    %first, remove the chain from the linear N axis
    FNaxisUse(BranchRgs(1,i),1)=FNaxisUse(BranchRgs(1,i),1)-1;
    %If the Rg of the current chain is closer to the higher linear Rg, add
    %the chain to that population.  if it is closer to the lower linear Rg,
    %add the chain there.  This is the "rounding".
    if ~isempty(find(BranchRgs(2,i)==linRg2b2,1))
        index=find(BranchRgs(2,i)==linRg2b2);
        FNaxisUse(index,1)=FNaxisUse(index,1)+1;
    elseif abs(BranchRgs(2,i)-linRg2b2(find(linRg2b2<BranchRgs(2,i),1,'last'),:)) > abs(BranchRgs(2,i)-linRg2b2(find(linRg2b2>BranchRgs(2,i),1,'first'),:))
        index=find(linRg2b2>BranchRgs(2,i),1,'first');
        FNaxisUse(index,1)=FNaxisUse(index,1)+1;
    elseif abs(BranchRgs(2,i)-linRg2b2(find(linRg2b2<BranchRgs(2,i),1,'last'),:)) < abs(BranchRgs(2,i)-linRg2b2(find(linRg2b2>BranchRgs(2,i),1,'first'),:))
        index=find(linRg2b2<BranchRgs(2,i),1,'last');
        FNaxisUse(index,1)=FNaxisUse(index,1)+1;
    elseif abs(BranchRgs(2,i)-linRg2b2(find(linRg2b2<BranchRgs(2,i),1,'last'),:)) == abs(BranchRgs(2,i)-linRg2b2(find(linRg2b2>BranchRgs(2,i),1,'first'),:))
        %if it is exactly equidistant between the two (unlikely, but
        %possible), add .5 to each
        index=find(linRg2b2<BranchRgs(2,i),1,'last');
        index2=find(linRg2b2>BranchRgs(2,i),1,'first');
        FNaxisUse(index,1)=FNaxisUse(index,1)+0.5;
        FNaxisUse(index2,1)=FNaxisUse(index2,1)+0.5;
    else
        disp('Radius large then linear axis')
    end
end
FNaxisRg=FNaxisUse;


function sse=RadiusCompare(param)
param
BranchRgs(2,:)=BranchRgs(2,:)+param;
%go through each column, remove a chain the appropriate N from the
%FNaxisUse, then locate that
%chain on the linear Rg axis and add a chain there.
for i=1:1:size(BranchRgs,2)
    %we will "round" to nearest Rg on the linear Rg axis.
    %first, remove the chain from the linear N axis
    FNaxisUse(BranchRgs(1,i),1)=FNaxisUse(BranchRgs(1,i),1)-1;
    %If the Rg of the current chain is closer to the higher linear Rg, add
    %the chain to that population.  if it is closer to the lower linear Rg,
    %add the chain there.  This is the "rounding".
    if ~isempty(find(BranchRgs(2,i)==linRg2b2,1))
        index=find(BranchRgs(2,i)==linRg2b2);
        FNaxisUse(index,1)=FNaxisUse(index,1)+1;
    elseif abs(BranchRgs(2,i)-linRg2b2(find(linRg2b2<BranchRgs(2,i),1,'last'),:)) > abs(BranchRgs(2,i)-linRg2b2(find(linRg2b2>BranchRgs(2,i),1,'first'),:))
        index=find(linRg2b2>BranchRgs(2,i),1,'first');
        FNaxisUse(index,1)=FNaxisUse(index,1)+1;
    elseif abs(BranchRgs(2,i)-linRg2b2(find(linRg2b2<BranchRgs(2,i),1,'last'),:)) < abs(BranchRgs(2,i)-linRg2b2(find(linRg2b2>BranchRgs(2,i),1,'first'),:))
        index=find(linRg2b2<BranchRgs(2,i),1,'last');
        FNaxisUse(index,1)=FNaxisUse(index,1)+1;
    elseif abs(BranchRgs(2,i)-linRg2b2(find(linRg2b2<BranchRgs(2,i),1,'last'),:)) == abs(BranchRgs(2,i)-linRg2b2(find(linRg2b2>BranchRgs(2,i),1,'first'),:))
        %if it is exactly equidistant between the two (unlikely, but
        %possible), add .5 to each
        index=find(linRg2b2<BranchRgs(2,i),1,'last');
        index2=find(linRg2b2>BranchRgs(2,i),1,'first');
        FNaxisUse(index,1)=FNaxisUse(index,1)+0.5;
        FNaxisUse(index2,1)=FNaxisUse(index2,1)+0.5;
    else
        disp('Radius large then linear axis')
    end
end

ErrorVector=FNaxisUse./repmat(sum(FNaxisUse,1),size(FNaxisUse,1),1)-FNcompare;
sse = sum(ErrorVector .^ 2);

end
end