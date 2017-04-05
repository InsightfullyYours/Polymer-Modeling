function [CombScisPureN] = CombinationScissionWithRg(NIn,FNIn,FNcompare,fraction)
%This function models the CombScis distribution.
%designed in conjunction with John Zhang

%define the distribution
% NIn=(1:1:10)';
% FNIn=[1;2;3;3;2;2;1;0;0;0];
% NIn=(1:1:100)';
% FNIn=round(10000./(10.*sqrt(2.*pi)).*exp(-(NIn-50).^2./(2.*10^2)));

FNIn=round(10000.*FNIn);

%Define the axes
Naxis=(1:1:max(NIn))';
Combaxis1=(1:1:max(Naxis))';
Combaxis=(1:1:max(Naxis))';
FNaxis=zeros(size(Naxis,1),1);

%populate the axes with F(N) distribution.  This allows distribution d
for i=1:1:size(NIn,1)
    FNaxis(Naxis==NIn(i,1))=FNIn(i);
end

%what percent of the bonds do we want to react?
%fraction=.10;       %fractional percent of all initial bonds that react
balance=0;       %distribution between pure combination (0) and pure scission (1)

%create a vector of every chain, by length.  The lengths are repeated once
%for every chain of that length.  So if there are 3 chains of length 1 and
%4 chains of length 2, the resulting vector will be
%  [1;1;1;2;2;2;2].
chains_static = cumsum(FNaxis);
%bond_static = cumsum(FNaxis.*Combaxis);

experiment_num=10;
%do the "experiment" experiment_num times and then we can average to get a better picture of
%the distribution.
CombScis_repeat=zeros(max(Combaxis),experiment_num);
CombScis_evolution=zeros(size(Combaxis,1),round(fraction.*max(chains_static)));
for j=1:1:experiment_num
    %the first column is the distribution; the second column counts the
    %number of single branched polymers in each part of the distribution;
    %third column counts double branched; four counts triple branched
    FNaxisUse=FNaxis;
    toobig=zeros(1,experiment_num);     %keep track of chains that get too big to be tracked
    toomanybranches=zeros(1,experiment_num);     %keep track of chains that have too many branches to be tracked, just in case
    %variables that track the combination locations of the branches.  These
    %are stored with each column a branched molecule.
    branchingarraychains=zeros(5,round(fraction.*max(chains_static)));
    branchingarraybonds=zeros(4,3.*round(fraction.*max(chains_static)));
    
    chain_counter=1;
    evolution_counter=1;
    while chain_counter <= round(fraction.*max(chains_static))
        random_reaction=rand(1);
        if random_reaction>balance
            
            bond=cumsum(FNaxisUse(:,1).*(Combaxis-1));
            %Where_It_Broke=[NaN,NaN,NaN];
            
            %choose a random bond
            bond_rand(1,1) = randi([1,max(bond)],1,1);
            %locate the length of chains selected by the bonds and remove a
            %chain
            index=find(bond<bond_rand(1,1),1,'last')+1;
            FNaxisUse(index,1)=FNaxisUse(index,1)-1;
            
            %locate the exact bond within that chain
            bond_number1=mod(bond(index,1)-bond_rand(1,1),Combaxis(index,1)-1);
            
            if bond_number1==0
                bond_number1=Combaxis(index,1)-1;
            end
            
            %we determine if this chain is branched
            %the branching array consists of 4 columns per branch, with 5 rows.  The
            %first column is the chain lengths in a branch, the second column the first
            %two bonds to react, the third the second two bonds to react,
            %and the
            %fourth the third two bonds to react.  The bottom row is the total length
            %of the chain. We assume the branched chains are at the beginnning of the
            %size.
            %first we make sure the chain lengths of the previously
            %branched molecules are correct
            branchingarraychains(5,:)=sum(branchingarraychains(1:4,:),1);
            %we know from above the length of the chosen chain; we search
            %the branching array for all chains of that length and compute
            %the number of bonds in them
            [row,column]=find(branchingarraychains(5,:)==index);
            branchbonds=sum(row,2).*(Combaxis(index,1)-1);        %calculate the number of bonds in branched molecules of that size
            chosenindex1=[];
            if branchbonds+bond(index-1,1)<bond_rand(1,1)    %if the randomly selected bond is higher than the number of bonds in the previous size plus the bonds in the branch for the current size, then it's not branched
                chain1is=0;  %no branching
                %find the first empty column in the branching matrix to
                %make a new branch
                chosenindex1=find(branchingarraychains(5,:)==0,1);
                branchingarraychains(1,chosenindex1)=index; %add the chain length
                branchingarraybonds(1,3.*chosenindex1-2)=bond_number1; %add the bond number
            elseif branchbonds+bond(index-1,1)>=bond_rand(1,1)       %if otherwise, it is branched
                %determine which of the branched molecules has been chosen,
                %and the bond within that molecule
                if rem(bond_rand(1,1)-bond(index-1,1),Combaxis(index,1)-1)==0
                    chosenchain=floor((bond_rand(1,1)-bond(index-1,1))/(Combaxis(index,1)-1));
                else
                    chosenchain=floor((bond_rand(1,1)-bond(index-1,1))/(Combaxis(index,1)-1))+1;
                end
                %get the chain itself from the column
                chosenindex1=column(:,chosenchain);
                %get info on chain-how many branches already?
                chain1is=find(branchingarraychains(1:4,chosenindex1)==0,1)-2;    %if 1, one branch.  if 2, two branch.  if three,we must eliminate it
                %allows infinite branching off the same bond
                if ~isempty(chain1is)   %if it has less than 3 branches, add the new reacting bond
                    %first determine what chain within the already branched chain reacted
                    if bond_number1<=branchingarraychains(1,chosenindex1);    %first chain
                        branchingarraybonds(1,3.*chosenindex1-2+chain1is)=bond_number1;   %add the bond that reacted to the bon array
                    elseif bond_number1>branchingarraychains(1,chosenindex1) && bond_number1<=sum(branchingarraychains(1:2,chosenindex1),1) %second chain
                        branchingarraybonds(2,3.*chosenindex1-2+chain1is)=bond_number1-branchingarraychains(1,chosenindex1);
                    elseif bond_number1>sum(branchingarraychains(1:2,chosenindex1),1) && bond_number1<=sum(branchingarraychains(1:3,chosenindex1),1)    %third chain
                        branchingarraybonds(3,3.*chosenindex1-2+chain1is)=bond_number1-sum(branchingarraychains(1:2,chosenindex1),1);
                    else
                        disp('This should never print CombinationScissionWithRg:1')
                    end
                elseif isempty(chain1is)
                    chain1is=3;
                end
            else
                disp('This should never print CombinationScissionWithRg:2')
            end
            
            %recount the bonds to account for the selected chain
            %(This sets interpolymer reactions only)
            bond=cumsum(FNaxisUse(:,1).*(Combaxis-1));
            %choose second random bond and remove the chain from the
            %distribution
            bond_rand(2,1) = randi([1,max(bond)],1,1);
            index2=find(bond<bond_rand(2,1),1,'last')+1;
            FNaxisUse(index2,1)=FNaxisUse(index2,1)-1;
            
            %locate the exact bond within that chain
            bond_number2=mod(bond(index2,1)-bond_rand(2,1),Combaxis(index2,1)-1);
            
            if bond_number2==0
                bond_number2=Combaxis(index2,1)-1;
            end
            
            
            %we determine again if this chain is branched
            %we know from above the length of the chosen chain; we search
            %the branching array for all chains of that length and compute
            %the number of bonds in them
            [row,column]=find(branchingarraychains(5,:)==index2);
            branchbonds=sum(row,2).*(Combaxis(index2,1)-1);        %calculate the number of bonds in branched molecules of that size
            chosenindex2=[];      %this insures there's not an error with the ~isempty(chosenindex2) below when checking for a second branched molecule.  Otherwise its value is held over from the previous iteration.
            if branchbonds+bond(index2-1,1)<bond_rand(2,1)    %if the randomly selected bond is higher than the number of bonds in the previous size plus the bonds in the branch for the current size, then it's not branched
                chain2is=0;  %no branching
            elseif branchbonds+bond(index2-1,1)>=bond_rand(2,1)       %if otherwise, it is branched
                if rem(bond_rand(2,1)-bond(index2-1,1),Combaxis(index2,1)-1)==0
                    chosenchain=floor((bond_rand(2,1)-bond(index2-1,1))/(Combaxis(index2,1)-1));
                else
                    chosenchain=floor((bond_rand(2,1)-bond(index2-1,1))/(Combaxis(index2,1)-1))+1;
                end
                %get the chain itself from the column
                chosenindex2=column(:,chosenchain);
                %get info on chain-how many branches already?
                chain2is=find(branchingarraychains(1:4,chosenindex2)==0,1)-2;    %if 1, one branch.  if 2, two branch.  if three,we must eliminate it
            end
            if isempty(chain2is)
                chain2is=3;
            end
            %we only allow up to 3 branches.  Check if there are more than
            %3 branches in the combination, and make the chain(s) dissapear
            %if there are
            if chain1is+chain2is+1>3
                branchingarraychains=[branchingarraychains(:,1:chosenindex1-1) branchingarraychains(:,chosenindex1+1:end)];
                branchingarraybonds=[branchingarraybonds(:,1:3.*(chosenindex1-1)) branchingarraybonds(:,3.*chosenindex1+1:end)];
                if ~isempty(chosenindex2)
                    branchingarraychains=[branchingarraychains(:,1:chosenindex2-1) branchingarraychains(:,chosenindex2+1:end)];
                    branchingarraybonds=[branchingarraybonds(:,1:3.*(chosenindex2-1)) branchingarraybonds(:,3.*chosenindex2+1:end)];
                end
                toomanybranches(1,experiment_num)=toomanybranches(1,experiment_num)+1;
                continue    %skip the rest of this loop
            end
            
            switch chain2is
                case 0 %if chain 2 is not branched; add it to the branching array
                    emptyrow=find(branchingarraychains(1:4,chosenindex1)==0,1);
                    branchingarraychains(emptyrow,chosenindex1)=index2; %add the chain length
                    branchingarraybonds(emptyrow,3.*chosenindex1-2+chain1is)=bond_number2; %add the bond number
                    
                case 1 %if it has 1 branch
                    %transplant the lengths directly to the first chain vector.
                    emptyrow=find(branchingarraychains(1:3,chosenindex1)==0,1);
                    emptycolumn=find(sum(branchingarraybonds(:,3.*chosenindex1-2:3.*chosenindex1),1)==0,1);  %should always be 2 for the same reason
                    branchingarraychains(emptyrow:emptyrow+1,chosenindex1)=branchingarraychains(1:2,chosenindex2);
                    branchingarraybonds(:,3.*chosenindex1-3+emptycolumn:3.*chosenindex1-2+emptycolumn)=circshift(branchingarraybonds(:,3.*chosenindex2-2:3.*chosenindex2-1),emptyrow-1);
                    %first determine what chain within the already branched
                    %chain reacted
                    if bond_number2<=branchingarraychains(1,chosenindex2);    %first chain
                        branchingarraybonds(emptyrow,3.*chosenindex1-2+chain1is)=bond_number2;   %add the bond that reacted to the bon array
                    elseif bond_number2>branchingarraychains(1,chosenindex2) && bond_number2<=sum(branchingarraychains(1:2,chosenindex2),1) %second chain
                        branchingarraybonds(emptyrow+1,3.*chosenindex1-2+chain1is)=bond_number2-branchingarraychains(1,chosenindex2);
                    else
                        disp('This should never print CombinationScissionWithRg:3')
                    end
                    
                    %remove chain 2 from the branching arrays
                    branchingarraychains=[branchingarraychains(:,1:chosenindex2-1) branchingarraychains(:,chosenindex2+1:end)];
                    branchingarraybonds=[branchingarraybonds(:,1:3.*(chosenindex2-1)) branchingarraybonds(:,3.*chosenindex2+1:end)];
                    
                case 2 %if it has 2 branches
                    %transplant the lengths directly to the first chain vector.
                    emptyrow=find(branchingarraychains(1:2,chosenindex1)==0,1); %should always be 2 as we just created this branch in the first chain section
                    emptycolumn=find(sum(branchingarraybonds(:,3.*chosenindex1-2:3.*chosenindex1),1)==0,1);  %should always be 2 for the same reason
                    branchingarraychains(emptyrow:emptyrow+2,chosenindex1)=branchingarraychains(1:3,chosenindex2);
                    branchingarraybonds(:,3.*chosenindex1-3+emptycolumn:3.*chosenindex1-2+emptycolumn)=circshift(branchingarraybonds(:,3.*chosenindex2-2:3.*chosenindex2-1),1);
                    %first determine what chain within the already branched chain reacted
                    if bond_number2<=branchingarraychains(1,chosenindex2);    %first chain
                        branchingarraybonds(emptyrow,3.*chosenindex1-4+emptycolumn)=bond_number2;   %add the bond that reacted to the bon array
                    elseif bond_number2>branchingarraychains(1,chosenindex2) && bond_number2<=sum(branchingarraychains(1:2,chosenindex2),1) %second chain
                        branchingarraybonds(emptyrow+1,3.*chosenindex1-4+emptycolumn)=bond_number2-branchingarraychains(1,chosenindex2);
                    elseif bond_number2>sum(branchingarraychains(1:2,chosenindex2),1) && bond_number2<=sum(branchingarraychains(1:3,chosenindex2),1)    %third chain
                        branchingarraybonds(emptyrow+2,3.*chosenindex1-4+emptycolumn)=bond_number2-sum(branchingarraychains(1:2,chosenindex2),1);
                    else
                        disp('This should never print CombinationScissionWithRg:4')
                    end
                    
                    %remove chain 2 from the branching arrays
                    branchingarraychains=[branchingarraychains(:,1:chosenindex2-1) branchingarraychains(:,chosenindex2+1:end)];
                    branchingarraybonds=[branchingarraybonds(:,1:3.*(chosenindex2-1)) branchingarraybonds(:,3.*chosenindex2+1:end)];
                    
                otherwise
                    disp('This should never print CombinationScissionWithRg:5')
            end
            
            
            
            %if the size of the combined chain is greater than 5 times the
            %maximum input array size (irrespective of if there are any chains at
            %that size; this is about just the array size) these chains
            %"dissapear" from the distribution and are counted as reacted,
            %though too big to be part of the output distribution.
            if index+index2>max(Combaxis)
                toobig(1,experiment_num)=toobig(1,experiment_num)+1;
                continue
            end
            
            
            %add the new, combined chain to the distribution.  This
            %variable is the pure linear length; radius of gyration changes
            %in output location are not considered.
            FNaxisUse(index+index2,1)=FNaxisUse(index+index2,1)+1;
            
            %now the distribution is adjusted.
            if evolution_counter==1  %if its the first experiment, the changes are just the distribution
                CombScis_evolution(:,evolution_counter)=FNaxisUse(:,1);
            else  %account for multiple experiments by averaging the new experiment with previous
                CombScis_evolution(:,evolution_counter)=mean([FNaxisUse(:,1), CombScis_evolution(:,evolution_counter)],2);
            end
            
            chain_counter=chain_counter+2;
            evolution_counter=evolution_counter+1;
            
        elseif random_reaction<balance
            
            bond=cumsum(FNaxisUse(:,1).*(Combaxis-1));
            %Where_It_Broke=[NaN,NaN,NaN];
            
            %choose a random bond
            bond_rand = randi([1,max(bond)],1,1);
            
            %locate the index (length of chain) of that bond
            index=find(bond<bond_rand,1,'last')+1;
            %locate the exact bond with the chain that broke
            bond_number=mod(bond(index,1)-bond_rand,Combaxis(index,1)-1);
            %The size of one piece is the bond_number (since a break of bond 2
            %produces a piece of size 2); the size of the other is the chain length
            %minus the bond number+1. Remove 1 chain from the original length and
            %add 1 chain each to each pieces' length.  A bond_number of zero
            %indicates the last bond in a chain; this is equivalent to
            %Naxis(index,1)-1
            if bond_number==0
                bond_number=Combaxis(index,1)-1;
            end
            FNaxisUse(index,1)=FNaxisUse(index,1)-1;
            FNaxisUse(bond_number,1)=FNaxisUse(bond_number,1)+1;
            FNaxisUse(Combaxis(index,1)-bond_number,1)=FNaxisUse(Combaxis(index,1)-bond_number,1)+1;
            
            %now the distribution is adjusted.
            if evolution_counter==1  %if its the first experiment, it's just the distribution
                CombScis_evolution(:,evolution_counter)=FNaxisUse(:,1);
            else  %account for multiple experiments by averaging the new experiment with previous
                CombScis_evolution(:,evolution_counter)=mean([FNaxisUse(:,1), CombScis_evolution(:,evolution_counter)],2);
            end
            chain_counter=chain_counter+1;
            evolution_counter=evolution_counter+1;
        elseif random_reaction == balance
            continue
        else
            disp('This should never print error CombinationScissionWithRg:6')
            keyboard
        end
    end
    CombScis_repeat(:,j)=FNaxisUse;
    branchingarraychains(5,:)=sum(branchingarraychains(1:4,:),1);
    [linRg2b2(:,j) FNaxisRg(:,j) singlebranchchains twobranchchains threebranchchains onebranchRg twobranchRg threebranchRg]=RadiusOfGyration(Combaxis,FNaxisUse,branchingarraychains,branchingarraybonds);

    
end

%----------------------------------------------------------------------
%output, theory, and display code

figure
CombScisPureN=mean(CombScis_repeat,2);
chains_static2 = sum(CombScisPureN,1);
plot(Combaxis,CombScisPureN,'og')
CombScisPureN=CombScisPureN./repmat(sum(CombScisPureN,1),size(CombScisPureN,1),1);  %scale to 1
title('Direct Chain Lengths')

linRg2b2Ave=mean(linRg2b2,2);
FNaxisRgAve=mean(FNaxisRg,2);
    BranchRgs=[onebranchRg twobranchRg threebranchRg];
 %   [FittedFNaxisRg(:,j) gestimates(:,j)]=RadiusOfGyrationCompare(linRg2b2Ave,BranchRgs,CombScisPureN,FNcompare);

figure
plot(flipud(Combaxis./6),FNcompare./repmat(sum(FNcompare,1),size(FNcompare,1),1))
hold
plot(linRg2b2Ave,FNaxisRgAve./repmat(sum(FNaxisRgAve,1),size(FNaxisRgAve,1),1),':r')
plot(Combaxis./6,CombScisPureN,':k')
title('Radius of Gyration/b^2')

%We can then use the input distribution to calculate the theory we've
%developed; convolution(F(N),F(N)). We do a few powers of N for comparison.
FNaxis=FNaxis(1:max(NIn),:);
oN4fitbeg=conv(FNaxis./Naxis.^4,FNaxis);
oN3fitbeg=conv(FNaxis./Naxis.^3,FNaxis);
oN2fitbeg=conv(FNaxis./Naxis.^2,FNaxis);
oNfitbeg=conv(FNaxis./Naxis,FNaxis);
fitbeg=conv(FNaxis,FNaxis);
Nfitbeg=conv(Naxis.*FNaxis,FNaxis);
N2fitbeg=conv(Naxis.^2.*FNaxis,FNaxis);
N3fitbeg=conv(Naxis.^3.*FNaxis,FNaxis);
N4fitbeg=conv(Naxis.^4.*FNaxis,FNaxis);

%convolution neglects the first value; the convolution of (1,1) and (1,1)
%is (2,1).  So we pad the begining with a zero.
oN4fit=[0;oN4fitbeg];
oN3fit=[0;oN3fitbeg];
oN2fit=[0;oN2fitbeg];
oNfit=[0;oNfitbeg];
fit1=[0;fitbeg];
Nfit=[0;Nfitbeg];
N2fit=[0;N2fitbeg];
N3fit=[0;N3fitbeg];
N4fit=[0;N4fitbeg];

%normalize by the magnitude.
FNaxisNorm=FNaxis./repmat(sum(FNaxis,1),size(FNaxis,1),1);
oN4fit=oN4fit./repmat(sum(oN4fit,1),size(oN4fit,1),1);
oN3fit=oN3fit./repmat(sum(oN3fit,1),size(oN3fit,1),1);
oN2fit=oN2fit./repmat(sum(oN2fit,1),size(oN2fit,1),1);
oNfit=oNfit./repmat(sum(oNfit,1),size(oNfit,1),1);
fit1=fit1./repmat(sum(fit1,1),size(fit1,1),1);
Nfit=Nfit./repmat(sum(Nfit,1),size(Nfit,1),1);
N2fit=N2fit./repmat(sum(N2fit,1),size(N2fit,1),1);
N3fit=N3fit./repmat(sum(N3fit,1),size(N3fit,1),1);
N4fit=N4fit./repmat(sum(N4fit,1),size(N4fit,1),1);

%We display the results
figure
plot(Naxis./6,FNaxisNorm,'bo')
hold
plot(Combaxis./6,CombScisPureN,'r^')
plot(Combaxis1./6,oN4fit,'-b')
plot(Combaxis1./6,oN3fit,'-r')
plot(Combaxis1./6,oN2fit,':b')
plot(Combaxis1./6,oNfit,':r')
plot(Combaxis1./6,fit1,'-.b')
plot(Combaxis1./6,Nfit,'-.r')
plot(Combaxis1./6,N2fit,'--b')
plot(Combaxis1./6,N3fit,'--r')
plot(Combaxis1./6,N4fit,'-k')
text(350,0.01,['Chains Before',num2str(max(chains_static))]);
text(350,0.005,['Chains After',num2str(max(chains_static2))]);
legend(['Input Distribution       ';'CombScis Product Dist.   ';'1/N^4                    ';'1/N^3                    ';'1/N^2                    ';'1/N                      ';'No N                     ';'N                        ';'N^2                      ';'N^3                      ';'N^4                      '],'Location','EastOutside')
xlabel('Number of Monomers in Chain')
ylabel('Number of Chains')
title('Comparison of Modelling with Theory of CombScis (CombScisModel)')
end

