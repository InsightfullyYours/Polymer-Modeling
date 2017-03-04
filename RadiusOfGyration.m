function [linRg2b2 FNaxisRg singlebranchchains twobranchchains threebranchchains onebranchRg twobranchRg threebranchRg]=RadiusOfGyration(Naxis,FNaxisUse,Chains,Bonds)
%This function uses the equations from
%Zimm,Stockmeyer, J.Chem.Physics,v17#12,1949,pg 1301
%To calculate from an input distribution and branching locations the
%relative radius of gyration for polymers with 0,1,2,3 branch units

%first calculate linear dimention.  These equations are in R^2/b^2 space to
%eliminate the rotational, bond angle, and solvent effects within a
%polymer.  In other words, all below equations are R^2/b^2=...

%check to insure the sums are correct
Chains(5,:)=sum(Chains(1:4,:),1);

%initial definitions in case the system does not have and branches of that
%size
onebranchRg=[];
twobranchRg=[];
threebranchRg=[];

linRg2b2=Naxis./6;        %linear Rg axis: this is merely R^2=b^2*N/6 changed to R^2/b^2=N/6

%one branch  R^2/b^2=1/N*sumv(Nv
%find the chains in Chains that have only 1 branch
[~,column1br]=find(Chains(1,:)~=0 & Chains(2,:)~=0 & Chains(3,:)==0 & Chains(4,:)==0 & Chains(5,:)~=0);   %because the third row of chains is only if there are 2 or more branches

if ~isempty(column1br)  %if there are any branches
    %we get the chains and bonds of just the single chain species.  To
    %calculate Rg, we need to count the number of monomers in each arm of the
    %single branch species.
    %this is an easy calculation.  The chains variable has the length of the
    %chain; the bonds variable the bond that's reacted
    %for each branched chain, create a 4 row array with the branch lengths.
    singlebranchchains=Chains(:,column1br);
    
    %we have to do a bit of manipulation to get the bond columns from the chain
    %one.
    singlebranchbonds=Bonds(:,3.*column1br-2);   %we only need these, as one branch only uses the first of the 3 available branch columns
    
    for i=1:1:size(singlebranchchains,2)
        branchesofsingle(:,i)=[singlebranchbonds(1:2,i);singlebranchchains(1:2,i)-singlebranchbonds(1:2,i)];
        onebranchRg(1,i)=singlebranchchains(5,i);
        onebranchRg(2,i)=(1./singlebranchchains(5,i)).*sum((branchesofsingle(:,i).^2./2-branchesofsingle(:,i).^3./(3.*singlebranchchains(5,i))),1);
        gfactor(:,i)=(6./singlebranchchains(5,i).^2).*sum((branchesofsingle(:,i).^2./2-branchesofsingle(:,i).^3./(3.*singlebranchchains(5,i))),1);
    end
end

%two branch
%find the chains in Chains that have 2 branches
[~,column2br]=find(Chains(3,:)~=0 & Chains(4,:)==0 & Chains(5,:)~=0);   %because the third row of chains is only if there are 2 or more branches

if ~isempty(column2br)      %if there are any 2 branch chains
    twobranchchains=Chains(:,column2br);
    
    %we have to do a bit of manipulation to get the bond columns from the chain
    %one.
    bonds1=3.*column2br-2;
    bonds2=3.*column2br-1;
    twobranchcolumns=sort([bonds1 bonds2]);
    
    twobranchbonds=Bonds(:,twobranchcolumns);   %we need the first two of the 3 available branch columns
    
    for i=1:1:size(twobranchchains,2)
        %we eliminate a bug, sofar of unknown source, whereby there is
        %sometimes a column of all zeros counted as a two branch species.
        if sum(twobranchbonds(:,2.*i-1),1)==0 || sum(twobranchbonds(:,2.*i),1)==0
            break
        end
        
        
        %we find the row in the bonds that has no zeros; this row is the
        %cross-chain that contains the interbranch ab chain. We make it a
        %boolean vector, with 1 the row and 0 the others
        row=twobranchbonds(:,2.*i-1)~=0 & twobranchbonds(:,2.*i)~=0;
        
        %we find the column of the row with the larger #; this helps
        %determine how to match up the ends.  If inequality is 0 (false), then
        %the first column is greater; if it's 1 (true) then the second is
        %greater.
        column=twobranchbonds(row,2.*i-1) < twobranchbonds(row,2.*i);
        %this allows the following definition:
        %   column index of the greater=2.*i-1+column;
        %   column index of the lesser2.*i-column;
        
        %NOTE: if the two values are equal, then the inequality is 0 and
        %everything is correct except the interbranch ab chain is length 0
        %and the equation reduces to the one branch equation, which is
        %valid for arbitrary functionality and therefore still valid in
        %this calculation (on branch of functionality 6).
        
        %we define branch a as the branch connected to the branch point
        %with the lesser of the two bonds connected to the interbranch ab
        %chain
        branchesoftwoa(:,i)=[
            twobranchbonds(xor(row,twobranchbonds(:,2.*i-column)),2.*i-column);     %xor finds the row that contains both a number, and not a number in the second column; the column is the lesser one as defined above
            twobranchchains(xor(row,twobranchbonds(:,2.*i-column)),i) - twobranchbonds(xor(row,twobranchbonds(:,2.*i-column)),2.*i-column);
            twobranchbonds(and(row,twobranchbonds(:,2.*i-column)),2.*i-column)                 %and finds the row that has reacted with the other column, and is of the lesser value.
            ];
        
        %branch b is the branch of the greater of the 2 connected bonds.
        %The only changes are the final row calculation, and the
        %substitution of 2.*i-1+column for the lesser code.
        branchesoftwob(:,i)=[
            twobranchbonds(xor(row,twobranchbonds(:,2.*i-1+column)),2.*i-1+column);     %xor finds the row that contains both a number, and not a number in the second column; the column is the greater one as defined above
            twobranchchains(xor(row,twobranchbonds(:,2.*i-1+column)),i) - twobranchbonds(xor(row,twobranchbonds(:,2.*i-1+column)),2.*i-1+column);
            twobranchchains(and(row,twobranchbonds(:,2.*i-1+column)),i)-twobranchbonds(and(row,twobranchbonds(:,2.*i-column)),2.*i-1+column)                 %and finds the row that has reacted with the other column, and is of the lesser value.
            ];
        
        %branch ab is the branch that connects branch a to b.  It's the
        %only branch with no ends.  It's length is the difference between
        %the greater and lesser value in the row that connects the columns.
        branchesoftwoab(:,i)=[
            twobranchbonds(and(row,twobranchbonds(:,2.*i-1+column)),2.*i-1+column)-twobranchbonds(and(row,twobranchbonds(:,2.*i-column)),2.*i-column)
            ];
        
        %to make the radius of gyration calculation easier, we put all the
        %branches into a single array for part of it.
        alltwobranches=[branchesoftwoa(:,i);branchesoftwob(:,i);branchesoftwoab(:,i)];
        
        twobranchRg(1,i)=twobranchchains(5,i);
        
        %first calculate the single branch equation using the single array
        parta=sum(alltwobranches.^2./2-alltwobranches.^3./(3.*twobranchchains(5,i)),1);
        %then the modification for two branches
        partb=branchesoftwoab(:,i)./twobranchchains(5,i).*sum(branchesoftwoa(:,i),1).*sum(branchesoftwob(:,i),1);
        %and finally put the two parts of the equation together
        twobranchRg(2,i)=(1./twobranchchains(5,i)).*(parta+partb);
        gfactor2(:,i)=(6./twobranchchains(5,i).^2).*(parta+partb);
    end
end
%three branch
%find the chains in Chains that have 3 branches
[~,column3br]=find(Chains(4,:)~=0 & Chains(5,:)~=0);   %because the fourth row of chains is only if there are 3 or more branches

if ~isempty(column3br) %if there are any three branch chains
    threebranchchains=Chains(:,column3br);
    
    %we have to do a bit of manipulation to get the bond columns from the
    %chain
    %one.
    bonds1=3.*column3br-2;
    bonds2=3.*column3br-1;
    bonds3=3.*column3br;
    threebranchcolumns=sort([bonds1 bonds2 bonds3]);
    
    threebranchbonds=Bonds(:,threebranchcolumns);   %we need all 3 of the available branch columns
    
    for i=1:1:size(threebranchchains,2)
        %we eliminate a bug, sofar of unknown source, whereby there is
        %sometimes a column of all zeros counted as a three branch species.
        if sum(threebranchbonds(:,3.*i-2),1)==0 || sum(threebranchbonds(:,3.*i-1),1)==0 || sum(threebranchbonds(:,3.*i),1)==0
            break
        end
        
        bonds=threebranchbonds(:,3.*i-2:3.*i);
        
        %this calculation is relatively simple because they treat the 3
        %branches as one "reference branch" which connects to both other
        %branches.  So we split this into 2 2-branch problems.  The
        %reference branch and a, and the reference branch and b.  This will
        %overcount the reference branch, but we can then remove the excess
        %branches.
        %first find the reference branch, which is the only one that has
        %both its reacted bonds connected to other branches
        refcolumn=sum((bonds>0 & repmat((bonds(:,1)~=0 & bonds(:,2)~=0) | (bonds(:,1)~=0 & bonds(:,3)~=0) | (bonds(:,2)~=0 & bonds(:,3)~=0),1,3)),1)>=2;  %we first boolean AND columns 1,2 and 2,3 and 1,3.  This gives us the matching row in each.  A boolean OR gives us the boolean signature of the reference column (the rows that aren't zero).  We repeat the signature into a 3 column matrix, then AND it with the bonds and sum that in the row direction.  The reference column will be the only one with sum==2
        
        %we have to account for the situation where one chain has all the
        %crosslinks (acts more like a comb polymer than a branch...but in
        %this code we still treat it like as a branch).  We know this has
        %happened when the sum of refcolumn is not 1 (more than 1
        %"refcolumn")
        if sum(refcolumn,2)~=1
            %get the row
            row=sum(bonds~=0,2)==3;
            %find the columns that correspond to the max and min
            [~,maxcol]=max(bonds(row,:));
            [~,mincol]=min(bonds(row,:));
            
            %we can directly get the b branches from this
            branchb(:,1)=[
                bonds(bonds(:,mincol)~=0,mincol);
                threebranchchains(xor(row,bonds(:,mincol)),i)-bonds(xor(row,bonds(:,mincol)),mincol)
                ];
            
            branchb(:,2)=[
                bonds(xor(row,bonds(:,maxcol)),maxcol);
                threebranchchains(xor(row,bonds(:,maxcol)),i)-bonds(xor(row,bonds(:,maxcol)),maxcol);
                threebranchchains(and(row,bonds(:,maxcol)),i)-bonds(and(row,bonds(:,maxcol)),maxcol)
                ];

            %we identify the middle column; this determines the reference
            %crosslink.  We do it by removing the max and min columns from
            %the bonds array
            if mincol<maxcol
                midcol=bonds(:,[1:mincol-1 mincol+1:maxcol-1 maxcol+1:end]);
            elseif mincol>maxcol
                midcol=bonds(:,[1:maxcol-1 maxcol+1:mincol-1 mincol+1:end]);
            end
            
            %we have the midcolumn thus
            brancha=[
                midcol(xor(row,midcol>0)); 
                threebranchchains(xor(row,midcol>0),i)-midcol(xor(row,midcol>0))
                ]';
            
            %last we calculate the connecting branches from the single chain. Use an and
            %to extract the three values and sort them lowest to highest.
            thechain=sort(bonds(and(repmat(row,1,3),bonds)));
            branchab=[thechain(2)-thechain(1)  thechain(3)-thechain(2)];
            
        elseif sum(refcolumn,2)==1
            %the negation of the refcolumn array gives us the nonreference bond
            %columns, while it itself gives us the reference bond column
            refbonds=bonds(:,refcolumn);
            nonrefbonds=bonds(:,~refcolumn);
            
            %we can then use the code for 2 branches on each pair,
            %nonreference1 and reference, and nonreference2 and reference.  We
            %know for a fact which of the 2 contains the reference branch; we
            %can therefore use that instead of the greater/lesser method fo
            %rpart of it.
            %we do is change the variable names and change the indexing from k
            %to i and from 2.*i-1+column to 1+column and 2.*i-column to 2-column
            for k=1:1:2
                twobonds=[refbonds nonrefbonds(:,k)];
                
                %we find the row in the bonds that has no zeros; this row is the
                %cross-chain that contains the interbranch ab chain. We make it a
                %boolean vector, with 1 the row and 0 the others
                row=twobonds(:,1)~=0 & twobonds(:,2)~=0;
                
                %we find the column of the row with the larger #; this helps
                %determine how to match up the ends.  If inequality is 0 (false), then
                %the first column is greater; if it's 1 (true) then the second is
                %greater.
                column=twobonds(row,1) < twobonds(row,2);
                %this allows the following definition:
                %   column index of the greater=1+column;
                %   column index of the lesser=2-column;
                
                %NOTE: if the two values are equal, then the inequality is 0
                %and
                %everything is correct except the interbranch ab chain is length 0
                %and the equation reduces to the one branch equation, which is
                %valid for arbitrary functionality and therefore still valid in
                %this calculation (on branch of functionality 6).
                
                %brancha is the 2 ends off the reference, a branch
                %branchb are both the other branchs
                %branchab is the distance between the branches and branch a
                %we know that the reference chain is in column 1;
                branchb(1:2,k)=[
                    twobonds(xor(row,twobonds(:,2)),2);     %xor finds the row that contains a number in the nonreference column but no number in the reference; this is one of the linear branch ends.  The number and the lenght of chain minus that number are 2 of the 3 arms.
                    threebranchchains(xor(row,twobonds(:,2)),i) - twobonds(xor(row,twobonds(:,2)),2)
                    ];
                %if the second column of the branch is the lesser, then that is
                %the value for the third arm.  If it's the greater, the value
                %is the chain length minus the value.  The same is true for
                %branchb; we can do both in one if loop
                if 1+column == 1
                    branchb(3,k)=twobonds(and(row,twobonds(:,2)),2);
                    brancha(1,k)=threebranchchains(and(row,twobonds(:,1)),i)-twobonds(and(row,twobonds(:,1)),1);
                elseif 1+column == 2
                    branchb(3,k)=threebranchchains(and(row,twobonds(:,2)),i)-twobonds(and(row,twobonds(:,2)),2);
                    brancha(1,k)=twobonds(and(row,twobonds(:,1)),1);
                else
                    disp('This should not print RadiusOfGyration:1')
                    keyboard
                end
                
                %branch ab is the branch that connects branch a to b.  It's the
                %only branch with no ends.  It's length is the difference between
                %the greater and lesser value in the row that connects the columns.
                branchab(:,k)=[
                    twobonds(and(row,twobonds(:,1+column)),1+column)-twobonds(and(row,twobonds(:,2-column)),2-column)
                    ];
            end
        else
            disp('This should not print RadiusOfGyration:2')
            keyboard
        end
        %to make the radius of gyration calculation easier, we put all the
        %branches into a single array for part of it.
        allthreebranches=[brancha;branchb;branchab];
        allthreebranches=reshape(allthreebranches,numel(allthreebranches),1);
        
        threebranchRg(1,i)=threebranchchains(5,i);
        
        %first the all chain part of the equation
        parta=sum(allthreebranches.^2./2-allthreebranches.^3/(3.*threebranchchains(5,i)),1);
        %then we deal with the branch a and branch b chain
        partb=fliplr(branchab)./repmat(threebranchchains(5,i),1,2).*(repmat(sum(brancha,2),1,2)+sum(branchb,1)+branchab).*sum(fliplr(branchb),1);
        %put it all together
        threebranchRg(2,i)=(1./threebranchchains(5,i)).*(parta+sum(partb,2));
        gfactor3(:,i)=6./(threebranchchains(5,i).^2).*(parta+sum(partb,2));
    end
end

%now we've determined the Rg of all the branched chains (up to 3 branches)
%we are going to combine this data with the input original distribution
%into a single Rg vs height

%one-,two-,and three-branchRg have the N in the top row and the Rg in the
%bottom.  We recombine them into one larger array; one branch then 2 branch
%then 3 branch.
BranchRgs=[onebranchRg twobranchRg threebranchRg];
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

FNaxisRg=FNaxisUse;
end





