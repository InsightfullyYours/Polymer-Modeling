function CombScis = CombinationScission(NIn,FNIn,fraction,balance,experiment_num,g,numChain)
%This function models the CombScis distribution.
%designed in conjunction with John Zhang
%close all
%define the distribution
%NIn=(1:1:100)';
%FNIn=[(1:1:50)';(50:-1:1)'];
% NIn=(1:1:100)';
% FNIn=round(50000./(10.*sqrt(2.*pi)).*exp(-(NIn-50).^2./(2.*10^2)));
FNIn=round(FNIn.*numChain);
FNIn(FNIn<0)=0;
Naxis=NIn;
FNaxis=FNIn;
g=round(abs(g));
%Define the axes
% Naxis=(1:1:max(NIn))';
%Combaxis1=(1:1:max(Naxis))';
Combaxis=(1:1:max(Naxis))';
% FNaxis=zeros(size(Naxis,1),1);


%populate the axes with F(N) distribution.  This allows distribution d
% for i=1:1:size(NIn,1)
%     FNaxis(Naxis==NIn(i,1))=FNIn(i);
% end

%what percent of the bonds do we want to react?
% fraction=0.1;
% balance=.5;

%create a vector of every chain, by length.  The lengths are repeated once
%for every chain of that length.  So if there are 3 chains of length 1 and
%4 chains of length 2, the resulting vector will be
%  [1;1;1;2;2;2;2].
chains_static = cumsum(FNaxis);
bond_static = cumsum(FNaxis.*Combaxis);
bond=cumsum(FNaxis.*(Combaxis-1));
%max(bond)
%experiment_num=10;
%do the "experiment" experiment_num times and then we can average to get a better picture of
%the distribution.
CombScis_repeat=zeros(max(Combaxis),experiment_num);
%CombScis_evolution=zeros(size(Combaxis,1),round(fraction.*max(chains_static)));
for j=1:1:experiment_num
    FNaxisUse=FNaxis;
    chain_counter=1;
    evolution_counter=1;
    while chain_counter <= round(fraction.*max(chains_static))
        random_reaction=rand(1);
        if random_reaction>balance
            
            bond=cumsum(FNaxisUse.*(Combaxis-1));
            %Where_It_Broke=[NaN,NaN,NaN];
            
            %choose a random bond
            bond_rand = randi([1,max(bond)],1,1);
            %locate the length of chains selected by the bonds
            index=find(bond<bond_rand,1,'last')+1;
            FNaxisUse(index,1)=FNaxisUse(index,1)-1;
            
            bond=cumsum(FNaxisUse.*(Combaxis-1));
            %choose 2 random bonds
            bond_rand = randi([1,max(bond)],1,1);
            index2=find(bond<bond_rand,1,'last')+1;
            
                        FNaxisUse(index2,1)=FNaxisUse(index2,1)-1;
            if index+index2>max(Combaxis)
                chain_counter=chain_counter+2;
                continue
            end
            
            

            
            indices=[index;index2];
            %locate the exact bonds with the chains that combine
            bond_numbers=mod(bond(indices,1)-bond_rand,Combaxis(indices,1)-1);
            
            %bond_number_test=[bond_number_test; bond_number];
            %The size of one piece is the bond_number (since a break of bond 2
            %produces a piece of size 2); the size of the other the chain length
            %minus the bond number+1. Remove 1 chain from the original length and
            %add 1 chain each to each pieces' length.  A bond_number of zero
            %indicates the last bond in a chain; this is equivalent to
            %Naxis(index,1)-1
            if bond_numbers(1,1)==0
                bond_numbers(1,1)=Combaxis(index,1)-1;
            end
            if bond_numbers(2,1)==0
                bond_numbers(2,1)=Combaxis(index2,1)-1;
            end
            
            FNaxisUse(index+index2-g,1)=FNaxisUse(index+index2-g,1)+1;
            
            %now the distribution is adjusted.  We track just the changes
            %separately just in case in another variable
            %Where_It_Broke=[Where_It_Broke; index , bond_number ,
            %Naxis(index,1)-bond_number];
            %CombScis_evolution(:,evolution_counter)=mean([FNaxisUse, CombScis_evolution(:,evolution_counter)],2);
            chain_counter=chain_counter+2;
            evolution_counter=evolution_counter+1;
        elseif random_reaction<balance
            
            bond=cumsum(FNaxisUse.*(Combaxis-1));
            %Where_It_Broke=[NaN,NaN,NaN];
            
            %choose a random bond
            bond_rand = randi([1,max(bond)],1,1);
            
            %locate the index (length of chain) of that bond
            index=find(bond<bond_rand,1,'last')+1;
            %locate the exact bond with the chain that broke
            bond_number=mod(bond(index,1)-bond_rand,Combaxis(index,1)-1);
            %The size of one piece is the bond_number (since a break of bond 2
            %produces a piece of size 2); the size of the other the chain length
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
            %now the distribution is adjusted.  We track just the changes
            %separately just in case in another variable
            %Where_It_Broke=[Where_It_Broke; index , bond_number , Naxis(index,1)-bond_number];
            %CombScis_evolution(:,evolution_counter)=mean([FNaxisUse, CombScis_evolution(:,evolution_counter)],2);
            chain_counter=chain_counter+1;
            evolution_counter=evolution_counter+1;
        elseif random_reaction == balance
            continue
        end
    end
    CombScis_repeat(:,j)=FNaxisUse;
    
end

%----------------------------------------------------------------------
%output, theory, and display code
%figure
CombScis=mean(CombScis_repeat,2);
%chains_static2 = cumsum(CombScis);
%plot(Combaxis,CombScis,'og')
CombScis=CombScis./repmat(sum(CombScis,1),size(CombScis,1),1);  %scale to 1

% %We can then use the input distribution to calculate the theory we've
% %developed; convolution(F(N),F(N)). We do a few powers of N for comparison.
% FNaxis=FNaxis(1:max(NIn),:);
% oN4fitbeg=conv(FNaxis./Naxis.^4,FNaxis);
% oN3fitbeg=conv(FNaxis./Naxis.^3,FNaxis);
% oN2fitbeg=conv(FNaxis./Naxis.^2,FNaxis);
% oNfitbeg=conv(FNaxis./Naxis,FNaxis);
% fitbeg=conv(FNaxis,FNaxis);
% Nfitbeg=conv(Naxis.*FNaxis,FNaxis);
% N2fitbeg=conv(Naxis.^2.*FNaxis,FNaxis);
% N3fitbeg=conv(Naxis.^3.*FNaxis,FNaxis);
% N4fitbeg=conv(Naxis.^4.*FNaxis,FNaxis);
%
% %convolution neglects the first value; the convolution of (1,1) and (1,1)
% %is (2,1).  So we pad the begining with a zero.
% oN4fit=[0;oN4fitbeg];
% oN3fit=[0;oN3fitbeg];
% oN2fit=[0;oN2fitbeg];
% oNfit=[0;oNfitbeg];
% fit1=[0;fitbeg];
% Nfit=[0;Nfitbeg];
% N2fit=[0;N2fitbeg];
% N3fit=[0;N3fitbeg];
% N4fit=[0;N4fitbeg];
%
% %normalize by the magnitude.
% FNaxisNorm=FNaxis./repmat(sum(FNaxis,1),size(FNaxis,1),1);
% oN4fit=oN4fit./repmat(sum(oN4fit,1),size(oN4fit,1),1);
% oN3fit=oN3fit./repmat(sum(oN3fit,1),size(oN3fit,1),1);
% oN2fit=oN2fit./repmat(sum(oN2fit,1),size(oN2fit,1),1);
% oNfit=oNfit./repmat(sum(oNfit,1),size(oNfit,1),1);
% fit1=fit1./repmat(sum(fit1,1),size(fit1,1),1);
% Nfit=Nfit./repmat(sum(Nfit,1),size(Nfit,1),1);
% N2fit=N2fit./repmat(sum(N2fit,1),size(N2fit,1),1);
% N3fit=N3fit./repmat(sum(N3fit,1),size(N3fit,1),1);
% N4fit=N4fit./repmat(sum(N4fit,1),size(N4fit,1),1);
%
% %We display the results
% figure
% plot(Naxis,FNaxisNorm,'bo')
% hold
% plot(Combaxis,CombScis,'r^')
% % plot(Combaxis1,oN4fit,'-b')
% % plot(Combaxis1,oN3fit,'-r')
% % plot(Combaxis1,oN2fit,':b')
% % plot(Combaxis1,oNfit,':r')
% % plot(Combaxis1,fit1,'-.b')
% % plot(Combaxis1,Nfit,'-.r')
% % plot(Combaxis1,N2fit,'--b')
% % plot(Combaxis1,N3fit,'--r')
% % plot(Combaxis1,N4fit,'-k')
% text(100,0.035,['Chains Before: ',num2str(max(chains_static))]);
% text(100,0.03,['Chains After: ',num2str(max(chains_static2))]);
% text(100,0.025,['Experiments: ',num2str(experiment_num)]);
% text(100,0.020,['Fraction: ',num2str(fraction)]);
% text(100,0.015,['Balance: ',num2str(balance)]);
% legend(['Input Distribution       ';'CombScis Product Dist.   ';'1/N^4                    ';'1/N^3                    ';'1/N^2                    ';'1/N                      ';'No N                     ';'N                        ';'N^2                      ';'N^3                      ';'N^4                      '],'Location','EastOutside')
% xlabel('Number of Monomers in Chain')
% ylabel('Number of Chains')
% title('Comparison of Modelling with Theory of CombScis (CombScisModel)')
end

