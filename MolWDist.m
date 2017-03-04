function [Mn Mw PDI]=MolWDist(MWaxis,Amplitude,range)
    
    %Mn calc
    %insure all positive amplitudes
    Amplitude=(Amplitude(range,:));
    %calculate the numerator (Number of chains*MW of chains)
    num=sum(Amplitude.*repmat(MWaxis(range,:),1,size(Amplitude,2)),1);
    %calculate the denominator (Number of chains)
    den=sum(Amplitude,1);
    
    Mn=(num./den)';
    
    %Mw calc
    %calculate the numerator (Number of chains*MW of chains^2)
    num2=sum(Amplitude.*(repmat(MWaxis(range,:),1,size(Amplitude,2))).^2,1);
    %calculate the denominator (Number of chains*MW of chains)
    den2=sum(Amplitude.*repmat(MWaxis(range,:),1,size(Amplitude,2)),1);
    
    Mw=(num2./den2)';
    
    %PDI Calc
    PDI=Mw./Mn;
end
    
    