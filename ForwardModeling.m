function [ResultingDistribution ControlDistwt DeadScission Deadendlinked Dead3arm Deadlong Dead4arm m]= ForwardModeling(uwt,f,InputXaxis,InputDistribution)
%For scission and combination calculations to work, the input InputXaxis
%and InputDistribution must go from low MW (index 1) to high MW (index max)
%outputs weight fractions

% InputXaxis=(1:1:1300)';
% %InputDistribution=round(1000000./(10.*sqrt(2.*pi)).*exp(-(InputXaxis-327).^2./(2.*66^2)));
% InputDistribution=[zeros(1,100) ones(1,200) zeros(1,1000)]';
% 
% uwt=.0005;
% f=.6;

%we insure this by testing for it and then flipping them both if not
if InputXaxis(1,1)>InputXaxis(end,1)
    InputXaxis=flipud(InputXaxis);
    InputDistribution=flipud(InputDistribution);
end

%define the input assuming it is a number fraction.
NIn=InputXaxis;
FNInnum=InputDistribution;

%insure that the distribution is properly scaled
FNInnum=FNInnum./repmat(sum(FNInnum,1),size(FNInnum,1),1);

%assume the input was a number fraction (as all data is) and convert it to
%a weight fraction.
ControlDist=InterchangeNumberandWeight(NIn,FNInnum,1);

%Calculate the distribution of chains that DID react (must be done to
%weight fraction)
ControlDistwt=ControlDist.*exp(-NIn.*uwt);  %calculate the distribution of chains that did NOT react
ChainsReactedwt=ControlDist-ControlDistwt;  %distribution of chains that DID react (formed radical)

%We know what DID react.  We introduce a new variable, f, to determine the
%fraction that DID react that underwent scission, and the corresponding 1-f
%that did NOT undergo scission. 0<=f<=1
%f=0.6;
ChainsWithScissionwt=ChainsReactedwt.*f;
ChainsWithRadicalwt=ChainsReactedwt.*(1-f);

%recalculate u for just scission from f.  This uses a fitting algoritm
start_point = uwt;
options=optimset('Display','none','MaxIter',10000,'MaxFunEvals',10000,'TolX',1E-80,'TolFun',1E-80);
uscission= fminsearch(@expfun,start_point,options);

%calculate the scission
%determine the largest and smalles chain to go to.
maxMonomerNum=max(NIn)-min(NIn)+1;
minMonomerNum=min(NIn);
%create variables
resultDist=zeros(size(FNInnum,1),1);

%use a for loop to cycle through each w(r,u) of equation 13 on the r axis
%this calculates Hamielec Scission
for r=minMonomerNum:1:maxMonomerNum
    % another for loop to cycle through the integral from r to inf in eq13
    count=1;
    integral=0;
    for s=r:1:maxMonomerNum
        integral=integral+(2+uscission.*(s-r))./s.*ChainsWithScissionwt(s,1);
        count=count+1;
    end
    %we divide by count to prevent from unintentionally multiplying by r,
    %as count=r, and the number of contributions to thesum will increase
    %with time. (at r=4, there will be 4 terms of thesum; at r=10 there
    %will be 10 terms of thesum; we divided by count to make a single term
    %and therefore remove the extra r)
    resultDist(r,1)=uscission.*r.*(integral./count).*exp(-uscission.*r);
end

%The areas are not all equal; they do not sum to 1.  We calculate the areas
%and set the area of the scission to that of the missing original
[~,AreaChainsReactedwt,~]=AreaNormalize(NIn,ChainsReactedwt,1:max(NIn)-min(NIn)+1);    %all chains that reacted
[~,AreaChainsWithScissionwt,~]=AreaNormalize(NIn,ChainsWithScissionwt,1:max(NIn)-min(NIn)+1);  %Chains that went to scission based on f
[~,AreaChainsWithRadicalwt,~]=AreaNormalize(NIn,ChainsWithRadicalwt,1:max(NIn)-min(NIn)+1);  %Chains that didn't scission based on 1-f

%give the scission peak the area 1
Hamielecwt=resultDist./repmat(sum(resultDist,1),size(resultDist,1),1);
[~,AreaHamielecwt,~]=AreaNormalize(NIn,Hamielecwt,1:max(NIn)-min(NIn)+1);

%------------------------------------------------------------------
%Now that we have all the distributions and have calculated the Hamielec
%scission, we can calculate the recombinations based on our model

%combination is based upon the number average.  Calculate those.

%convert the weight fractions to number fractions and make sure areas are
%appropriate (f and (1-f) for scission and radicals, respectively)

%Convert results to number fractions from weight fractions
ChainsReactednum=InterchangeNumberandWeight(NIn,ChainsReactedwt,2);
Hamielecnum=InterchangeNumberandWeight(NIn,Hamielecwt,2);
ChainsWithRadicalnum=InterchangeNumberandWeight(NIn,ChainsWithRadicalwt,2);
ControlDistnum=InterchangeNumberandWeight(NIn,ControlDistwt,2);

%recalculate areas to confirm the same
[~,AreaChainsReactednum,~]=AreaNormalize(NIn,ChainsReactednum,1:max(NIn)-min(NIn)+1);
[~,AreaHamielecnum,~]=AreaNormalize(NIn,Hamielecnum,1:max(NIn)-min(NIn)+1);
[~,AreaChainsWithRadicalnum,~]=AreaNormalize(NIn,ChainsWithRadicalnum,1:max(NIn)-min(NIn)+1);
[~,AreaControlDistwt,~]=AreaNormalize(NIn,ControlDistwt,1:max(NIn)-min(NIn)+1);

%A check vector is output.
[AreaChainsReactedwt AreaChainsReactednum AreaChainsWithScissionwt AreaChainsWithScissionwt AreaChainsWithRadicalwt AreaChainsWithRadicalnum AreaHamielecwt AreaHamielecnum];

%---------------------------------------------------------------------
%Calculate the Combination
%we employ another m-file, as there will be multiple combination
%calculations

%calculate the recombination of the control with itself
%(4-arm star)
ContContCombnum=MacroRadicalCombination(NIn,ChainsWithRadicalnum,ChainsWithRadicalnum);

%calculate the recombination of the Hamielec scission with itself
%(end-linking)
HamHamCombnum=MacroRadicalCombination(NIn,Hamielecnum,Hamielecnum);

%calculate the recombination of the control with the Hamielec scission
%(3-arm star)
HamContCombnum=MacroRadicalCombination(NIn,Hamielecnum,ChainsWithRadicalnum);

%-------------------------------------------------------------------
%check all important areas are 1 (of the number fractions)

[~,AreaHHn,~]=AreaNormalize(NIn,HamHamCombnum,1:max(NIn)-min(NIn)+1);
[~,AreaHCn,~]=AreaNormalize(NIn,HamContCombnum,1:max(NIn)-min(NIn)+1);
[~,AreaCCn,~]=AreaNormalize(NIn,ContContCombnum,1:max(NIn)-min(NIn)+1);
[~,AreaFNInn,~]=AreaNormalize(NIn,FNInnum,1:max(NIn)-min(NIn)+1);
[~,AreaHamn,~]=AreaNormalize(NIn,Hamielecnum,1:max(NIn)-min(NIn)+1);

[AreaHHn AreaHCn AreaCCn AreaFNInn AreaHamn];

%------------------------------------------------------------------
%now apply the probabilities.  We know m, which is the area of the chains
%that have reacted (dep) and f (f).  We use these
%to scale the signals.  All signals at this pt have an area of 1 s
%confirmed above.
%AreaChainsReacted is m (dep).  This is value derived from u.

m=-AreaChainsReactedwt;

Pdead=m.*f;
Pendlinked=m.*f.^2;
Pdead2=m.*f;
P3arm=m.*f.*(1-f);
Pdeadlong=m.*(1-f);
P4arm=m.*(1-f).^2;

DeadScission=(Pdead+Pdead2).*Hamielecnum;
Deadendlinked=Pendlinked.*HamHamCombnum;
Dead3arm=P3arm.*HamContCombnum;
Deadlong=Pdeadlong.*FNInnum;
Dead4arm=P4arm.*ContContCombnum;

%this is a number fraction
ResultingDistribution=ControlDistwt+DeadScission+Dead3arm+Deadendlinked+Deadlong+Dead4arm;
ResultingDistribution=ResultingDistribution./repmat(sum(ResultingDistribution,1),size(ResultingDistribution,1),1);

%all the data is usually from high MW to low MW.  This m-file assumed the
%opposite and forced data to low MW to high MW.  We therefore reverse this.
% ResultingDistribution=flipud(ResultingDistribution);
% ControlDistwt=flipud(ControlDistwt);
% DeadScission=flipud(DeadScission);
% Dead3arm=flipud(Dead3arm);
% Deadendlinked=flipud(Deadendlinked);
% Deadlong=flipud(Deadlong);
% Dead4arm=flipud(Dead4arm);

% figure
% hold
% plot(NIn,ControlDist,'b')
% plot(NIn,ControlDistwt,'b')
% plot(NIn,DeadScissionwt,'r')
% plot(NIn,Deadendlinked,':b')
% plot(NIn,Dead3arm,'-.b')
% plot(NIn,Deadlong,'--b')
% plot(NIn,Dead4arm,':r')
% plot(NIn,ControlDistwt+DeadScissionwt+Dead3arm+Deadendlinked+Deadlong+Dead4arm,':k')
%
% text(500,.007,['u = ',num2str(u)])
% text(500,.006,['f = ',num2str(f)])
% text(500,.005,[num2str(m.*100),'% Original Chains Reacted'])
%
% figure
% hold
% plot(NIn,ControlDist,'k')
% plot(NIn,ControlDistnum,'b')
% plot(NIn,Hamielecnum,'r')
% plot(NIn,ChainsWithRadicalnum,':b')
% plot(NIn,ContContCombnum,':b')
% plot(NIn,HamHamCombnum,':b')
% plot(NIn,HamContCombnum,':b')

    function sse = expfun(params)
        u = params(1);
        FittedCurve = ControlDist-ControlDist.*exp(-NIn.*u);
        ErrorVector = FittedCurve - ChainsWithScissionwt;
        sse = sum(ErrorVector .^ 2);
    end

end
