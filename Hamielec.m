function [NumberFractions WeightFractions uwt unum ufromsubtraction]=Hamielec(xaxis,DataDist,ControlDist,range,g1,g2,g3)
%This function models random chain scission from an input distribution.
%
%The code is based upon equations 10 and 13 from
%Triacca, V.J. et al; Polymer Engineering and Science; 33(8) 1993 pg445.
%
%xaxis: the number of monomers per chain
%DataDist: the distribution of polymer after degradation
%ControlDist: the distribution of polymer before degradation
%range: the full width of the distribution for normalization
%g1 to g3: the g-factors for end-to-end, 3-arm, and 4-arm branches.  make
%them all 1 to have no effect on anything.

%we don't want to incorporate species too small or too large into the
%calculations.  We want the distributions to reflect the GPC range.  We
%therefore zero all the values of the inputs outside of the range
DataDist([1:range(1), range(end):size(DataDist,1)]',:)=0;
ControlDist([1:range(1), range(end):size(ControlDist,1)]',:)=0;

%we assume the input distributions are the number fraction.
%to make sure it is normalized correctly, we renormalize the number
%fraction distributions within the range of the GPC columns.  If it's
%already normalized correctly this won't do anything; if not this will
%correct it.  This will also normalize just the range we're looking at
DataDistnum=DataDist./repmat(sum(DataDist,1),size(DataDist,1),1);
ControlDistnum=ControlDist./repmat(sum(ControlDist,1),size(ControlDist,1),1);

%after normalization the sums should be 1
disp('The number fraction normalization; data and control')
[sum(DataDistnum,1) sum(ControlDistnum,1)]

%we then convert to weight fraction
%numerator
DataDistwt=DataDistnum.*xaxis.*repmat(sum(DataDistnum,1),size(DataDistnum,1),1);
ControlDistwt=ControlDistnum.*xaxis.*repmat(sum(ControlDistnum,1),size(ControlDistnum,1),1);
%Area normalize to produce fraction on the same interval.
%This calculates and scales by the denominator
DataDistwt=DataDistwt./repmat(sum(DataDistwt,1),size(DataDistwt,1),1);
ControlDistwt=ControlDistwt./repmat(sum(ControlDistwt,1),size(ControlDistwt,1),1);

%after normalization the sums should be 1
disp('The weight fraction normalization; data and control')
[sum(DataDistwt,1) sum(ControlDistwt,1)]

%we now have the scaled number and weight fraction of the input
%distributions.  We need to determine u.

%first we calculate the u based upon the weight fraction, with a range that
%is chosen to be a 100 points around the peak of the distributions.  This
%guarantees a good fit around the peak, as there fewest changes are there
%due to the large number of original MW species remaining.

%find the index of the max value of the control within the range already provided
[~,indexwt]=max(ControlDistwt(range,:));
[~,indexnum]=max(ControlDistnum(range,:));
%use the index of the max as the center and calculate the fit to get u and
%the percent depletion.
[fittedwt uwt depwt]=FitShape(xaxis,DataDistwt,ControlDistwt,(range(1)+indexwt-50):(range(1)+indexwt+50));
[fittednum unum depnum]=FitShape(xaxis,DataDistnum,ControlDistnum,(range(1)+indexnum-50):(range(1)+indexnum+50));

%the depletion is the constant we multiply the control by to fit the max
%peak of the data;  we subtract from 1 to get the fractional percent that
%actually reacted.
depwt=1-depwt;
depnum=1-depnum;

%estimates is the u; fitted the fitted control peaks.
%we subtract the fitted controls from the data to get the subtraction dist
Subtractionwt=DataDistwt-fittedwt;
Subtractionnum=DataDistnum-fittednum;

%We assume everything bigger than the control peak location is combination
%and everything smaller is scission.  We put them in their own variables.
%NOTE: The data files are stored from the largest MW (first) to the smallest MW (last).
SubtractionwtComb=zeros(size(Subtractionwt,1),1);
SubtractionwtScis=zeros(size(Subtractionwt,1),1);
SubtractionnumComb=zeros(size(Subtractionnum,1),1);
SubtractionnumScis=zeros(size(Subtractionnum,1),1);

SubtractionwtComb(range(1):range(1)+indexwt,:)=Subtractionwt(range(1):range(1)+indexwt,:);
SubtractionwtScis(range(1)+indexwt:range(end),:)=Subtractionwt(range(1)+indexwt:range(end),:);
SubtractionnumComb(range(1):range(1)+indexnum,:)=Subtractionnum(range(1):range(1)+indexnum,:);
SubtractionnumScis(range(1)+indexnum:range(end),:)=Subtractionnum(range(1)+indexnum:range(end),:);

SubtractionwtComb=SubtractionwtComb./repmat(sum(SubtractionwtComb,1),size(SubtractionwtComb,1),1);
SubtractionwtScis=SubtractionwtScis./repmat(sum(SubtractionwtScis,1),size(SubtractionwtScis,1),1);
SubtractionnumComb=SubtractionnumComb./repmat(sum(SubtractionnumComb,1),size(SubtractionnumComb,1),1);
SubtractionnumScis=SubtractionnumScis./repmat(sum(SubtractionnumScis,1),size(SubtractionnumScis,1),1);

%we now have the weight average and number average subtraction plots.  
%As a check, we calculate u from the number average scission

%in control
[Mncontrol,~,~]=MolWDist(xaxis,ControlDistnum,range);
%in scission
[Mnscis,~,~]=MolWDist(xaxis,SubtractionnumScis,range);

ufromsubtraction=(1./Mnscis)-(1./Mncontrol);

%Now that we have both the weight and number fraction, and we have also
%isolated the scission area of the distribution, we use the equations of
%Hamielec and Saito to predict the scission from the control weight
%fraction

%the data is usually input as MW high-low (the first value in the array is
%the highest MW, the last the lowest.  We check if this is true, and if so,
%flip the arrays to make it increase for loop work.
if xaxis(end,1)==1
    xaxisflip=flipud(xaxis);
    ControlDistwtflip=flipud(ControlDistwt);
    ControlDistnumflip=flipud(ControlDistnum);
end

%---------------------------------------------------------
%Hamielec Scission

%determine the largest chain to go to.
maxMonomerNum=xaxis(min(range));
%create variables
resultDist=zeros(size(ControlDistwtflip,1),1);

%use a for loop to cycle through each w(r,u) of equation 13 on the r axis
for r=1:1:maxMonomerNum
    % another for loop to cycle through the integral from r to inf in eq13
    count=1;
    integral=0;
    for s=r:1:maxMonomerNum
        integral=integral+(2+uwt.*(s-r))./s.*ControlDistwtflip(s,1);
        count=count+1;
    end
    resultDist(r,1)=uwt.*r.*(integral./count).*exp(-uwt.*r);
end

%zero the resulting scission distribution everywhere except where the data
%exists.
%flip it back over so it again has the highest MWs first
resultDist=flipud(resultDist);
resultDist([1:range(1), range(end):size(resultDist,1)]',:)=0;

%then calculate the weight and number fraction
resultDistwt=resultDist./repmat(sum(resultDist,1),size(resultDist,1),1); 
resultDistnum=resultDistwt./xaxisflip.*repmat(sum(resultDistwt,1),size(resultDistwt,1),1);
%make the number distribution a number fraction
resultDistnum=resultDistnum./repmat(sum(resultDistnum,1),size(resultDistnum,1),1); 

%------------------------------------------------------------------
%Now that we have all the distributions and have calculated the Hamielec
%scission, we can calculate the recombinations based on our model

%we employ another m-file, as there will be multiple combination
%calculations

%calculate the recombination of the control with itself
%(4-arm star)
ContContCombnum=MacroRadicalCombination(xaxis,ControlDistnumflip,ControlDistnumflip,range);

%calculate the recombination of the Hamielec scission with itself
%(end-linking)
HamHamCombnum=MacroRadicalCombination(xaxis,flipud(resultDistwt),flipud(resultDistwt),range);

%calculate the recombination of the control with the Hamielec scission
%(3-arm star)
HamContCombnum=MacroRadicalCombination(xaxis,ControlDistnumflip,flipud(resultDistwt),range);

%convert the combination number fractions to weight fractions
ContContCombwt=ContContCombnum.*xaxis.*repmat(sum(ContContCombnum,1),size(ContContCombnum,1),1);
ContContCombwt=ContContCombwt./repmat(sum(ContContCombwt,1),size(ContContCombwt,1),1);

HamHamCombwt=HamHamCombnum.*xaxis.*repmat(sum(HamHamCombnum,1),size(HamHamCombnum,1),1);
HamHamCombwt=HamHamCombwt./repmat(sum(HamHamCombwt,1),size(HamHamCombwt,1),1);

HamContCombwt=HamContCombnum.*xaxis.*repmat(sum(HamContCombnum,1),size(HamContCombnum,1),1);
HamContCombwt=HamContCombwt./repmat(sum(HamContCombwt,1),size(HamContCombwt,1),1);

%-----------------------------------------------------
%we have the scission and the combination distributions.  We here fit them
%to the data in amplitude.  This destroys the weight or number fraction,
%but allows us to see how the data overlays

%to set no end-linking, set noendlink=1.  With end-linking set noendlink=0
noendlink=0;
%Full fit
[peakfits, peakestimates]=FitPeaks(xaxis,ControlDistwt,DataDistwt,resultDistwt,HamHamCombwt,HamContCombwt,ContContCombwt,g1,g2,g3,range,noendlink);
%Subtraction fit
%[peakfits, peakestimates]=FitPeaks(xaxis,Subtractionwt,Subtractionwt,resultDistwt,HamHamCombwt,HamContCombwt,ContContCombwt,g1,g2,g3,range,noendlink);


%the peaks are fitted to the subtraction.  We calculate the area of the
%best fit combination and scission peaks.  We use the data subtraction,
%only the >0 parts, as we will be using this to calculate probabilities and
%the data subtraction is obviously the data fit.  This reduces the error.
%The recombination areas are from the best fit peaks.  Every input must be
%large MW first, as the area calc uses a trapezoidal riemann sum using the
%xaxis.

[ScissionArea EndLinkingArea ThreeArmArea FourArmArea ControlArea]=AreaofPeaks(xaxis(:,1),peakestimates(4).*resultDistwt,peakestimates(1).*HamHamCombwt,peakestimates(2).*HamContCombwt,peakestimates(3).*ContContCombwt,peakestimates(5).*ControlDistwt);
%Subtraction:
%[ScissionArea EndLinkingArea ThreeArmArea FourArmArea]=AreaofPeaks(xaxis(:,1),Subtractionwt(range(1)+indexwt:range(end),:),peakestimates(1).*HamHamCombwt,peakestimates(2).*HamContCombwt,peakestimates(3).*ContContCombwt);

%now that we have the areas, we use the probability arguments to determine
%f, the fraction scission.  We assume that depwt, the wt percent fractional 
%depletion that is the equivalent of exp(-x*u), is all the reacted species 
%and is equivalent to m.
[fscissionm, fendlinkm, fthreearmm, ffourarmm, fscissionb, fendlinkb, fthreearmb, ffourarmb, f, Pmono, Pbis]=Probabilities(depwt,ScissionArea, EndLinkingArea, ThreeArmArea, FourArmArea, ControlArea);



%---------------------------------------------------
%plot the data

figure
plot(xaxis,DataDistwt,'.b')
hold
plot(xaxis,peakfits,'ob')
plot(xaxis,peakestimates(5).*ControlDistwt,':r')
plot(xaxis,peakestimates(4).*resultDistwt,':b')
plot(xaxis.*g1,(1-noendlink).*peakestimates(1).*HamHamCombwt,'r')
plot(xaxis.*g2,peakestimates(2).*HamContCombwt,'g')
plot(xaxis.*g3,peakestimates(3).*ContContCombwt,'k')
legend({'Data';'Fit';'Scaled Control';'ScaledScission';'End linked';'3-arm';'4-arm'})
text(100,2E-4,{'End-link g-factor= ',num2str(g1),'  3-arm g-factor= ',num2str(g2),'4-arm g-factor= ',num2str(g3)})
axis([100 1400 -1E-5 3.5E-4])

%-----------------------------------------------------
%Prepare data for output

NumberFractions=[ControlDistnum DataDistnum fittednum Subtractionnum SubtractionnumComb SubtractionnumScis resultDistnum ContContCombnum HamHamCombnum HamContCombnum];
WeightFractions=[ControlDistwt DataDistwt fittedwt Subtractionwt SubtractionwtComb SubtractionwtScis resultDistwt ContContCombwt HamHamCombwt HamContCombwt];


end




