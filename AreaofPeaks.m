
function [ScissionArea EndLinkingArea ThreeArmArea FourArmArea ControlArea]=AreaofPeaks(xaxis,ScissionPeak,EndLinkingPeak,ThreeArmPeak,FourArmPeak,ControlPeak)
%this function calculates the areas of the input peaks.  It calculates the
%ENTIRE AREA input; this is feasible because the combination peaks are zero
%everywhere except at the peaks.  We artificially impose this on the
%ScissionPeak, which is data, by making all data<0=0.  We also only use the
%scission section, but that is set through the input.

%This calculates area by doing a trapezoidal reiman sum of the
%form:  (1/2)*Q*[f(a)+2f(a+Q)+2f(a+2Q)+2f(a+3Q)+...+f(b)]
%where a is the begining of the xaxis, b the end, and Q the step
%size.  So it's .5*rate*[y(1)+2y(2)+2y(3)+...+2y(b-1)+y(b)]

%Pad the scisison data with zeros at the high MW so it's the same size
%vector as the xaxis and then set the data scission peak to just the 
%positive scission for the area calc
% ScissionPeak=[zeros(size(xaxis,1)-size(ScissionPeak,1),1); ScissionPeak];
% ScissionPeak(ScissionPeak<0,:)=0;

%all xaxis are negative to offset the highMW to lowMW
ScissionArea=trapz(-xaxis,ScissionPeak,1);
EndLinkingArea=trapz(-xaxis,EndLinkingPeak,1);
ThreeArmArea=trapz(-xaxis,ThreeArmPeak,1);
FourArmArea=trapz(-xaxis,FourArmPeak,1);
ControlArea=trapz(-xaxis,ControlPeak,1);

%if any of the areas are near zero, set the area to zero.  That way the
%floating point calculation errors don't propogate
if ScissionArea<1E-10
    ScissionArea=0;
end
if EndLinkingArea<1E-10
    EndLinkingArea=0;
end
if ThreeArmArea<1E-10
    ThreeArmArea=0;
end
if FourArmArea<1E-10
    FourArmArea=0;
end
if ControlArea<1E-10
    ControlArea=0;
end
end