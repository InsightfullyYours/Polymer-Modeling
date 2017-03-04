function [fscissionm, fendlinkm, fthreearmm, ffourarmm, fscissionb, fendlinkb, fthreearmb, ffourarmb, f, Pmono, Pbis]=Probabilities(m,ScissionArea, EndLinkingArea, ThreeArmArea, FourArmArea, ControlArea)
%This function uses the probabilities and the depletion percent (in
%decimal form) to calculate, f, the balance of chains that underwent
%scission vs those that didn't.

%the equations used are shown in the thesis

%we assume that the scission area reflects the dead chains that exist
%after scission events

%we are assuming that the depletion, exp(-ru), is equal to the sum of all
%the processes.  IE mass conservation. So:
%depletion=dead polymer (scission) + endlinking + 3arm stars + 4arm stars
%we calculate this so that each area is a fraction of the depletion

%all areas should sum to 1
totalarea=ScissionArea+EndLinkingArea+ThreeArmArea+FourArmArea+ControlArea
ScissionArea=ScissionArea./totalarea;%.*m;
EndLinkingArea=EndLinkingArea./totalarea;%.*m;
ThreeArmArea=ThreeArmArea./totalarea;%.*m;
FourArmArea=FourArmArea./totalarea;%.*m;
ControlArea=ControlArea./totalarea;%.*m;

%for mono-functional
fscissionm=ScissionArea./(m+m.^2);
fendlinkm=sqrt(EndLinkingArea./m.^2);
fthreearmm=roots([m.^2;-m.^2;ThreeArmArea]);
ffourarmm=roots([m.^2;-2.*m.^2;m.^2-FourArmArea]);

%for bi-functional
fscissionb=ScissionArea./(m+m.^2+m.^4);
fendlinkb=sqrt(EndLinkingArea./(m.^2+m.^4));
fthreearmb=roots([m.^4+m.^2;-1.*(m.^4+m.^2);ThreeArmArea]);
ffourarmb=roots([m.^4+m.^2;-2.*(m.^4+m.^2);m.^4+m.^2-FourArmArea]);

% %choose the root between 0 and 1, if it exists
% for i=1:1:2
%     if (0<=fthreearmm(i) && 1>=fthreearmm(i))
%         fthreearmm=fthreearmm(i);
%     end
%     if (0<=ffourarmm(i) && 1>=ffourarmm(i))
%         ffourarmm=ffourarmm(i);
%     end
%     if (0<=fthreearmb(i) && 1>=fthreearmb(i))
%         fthreearmb=fthreearmb(i);
%     end
%     if (0<=ffourarmb(i) && 1>=ffourarmb(i))
%         ffourarmb=ffourarmb(i);
%     end
% end



%create vectors of the values of the probabilities at all f for display
%purposes
f=(0:.01:1)';

%for mono-functional
Pmono(:,1)=f.*(m+m.^2);
Pmono(:,2)=f.^2.*m.^2;
Pmono(:,3)=polyval([m.^2;-m.^2;0],f);
Pmono(:,4)=polyval([m.^2;-2.*m.^2;m.^2],f);

%for bi-functional
Pbis(:,1)=f.*(m+m.^2+m.^4);
Pbis(:,2)=(m.^2+m.^4).*f.^2;
Pbis(:,3)=-polyval([m.^4+m.^2;-1.*(m.^4+m.^2);0],f);
Pbis(:,4)=polyval([m.^4+m.^2;-2.*(m.^4+m.^2);m.^4+m.^2],f);

FourArmArea./ThreeArmArea
FourArmArea./EndLinkingArea
ThreeArmArea./EndLinkingArea


end

