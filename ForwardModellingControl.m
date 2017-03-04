function ForwardModellingControl
%This function allows the iterative calculations with
%ForwardModelling by using for loops.

%define the input (for pure gaussian), assuming this is weight fraction
NIn=(1300:-1:1)';
FNIn=round(1000000./(10.*sqrt(2.*pi)).*exp(-(NIn-227).^2./(2.*66^2)));
FNIn=FNIn./repmat(sum(FNIn,1),size(FNIn,1),1);

ResultingDistribution=FNIn;
Scaled2=FNIn;

%for f=.1:.1:.9

%for direct
uwts=0.0001:0.0001:0.0025;
f=0.5;
for i=1:1:size(uwts,2)
    uwt=uwts(i)
    %for direct calculation. It outputs number fraction, with data from
    %high MW to low MW
    [ResultingDistribution(:,i) Unreacted DeadScissionwt Deadendlinked Dead3arm Deadlong Dead4arm] = ForwardModeling(uwt,f,NIn,FNIn);
    
    [~,AreaResultingDistribution,ppp]=AreaNormalize(NIn,ResultingDistribution(:,i),1:max(NIn));
    
    [~,AreaDendlinked(i),~]=AreaNormalize(NIn,Deadendlinked./AreaResultingDistribution,1:max(NIn));
    [~,AreaD3arm(i),~]=AreaNormalize(NIn,Dead3arm./AreaResultingDistribution,1:max(NIn));
    [~,AreaD4arm(i),~]=AreaNormalize(NIn,Dead4arm./AreaResultingDistribution,1:max(NIn));
    [~,AreaDScissionwt(i),~]=AreaNormalize(NIn,DeadScissionwt./AreaResultingDistribution,1:max(NIn));
    [~,AreaDlong(i),~]=AreaNormalize(NIn,Deadlong./AreaResultingDistribution,1:max(NIn));
    [~,AreaUnreacted(i),~]=AreaNormalize(NIn,Unreacted./AreaResultingDistribution,1:max(NIn));
end
%for recursion
% %for ppppppp=1:1:uwt./uwt2
%     ppppppp
%     %recursive
%     [ResultingDistribution2 Unreacted2 DeadScissionwt2 Deadendlinked2 Dead3arm2 Deadlong2 Dead4arm2] = ForwardModeling(uwt2,f,NIn,Scaled2);
%
%     [~,AreaResultingDistribution2,ppp]=AreaNormalize(flipud(NIn),ResultingDistribution2,1:max(NIn));
%     Scaled2=ResultingDistribution2./repmat(sum(ResultingDistribution2,1),size(ResultingDistribution2,1),1);
%
%     [~,AreaDendlinked2(i),~]=AreaNormalize(flipud(NIn),Deadendlinked2./AreaResultingDistribution2,1:max(NIn));
%     [~,AreaD3arm2(i),~]=AreaNormalize(flipud(NIn),Dead3arm2./AreaResultingDistribution2,1:max(NIn));
%     [~,AreaD4arm2(i),~]=AreaNormalize(flipud(NIn),Dead4arm2./AreaResultingDistribution2,1:max(NIn));
%     [~,AreaDScissionwt2(i),~]=AreaNormalize(flipud(NIn),DeadScissionwt2./AreaResultingDistribution2,1:max(NIn));
%     [~,AreaDlong2(i),~]=AreaNormalize(flipud(NIn),Deadlong2./AreaResultingDistribution2,1:max(NIn));
%     [~,AreaUnreacted2(i),~]=AreaNormalize(flipud(NIn),Unreacted2./AreaResultingDistribution2,1:max(NIn));
%
%
%     plot(NIn(20:end), FNIn(20:end),'k')
%     %         plot(NIn(20:end), ResultingDistribution(20:end),'b')
%     plot(NIn(20:end), Scaled(20:end),':b')
%     plot(NIn(20:end), Scaled2(20:end),':r')
%     %         plot(NIn(20:end), Deadendlinked(20:end)./AreaResultingDistribution,'r:')
%     %         plot(NIn(20:end), Dead3arm(20:end)./AreaResultingDistribution,'--r')
%     %         plot(NIn(20:end),
%     %         Dead4arm(20:end)./AreaResultingDistribution,'-.r')
%     theuwt(i)=uwt
%     quickly=quickly+uwt2
%     i=i+1
% end
% figure
% hold
% plot(theuwt,AreaD3arm,'r')
% plot(theuwt,AreaDendlinked,'b')
% plot(theuwt,AreaD4arm,'k')
% 

% %end
% figure
% plot(theuwt,AreaD3arm+AreaDendlinked+AreaD4arm+AreaDlong+AreaDScissionwt)
% hold
% plot(theuwt,AreaD3arm+AreaDendlinked+AreaD4arm+AreaDlong,':')
% plot(theuwt,AreaUnreacted+AreaD3arm+AreaDendlinked+AreaD4arm+AreaDlong+AreaDScissionwt,'r')
% plot(theuwt,AreaUnreacted+AreaD3arm+AreaDendlinked+AreaD4arm+AreaDlong,':r')
end