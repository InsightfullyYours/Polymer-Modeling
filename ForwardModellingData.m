function ForwardModellingData(Naxis,yData,range)
%This function calculates the best fits for data.  It firsts calculates the
%u value for each data pair (assuming the first of the yData columns is the
%control), and then fits the f to minimize the least squares.  Naxis is
%assumed to be the same for all yData.

%we need as much processing as possible.  Use the matlabpool parallel
%processing.
if matlabpool('size') > 0
    matlabpool close
end
matlabpool open


%we extract just the data within the range
xaxis=Naxis(range,:);
yaxis=yData(range,:);

%we flip the data and the xaxis so they go from low MW to high MW
if Naxis(1,1)>Naxis(end,1)
    xaxisall=flipud(xaxis);
    yaxisall=flipud(yaxis);
end

%then we insure that it is normalized to an area of 1
[yaxis,~,AreaAfter]=AreaNormalize(xaxis,yaxis,1:size(yaxis,1));
AreaAfter      %output control variable: should all be 1

%now that the data has been prepared, we need to determine the u value for
%each column of the yaxis.  We assume the first column is the control data
%(minimume 2 columns, in other words).
xaxis=xaxisall(:,1);
yaxisnumcon=yaxisall(:,1);
yaxiswtcon=InterchangeNumberandWeight(xaxis,yaxisall(:,1),1);
window=100; %used in the fitting functions at the very bottom.
for i=2:1:size(yData,2)
    
    %assume the input was a number fraction (as all data is) and convert it to
    %a weight fraction.  Also define the working variables
    yaxisnum=yaxisall(:,i);
    yaxiswt=InterchangeNumberandWeight(xaxis,yaxisnum,1);
    
    %calculate u for the data (compare the control to the data, calculate
    %the u that shifts it to that best.  This is the total degradation u,
    %not the scission degradation u.
    start_point = .0005;
    options=optimset('Display','none','MaxIter',10000,'MaxFunEvals',10000,'TolX',1E-80,'TolFun',1E-80);
    uwt(:,i)= fminsearch(@expfun,start_point,options);
    
    ControlDist=InterchangeNumberandWeight(xaxis,yaxisnumcon,1);
    ControlDistwt=ControlDist.*exp(-xaxis.*uwt(:,i));  %calculate the distribution of chains that did NOT react
    ChainsReactedwt=ControlDist-ControlDistwt;
    
    [~,AreaChainsReactedwt,~]=AreaNormalize(xaxis,ChainsReactedwt,1:max(xaxis)-min(xaxis)+1)
    
    %now we fit f, the fraction of scission vs combination (f=0 means 100%
    %combination).  The u is COMPLETELY determined by comparing the data
    %and control.
    start_point = .5;
    options=optimset('Display','iter','MaxIter',50,'MaxFunEvals',50,'TolX',1E-40,'TolFun',1E-40);
    fscis(:,i)= fminsearch(@expfun2,start_point,options);
    
    [ResultingDistribution(:,i) Unreacted DeadScission Deadendlinked Dead3arm Deadlong Dead4arm m] = ForwardModeling(uwt(:,i),fscis(:,i),xaxis,yaxisnumcon);
    [~,AreaResultingDistribution,ppp]=AreaNormalize(xaxis,ResultingDistribution(:,i),1:max(xaxis)-min(xaxis)+1);
    [AreaResultingDistribution ppp]
    
    ControlDistwt=yaxisnumcon.*exp(-xaxis.*uwt(:,i));
    
    figure
    hold
    %plot(xaxis,Unreacted,'.k')
    plot(xaxis,yaxisnum,'b')
    plot(xaxis,yaxisnumcon,'r')
    %plot(xaxis,DeadScissionwt,'r')
    %plot(xaxis,Deadendlinked,':b')
    %plot(xaxis,Dead3arm,'-.b')
    %plot(xaxis,Deadlong,'--b')
    %plot(xaxis,Dead4arm,':r')
    plot(xaxis,ControlDistwt+DeadScission+Dead3arm+Deadendlinked+Deadlong+Dead4arm,'k')
    
    text(800,.007,['u = ',num2str(uwt(:,i))])
    text(800,.006,['f = ',num2str(fscis(:,i))])
    text(800,.005,[num2str(m.*100),'% Original Chains Reacted'])
    
    
    
    
    
    %
    %     %for direct calculation
    %     [ResultingDistribution(:,i) Unreacted DeadScissionwt Deadendlinked Dead3arm Deadlong Dead4arm] = ForwardModeling(uwt(:,i),f(:,i),xaxis,yaxis);
    %
    %     [~,AreaResultingDistribution,ppp]=AreaNormalize(flipud(NIn),ResultingDistribution(:,i),1:max(NIn));
    %     Scaled(:,i)=ResultingDistribution(:,i)./repmat(sum(ResultingDistribution(:,i),1),size(ResultingDistribution(:,i),1),1);
    %
    %     [~,AreaDendlinked(i),~]=AreaNormalize(flipud(NIn),Deadendlinked./AreaResultingDistribution,1:max(NIn));
    %     [~,AreaD3arm(i),~]=AreaNormalize(flipud(NIn),Dead3arm./AreaResultingDistribution,1:max(NIn));
    %     [~,AreaD4arm(i),~]=AreaNormalize(flipud(NIn),Dead4arm./AreaResultingDistribution,1:max(NIn));
    %     [~,AreaDScissionwt(i),~]=AreaNormalize(flipud(NIn),DeadScissionwt./AreaResultingDistribution,1:max(NIn));
    %     [~,AreaDlong(i),~]=AreaNormalize(flipud(NIn),Deadlong./AreaResultingDistribution,1:max(NIn));
    %     [~,AreaUnreacted(i),~]=AreaNormalize(flipud(NIn),Unreacted./AreaResultingDistribution,1:max(NIn));
    %
    %     %for recursion
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

    function sse = expfun(params)
        %we use the window to make the fit better by focusing on the area
        %of the peak
        u = params(1);
        FittedCurve = yaxiswtcon.*exp(-xaxis.*u);
        [~, index]=max(FittedCurve);
        ErrorVector = FittedCurve([index-window:index+window],:) - yaxiswt([index-window:index+window],:);
        sse = sum(ErrorVector .^ 2);
    end

    function sse = expfun2(params)
        %we use the window to make the fit better by focusing on the area
        %not around the peak; just the areas that have changed.
        f = params(1);
        if f>1 || f<0
            f=.5;
        end
        Output = ForwardModeling(uwt(:,i),f,xaxis,yaxisnumcon);
        FittedCurve = flipud(Output);
        [~, index]=max(yaxisnumcon);
        ErrorVector = FittedCurve([1:index-window, index+window:end],:) - yaxisnum([1:index-window, index+window:end],:);
        sse = sum(ErrorVector .^ 2);
    end
matlabpool close
end