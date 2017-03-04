function fitstruct= KineticFits2(xData,yData)
%this computes the linear (zeroeth order) and exponential (first order) fit
%to kinetic data and reports the values, errors, and residuals
xDataOut=(0:0.01:max(xData,[],1))';
if size(xData,2)==1
    xData=repmat(xData,1,size(yData,2));
end
xDataOut=repmat(xDataOut,1,size(xData,2));

%change it to fractions (Ca/C0) st. the sum of the area of all three
%regions is 1.
yDataratio=yData./repmat(sum(yData,2),1,size(yData,2));

for i=1:size(yData,2)
    imaginarynumber=1;
    while imaginarynumber==1
        y=yData(:,i);
        yrat=yDataratio(:,i);
        x=xData(:,i);
        start_point1 = rand(1,2);
        signs=[-1 1 -1];
        options=optimset('Display','none','MaxIter',10000,'MaxFunEvals',10000,'TolX',1E-40,'TolFun',1E-40);
        estimates(i,:) = fminsearch(@kinfun, start_point1,options);
        Man=estimates(i,1).*y(1,:).^(estimates(i,2)-1).*xDataOut(:,i);
        Man2=estimates(i,1).*y(1,:).^(estimates(i,2)-1).*xData(:,i);
        if estimates(i,2)~=1
            fitted(:,i)=y(1,:)./sum(yData(1,:),2).*(1+signs(i).*(estimates(i,2)-1).*Man).^(1./(1-estimates(i,2)));
        elseif estimates(i,2)==1
            fitted(:,i)=y(1,:)./sum(yData(1,:),2).*exp(-signs(i).*estimates(i,1).*xDataOut(:,i));
        end
        if estimates(i,2)~=1
            linfitted=y(1,:)./sum(yData(1,:),2).*(1+signs(i).*(estimates(i,2)-1).*Man2).^(1./(1-estimates(i,2)));
        elseif estimates(i,2)==1
            linfitted=y(1,:)./sum(yData(1,:),2).*exp(-signs(i).*estimates(i,1).*xData(:,i));
        end
        sse(:,i)=sum((linfitted-yrat).^2);
        residuals(:,i)=yrat-linfitted;
        if isreal(sse(:,i))
            imaginarynumber=0;
        else
            imaginarynumber=1;
        end
    end
end
fitstruct.xDataOut=xDataOut;
fitstruct.fitted=fitted;
fitstruct.estimates=estimates;
fitstruct.sse=sse;
fitstruct.residuals=residuals;
fitstruct.yData=yData;
fitstruct.yDataratio=yDataratio;
fitstruct.xData=xData;


% figure
% subplot(2,1,1)
% errorbar(xData,yDataratio,yDataratio.*.2,'o')
% text(1,1,['ks=',num2str(estimates(:,1)')])
% text(1,.5,['ns=',num2str(estimates(:,2)')])
% hold
% plot(xDataOut,fitted,'-')
% subplot(2,1,2)
% plot(xData,residuals,'o')
% text(1,.01,['SSE=',num2str(sse)])
% subplot(2,1,1)

    function sse = kinfun(params)
        ka = params(1);
        n = params(2);
        Man=ka.*y(1,:).^(n-1).*x;
        if n~=1
            line = y(1,:)./sum(yData(1,:),2).*(1+signs(i).*(n-1).*Man).^(1./(1-n));
        elseif n==1
            line = y(1,:)./sum(yData(1,:),2).*exp(-signs(i).*ka.*x);
        end
        FittedCurve=line;
        ErrorVector = FittedCurve - yrat;
        sse = sum(ErrorVector .^ 2);
    end

end