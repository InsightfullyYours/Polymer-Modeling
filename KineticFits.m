function [linstruct expstruct expstruct2 expstruct3 xDataOut] = KineticFits(xData,yData)
%this computes the linear (zeroeth order) and exponential (first order) fit
%to kinetic data and reports the values, errors, and residuals
xDataOut=(0:0.01:max(xData,[],1))';
if size(xData,2)==1
    xData=repmat(xData,1,size(yData,2));
end
xDataOut=repmat(xDataOut,1,size(xData,2));



%yData=abs(yData);
for i=1:size(yData,2)
    y=yData(:,i);
    x=xData(:,i);
    start_point1 = rand(1,2);
    start_point2 = rand(1,3);
    options=optimset('Display','none','MaxIter',10000,'MaxFunEvals',10000,'TolX',1E-40,'TolFun',1E-40);
    linestimates(i,:) = fminsearch(@linfun, start_point1,options);
    expestimates(i,:) = fminsearch(@expfun, start_point2,options);
    expestimates2(i,:) = fminsearch(@expfun2, start_point2,options);
    expestimates3(i,:) = fminsearch(@expfun3, start_point2,options);
    linfitted(:,i)=linestimates(i,1).*xData(:,i)+linestimates(i,2);
    linsse(:,i)=sum((linfitted(:,i)-y).^2);
    linresiduals(:,i)=y-linfitted(:,i);
    expfitted(:,i)=1-expestimates(i,1).*exp(expestimates(i,2).*xDataOut(:,i))+expestimates(i,3);
    expfitted5(:,i)=1-expestimates(i,1).*exp(expestimates(i,2).*xData(:,i))+expestimates(i,3);
    expsse(:,i)=sum((expfitted5(:,i)-y).^2);
    expresiduals(:,i)=y-expfitted5(:,i);
    expfitted2(:,i)=expestimates2(i,1).*exp(expestimates2(i,2).*xDataOut(:,i))+expestimates2(i,3);
    expfitted5(:,i)=expestimates2(i,1).*exp(expestimates2(i,2).*xData(:,i))+expestimates2(i,3);
    expsse2(:,i)=sum((expfitted5(:,i)-y).^2);
    expresiduals2(:,i)=y-expfitted5(:,i);
    expfitted3(:,i)=-expestimates3(i,1).*exp(expestimates3(i,2).*xDataOut(:,i))+expestimates3(i,3);
    expfitted5(:,i)=-expestimates3(i,1).*exp(expestimates3(i,2).*xData(:,i))+expestimates3(i,3);
    expsse3(:,i)=sum((expfitted5(:,i)-y).^2);
    expresiduals3(:,i)=y-expfitted5(:,i);
end
linstruct.fitted=linfitted;
linstruct.estimates=linestimates;
linstruct.sse=linsse;
linstruct.residuals=linresiduals;
expstruct.fitted=expfitted;
expstruct.estimates=expestimates;
expstruct.sse=expsse;
expstruct.residuals=expresiduals;
expstruct2.fitted=expfitted2;
expstruct2.estimates=expestimates2;
expstruct2.sse=expsse2;
expstruct2.residuals=expresiduals2;
expstruct3.fitted=expfitted3;
expstruct3.estimates=expestimates3;
expstruct3.sse=expsse3;
expstruct3.residuals=expresiduals3;

figure
subplot(2,1,1)
plot(xData,yData,'o')
hold
plot(xData,linfitted,'-')
plot(xDataOut,expfitted,':')
plot(xDataOut,expfitted2,'.-')
plot(xDataOut,expfitted3,'--')
subplot(2,1,2)
plot(xData,linresiduals,'o')
hold
plot(xData,expresiduals,'^')
text(1,.01,['lin SSE=',num2str(linsse(:,i))])
text(1,.01,['exp SSE=',num2str(expsse(:,i))])



    function sse = linfun(params)
        m = params(1);
        b = params(2);
        line = m.*x+b;
        FittedCurve=line;
        ErrorVector = FittedCurve - y;
        sse = sum(ErrorVector .^ 2);
    end
    function sse = expfun(params)
        m = params(1);
        b = params(2);
        c = params(3);
        line = 1-m.*exp(b.*x)+c;
        FittedCurve=line;
        ErrorVector = FittedCurve - y;
        sse = sum(ErrorVector .^ 2);
    end
    function sse = expfun2(params)
        m = params(1);
        b = params(2);
        c = params(3);
        line = m.*exp(b.*x)+c;
        FittedCurve=line;
        ErrorVector = FittedCurve - y;
        sse = sum(ErrorVector .^ 2);
    end
    function sse = expfun3(params)
        m = params(1);
        b = params(2);
        c = params(3);
        line = -m.*exp(b.*x)+c;
        FittedCurve=line;
        ErrorVector = FittedCurve - y;
        sse = sum(ErrorVector .^ 2);
    end

end