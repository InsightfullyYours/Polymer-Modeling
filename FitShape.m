function [fitteddep uestimates depestimates]=FitShape(xData,yData,referenceData,range)
%this function fits the control distribution to the data by multiplying by
%a function of the scission extent.  It also calculates the depletion,
%defined as the fraction of the reference that best matches the data.  This
%is a single number equivalent of the exp(-x*u) used to fit u, the scission
%extent.  The depletion, in other words, has no x-dependence.
if size(xData,2)==1
    xData=repmat(xData,1,size(yData,2));
end
if size(referenceData,2)==1
    referenceData=repmat(referenceData,1,size(yData,2));
end

%this for loop cycles through all data and fits to both the exponential and
%to a pure multiple.  The pure multiple is the %depletion: what percent of
%the original chains have reacted
for i=1:size(yData,2)
    r=referenceData(range,i);
    y=yData(range,i);
    x=xData(range,i);
    
    %exponential fit
    start_point = 0.0001;
    options=optimset('Display','none','MaxIter',10000,'MaxFunEvals',10000,'TolX',1E-40,'TolFun',1E-40);
    uestimates(i,:) = fminsearch(@expfun, start_point,options);
    fitted(:,i)=referenceData(:,i).*exp(-xData(:,i).*uestimates(i,1));
    
    %depletion fit
    start_point = rand(1,1);
    options=optimset('Display','none','MaxIter',10000,'MaxFunEvals',10000,'TolX',1E-40,'TolFun',1E-40);
    depestimates(i,:) = fminsearch(@depletionfun, start_point,options);
    fitteddep(:,i)=referenceData(:,i).*depestimates(i,1);
    
end
%function that calculates the sum of squared errors for the exponential
    function sse = expfun(params)
        m = params(1);
        line = r.*exp(-x.*m);
        FittedCurve=line;
        ErrorVector = y - FittedCurve;
        sse = sum(ErrorVector .^ 2);
    end

%function that calculates sum of squared errors for the depletion
    function sse = depletionfun(params)
        m = params(1);
        line = r.*m;
        FittedCurve=line;
        ErrorVector = y - FittedCurve;
        sse = sum(ErrorVector .^ 2);
    end


end