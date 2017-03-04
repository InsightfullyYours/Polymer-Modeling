function [fitted xaxis estimates Man yData yDataRatio xData]=structstrip(struct)
fitted=struct.fitted;
xaxis=struct.xDataOut;
estimates=struct.estimates;
yData=struct.yData;
yDataRatio=struct.yDataratio;
Man=repmat((estimates(:,1).*yData(1,:)'.^(estimates(:,2)-1))',size(xaxis,1),1).*xaxis;
xData=struct.xData;
end