function [fitted estimates]=FitPeaks(xData,controlData,yData,Scission,convolution1,convolution2,convolution3,g1,g2,g3,range,noendlink)
global j
global sensitivity
global qqqqq
if size(xData,2)==1
    xData=repmat(xData,1,size(yData,2));
end

%first we change the convolution in the xaxis using the g-factors
%find the index of the max
[~,index1]=max(convolution1);
%we find the MW at g1
mw1=xData(index1,1);
%we scale the MW by g1 and round to the nearest whole number
newmw1=round(mw1.*g1);
%we find the index of the new mw
newindex1=find(xData(:,1)==newmw1,1);
%calculate the shift
shift1=newindex1-index1;
%add that shift in zeros to the beginning of the convolution
convolution1=[zeros(shift1,1);convolution1(1:end-shift1,:)];

%we do it again for the second g-factor
%find the index of the max
[~,index2]=max(convolution2);
%we find the MW at g1
mw2=xData(index2,1);
%we scale the MW by g1 and round to the nearest whole number
newmw2=round(mw2.*g2);
%we find the index of the new mw
newindex2=find(xData(:,1)==newmw2,1);
%calculate the shift
shift2=newindex2-index2;
%add that shift in zeros to the beginning of the convolution
convolution2=[zeros(shift2,1);convolution2(1:end-shift2,:)];

%we do it again for the second g-factor
%find the index of the max
[~,index3]=max(convolution3);
%we find the MW at g1
mw3=xData(index3,1);
%we scale the MW by g1 and round to the nearest whole number
newmw3=round(mw3.*g3);
%we find the index of the new mw
newindex3=find(xData(:,1)==newmw3,1);
%calculate the shift
shift3=newindex3-index3;
%add that shift in zeros to the beginning of the convolution
convolution3=[zeros(shift3,1);convolution3(1:end-shift3,:)];

start_point_testing=rand(1,5);
% [
% %     0.3171    0.6463    0.3404    0.1386    0.2511;
% %     0.9502    0.7094    0.5853    0.1493    0.6160;
% %     0.0344    0.7547    0.2238    0.2575    0.4733;
% %     0.4387    0.2760    0.7513    0.8407    0.3517;
% %     0 0 0 0 0;
% %     0 0 1 0 0;
% %     1 1 0 0 0;
% %     0 1 1 0 0;
% %     1 0 1 1 0;
% %     1 1 1 0 0;
% %     0 1 1 1 0;
% %     1 1 0 1 1;
% %     1 1 1 1 1;
% %     100 100 100 100 100;
% %     100000 100000 100000 100000 100000;
% %     %0.7537    0.0052	0.0028	0.0020	0.9892
% %     0 0 0 0 0
%     ];
    
    

for i=1:size(yData,2)
    for startpointtestingindex=1:1:size(start_point_testing,1)
        startpointtestingindex
        
        [~,indexwt]=max(controlData(range,:));
        fitted=FitShape(xData,yData,controlData,(range(1)+indexwt-50):(range(1)+indexwt+50));
        Subtractionwt=yData-fitted;
        
        
        %range(1):6076
        y2=yData(range(1):6076,i);
        x2=xData(range(1):6076,i);
        start_point2 = start_point_testing(startpointtestingindex,:);
        start_point2 = start_point2(1,1:3)';
        options=optimset('Display','none','MaxIter',10000,'MaxFunEvals',10000,'TolX',1E-80,'TolFun',1E-80);
        [estimates2(i,:) sse2(i,:)]= fminsearch(@expfun2,start_point2,options);
        estimates2=abs(estimates2);
        startpointtestingoutput(startpointtestingindex,:)=[estimates2(i,:) 0 0];
        sses2(startpointtestingindex,:)=sse2(i,:);
        fitted2(:,i)=(1-noendlink).*estimates2(i,1).*convolution1+estimates2(i,2).*convolution2+estimates2(i,3).*convolution3;%+0.*estimates(4).*Scission+estimates(4).*controlData;
        estimates2(i,1)=estimates2(i,1).*(1-noendlink);
        
                %range(1):6076
        y3=yData(6476:range(end),i);
        x3=xData(6476:range(end),i);
        start_point3 = start_point_testing(startpointtestingindex,:);
        start_point3 = start_point3(1,4:5)';
        options=optimset('Display','none','MaxIter',10000,'MaxFunEvals',10000,'TolX',1E-80,'TolFun',1E-80);
        [estimates3(i,:) sse3(i,:)]= fminsearch(@expfun3,start_point3,options);
        estimates3=abs(estimates3);
        startpointtestingoutput(startpointtestingindex,:)=[estimates2(i,:) estimates3(i,:)];
        sses3(startpointtestingindex,:)=sse3(i,:);
        fitted3(:,i)=estimates3(1).*Scission+estimates3(2).*controlData;
        
        y4=yData(6176:6476,i);
        x4=xData(6176:6476,i);
        start_point4 = start_point_testing(startpointtestingindex,:);
        start_point4 = start_point4(1,4:5)';
        options=optimset('Display','none','MaxIter',10000,'MaxFunEvals',10000,'TolX',1E-80,'TolFun',1E-80);
        [estimates4(i,:) sse4(i,:)]= fminsearch(@expfun4,start_point4,options);
        estimates4=abs(estimates3);
        startpointtestingoutput(startpointtestingindex,:)=[estimates2(i,:) estimates3(i,1) estimates4(2)];
        sses4(startpointtestingindex,:)=sse3(i,:);
        fitted4(:,i)=estimates4(1).*Scission+estimates4(2).*controlData;
        
        
        %range(1):6076
        y=yData(range,i);
        x=xData(range,i);
        %start_point = [estimates2 start_point_testing(startpointtestingindex,4:5)]';  
        %start_point = start_point_testing(startpointtestingindex,:);              %rand(1, 3)
        start_point = startpointtestingoutput(startpointtestingindex,:);
        options=optimset('Display','none','MaxIter',10000,'MaxFunEvals',10000,'TolX',1E-80,'TolFun',1E-80);
        [estimates(i,:) sse(i,:)]= fminsearch(@expfun,start_point,options);
        estimates=abs(estimates);
        startpointtestingoutput(startpointtestingindex,:)=estimates(i,:);
        sses(startpointtestingindex,:)=sse(i,:);
        fitted(:,i)=(1-noendlink).*estimates(i,1).*convolution1+estimates(i,2).*convolution2+estimates(i,3).*convolution3+estimates(4).*Scission+estimates(5).*controlData;
        estimates(i,1)=estimates(i,1).*(1-noendlink);
        
    end
    j
    crazy=[start_point_testing startpointtestingoutput sses.*1000000000];
    sensitivity{qqqqq}=crazy;
    startpointtestingoutput
end
disp('done')




%function that calculates the sum of squared errors
    function sse = expfun(params)
        m = abs(params(1));
        b = abs(params(2));
        q = abs(params(3));
        w = abs(params(4));
        r = abs(params(5));
        line = m.*convolution1.*(1-noendlink)+b.*convolution2+q.*convolution3+w.*Scission+r.*controlData;
        %range(1):6076
        %FittedCurve=line([range(1):6176,6476:range(end)],1);
        %ErrorVector = FittedCurve - yData([range(1):6176,6476:range(end)],:);
        FittedCurve=line(range,1);
        ErrorVector = FittedCurve - yData(range,:);
        sse = sum(ErrorVector .^ 2);
    end

    function sse2 = expfun2(params)
        m = abs(params(1));
        b = abs(params(2));
        q = abs(params(3));
     %   w = abs(params(4));
     %   r = abs(params(4));
        line = m.*convolution1.*(1-noendlink)+b.*convolution2+q.*convolution3;%+0.*w.*Scission+0.*r.*controlData;
        FittedCurve=line(range(1):6076,1);
        ErrorVector = FittedCurve - y2;
        sse2 = sum(ErrorVector .^ 2);
    end

    function sse3 = expfun3(params)
%         m = abs(params(1));
%         b = abs(params(2));
%         q = abs(params(3));
        w = abs(params(1));
        r = abs(params(2));
        line = w.*Scission+r.*controlData;
        FittedCurve=line(6476:range(end),1);
        ErrorVector = FittedCurve - y3;
        sse3 = sum(ErrorVector .^ 2);
    end

    function sse4 = expfun4(params)
%         m = abs(params(1));
%         b = abs(params(2));
%         q = abs(params(3));
        w = abs(params(1));
        r = abs(params(2));
        line = w.*Scission+r.*controlData;
        FittedCurve=line(6176:6476,1);
        ErrorVector = FittedCurve - y4;
        sse4 = sum(ErrorVector .^ 2);
    end
end