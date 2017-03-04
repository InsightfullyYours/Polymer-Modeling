function PolymerBlend

%Polymer Blend Model
%   written by Nicholas Carbone, Columbia University, 7/06

%Define all internal constants
lam0=2/3;       %probability of staying in same layer
lam1=1/6;       %probability of moving +- 1 layer

%User Inputs and defaults
%# of layers in the model
n=input(['\nHow many layers do you wish the model to have?\n',...
    'Default is 100.\n']);
if isempty(n)
    n=100;
end

%# segments in functional polymer
r=input(['\nHow many segments in the functionalized polymer?\n',...
    'Default is 20.\n']);
if isempty(r)
    r=20;
end

%# segments in nonfunctional polymer
rprime=input(['\nHow many segments in the nonfunctionalized polymer?\n',...
    'Default is 20.\n']);
if isempty(rprime)
    rprime=20;
end

%#initial guesses
x=input(['\nPlease enter the initial guesses for free segment probability\n',...
    'in MATLAB vector format.\n','Default is [1.2 1 ... 1].\n']);
%check if x vector has one value per layer; if not, see if the transpose
%does; if still not, display message and switch to default
check=size(x);
if isempty(x)
    x=ones(n,1);
    x(1)=1.2;
elseif check(1)==n & check(2)==1
elseif check(1)==1 & check(2)==n
    x=x';
else
    disp(sprintf('x must be an (n x 1) or (1 x n) matrix.  Default will be used.'));
    x=ones(n,1);
    x(1)=1.2;
end

%initial chi(s) values for the two end groups
chi=input(['\nPlease enter the chi(s) values for the first and second\n',...
    'functional groups as a column vector.  Positive values repel from\n',...
    'the surface.  Default is [-2 -2].\n']);
if isempty(chi)
    chi=[2 2];
else
    chi=-chi; %this sign is reversed to mesh with Theodorou's conventions (see ref)
end

%volume fraction of components
phi=input(['\nWhat is the volume fraction of the functional polymer?\n',...
    'Default is 0.50.\n']);
phip=1-phi;      %defines volume fraction of nonfunc polymer
if isempty(phi)
    phi=0.5;
    phip=0.5;
elseif phi <= 1 & phi >= 0
else
    disp(sprintf('phi must be between 0 and 1.  Default will be used'));
    phi=0.5;
    phip=0.5;
end

%all the values are entered. Now, let us run a nonlinear least squares
%optimization with the initial probability guess (x) as the initial values.
%to determine the correct values.
options=optimset('TolX',1e-12,'TolFun',1e-12,'MaxIter',1000);
xmin=lsqnonlin(@blend,x,[],[],options); %outputs the optimum values of x

%Because the non-linear least squares calculation above only returns the
%optimum x values, we have to rerun the calculation with those values in
%order to extract the data we want: segment probabilities in each layer.
%minout is the value of the function that was minimized by lsqnonlin. It 
%should be driven close to zero.  func has the segment probabilities of 
%the functional polymer, nonfunc has the probabilities of the 
%nonfunctional, and intermed is a giant matrix with the final values of 
%most of the intermediate matrices used in the calculation.
[minout func nonfunc intermed]=blend(xmin);

vprob=[phi.*(func./r) phip.*(nonfunc./rprime)]; %probabilities accounting for volume fraction of functional and nonfunctional polymer
vsum=sum(vprob,2); %sum of all probabilities

FinalProb=[(1:1:n)' vprob./repmat(vsum,1,4)]; %normalized probabilities for all polymers.

[FinalProb sum(FinalProb(:,2:5),2)] %this shows FinalProb and a extra column that sums all the probs to confirm they equal 1



function [minimization par parnf intermediates]=blend(x)
    %We initialize the matrices.  pplus corresponds to p+ in Theodorou, 
    %pminus is p-, w is w+ and w-.  
    %See Theodorou, Macromolecules, 21(5) 1988 pg. 1411(1415)
    pplus=[x zeros(n,r-1)];   %first column is guess for layer probabilities, all others initialized zero
    pminus=[x zeros(n,r-1)];

    %We must take the surface interaction into account.  This is only relevent
    %for the functional group locations when the functional group is at the
    %surface.  Functional groups are at beinning and end of polymer, so first
    %and last column, row one (surface) must be modified.  Since the last
    %column is not yet present, we do not edit it here.
    pplus(1,1)=pplus(1,1).*exp(chi(1));
    pminus(1,1)=pminus(1,1).*exp(chi(2));

    %w is the matrix of layer changing probablilities.  It corresponds to w+ and
    %w- in Theodorou.  This can be done because it does not actually changein
    %the forward or reverse direction because we consider the polymer chain
    %flipped.  In other words, we consider abcde for forward and edcba for
    %backward, while Theodorou looks at abcde only, and varies direction.
    w=diag(x*lam1,1)+diag(x*lam1,-1);  %creates an off-center diagonal matrix of layer change probabilities
    w=w(1:length(w)-1,1:length(w)-1);  %cuts off last row and column to create correct size (n x n)
    w=diag(x*lam0)+w;              %generates the final matrix with 3 diagonals and the rest zeroes


    %Now we calculate the layer probabilities while moving along the chain in
    %both a positive and negative direction.  Each column of pplus and pminus
    %is a segment of the chain.  From the initial guess of the layer the first
    %segment will be in (x) we multiply by the matrix of probabilities that the
    %segment is in the same layer or moved (w) to find the probability the
    %second segment is in each layer (pplus(:,2)).  We do this to the end of
    %the chain (r).  The last segment has a interaction energy as well, so
    %after computing the pplus and pminus matrices, we modify the surface
    %(first row) probability by multiplying by the interaction parameter.
    %Below calculates equations (60) and (62) in Theodorou
    for i=1:r-1
        if i==r-1 %if it is the last segment, we must account for last segment surface interation
            pplus(1:2,i)=pplus(1:2,i).*exp(chi(2));%these two add the surface interaction in the second to last segment so when the equations are iterated the seg. interaction will show up correctly in the last segment
            pminus(1:2,i)=pminus(1:2,i).*exp(chi(1));
            pplus(:,i+1)=w*pplus(:,i);   %calculate probability of last segment from modified second to last segment
            pminus(:,i+1)=w*pminus(:,i);
            pplus(1:2,i)=pplus(1:2,i).*exp(-chi(2));%restores second to last segment by multiplying by the inverse surface interaction
            pminus(1:2,i)=pminus(1:2,i).*exp(-chi(1));
        else
        pplus(:,i+1)=w*pplus(:,i);   %calculate probability of next segment from previous for all but last segment
        pminus(:,i+1)=w*pminus(:,i);
        end
    end
    

    %We want to look at 3 things: the probability functional group one is in
    %each layer, the probability a bulk segment is in each layer, and the
    %probability functional group two is in each layer.  
    %We manipulate pminus a bit first so that we can use matrix multiplication
    %to add up our overall bulk probabilities in each layer for us.
    pminusflip=fliplr(pminus); %flips pminus around a vertical axis
    
    %now we create the matrix of probability sums for the functional polymer.
    %We take the first segment probabilities in all layers in the positive
    %direction multiplied by the same in the negative direction:
    %(pplus(:,1).*pminusflip(:,1))...the fliplr command makes segment 1 column
    %one in pminusflip.  This is the first column.  The rest of the columns are
    %the diagonal of the matrix multiplication of the bulk parts of the
    %positive (pplus) and negative (pminusflip) directions.  We only want the
    %diagonal because that is where the same segments are multiplied (we want
    %B1*b1+C1*c1+...+Y1*y1, not B1*b2+...).  The last column is the last
    %segment probability.  The first segment and bulk are all normalized with the
    %first segment probability, while the last segment probability is
    %normalized with the last segment probability
    par=[pplus(:,1).*pminusflip(:,1) ... %First segment probability
        diag(pplus(:,2:r-1)*pminusflip(:,2:r-1)') ... %bulk probability
        pplus(:,r).*pminusflip(:,r)]; %last segment probability
    par=[par(:,1:2)./repmat(pplus(:,1),1,2) par(:,3)./pminus(:,1)]; %normalize as discussed above

    %Now we look at the nonfunctional polymer.  We follow the same proceedure,
    %only without using surface interaction.  We also use the same initial
    %guesses for the first segment.
    nfplus=[x zeros(n,rprime-1)]; %initialize probability matrices
    nfminus=[x zeros(n,rprime-1)];
    for i=1:rprime-1
        nfplus(:,i+1)=w*pplus(:,i);   %calculate probability of next segment from previous for all but last segment
        nfminus(:,i+1)=w*pminus(:,i);
    end
    nfminusflip=fliplr(nfminus); %flip nfminus for matrix multiplication
    
    parnf=diag(nfplus*nfminusflip')./nfplus(:,1); %probabilities of nonfunctional polymer
    
    minimization=1-phi.*(sum(par,2)./r)-phip.*(sum(parnf,2)./r); %this is based on phi+phip=1...minimization=0 is optimal
    intermediates=[x pplus pminus nfplus nfminus w]; %outputs intermediate matrices for debugging
    end
end
