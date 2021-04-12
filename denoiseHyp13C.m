function [ denoised, coredim, nT, detailmap ] = denoiseHyp13C( kspace, nsmap, target_value )
%Script for determining a satisfactory denoising that satisfies the
%inequality in expression [3] from the paper while maintaining low rank in
%all dimensions.
%This code matches the pseudocode and description given in supplementary
%methods 3. 

%This function calls the denoising and sum of SSnoise calculation function:
% - tridenoise
%   - ftingrti
%
%INPUTS
%kspace - a three dimensional matrix of dimensions 32x32x128 consisting of
%raw acquired MRSI data in kspace and time dimensions
%nsmap - a 32x32 matrix indicating the kspace points that have significant
%signal
%target_value - the expected sum of SSnoise for all the kspace points in
%nsmap - any found denoising must have a summed residual sum of squares
%smaller than this value.
%
%OUTPUTS
%denoised - the final denoised kspace data
%coredim - a vector giving the dimensions of the reduced rank used in TRI
%nT - th final sum of SSnoise value
%detailmap - a matrix of sum of SSnoise values at different core dimensions
%checked in the final analysis



    nlist=cat(2,[1:32]',[1:32]',[1:32]');%%%% incremental lists of reduced
    nlist2=cat(2,[1:32]',[1:32]',[3:3:96]');% rank to cycle throught
    nTv=zeros(32,2);
    zpoint=1;
    zblock=0;
    for n=1:32
        coredim=nlist(n,:);
        [denoised,nT]=tridenoise(kspace,nsmap,coredim); %performing the denoising
        nTv(n,1)=nT;%summed for all significant kspace points
        if (zblock==0)&&(nTv(n,1)<target_value)
            zpoint=coredim;
            zblock=1;
        end
        clear coredim nT 
    end
    topright=zpoint;
    clear zpoint;
    
    %do the denoising run again with the time dimension incremented in 3s.
    zpoint=3;
    zblock=0;
    for n=1:32
        coredim=nlist2(n,:);
        [denoised,nT]=tridenoise(kspace,nsmap,coredim);
        nTv(n,2)=nT;%summed for all significant kspace points
        clear nT 
        if (zblock==0)&&(nTv(n,2)<target_value)
            zpoint=coredim;
            zblock=1;
        end
        clear coredim
    end
    bottomlefts=zpoint;
    clear zpoint;    
    
   %This section of code checks to see if the best solution of both runs
   %has the same spatial dimension, if it does the smaller time dimension
   %is selected as the final reduced rank. If not, a detailed grid search
   %is performed on the range of possible spatial and time dimension ranks.
    
   if bottomlefts(1,1)<topright(1,1);
       AA=[bottomlefts(1,1):1:topright(1,1)];
       BB=AA;
       if bottomlefts(1,3)>topright(1,3)
           CC=[topright(1,3):bottomlefts(1,3)];
           detailmap=zeros(AA(1,end)-AA(1,1)+1,BB(1,end)-BB(1,1)+1,CC(1,end)-CC(1,1)+1);
           %construct the 3D search grid
           for a=AA
               for b=BB
                   for c=CC
                       coredim=[a,b,c];
                       [denoised,nT]=tridenoise(kspace,nsmap,coredim);           
                       detailmap(a+1-AA(1,1),b+1-BB(1,1),c+1-CC(1,1))=nT;
                   end
               end
           end
           %find the highest below target value
           cp=detailmap-target_value;
           Tz=sign(cp+abs(cp)); %values below target_value are now zeros in Tz
           dp=cp+target_value*Tz; %values above target_value are now all >Tz away
           ep=abs(dp); %matrix to find the minimum index
           dops=size(ep);
           U=zeros(3,dops(3));
           for c=1:dops(3);
               Sh=ep(:,:,c);
               [M I]=min(Sh,[],1);
               V(1,1:dops(2))=M;
               V(2,1:dops(2))=I; %this indicates the row where the minimum happened
               clear M I
               [M I]=min(V(1,1:dops(2))); %this finds the column of the minimal minimum row
               U(1,c)=M;
               U(2,c)=V(2,I);
               U(3,c)=I;
           end
           [M I]=min(U(1,1:dops(3)));
           coredim=[AA(U(2,I)),BB(U(3,I)),CC(I)];
           [denoised,nT]=tridenoise(kspace,nsmap,coredim);
       elseif bottomlefts(1,3)<=topright(1,3)
           coredim=bottomlefts;
           AA=coredim(1);BB=coredim(2);CC=coredim(3);
           [denoised,nT]=tridenoise(kspace,nsmap,coredim);
           detailmap=nT;
       end
   elseif bottomlefts(1,1)>=topright(1,1)
       coredim=topright;
       [denoised,nT]=tridenoise(kspace,nsmap,coredim);
       detailmap=nT;
   end






end

