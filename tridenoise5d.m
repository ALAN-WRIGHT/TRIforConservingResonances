function [denoised,rrss]=tridenoise5d(kspace,nsmap,coredims,npo);
%Script to perform denoising by tucker tensor decocmposition
%calculates the sum of residual over defined kspace voxels
%this is the 5 dimensional version where the dimensions are, in order,
%1) 2) 3) spatial/ inverse spatial
%4) time(fid)/ frequency(spectrum)
%5) dynamic time series - this dimension is not Fourier transformed
%np - this specifies the number of points to evaluate (from the fid start)
%for the purposes of calculating the rrss number

%functions called
%- ftingtri

%INPUTS
%kspace - a three dimensional matrix of dimensions 32x32x128 consisting of
%raw acquired MRSI data in kspace and time dimensions
%nsmap - a 32x32 matrix indicating the kspace points that have significant
%signal
%coredims - a vector of the rank order to reduce the data to in the denoising

%OUTPUTS
%denoised - the kspace data denoised
%rrss - the sum of the sum of rssFID over the nsmap points.

        dims=size(kspace);
        nr=dims(5);
if npo<=dims(4);
        sspace=zeros(dims);
        for n=1:nr
            sspace(:,:,:,:,n)=ftingtri(kspace(:,:,:,:,n),2); 
        end
        tensorr=tensor(real(sspace));           
        tensori=tensor(imag(sspace));           
        realT= tucker_als(tensorr,coredims ,'printitn',100);
        imagT= tucker_als(tensori,coredims,'printitn',100 );    
        rT=double(realT);
        iT=double(imagT);
        sc=complex(rT,iT);
        kc=zeros(dims);
        for n=1:nr
            kc(:,:,:,:,n)=ftingtri(sc(:,:,:,:,n),5); 
        end        
        resmat=(kc-kspace);          %the residual
        resmat2=resmat(:,:,:,1:npo,:);
        rtwo=resmat2.*conj(resmat2);   %as an absolute square
        rflat=sum(rtwo,4);            %sum over the whole fid
        rflat2=permute(rflat,[1 2 3 5 4]); %rearange so dimensions the same as nsmap
        rrss=sum(sum(sum(sum(rflat2.*nsmap))));%summed for all significant kspace points
        denoised=kc;
else
    error='dimensions inconsistent with fid test length'
end

