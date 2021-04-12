function [denoised,rrss]=tridenoise(kspace,nsmap,coredims);
%Script to perform denoising by tucker tensor decocmposition
%calculates the sum of residual over defined kspace voxels

%functions called
%- ftingtri
%- funtctions of the tensor toolbox

%INPUTS
%kspace - a three dimensional matrix of dimensions 32x32x128 consisting of
%raw acquired MRSI data in kspace and time dimensions
%nsmap - a 32x32 matrix indicating the kspace points that have significant
%signal
%coredims - a vector of the rank order to reduce the data to in the denoising

%OUTPUTS
%denoised - the kspace data denoised
%rrss - the sum of the sum of rssFID over the nsmap points.

        sspace=ftingtri(kspace,2);              
        tensorr=tensor(real(sspace));           
        tensori=tensor(imag(sspace));           
        realT= tucker_als(tensorr,coredims ,'printitn',100);
        imagT= tucker_als(tensori,coredims,'printitn',100 );    
        rT=double(realT);
        iT=double(imagT);
        sc=complex(rT,iT);
        kc=ftingtri(sc,5);
        resmat=(kc-kspace);          %the residual
        rtwo=resmat.*conj(resmat);   %as an absolute square
        rflat=sum(rtwo,3);            %sum over the whole fid
        rrss=sum(sum(rflat.*nsmap));%summed for all significant kspace points
        denoised=kc;
end

