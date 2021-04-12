%This is a script to give an example of how to do the tri processing
%It refers to sections in the paper and supplementary methods of 
%"Tensor rank truncation for image enhancement of spectroscopic imaging
% data while conserving low amplitude resonances"
%Wright et al. 

%The first step is to load in an example kspace data set of 32x32
%dimensions with 128 fid points. The centre of kspace is at row 17 and
%column 17.
load('example_dataset.mat');
kspace=ex_kspace_1;

%The dataset is processed to determine the points in kspace that have
%significant signal in their fid
[nsmap,fsmap]=kpoints(kspace);

%The expected noise of these points needs to be calculated
%1) Count the number of points that have no significant signal 
fno=sum(sum(fsmap)); %number of excluded kspace points
%2)reduce the time dimension of the data to a residual sum of squares
ksrss=sum(kspace.*conj(kspace),3); %2D kspace of rss values
%3)Calculate the average rss value per kspace point by summing over all
%signal free kspace points and dividing by fno
fsa=sum(sum(ksrss.*fsmap))/fno;  %mean point summed rss
%4)Count the number of points that have significant signal
nno=sum(sum(nsmap)); %number of included kspace points
%5)Calculate ΣkSSnoise for use in [3] (the expected value of the summed rss
%noise in all kspace points with significant signal) by multiplying the mean-point
%summed rss by the number of kspace points with signal.
tno=fsa*nno;


%The significant positions in kspace are now used to determine at what
%dimension the dataset can be denoised without losing data. 
%Cycle through potential rank orders of tri and then pick the largest
%denoising that still satisfies the inequality in expression [3] from the
%paper. The inputs are the kspace data, the position of kspace points with
%significant signal and the expected sum of rss noise over these points
[denoised, core,snoisep,details] = denoiseHyp13C(kspace, nsmap, tno);

%Fourier transforming to spectral-spatial data
spectra=ftingtri(kspace,2);
dspectra=ftingtri(denoised,2);