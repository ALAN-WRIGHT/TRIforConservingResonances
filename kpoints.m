function [ nsmap,fsmap ] = kpoints( kspace )
%Function for defining the kspace points that contain signal and those that
%just contain noise for Hyperpolarised [1-13C]pytuvate datasets as
%described in Supplementary methods 2.
%
%INPUTS
%kspace a three dimensional matrix containing kspace and time MRSI data

%OUTPUTS
%nsmap: 2D map (1s in 0s) of the kspace points that contain significant
%signal
%fsmap: 2D map (1s in 0s) of the kspace points that contain no significant
%signal

    %setup%%%%%%%%%%%%%%%%
    dims=size(kspace);%%%%
    np=dims(3);%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%
    
    %simulate fid noise and calculate the sum of its autocorrelation function%%%%
    sumsim=zeros(10000,1);
    for a=1:10000
        A=randn([1 np]);
        simline=xcorr(A,'coeff');
        sumsim(a,1)=sum(abs(simline));
        clear simline A
    end
    acav=mean(sumsim(:,1));
    acs=std(sumsim(:,1));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %define signal and noise fids%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ns=0
    nf=0;
    fsmap=zeros(dims(1),dims(2));
    nsmap=zeros(dims(1),dims(2));
    for a=1:dims(1)
        for b=1:dims(2)
            fidline(1,1:np)=kspace(a,b,1:np);
            fidacr=sum(abs(xcorr(real(fidline),'coeff')));
            fidaci=sum(abs(xcorr(imag(fidline),'coeff')));
            if ((fidacr<=acav)&&(fidaci<(acav+2*acs)))||((fidaci<=acav)&&(fidacr<(acav+2*acs)))%condition for defining kspace points as noise
                nf=nf+1;
                fsmap(a,b)=1;
            end
            
            if ((fidacr>acav)&&(fidaci>(acav+2*acs)))||((fidaci>acav)&&(fidacr>(acav+2*acs)))%condition for defining kspace points as signal
                ns=ns+1;
                nsmap(a,b)=1;
            end
            
            
            
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end

