function [ matrixout ] = ftingtri( matrixin, argument )
%Function for doing a fourier transform of multidimensional data
%Works for 2 spatial and one spectral dimension or 3 spatial and one
%spectral dimension data (2D or 3D MRSImages)
dims=size(matrixin);
m1=zeros(dims);
m2=m1;
m3=m1;
if length(dims)==3
    %2D MRSImage
   rr=dims(1);
   cc=dims(2);
   np=dims(3);
   if argument==1
       %Fourier transform just the inverse spatial dimensions to spatial
       %dimensions leave (time) dimension alone
       
       %Fourier transform along columns
        for r=1:rr
                for n=1:np
                    m1(r,1:cc,n)=fftshift(fft(ifftshift(matrixin(r,1:cc,n))));
                end
        end

        %Fourier transform along rows
        for c=1:cc
                for n=1:np
                    m2(1:rr,c,n)=fftshift(fft(ifftshift(m1(1:rr,c,n))));
                end
        end
       
        matrixout=m2;
       
   elseif argument==2
       %Fourier transform full kspace to full spectral spatial
       
       %Fourier transform along columns
        for r=1:rr
                for n=1:np
                    m1(r,1:cc,n)=fftshift(fft(ifftshift(matrixin(r,1:cc,n))));
                end
        end

        %Fourier transform along rows
        for c=1:cc
                for n=1:np
                    m2(1:rr,c,n)=fftshift(fft(ifftshift(m1(1:rr,c,n))));
                end
        end
        
        %Fourier transform along time
        for r=1:rr
                for c=1:cc
                    m3(r,c,1:np)=fftshift(fft(m2(r,c,1:np)));
                end
        end        
        
       
        matrixout=m3;
       
   elseif argument==3
        %just Fourier transform the time dimension to spectral
        for r=1:rr
                for c=1:cc
                    m1(r,c,1:np)=fftshift(fft(matrixin(r,c,1:np)));
                end
        end
        matrixout=m1;
   elseif argument==4
        %Fourier transform just the inverse spatial dimensions to k-spatial
       %dimensions leave (time) dimension alone
       
       %Fourier transform along columns
        for r=1:rr
                for n=1:np
                    m1(r,1:cc,n)=fftshift(ifft(ifftshift(matrixin(r,1:cc,n))));
                end
        end

        %Fourier transform along rows
        for c=1:cc
                for n=1:np
                    m2(1:rr,c,n)=fftshift(ifft(ifftshift(m1(1:rr,c,n))));
                end
        end
       
        matrixout=m2;
       
       
   elseif argument==5          
   %Fourier transform along columns
        for r=1:rr
                for n=1:np
                    m1(r,1:cc,n)=fftshift(ifft(ifftshift(matrixin(r,1:cc,n))));
                end
        end

        %Fourier transform along rows
        for c=1:cc
                for n=1:np
                    m2(1:rr,c,n)=fftshift(ifft(ifftshift(m1(1:rr,c,n))));
                end
        end
        
        %Fourier transform along time
        for r=1:rr
                for c=1:cc
                    m3(r,c,1:np)=ifft(ifftshift(m2(r,c,1:np)));
                end
        end  
        
        matrixout=m3;
        
   elseif argument==6     
        %Fourier transform along time
        for r=1:rr
                for c=1:cc
                    m1(r,c,1:np)=ifft(ifftshift(matrixin(r,c,1:np)));
                end
        end  
        
        matrixout=m1;        
   end
        
elseif length(dims)==4
    %3D MRSImage
   m4=m1;
   rr=dims(1);
   cc=dims(2);
   ss=dims(3);
   np=dims(4);
   if argument==1
       %Fourier transform just the inverse spatial dimensions to spatial
       %dimensions leave (time) dimension alone
       
       %Fourier transform along slices
        for r=1:rr
            for c=1:cc
                for n=1:np
                    m1(r,c,1:ss,n)=fftshift(fft(ifftshift(matrixin(r,c,1:ss,n))));
                end
            end
        end
       
       %Fourier transform along columns
        for r=1:rr
            for s=1:ss
                for n=1:np
                    m2(r,1:cc,s,n)=fftshift(fft(ifftshift(m1(r,1:cc,s,n))));
                end
            end
        end

        %Fourier transform along rows
        for c=1:cc
            for s=1:ss
                for n=1:np
                    m3(1:rr,c,s,n)=fftshift(fft(ifftshift(m2(1:rr,c,s,n))));
                end
            end
        end
       
        matrixout=m3;
       
   elseif argument==2
       %Fourier transform full kspace to full spectral spatial
       
       %Fourier transform just the inverse spatial dimensions to spatial
       %dimensions leave (time) dimension alone
       
       %Fourier transform along slices
        for r=1:rr
            for c=1:cc
                for n=1:np
                    m1(r,c,1:ss,n)=fftshift(fft(ifftshift(matrixin(r,c,1:ss,n))));
                end
            end
        end
       
       %Fourier transform along columns
        for r=1:rr
            for s=1:ss
                for n=1:np
                    m2(r,1:cc,s,n)=fftshift(fft(ifftshift(m1(r,1:cc,s,n))));
                end
            end
        end

        %Fourier transform along rows
        for c=1:cc
            for s=1:ss
                for n=1:np
                    m3(1:rr,c,s,n)=fftshift(fft(ifftshift(m2(1:rr,c,s,n))));
                end
            end
        end
        
        %Fourier transform along time
        for r=1:rr
            for c=1:cc
                for s=1:ss
                    m4(r,c,s,1:np)=fftshift(fft(m3(r,c,s,1:np)));
                end
            end
        end        
        
       
        matrixout=m4;
       
   elseif argument==3
        %Fourier transform along time
        for r=1:rr
            for c=1:cc
                for s=1:ss
                    m1(r,c,s,1:np)=fftshift(fft(matrixin(r,c,s,1:np)));
                end
            end
        end        
        
       
    
        matrixout=m1;
   elseif argument==4
        %Fourier transform just the inverse spatial dimensions to spatial
       %dimensions leave (time) dimension alone
       
        %Fourier transform along slices
        for r=1:rr
            for c=1:cc
                for n=1:np
                    m1(r,c,1:ss,n)=fftshift(ifft(ifftshift(matrixin(r,c,1:ss,n))));
                end
            end
        end
       
       %Fourier transform along columns
        for r=1:rr
            for s=1:ss
                for n=1:np
                    m2(r,1:cc,s,n)=fftshift(ifft(ifftshift(m1(r,1:cc,s,n))));
                end
            end
        end

        %Fourier transform along rows
        for c=1:cc
            for s=1:ss
                for n=1:np
                    m3(1:rr,c,s,n)=fftshift(ifft(ifftshift(m2(1:rr,c,s,n))));
                end
            end
        end
       
        matrixout=m3;
       
       
   elseif argument==5          

        %Fourier transform along slices
        for r=1:rr
            for c=1:cc
                for n=1:np
                    m1(r,c,1:ss,n)=fftshift(ifft(ifftshift(matrixin(r,c,1:ss,n))));
                end
            end
        end
       
       %Fourier transform along columns
        for r=1:rr
            for s=1:ss
                for n=1:np
                    m2(r,1:cc,s,n)=fftshift(ifft(ifftshift(m1(r,1:cc,s,n))));
                end
            end
        end

        %Fourier transform along rows
        for c=1:cc
            for s=1:ss
                for n=1:np
                    m3(1:rr,c,s,n)=fftshift(ifft(ifftshift(m2(1:rr,c,s,n))));
                end
            end
        end
        
    %Fourier transform along time
        for r=1:rr
            for c=1:cc
                for s=1:ss
                    m4(r,c,s,1:np)=ifft(ifftshift(m3(r,c,s,1:np)));
                end
            end
        end   
        
        matrixout=m4;
        
   elseif argument==6     
    %Fourier transform along time
        for r=1:rr
            for c=1:cc
                for s=1:ss
                    m1(r,c,s,1:np)=ifft(ifftshift(matrixin(r,c,s,1:np)));
                end
            end
        end   
        
        matrixout=m1;        
   end
 
    
    
    
    
end



end

