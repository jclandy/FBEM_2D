function [xk] = dft(xn,a,b)

ln=length(xn); %find the length of the sequence
xk=zeros(1,ln); %initialize an array of same size as that of input sequence

%code block to find the DFT of the sequence

i=sqrt(-1);
%-----------------------------------------------------------
for k=0:ln-1
    for n=0:ln-1
        xk(k+1)=xk(k+1)+(xn(n+1)*(1/(ln^((1-a)/2)))*exp(b*i*2*pi*k*n/ln));
    end
end
%------------------------------------------------------------

end

