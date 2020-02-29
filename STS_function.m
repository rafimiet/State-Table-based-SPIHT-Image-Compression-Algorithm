function [Encoder_time,Decoder_time,RIm] = STS_function(Imge, rate, fg)
% Acquisition
    Orig_I = double(Imge);
    [R,C] = size(Orig_I);
    if R > 512
        L = 8;
    else
        if R > 256
            L = 7;
        else
            if R > 128
                L = 6;
            else
                L = 5;
            end
        end
    end
% DWT Transform
    type = 'bior4.4';
    [Lo_D,Hi_D,Lo_R,Hi_R] = wfilters(type);
    [S, S_m] = func_DWT(Orig_I, L, Lo_D, Hi_D);
    M = max(max(S));
    n = floor(log(M)/log(2));
    Tint = 2^n; % THRESHOLD
% STS Encoder
tic
EIm = STS_Encoder(S,rate,L);
Encoder_time = toc; 
len = 0;
for i = 1:size(EIm,2)
    len = len + length(EIm{i});
end
actual_bpp = len/numel(Imge);
CR = 8/actual_bpp;
% STS Decoder
tic 
[DIm] = STS_Decoder(EIm,size(S,1),size(S,2),L,Tint);
if DIm<0, clear, return; end;
Decoder_time = toc; 

% Step 6: Reverse wavelet transform; 

RIm = func_InvDWT(DIm, S_m, Lo_R, Hi_R, L);
return

