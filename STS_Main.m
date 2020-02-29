% Main body of STS Algorithm
% Please cite this paper as:
% Lone, Mohd Rafi & HAKIM, Najeeb-ud-Din. (2018). FPGA implementation of a low-power and area-efficient state-table-based compression algorithm for DSLR cameras. TURKISH JOURNAL OF ELECTRICAL ENGINEERING & COMPUTER SCIENCES. 26. 2928-2943. 10.3906/elk-1804-208. 
clear all;clc;
%filename = 'lena512.bmp';
%Orig_I = imread(filename);
Orig_I = Read_Raw('lena512.raw',512,512);
rates = [1];results = [];
for i = 1:size(rates,2)
    [Encoder_time,Decoder_time,RIm] = STS_function(Orig_I, rates(i));
    % Parameter Evaluation
    Q = 255;[R,C] = size(Orig_I);
    MSE = sum(sum((RIm-double(Orig_I)).^2))/R / C;
    bt = 0;
    peak_SNR = psnr(double(Orig_I),RIm,255);
    PAE = max(max(abs(double(Orig_I)-RIm)));
    MAE = sum(sum(abs(double(Orig_I)-RIm)))/numel(Orig_I);
    results = [results;peak_SNR,Encoder_time,Decoder_time,PAE,MAE,MSE];
end
imshow(Orig_I), figure, imshow(RIm,[])
clearvars -except results