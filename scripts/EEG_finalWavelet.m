%NEURO MARKETING
clear all 
close all
clc
input_data = csvread('jai surya_EPOCX_217920_2024.04.10T13.31.18+05.30.md.csv',2,0); % This excludes the first two rows
size(input_data);
eegcols = 5:18; % EEG Columns.
eeg_data = input_data(:, eegcols);
% plot(eeg_data);
size(eeg_data);
N=4;
Fc1 = 0.16; % First Cutoff Frequency
Fc2 = 64;  % Second Cutoff Frequency
fs = 256;  % Sampling frequency (update as per your data)

% Design bandpass filter
h_bandpass = fdesign.bandpass('N,F3dB1,F3dB2', N, Fc1, Fc2, fs);
filter_bandpass = design(h_bandpass, 'butter');

% Apply bandpass filter
butter_filtered = filter(filter_bandpass, eeg_data);


notch = designfilt('bandstopiir','FilterOrder',2, ...
    'HalfPowerFrequency1',49,'HalfPowerFrequency2',51, ...
    'DesignMethod','butter','SampleRate',fs);    %notch filter to remove power line interference at 50 Hz
notchfiltered=filter(notch,butter_filtered);
eeg_filtered = notchfiltered;
size(eeg_filtered)
N=length(eeg_filtered);
[r,c]=size(eeg_filtered)
ST = input('enter starting sample of events\n');%entering start of any segment
ET = input('enter ending sample of events\n');%end sample no pf that particular segment
epoch = input('enter the epoch (sec)\n'); % enter epoch (window period) in sec 
TotalTime=(ET-ST);% total time of corresponding segment ex 5:90000
j=0;
Tflag=round((TotalTime))% roundinf off of total time ex 90
 for flag=1:epoch:(Tflag)% finding how many no of segments to run 1:10:90 ( 9 times loop to run) This loop iterates over each segment defined by the epoch length until the total time duration is reached. It increments j for each iteration, which represents the number of segment
     j=j+1
    for i= 1:c % picking column of filter signal from y9 and ( ex 1 to 10) This loop iterates over each column of the filtered EEG data. It extracts a column (z1) and transposes it to make it a row vector.
       z1=eeg_filtered(:,i);%  
z1=z1';% first row transpose  

size(z1);
x22=z1(ST:(ST+(epoch*fs)));% This line picks a segment of data from z1 starting from the ST sample and ending at ST + (epoch * fs). This segment represents the epoch window for analysis.
  
N=length(x22);
 

waveletFunction = 'db8'; %This line performs a wavelet decomposition of the segment x22 using the Daubechies wavelet ('db8') up to level 5. It returns the wavelet coefficients (C) and the bookkeeping matrix (L).

% The following lines calculate various features from the wavelet coefficients:
              
[C,L] = wavedec(x22,5,waveletFunction);
       
                
cD1 = detcoef(C,L,1);
cD2 = detcoef(C,L,2);
cD3 = detcoef(C,L,3);
cD4 = detcoef(C,L,4);
cD5 = detcoef(C,L,5);
               
cA8 = appcoef(C,L,waveletFunction,5); 
D1 = abs(wrcoef('d',C,L,waveletFunction,1));
D2 = abs(wrcoef('d',C,L,waveletFunction,2));
D3 = abs(wrcoef('d',C,L,waveletFunction,3));
D4 = abs(wrcoef('d',C,L,waveletFunction,4));
D5 = abs(wrcoef('d',C,L,waveletFunction,5)); 

A1 = wrcoef('a',C,L,waveletFunction,1);
A2 = wrcoef('a',C,L,waveletFunction,2);
A3 = wrcoef('a',C,L,waveletFunction,3);
A4 = wrcoef('a',C,L,waveletFunction,4);
A5 = abs(wrcoef('a',C,L,waveletFunction,5)); 

abpDelta=mean(A5)  %These lines calculate the mean values of the wavelet coefficients at different frequency bands, which are typically associated with different EEG rhythms.
abpTheta=mean(D5)
abpAlpha=mean(D4)
abpBeta=mean(D3)
abpGamma=mean(D2)

TBR = abpTheta ./ abpBeta;
HR=abpTheta ./abpAlpha; 
AI=abpBeta ./abpAlpha;
NA=abpBeta ./abpTheta;
Sync=abpDelta ./abpTheta;                     % Synchronisation
BP=abpAlpha ./abpDelta;                       %Brain Perfusion
CNSA=abpTheta ./abpBeta;                      % CNS Arousal
Dsync=abpAlpha ./abpBeta;                     %Desynchronisation
CPI=abpBeta ./(abpAlpha+abpTheta);            %Cognitive performance index
ELI=(abpDelta+abpTheta) ./abpAlpha;           %executive load index
PEI=abpAlpha ./abpTheta;                      %performance enhancement index
LF=(abpTheta+abpDelta) ./(abpAlpha+abpBeta);  %LF to HF ratio
LTI=abpTheta ./abpAlpha;                      %task load index
VI=(abpTheta+abpAlpha) ./abpBeta;             %Vigilence index
AGR=abpAlpha ./abpGamma
 

abpDelta1(:,i)=abpDelta;
abpTheta1(:,i)=abpTheta;
abpAlpha1(:,i)=abpAlpha;
abpBeta1(:,i)=abpBeta;
abpGamma1(:,i)=abpGamma;
HR1(:,i)=HR;                       % Heart rate index
AI1(:,i)=AI;                        %Arousal index
NA1(:,i)=NA;                        % Neuronal Activity
Sync1(:,i)=Sync;                     % Synchronisation
BP1(:,i)=BP;                      %Brain Perfusion
CNSA1(:,i)=CNSA;                      % CNS Arousal
Dsync1(:,i)=Dsync;                    %Desynchronisation
CPI1(:,i)=CPI;                     %Cognitive performance index
ELI1(:,i)=ELI;                      %executive load index
PEI1(:,i)=PEI;                     %performance enhancement index
LF1(:,i)=LF;                      %LF to HF ratio
LTI1(:,i)=LTI;                     %task load index
VI1(:,i)=VI;                        %Vigilence index
AGR1(:,i)=AGR;                 
 end
    
ST=ST+(epoch*fs);
abpDelta2(:,j)=abpDelta1';
abpTheta2(:,j)=abpTheta1';
abpAlpha2(:,j)=abpAlpha1';
abpBeta2(:,j)=abpBeta1';
abpGamma2(:,j)=abpGamma1';
HR2(:,j)=HR1';                         % Heart rate index
AI2(:,j)=AI1';                         %Arousal index
NA2(:,j)=NA1';                         % Neuronal Activity
Sync2(:,j)=Sync1';                     % Synchronisation
BP2(:,j)=BP1';                         %Brain Perfusion
CNSA2(:,j)=CNSA1';                     % CNS Arousal
Dsync2(:,j)=Dsync1';                   %Desynchronisation
CPI2(:,j)=CPI1';                       %Cognitive performance index
ELI2(:,j)=ELI1';                       %executive load index
PEI2(:,j)=PEI1';                       %performance enhancement index
LF2(:,j)=LF1';                         %LF to HF ratio
LTI2(:,j)=LTI1';                       %task load index
VI2(:,j)=VI1';                       %Vigilence index
AGR2(:,j)=AGR1'; 
      end

result=[abpDelta2' abpTheta2' abpAlpha2' abpBeta2' abpGamma2' HR2' AI2' NA2' Sync2' BP2' CNSA2' Dsync2' CPI2' ELI2' PEI2' LF2' LTI2' VI2' AGR2'];




