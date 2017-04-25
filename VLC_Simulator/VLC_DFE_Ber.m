%%
%Volterra NLMS Equalyzer using PAM symbols, different evaluation of SNR and of
%the nonlinearity



clear;
clc;
close all;

addpath(['..' filesep 'Equalizer' filesep 'resultsMSE']);
addpath(['.' filesep 'LED Parameters']);

load whiteLED_334-15.mat;

load testDFEEq.mat;

%-------------------------Adaptive Filtering Parameters--------------------

M = 4;
numberOfSymbols = 40000;
numberOfBits = log2(M);

blockLength = numberOfSymbols*numberOfBits;
monteCarloLoops = 1000;

feedforwardLength = 12;
feedbackLength = 12;

N = feedforwardLength;

adaptfiltFF = feedforwardLength;
adaptfiltFB = feedbackLength;

adapFiltLength = adaptfiltFF + adaptfiltFB;

%-------------------------Adaptive Filtering Parameters--------------------



%-------------------------LED Parameters-----------------------------------

Poptical = @(ledLuminousEfficacy,electricalPower,k) (ledLuminousEfficacy.*electricalPower)./((1 + (ledLuminousEfficacy.*electricalPower./(maxLuminousIntensityLED/1000)).^(2*k)).^(1/(2*k)));

%-------------------------LED Parameters-----------------------------------

%-------------------------Photodiode Parameters----------------------------

A = 1e-4; %photodiode area (cm)
d = 10e-2; %distance between LED and photodiode (cm)
R = 0.5;
FOV = deg2rad(25);

%-------------------------Photodiode Parameters----------------------------


%-------------------------Pre Amplifier Parameters-------------------------

%-------------------------Pre Amplifier Parameters-------------------------


%-------------------------Transmission Parameters--------------------------

kNonLinearity = 2;

fs = 2e6;
LEDfreqRespPoints = 1000;

theta = 0;
phi = 0;

H_0 = A/d^2 * (n+1)/(2*pi) * cos(phi)^n * cos(theta) * rectangularPulse(-1,1,theta/FOV);

VDC = 3.25; 
maxAbsoluteValueModulation = 3;

maxModulationIndex = (maxLEDVoltage - VDC)/VDC;

modulationIndexVector = [0.05 0.075 0.1];

%-------------------------Transmission Parameters--------------------------

SNR = 0:5:30;

berAux = zeros(monteCarloLoops,1);
ber = zeros(length(SNR),length(modulationIndexVector));


for index = 1:length(modulationIndexVector)
    
    modulationIndex = modulationIndexVector(index);
    
    if modulationIndex > maxModulationIndex
        warning('Modulation Index may cause undesired nonlinear effects')
    end

    maxVoltage = VDC*(1+modulationIndex);
    deltaV = maxVoltage - VDC;

    equalyzerFilter = wFinal(index,:).';

    for SNRIndex = 1:length(SNR)


        for j = 1:monteCarloLoops
                j
            equalyzedSignal = zeros(numberOfSymbols,1);

            binaryInputData = randi([0,1],blockLength+100,1);
            binaryInputData = reshape(binaryInputData,[],numberOfBits);
            deciInputData = bi2de(binaryInputData);    
            pilot = real(pammod(deciInputData,2^numberOfBits,0,'gray'));

            convLength = length(pilot) + LEDfreqRespPoints -1;
            NFFT = 2^nextpow2(convLength);

            pilotFreq = fft(pilot,NFFT);

            f = fs/2*linspace(0,1,NFFT/2 + 1)*2*pi;

            w = [-fliplr(f(2:end-1)) f];

            LEDResp = freqRespLED(w);

            filteredVinAux = real(ifft(pilotFreq.*fftshift(LEDResp))); 

            filteredVin = filteredVinAux(1:length(pilot));

            VoltageConstant = modulationIndex*maxVoltage/((1+modulationIndex)*max(filteredVin));

            filteredVin = filteredVin*VoltageConstant + VDC;

            iLEDOutput = I_V_Fun(filteredVin,VT,nLED,ISat);

            eletricalPowerOutput = filteredVin.*iLEDOutput;

            opticalPowerOutput = Poptical(ledLuminousEfficacy,eletricalPowerOutput,kNonLinearity);

            opticalPowerOutputConvolved = opticalPowerOutput*H_0;

            n = randn(length(opticalPowerOutputConvolved),1); %noise signal

            receivedCurrentSignal = opticalPowerOutputConvolved*R*A;
            receivedCurrentSignalAC = receivedCurrentSignal - mean(receivedCurrentSignal);
            receivedCurrentSignalPower = receivedCurrentSignalAC'*receivedCurrentSignalAC/length(receivedCurrentSignal);

            powerNoiseAux = n'*n/(length(n));
            powerNoise = (receivedCurrentSignalPower/db2pow(SNR(SNRIndex)));
            n = n.*sqrt(powerNoise/powerNoiseAux);

            receivedVoltageSignalAux = (receivedCurrentSignal + n);
            receivedVoltageSignalAux = receivedVoltageSignalAux - mean(receivedVoltageSignalAux);
            receivedVoltageSignal =  receivedVoltageSignalAux*sqrt(var(pilot)/var(receivedVoltageSignalAux));
            
            xAux = [zeros(N-1,1);receivedVoltageSignal];
            
            outputFF = zeros(length(pilot),1);
            outputFB = zeros(length(pilot),1);
            
            for k = feedforwardLength:length(pilot)

                xFlip = xAux(k:-1:k-feedforwardLength+1);

                outputFF(k) =  (equalyzerFilter(1:adaptfiltFF))'*xFlip;
                equalyzedSignal(k,1) = pamHardThreshold(outputFF(k) + outputFB(k-1));
                inputFB = equalyzedSignal(k:-1:k-feedbackLength + 1); 

                outputFB(k) = (equalyzerFilter(adaptfiltFF+1:end))'*inputFB;

            end
            equalyzedSignal2 = (equalyzedSignal - VDC)/VoltageConstant;

           [corr,lags] = xcorr(equalyzedSignal,xAux(N:end));
           [~,idx] = max(abs(corr));
           delay = abs(lags(idx));

           decDemodSignal = pamdemod(equalyzedSignal,2^numberOfBits,0,'gray');

           binaryOutputData = de2bi(decDemodSignal,numberOfBits);

           berAux(j) = sum(sum(abs(binaryOutputData(delay+1:(blockLength/2) + delay,:) - binaryInputData(1:(blockLength/2),:))))./blockLength;
               
        end

        ber(SNRIndex,index) = mean(berAux);   

    end

end


save(['.' filesep 'resultsBER' filesep 'resultsBERDFE.mat'],'ber','SNR');

rmpath(['..' filesep 'Equalizer' filesep 'resultsMSE']);
rmpath(['.' filesep 'LED Parameters']);



