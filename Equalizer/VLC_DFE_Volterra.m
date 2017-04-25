%%
%Volterra NLMS Equalyzer using PAM symbols, different evaluation of SNR and of
%the nonlinearity



clear;
clc;
close all;

addpath(['..' filesep 'VLC_Simulator' filesep]);
addpath(['..' filesep 'VLC_Simulator' filesep 'LED Parameters']);

load whiteLED_334-15.mat;


%-------------------------Adaptive Filtering Parameters--------------------
numberOfBits = 2;
maxRuns = 15000;
maxIt = 1000;
gamma = 1e-12;
SNR = 30;
mu = 0.8;

volterraFFFlag = 1;
volterraFBFlag = 0;

feedforwardLength = 12;
feedbackLength = 12;

N = feedforwardLength;

adaptfiltFF = (feedforwardLength^2+feedforwardLength)/2 + feedforwardLength;
adaptfiltFB = (feedbackLength^2+feedbackLength)/2 + feedbackLength;

adaptfilt = adaptfiltFF + adaptfiltFB;

auxMatrix = triu(ones(feedforwardLength));
[l1FF,l2FF] = find(auxMatrix);

auxMatrix = triu(ones(feedbackLength));
[l1FB,l2FB] = find(auxMatrix);

delayinSamples = 14;


if ~volterraFFFlag
    adaptfiltFF = feedforwardLength;
end

if ~volterraFBFlag
    adaptfiltFB = feedbackLength;
end

    
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


wFinal = zeros(length(modulationIndexVector),adapFiltLength);
e3 = zeros(length(modulationIndexVector),maxRuns + delayinSamples + 1 + adapFiltLength);

for index = 1:length(modulationIndexVector)
    modulationIndex = modulationIndexVector(index);
    
    if modulationIndex > maxModulationIndex
        warning('Modulation Index may cause undesired nonlinear effects')
    end

    maxVoltage = VDC*(1+modulationIndex);
    deltaV = maxVoltage - VDC;


    w2 = zeros(adapFiltLength,maxRuns,maxIt);
    e2 = zeros(maxRuns + delayinSamples + 1 + adapFiltLength,maxIt);
    
    for j = 1:maxIt
        j


        input = randi([0,2^numberOfBits-1],maxRuns*2,1);
        pilot = real(pammod(input,2^numberOfBits,0,'gray'));

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
        powerNoise = (receivedCurrentSignalPower/db2pow(SNR));
        n = n.*sqrt(powerNoise/powerNoiseAux);

        receivedVoltageSignalAux = (receivedCurrentSignal + n);
        receivedVoltageSignalAux = receivedVoltageSignalAux - mean(receivedVoltageSignalAux);
        receivedVoltageSignal =  receivedVoltageSignalAux*sqrt(var(pilot)/var(receivedVoltageSignalAux));

        xAux = [zeros(N-1,1);receivedVoltageSignal];

        w = zeros(adapFiltLength,maxRuns);

        d = zeros(maxRuns + delayinSamples + 1,1);
        e = zeros(maxRuns + delayinSamples + 1,1);
        
        for k = delayinSamples + adapFiltLength + 1:maxRuns + delayinSamples + 1 + adapFiltLength

            x = xAux(k:-1:k-feedforwardLength+1);

            yHat = (pilot(-delayinSamples + k + 1 -1:-1:-delayinSamples + k + 1 - feedbackLength - 1 + 1));

            if volterraFFFlag

                aux = zeros((feedforwardLength^2+feedforwardLength)/2,1);

                for lIndex = 1:length(l1FF)
                    aux(lIndex,1) = x(l1FF(lIndex),1)*(x(l2FF(lIndex),1));
                end
                xConc = [x;aux];
            else
                xConc = x;
            end


            if volterraFBFlag
                aux = zeros((feedbackLength^2+feedbackLength)/2,1);
                for lIndex = 1:length(l1FB)
                    aux(lIndex,1) = yHat(l1FB(lIndex),1)*(yHat(l2FB(lIndex),k1));
                end

                yHatConc = [yHat;aux];
            else
                yHatConc = yHat;
            end

            if ~volterraFFFlag && ~volterraFBFlag 
                xConc = x;
                yHatConc = yHat;
            end

            z = [xConc;yHatConc];

            d(k) = (pilot(-delayinSamples + k + 1)); 

            e(k) = d(k) - w(:,k)'*z;

            w(:,k+1) = w(:,k) + mu*z*((z'*z+gamma*eye(1))\eye(1))*conj(e(k));

        end
        w2(:,:,j) = (w(:,1:maxRuns));
        e2(:,j) = abs(e).^2;
    end


    w3 = mean(w2,3);
    wFinal(index,:) = w3(:,end);

    e3(index,:) = mean(e2,2);

end


save(['.' filesep 'resultsMSE' filesep 'testDFEVolterraEq.mat'],'wFinal','e3');


rmpath(['..' filesep 'VLC_Simulator' filesep]);
rmpath(['..' filesep 'VLC_Simulator' filesep 'LED Parameters']);

