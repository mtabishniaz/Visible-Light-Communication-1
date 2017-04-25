clear;
clc;
close all;




addpath(['.' filesep 'resultsBER']);



load resultsBERLinEq.mat;
figure
for i = 1:size(ber,2)
    semilogy(SNR,ber(:,i))
    hold on
end
xlabel('SNR [dB]','interpreter','latex');
ylabel('BER','interpreter','latex');

H = legend('MI = 0.05','MI = 0.075','MI = 0.1');

set(H,'interpreter','latex','location','SouthWest')
xlim([min(SNR) max(SNR)]);
ylim([1e-8 1]);


load resultsBERVolterraEq.mat;
figure
for i = 1:size(ber,2)
    semilogy(SNR,ber(:,i))
    hold on
end
xlabel('SNR [dB]','interpreter','latex');
ylabel('BER','interpreter','latex');

H = legend('MI = 0.05','MI = 0.075','MI = 0.1');

set(H,'interpreter','latex','location','SouthWest')
xlim([min(SNR) max(SNR)]);
ylim([1e-8 1]);


load resultsBERDFE.mat;
figure
for i = 1:size(ber,2)
    semilogy(SNR,ber(:,i))
    hold on
end
xlabel('SNR [dB]','interpreter','latex');
ylabel('BER','interpreter','latex');

H = legend('MI = 0.05','MI = 0.075','MI = 0.1');

set(H,'interpreter','latex','location','SouthWest')
xlim([min(SNR) max(SNR)]);
ylim([1e-8 1]);


load resultsBERVolterraDFE.mat;
figure
for i = 1:size(ber,2)
    semilogy(SNR,ber(:,i))
    hold on
end
xlabel('SNR [dB]','interpreter','latex');
ylabel('BER','interpreter','latex');

H = legend('MI = 0.05','MI = 0.075','MI = 0.1');

set(H,'interpreter','latex','location','SouthWest')
xlim([min(SNR) max(SNR)]);
ylim([1e-8 1]);



