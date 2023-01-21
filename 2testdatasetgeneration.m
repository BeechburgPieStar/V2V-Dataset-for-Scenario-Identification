clc
close all
clear all
% Link parameters
mcs = 1;       % 
psduLen = 500; % PSDU length in bytes
channeltype=["Rural LOS","Urban approaching LOS","Urban NLOS","Highway LOS","Highway NLOS"];
NC=10000;
H = [];
for snr = 15:40
    name = ['V2V_H_SNR=',num2str(snr),'.mat'];
for j=1:5
    for i = 1:NC
    % Create a format configuration object for an 802.11p transmission
    cfgNHT = wlanNonHTConfig;
    cfgNHT.ChannelBandwidth = 'CBW10';
    cfgNHT.PSDULength = psduLen;
    cfgNHT.MCS = mcs;
    % Create and configure the channel
    fs = wlanSampleRate(cfgNHT); % Baseband sampling rate for 10 MHz
    
    chan = V2VChannel;
    chan.SampleRate = fs;
    chan.DelayProfile = channeltype(j);
    temp = v2vChanEstSimulator(cfgNHT, chan, snr);
    H(i+NC*(j-1),:,1)= real(temp);
    H(i+NC*(j-1),:,2)= imag(temp);
    end
end
   save(name,'H')
end