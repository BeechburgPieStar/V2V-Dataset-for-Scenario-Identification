function chanEst = v2vChanEstSimulator(cfgNHT, chan, snr)
%V2VPERSimulator Simulates the non-HT transmit-receive link over a
%Vehicle-to-Vehicle (V2V) fading channel

%   Copyright 2018-2019 The MathWorks, Inc.

% Waveform generation parameters
idleTime = 0;
numPkts = 1;
winTransTime = 0; % No windowing

fs = wlanSampleRate(cfgNHT); % Baseband sampling rate

% Indices for accessing each field within the time-domain packet
ind = wlanFieldIndices(cfgNHT);

% Get the number of occupied subcarriers and FFT length
ofdmInfo = wlan.internal.wlanGetOFDMConfig(cfgNHT.ChannelBandwidth, 'Long', 'Legacy');
Nst = numel(ofdmInfo.DataIndices)+numel(ofdmInfo.PilotIndices); % Number of occupied subcarriers

% Create an instance of the AWGN channel per SNR point simulated
awgnChannel = comm.AWGNChannel;
awgnChannel.NoiseMethod = 'Signal to noise ratio (SNR)';
awgnChannel.SignalPower = 1; % Unit power
awgnChannel.SNR = snr-10*log10(ofdmInfo.FFTLength/Nst); % Account for energy in nulls

chDelay = 100; % Arbitrary delay to account for V2V channel scenarios

% Loop to simulate multiple packets
numPacketErrors = 1;
numPkt = 1;          % Index of packet transmitted
maxNumErrors = 10;
maxNumPackets = 100;

while numPacketErrors<=maxNumErrors && numPkt<=maxNumPackets
    % Generate a packet waveform
    inpPSDU = randi([0 1], cfgNHT.PSDULength*8, 1); % PSDULength in bytes
    
    tx = wlanWaveformGenerator(inpPSDU, cfgNHT, 'IdleTime', idleTime, ...
        'NumPackets', numPkts, 'WindowTransitionTime', winTransTime);
    
    % Add trailing zeros to allow for channel delay
    padTx = [tx; zeros(chDelay, 1)];
    
    % Pass through V2V channel model
    rx = chan(padTx);
    reset(chan); % Reset channel to create different realizations
    
    % Add noise
    rx = awgnChannel(rx);
    
    % Packet detect and determine coarse packet offset
    coarsePktOffset = wlanPacketDetect(rx, cfgNHT.ChannelBandwidth);
    if isempty(coarsePktOffset) % If empty no L-STF detected; packet error
        numPacketErrors = numPacketErrors+1;
        numPkt = numPkt+1;
        continue; % Go to next loop iteration
    end
    
    % Extract L-STF and perform coarse frequency offset correction
    lstf = rx(coarsePktOffset+(ind.LSTF(1):ind.LSTF(2)), :);
    coarseFreqOff = wlanCoarseCFOEstimate(lstf, cfgNHT.ChannelBandwidth);
    rx = helperFrequencyOffset(rx, fs, -coarseFreqOff);
    
    % Extract the non-HT fields and determine fine packet offset
    nonhtfields = rx(coarsePktOffset+(ind.LSTF(1):ind.LSIG(2)), :);
    finePktOffset = wlanSymbolTimingEstimate(nonhtfields, cfgNHT.ChannelBandwidth);
    
    % Determine final packet offset
    pktOffset = coarsePktOffset+finePktOffset;
    
    % If packet detected outside the range of expected delays from the
    % channel modeling; packet error
    if pktOffset>chDelay
        numPacketErrors = numPacketErrors+1;
        numPkt = numPkt+1;
        continue; % Go to next loop iteration
    end

    
    % Extract L-LTF and perform fine frequency offset correction
    lltf = rx(pktOffset+(ind.LLTF(1):ind.LLTF(2)), :);
    fineFreqOff = wlanFineCFOEstimate(lltf, cfgNHT.ChannelBandwidth);
    rx = helperFrequencyOffset(rx, fs, -fineFreqOff);

    % Extract L-LTF samples from the waveform, demodulate and perform
    % channel estimation
    lltf = rx(pktOffset+(ind.LLTF(1):ind.LLTF(2)), :);
    lltfDemod = wlanLLTFDemodulate(lltf, cfgNHT, 1);
    chanEst = wlanLLTFChannelEstimate(lltfDemod, cfgNHT);
    return
    
end