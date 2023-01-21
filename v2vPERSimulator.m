function packetErrorRate = v2vPERSimulator(cfgNHT, chan, snr, ...
    maxNumErrors, maxNumPackets, enableChanTracking)
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
numPacketErrors = 0;
numPkt = 0; % Index of packet transmitted
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
    
    % Get estimate of the noise power from L-LTF
    nVar = helperNoiseEstimate(lltfDemod);
    % Extract non-HT Data samples from the waveform and recover the PSDU
    nhtdata = rx(pktOffset+(ind.NonHTData(1):ind.NonHTData(2)), :);
    if enableChanTracking
        rxPSDU = V2VDataRecover(nhtdata, chanEst, nVar, cfgNHT);
    else
        rxPSDU = wlanNonHTDataRecover(nhtdata, chanEst, nVar, cfgNHT);
    end
    % Determine if any bits are in error, i.e. a packet error
    packetError = any(biterr(inpPSDU, rxPSDU));
    numPacketErrors = numPacketErrors+packetError;
    numPkt = numPkt+1;
end

% Calculate packet error rate (PER) at SNR point
if enableChanTracking
    trackingStr = 'with channel tracking';
else
    trackingStr = 'without channel tracking';
end
packetErrorRate = numPacketErrors/numPkt;
disp(['SNR ' num2str(snr) ' dB ' trackingStr ' completed after ' ...
    num2str(numPkt) ' packets, PER: ' num2str(packetErrorRate)]);

end


function [bits, eqDataSym] = V2VDataRecover( ...
    rxNonHTData, chanEst, noiseVarEst, cfgNonHT)
% V2VDataRecover Recover information bits from non-HT Data field

s = validateConfig(cfgNonHT);
numOFDMSym = s.NumDataSymbols;

% Data recovery configuration Parameters
symOffset = 0.75;
eqMethod  = 'MMSE';

mcsTable = wlan.internal.getRateTable(cfgNonHT);

% Get OFDM configuration
cfgOFDM = wlan.internal.wlanGetOFDMConfig(cfgNonHT.ChannelBandwidth, 'Long', 'Legacy');

minInputLen = numOFDMSym*(cfgOFDM.FFTLength+cfgOFDM.CyclicPrefixLength);

% Processing OFDM demodulation
[ofdmDemodData, ofdmDemodPilots] = wlan.internal.wlanOFDMDemodulate(rxNonHTData(1:minInputLen, :), cfgOFDM, symOffset);

% CPE correction and channel tracking
[hEst, ofdmDemodData] = V2VChanEst(ofdmDemodData, ofdmDemodPilots, chanEst, cfgNonHT, noiseVarEst);

% Equalization
csiData = zeros(size(ofdmDemodData));
eqDataSym = zeros(size(ofdmDemodData));

for cntr = 1:numOFDMSym
    [eqDataSym(:, cntr), csiData(:, cntr)] = wlan.internal.wlanEqualize(ofdmDemodData(:, cntr), hEst(:, cntr), eqMethod, noiseVarEst);
end
csiData = reshape(csiData, 1, [], numOFDMSym);

% Constellation demapping
qamDemodOut = wlanConstellationDemap(eqDataSym, noiseVarEst, mcsTable.NBPSCS);

% Apply bit-wise CSI and concatenate OFDM symbols in the first dimension
qamDemodOut = bsxfun(@times, ...
    reshape(qamDemodOut, mcsTable.NBPSCS, [], numOFDMSym), csiData); % [Nbpscs Nsd Nsym]
qamDemodOut = reshape(qamDemodOut, [], 1);

% Deinterleave
deintlvrOut = wlanBCCDeinterleave(qamDemodOut, 'Non-HT', mcsTable.NCBPS);

% Channel decoding
decBits = wlanBCCDecode(deintlvrOut, mcsTable.Rate);

% Derive initial state of the scrambler
scramInit = wlan.internal.scramblerInitialState(decBits(1:7));

% Remove pad and tail bits, and descramble
if all(scramInit==0)
    % Scrambler initialization invalid (0), therefore do not descramble
    descramDataOut = decBits(1:(16+8*cfgNonHT.PSDULength));
else
    descramDataOut = wlanScramble(decBits(1:(16+8*cfgNonHT.PSDULength)), scramInit);
end

% Remove the 16 service bits
bits = descramDataOut(17:end);

end

function [hEst, ofdmDemodDataChanEst] = V2VChanEst(ofdmDemodData, ofdmDemodPilots, chanEst, cfgNonHT, noiseVarEst)
%V2VChanEst Channel estimation at Data Subcarriers of Non-HT Data field as
%presented in J. A. Fernandez et al [1].
%  [1] J. A. Fernandez, K. Borries, L. Cheng, B. V. K.
%  Vijaya Kumar, D. D. Stancil and F. Bai, "Performance of the 802.11p
%  Physical Layer in Vehicle-to-Vehicle Environments," in IEEE Transactions
%  on Vehicular Technology, vol. 61, no. 1, pp. 3-14, Jan. 2012.

numOFDMSym = size(ofdmDemodData, 2);
mcsTable = wlan.internal.getRateTable(cfgNonHT);

% Get OFDM configuration
[ofdmInfo, dataInd, pilotInd] = wlan.internal.wlanGetOFDMConfig(cfgNonHT.ChannelBandwidth, 'Long', 'Legacy');

% Equalization method
eqMethod = 'MMSE';

% Extract data and pilot subcarriers from LTF channel estimate
chanEstData = chanEst(dataInd, :, :);
chanEstPilotsLTF = chanEst(pilotInd, :, :);

%% Common phase error
chanEstPilots = zeros(size(ofdmDemodPilots));
chanEstPilots(:, 1) = chanEstPilotsLTF;

ofdmDemodDataChanEst = ofdmDemodData;
% Pilot phase tracking, get reference pilots, from IEEE Std 802.11-2012, Eqn
% 18-22
z = 1; % Offset by 1 to account for L-SIG pilot symbol
refPilots = wlan.internal.nonHTPilots(numOFDMSym, z);
cpe = zeros(1, numOFDMSym);
for ctr = 1:numOFDMSym
    % Estimate CPE and phase correct symbols
    cpe(ctr) = wlan.internal.commonPhaseErrorEstimate(ofdmDemodPilots(:, ctr), chanEstPilots(:, ctr), refPilots(:, ctr));
    chanEstPilots(:, ctr+1) = ofdmDemodPilots(:, ctr)./refPilots(:, ctr);
    ofdmDemodDataChanEst(:, ctr) = wlan.internal.commonPhaseErrorCorrect(ofdmDemodData(:, ctr), cpe(ctr));
end

% Spectral Temporal Averaging (STA) channel Estimation
hEst = zeros(size(ofdmDemodDataChanEst));
chanEstDC = nan(ofdmInfo.FFTLength, 1); % Channel estimate with DC nan

alpha = 2; % Updating parameter
beta = 3; % Number of subcarriers used for frequency smoothing

for cntr = 1:numOFDMSym
    % Equalization
    eqDataSymTemp = wlan.internal.wlanEqualize(ofdmDemodDataChanEst(:, cntr), chanEstData, eqMethod, noiseVarEst);
    
    % Constellation demapping
    qamDemodOut = wlanConstellationDemap(eqDataSymTemp, noiseVarEst, mcsTable.NBPSCS, 'hard');
    
    % Constellation mapping
    mappedSymb = wlanConstellationMap(qamDemodOut, mcsTable.NBPSCS);
    
    % Channel estimation
    hEst(:, cntr) = ofdmDemodDataChanEst(:, cntr)./mappedSymb;
    
    % Combining subcarriers for Frequency smoothing
    chanEstDC(ofdmInfo.DataIndices) = hEst(:, cntr);
    chanEstDC(ofdmInfo.PilotIndices) = chanEstPilots(:, cntr);
    
    % Frequency smoothing
    hUpdate = movmean(chanEstDC, beta, 'omitnan');
    
    % Time averaging
    hEst(:, cntr) = (1/alpha)*hUpdate(ofdmInfo.DataIndices)+(1-1/alpha)*chanEstData;
    chanEstData = hEst(:, cntr);
end
end