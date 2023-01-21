%% 802.11p Packet Error Rate Simulation for a Vehicular Channel
%
% This example shows how to measure the packet error rate (PER) of an
% IEEE(R) 802.11p(TM) link using an end-to-end simulation with a
% Vehicle-to-Vehicle (V2V) fading channel and additive white Gaussian
% noise. The PER performance of a receiver with and without channel
% tracking is compared. In a vehicular environment (high Doppler), a
% receiver with channel tracking performs better.

% Copyright 2018 The MathWorks, Inc.

%% Introduction
% IEEE 802.11p [ <#10 1> ] is an approved amendment to the IEEE 802.11(TM)
% standard to enable support for wireless access in vehicular environments
% (WAVE). Using the half-clocked mode with a 10 MHz channel bandwidth, it
% operates in 5.85-5.925 GHz bands to support applications for Intelligent
% Transportation Systems (ITS) [ <#10 2> ].
%
% In this example, an end-to-end simulation is used to determine the packet
% error rate for an 802.11p [ <#10 1> ] link with a fading channel at a
% selection of SNR points with and without channel tracking. For each SNR
% point, multiple packets are transmitted through a V2V channel,
% demodulated and the PSDUs are recovered. The PSDUs are compared to those
% transmitted to determine the number of packet errors. For each packet,
% packet detection, timing synchronization, carrier frequency offset
% correction and phase tracking are performed at the receiver. For channel
% tracking, decision directed channel estimation [ <#10 3> ] is used to
% compensate for the high Doppler spread. The figure below shows the
% processing chain with channel tracking.
%
% <<../nonHTPERSchematicChanTracking.png>>

%% Waveform Configuration
% An 802.11p non-HT format transmission is simulated in this example. A
% non-HT format configuration object contains the format specific
% configuration of the transmission. This object is created using the
% <matlab:doc('wlanNonHTConfig') wlanNonHTConfig> function. In this
% example, the object is configured for a 10 MHz channel bandwidth and QPSK
% rate 1/2 (MCS 2) operation.
    
% Link parameters
mcs = 2;       % QPSK rate 1/2
psduLen = 500; % PSDU length in bytes

% Create a format configuration object for an 802.11p transmission
cfgNHT = wlanNonHTConfig;
cfgNHT.ChannelBandwidth = 'CBW10';
cfgNHT.PSDULength = psduLen;              
cfgNHT.MCS = mcs;                     

%% Channel Configuration
% The V2V radio channel model defines five scenarios to represent fading
% conditions within a vehicular environment. In this example, 'Urban NLOS'
% [ <#10 4> ] scenario is used. This corresponds to a scenario with two
% vehicles crossing each other at an urban blind intersection with building
% and fences present on the corners.

% Create and configure the channel
fs = wlanSampleRate(cfgNHT); % Baseband sampling rate for 10 MHz

chan = V2VChannel;
chan.SampleRate = fs;
chan.DelayProfile = 'Urban NLOS';

%% Simulation Parameters
% For each SNR (dB) point in the vector |snr| a number of packets are
% generated, passed through a channel and demodulated to determine the
% packet error rate.

snr = 15:5:30;

%%
% The number of packets tested at each SNR point is controlled by two
% parameters:
%
% # |maxNumErrors| is the maximum number of packet errors simulated at each
% SNR point. When the number of packet errors reaches this limit, the
% simulation at this SNR point is complete.
% # |maxNumPackets| is the maximum number of packets simulated at each SNR
% point. It limits the length of the simulation if the packet error limit
% is not reached.
%
% The numbers chosen in this example lead to a short simulation. For
% statistical meaningful results these numbers should be increased.

maxNumErrors = 20;   % The maximum number of packet errors at an SNR point
maxNumPackets = 200; % Maximum number of packets at an SNR point

% Set random stream for repeatability of results
s = rng(98);

%% Processing SNR Points
% For each SNR point, a number of packets are tested and the packet error
% rate is calculated. For each packet the following processing steps occur:
%
% # A PSDU is created and encoded to create a single packet waveform.
% # The waveform is passed through the channel. Different channel
% realizations are used for each transmitted packet.
% # AWGN is added to the received waveform to create the desired average
% SNR per subcarrier after OFDM demodulation. comm.AWGNChannel is
% configured to provide the correct SNR. The configuration accounts for
% normalization within the channel by the number of receive antennas, and
% the noise energy in unused subcarriers which are removed during OFDM
% demodulation.
% # The per-packet processing includes packet detection, coarse carrier
% frequency offset estimation and correction, symbol timing and fine
% carrier frequency offset estimation and correction.
% # The L-LTF is extracted from the synchronized received waveform. The
% L-LTF is OFDM demodulated and initial channel estimates are obtained.
% # Channel tracking can be enabled using the switch |enableChanTracking|.
% If enabled, the channel estimates obtained from L-LTF are updated per
% symbol using decision directed channel tracking as presented in J. A.
% Fernandez et al in [ <#10 3> ]. If disabled, the initial channel
% estimates from L-LTF are used for the entire packet duration.
% # The non-HT Data field is extracted from the synchronized received
% waveform. The PSDU is recovered using the extracted data field and the
% channel estimates and noise power estimate.

% Set up a figure for visualizing PER results
h = figure;
grid on;
hold on;
ax = gca;
ax.YScale = 'log';
xlim([snr(1), snr(end)]);
ylim([1e-3 1]);
xlabel('SNR (dB)');
ylabel('PER');
h.NumberTitle = 'off';
h.Name = '802.11p ';
title(['MCS ' num2str(mcs) ', V2V channel - ' chan.DelayProfile ' profile']);

% Simulation loop for 802.11p link
S = numel(snr);
per_LS = zeros(S,1);
per_STA = per_LS;
for i = 1:S
    enableChanTracking = true;
    % 802.11p link with channel tracking
    per_STA(i) = v2vPERSimulator(cfgNHT, chan, snr(i), ...
        maxNumErrors, maxNumPackets, enableChanTracking);
    
    enableChanTracking = false;
    % 802.11p link without channel tracking
    per_LS(i) = v2vPERSimulator(cfgNHT, chan, snr(i), ...
        maxNumErrors, maxNumPackets, enableChanTracking);
    
    semilogy(snr, per_STA, 'bd-');
    semilogy(snr, per_LS, 'ro--');
    legend('with Channel Tracking','without Channel Tracking')
    drawnow;
end

axis([10 35 1e-3 1])
hold off;

% Restore default stream
rng(s);

%% 
% For meaningful results |maxNumErrors|, |maxNumPackets| should be
% increased. The below plot provides results for |maxNumErrors: 1000| and
% |maxNumPackets: 10000|.
%
% <<../V2VPERresultsUrbanNLOS.png>>

%% Further Exploration
% Try changing the channel delay profile, the length of the packet or the data
% rate ( |mcs| values ) and observe the performance of channel tracking.
% For some configurations channel tracking provides little performance
% improvement. For a small number of OFDM symbols (small PSDU length or
% high MCS), temporal averaging performed during decision directed channel
% tracking may not be effective. The characteristics of the channel may
% also limit the performance for higher order modulation schemes ( |mcs| >
% 5 ).
%% Appendix
% This example uses the following helper functions and objects:
%
% * <matlab:edit('v2vPERSimulator') v2vPERSimulator.m>
% * <matlab:edit('V2VChannel') V2VChannel.m>

%% Selected Bibliography
% # IEEE Std 802.11p-2010: IEEE Standard for Information technology -
% Telecommunications and information exchange between systems - Local and
% metropolitan area networks - Specific requirements, Part 11: Wireless LAN
% Medium Access Control (MAC) and Physical Layer (PHY) Specifications,
% Amendment 6: Wireless Access in Vehicular Environments, IEEE, New York,
% NY, USA, 2010.
% # ETSI,
% https://www.etsi.org/technologies-clusters/technologies/automotive-intelligent-transport.
% # J. A. Fernandez, D. D. Stancil and F. Bai, "Dynamic channel
% equalization for IEEE 802.11p waveforms in the vehicle-to-vehicle
% channel," 2010 48th Annual Allerton Conference on Communication, Control,
% and Computing (Allerton), Allerton, IL, 2010, pp. 542-551. doi:
% 10.1109/ALLERTON.2010.5706954
% # P. Alexander, D. Haley and A. Grant, "Cooperative Intelligent Transport
% Systems: 5.9-GHz Field Trials," in Proceedings of the IEEE, vol. 99, no.
% 7, pp. 1213-1235, July 2011.
