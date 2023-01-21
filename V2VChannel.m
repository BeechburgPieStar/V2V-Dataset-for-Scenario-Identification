classdef V2VChannel < matlab.System
%V2VChannel Filter input signal through an 802.11p fading channel
%   chan = V2VChannel creates a System Object, chan, for the time and
%   frequency selective multipath Rayleigh fading channel as specified by
%   Alexander et al in [1]. This object filters a real or complex input
%   signal through the multipath, vehicle-to-vehicle (V2V) channel to
%   obtain the channel impaired signal.
%
%   chan = V2VChannel(Name,Value) creates a V2V channel object, chan,
%   with the specified property Name set to the specified Value. You can
%   specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
% 
%   Step method syntax:
% 
%   Y = step(chan,X) filters input signal X through a V2V channel and
%   returns the result in Y. The input X must be a double precision, real
%   or complex matrix of size Ns-by-1 where Ns is the number of samples and
%   the singleton dimension corresponds to a single transmit antenna in
%   802.11p. Y is the output signal of size Ns-by-1, where the singleton
%   dimension corresponds to a single receive antenna at the receiver. Y is
%   of double precision data type with complex values.
%   
%   [Y,PATHGAINS] = step(chan,X) returns the channel path gains of the
%   Rayleigh fading process in PATHGAINS. PATHGAINS is of size Ns-by-Np,
%   where Np is the number of resolvable paths, that is, the number of
%   paths defined for the case specified by the DelayProfile property.
%   PATHGAINS is of double precision data type with complex values.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(obj, x) and y = obj(x) are
%   equivalent.
% 
%   V2VChannel methods:
% 
%   step     - Filter the input signal through a Vehicular channel (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create V2V channel object with same property values
%   isLocked - Locked status (logical)
%   <a href="matlab:help matlab.System/reset   ">reset</a>    - Reset states of filters, and random stream if the
%              RandomStream property is set to 'mt19937ar with seed'
% 
%   V2VChannel properties:
%
%   SampleRate              - Input signal sample rate (Hz)
%   DelayProfile            - Vehicular channel delay profile
%   NormalizePathGains      - Normalize path gains (logical)
%   RandomStream            - Source of random number stream
%   Seed                    - Initial seed of mt19937ar random number stream
%   
%   % References
%   % [1] P.Alexander, D. Haley and A. Grant, "Cooperative
%   % Intelligent Transport Systems: 5.9-GHz Field Trials," in Proceedings of
%   % the IEEE, vol. 99, no. 7, pp. 1213-1235, July 2011.

% Copyright 2018 The MathWorks, Inc.
    
properties (Nontunable)
    %SampleRate Sample rate (Hz)
    %   Specify the sample rate of the input signal in Hz as a double
    %   precision, real, positive scalar. The default is 10e6 Hz.
    SampleRate = 1e7;
    %Delay Profile Channel Profile
    %   Specify Channel profile as one of 'Rural LOS'|'Urban approaching
    %   LOS'|'Urban NLOS'|'Highway LOS'| 'Highway NLOS'.The default is
    %   'Rural LOS'
    DelayProfile = 'Rural LOS';
end

properties (Nontunable, Logical)
    %NormalizePathGains Normalize average path gains to 0 dB
    %   Set this property to true to normalize the fading processes such
    %   that the total power of the path gains, averaged over time, is 0dB.
    %   The default is true.
    NormalizePathGains = true;
end

properties (Nontunable)
    %RandomStream Source of random number stream
    %   Specify the source of random number stream as one of 'Global
    %   stream' | 'mt19937ar with seed'. If RandomStream is set to 'Global
    %   stream', the current global random number stream is used for
    %   normally distributed random number generation, in which case the
    %   reset method only resets the filters. If RandomStream is set to
    %   'mt19937ar with seed', the mt19937ar algorithm is used for normally
    %   distributed random number generation, in which case the reset
    %   method not only resets the filters but also re-initializes the
    %   random number stream to the value of the Seed property. The default
    %   value of this property is 'Global stream'.
    RandomStream = 'Global stream';
    %Seed Initial seed
    %   Specify the initial seed of a mt19937ar random number generator
    %   algorithm as a double precision, real, nonnegative integer scalar.
    %   This property applies when you set the RandomStream property to
    %   'mt19937ar with seed'. The Seed is to re-initialize the mt19937ar
    %   random number stream in the reset method. The default value of this
    %   property is 73.
    Seed = 73;
end

properties(Constant, Hidden)
    DelayProfileSet = matlab.system.StringSet({'Rural LOS', ...
        'Urban approaching LOS', 'Urban NLOS', 'Highway LOS', ...
        'Highway NLOS'});
end

properties (Nontunable, Access = private)
    pChannel1;
    pChannel2;   
    pPathGainNormFactor;
end

methods
  function obj = V2VChannel(varargin) % Constructor
    setProperties(obj, nargin, varargin{:});
  end

  function set.SampleRate(obj, Rs)
    propName = 'SampleRate';
    validateattributes(Rs, {'double'}, ...
        {'real','scalar','positive','finite'}, ...
        [class(obj) '.' propName], propName);  

    obj.SampleRate = Rs;
  end
  
  function set.Seed(obj, seed)
    propName = 'Seed';
    validateattributes(seed, {'double'}, ...
        {'real','scalar','integer','nonnegative','finite'}, ...
        [class(obj) '.' propName], propName);

    obj.Seed = seed;
  end
end

methods(Access = protected)
  function validateInputsImpl(obj, x, varargin)
    validateattributes(x, {'double','single'}, {'column','finite'}, ...
        class(obj), 'signal input');
  end
  
  function setupImpl(obj)
    switch obj.DelayProfile
      case 'Rural LOS'
        pathDelays = [0 83 183]*1e-9;
        avgPathGains = [0 -14 -17];
        dopplerShifts = [0 492 -295];
      case 'Urban approaching LOS'
        pathDelays = [0 117 183 333]*1e-9;
        avgPathGains = [0 -8 -10 -15];
        dopplerShifts = [0 236 -157 492];
      case 'Urban NLOS'
        pathDelays = [0 267 400 533]*1e-9;
        avgPathGains = [0 -3 -5 -10];
        dopplerShifts = [0 295 -98 591];
      case 'Highway LOS'
        pathDelays = [0 100 167 500]*1e-9;
        avgPathGains = [0 -10 -15 -20];
        dopplerShifts = [0 689 -492 886];
      case 'Highway NLOS'
        pathDelays = [0 200 433 700]*1e-9;
        avgPathGains = [0 -2 -5 -7];
        dopplerShifts = [0 689 -492 886];
    end
    maxDopplerShift = max(abs(dopplerShifts));
    numPaths = length(pathDelays);
    dopplerSpectrum = cellfun(@(x)doppler('Asymmetric Jakes', sort([0 x/maxDopplerShift])), ...
        num2cell(dopplerShifts(2:end)), 'UniformOutput', false);
    
    % Configure a channel for the 1st static tap
    obj.pChannel1 = comm.RayleighChannel(...
        'SampleRate',          obj.SampleRate, ...
        'PathDelays',          pathDelays, ...
        'AveragePathGains',    avgPathGains - [0 200*ones(1, numPaths-1)], ...
        'NormalizePathGains',  false, ...
        'MaximumDopplerShift', 0, ...
        'PathGainsOutputPort', true, ...
        'RandomStream',        obj.RandomStream);

    % Configure a channel for all the rest taps
    obj.pChannel2 = comm.RayleighChannel(...
        'SampleRate',          obj.SampleRate, ...
        'PathDelays',          pathDelays, ...
        'AveragePathGains',    avgPathGains - [200 zeros(1, numPaths-1)], ...
        'NormalizePathGains',  false, ...
        'MaximumDopplerShift', maxDopplerShift, ...
        'DopplerSpectrum',     [{doppler('Jakes')}, dopplerSpectrum], ...
        'PathGainsOutputPort', true, ...
        'RandomStream',        obj.RandomStream);
    
    if ~strcmp(obj.RandomStream, 'Global stream')
        obj.pChannel1.Seed = obj.Seed;
        obj.pChannel2.Seed = obj.Seed;
    end
    
    % Set up normalization factor
    if obj.NormalizePathGains
        obj.pPathGainNormFactor = sum((10.^(avgPathGains/10)));
    end
  end
  
  function resetImpl(obj)
    reset(obj.pChannel1);
    reset(obj.pChannel2);
  end
  
  function [y, g] = stepImpl(obj, x)
     [y1, g1] = step(obj.pChannel1, x);
     [y2, g2] = step(obj.pChannel2, x);
     
     g = [g1(:,1), g2(:,2:end)];
     y = y1 + y2;
     
     if obj.NormalizePathGains
         g = g/obj.pPathGainNormFactor;
         y = y/obj.pPathGainNormFactor;
     end
  end
  
  function releaseImpl(obj)
    release(obj.pChannel1);
    release(obj.pChannel2);
  end
  
  function flag = isInactivePropertyImpl(obj, prop)
    flag = strcmp(prop, 'Seed') && strcmp(obj.RandomStream, 'Global stream');
  end
  
  function s = saveObjectImpl(obj)
    s = saveObjectImpl@matlab.System(obj);
    if isLocked(obj)
        s.pChannel1 = matlab.System.saveObject(obj.pChannel1);
        s.pChannel2 = matlab.System.saveObject(obj.pChannel2);
        s.pPathGainNormFactor = obj.pPathGainNormFactor;
    end
  end
  
  function loadObjectImpl(obj, s, wasLocked)
    if wasLocked
        obj.pChannel1 = matlab.System.loadObject(s.pChannel1);
        obj.pChannel2 = matlab.System.loadObject(s.pChannel2);
        obj.pPathGainNormFactor = s.pPathGainNormFactor;
    end
    loadObjectImpl@matlab.System(obj, s);    
  end
end

end
