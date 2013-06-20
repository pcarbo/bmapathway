% INITRNG(SEED) initializes the random number generator. While there is a
% simpler way to do this in MATLAB 7.14 (R2012a) using the function RNG, I
% set the random number generator seed in this way here because these
% function calls are compatible with some earlier versions of MATLAB. 
function initrng (seed)

  % Here, 'mt19937ar' refers to the Mersenne Twister, which is the default
  % random number generator in MATLAB 7.14.
  warning('off','MATLAB:RandStream:SetDefaultStream');
  s = RandStream('mt19937ar','Seed',seed);
  RandStream.setDefaultStream(s);
