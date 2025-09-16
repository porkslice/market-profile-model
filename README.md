# market-profile-model
A Bayesian, real-time Market-Profile that models a session as a scale-mixture density over price. 
Scales come from a wave (MODWT) decomposition of the profile. 
The density is smoothed by convolution (σ adapts to noise). 
Parameters (scale weights w and smoothing σ) are fit with Metropolis–Hastings (MH) to get uncertainty bands. 
Two context modules—Opening Type and Day Type—output probabilities (not hard labels) that steer priors, smoothing, and a real-time covariance panel. 
