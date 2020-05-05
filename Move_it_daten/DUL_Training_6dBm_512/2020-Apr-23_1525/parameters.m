% Parameter File for Model DUL_Training 
% Modified 23-Apr-2020 

% Common Simulation Parameters 

common.loopmax=1;                                  %number of iterations
common.numsym=32768;                              
common.dimension=6;                               
common.f_sym=64000000000;                         
common.f_sim=32*common.f_sym;                     
common.K=2;                                       
common.f_adc=common.K*common.f_sym;               
common.alpha=0.3;                                 
common.trunc=256;                                 
common.plaunch=0;                                 
common.num_span=1;                                
common.Vpp=2.2;                                   
common.volt_train=512;                            
common.taps_NL2=1;                                
common.taps_NL3=4;                                
common.NL2_on=1;                                  
common.buffer=100;                                
                         



% prms

prms_para.bl=common.numsym;                        %Block length (Length of output vector)
prms_para.rand=1;                                  %Random (True) or pseudo-random (False)
prms_para.dimension=common.dimension;              %dimension (i.e. number of parallel bit streams)
prms_para.order=10;                                %order of PRMS (needs to be adapted to memory of DUT)
prms_para.skip=0;                                  %Skip the first n Bits of sequence
prms_para.bruijn=1;                                %extends the PRBS sequence of length 2^n-1 to a deBruijn sequence of 2^n
prms_para.reset_prms=0;                            %reset the sequence generator every loop
prms_para.method=1;                                %false= old method , true = new method

% digi_mod

digi_mod_para.format=0;                            %real ASK Modulation (True) or complex QAM modulation (False)
digi_mod_para.ask_unipolar=0;                      %true: unipolar ASK constellation. false: bipolar ASK constellation
digi_mod_para.PDM=0;                               %true: PDM. false: no PDM
digi_mod_para.order=common.dimension;              %modulation order
digi_mod_para.scale_mode=2;                        %scale constellation power to 1

% pulsef

pulsef_para.enabled=1;                             %enable/disable
pulsef_para.f_sym=common.f_sym;                    %symbol frequency (input) 
pulsef_para.fs=common.f_sim;                       %sampling frequency
                                                   %needs to be a multiple of f_sym
pulsef_para.pulse=2;                               %1: Raised-cosine
                                                   %2: square root of raised-cosine
                                                   %3: Rectangle
                                                   %4: Sine Trapezoid
                                                   %5: Trapezoid
                                                   %6: Half sine
                                                   %7: Gauss
                                                   %8: External pulse from data file (time domain)
                                                   %9: External pulse from data file (frequency domain)
pulsef_para.conv=1;                                %0: linear convolution with overlap-add; 1: cyclic convolution
pulsef_para.method=1;                              %fast convolution, optimized fast convolution via fft
                                                   %or conventional convolution in the time domain
                                                   %ignored in case of cyclic convolution
pulsef_para.window=0;                              %1: Windowing of pulse,  0: No windowing  
pulsef_para.non_int=0;                             %non integer resampling for racos 
pulsef_para.amp=-1;                                %Amplitude of pulse in V or A
pulsef_para.jitoff=0;                              %Jitter offset in [s]
                                                   %for linear convolution no negative jitter is possible
                                                   %for convolution in time domain jitter is ignored
pulsef_para.jitvar=0;                              %Jitter (standard deviation, normally distributed) in [s]
                                                   %for linear convolution no negative jitter is possible
                                                   %for convolution in time domain jitter will be ignored
pulsef_para.alpharacos=common. alpha;              %Roll-off factor (alpha) of raised-cosine (racos) filter
pulsef_para.truncracos=common. trunc;              %truncate the impulse response to a number of symbol periods according to this parameter
pulsef_para.recwidth=1;                            %normalized width of rectangle pulse
pulsef_para.risetrapez=0.7;                        %normalized rise time
pulsef_para.durtrapez=0.3;                         %normalized duration time
pulsef_para.falltrapez=0.7;                        %normalized fall time
pulsef_para.halfsine=0.6;                          %normalized width of halfsine
pulsef_para.widthgauss=0.2;                        %FWHM of Gaussian pulse
pulsef_para.truncgauss=8;                          %truncate the impulse response to a multiple of FWHM according to this parameter
pulsef_para.fileextern='';                         %File name (Matlab '.mat'-file)
                                                   %File MUST contain a vector named 'pulse'.
pulsef_para.tfextern=0.1;                          %Pulse duration    [s]: time domain       OR
                                                   %Max. frequency   [Hz]: frequency domain

% user_function

user_function_para.enable=1;                       %enable?
user_function_para.func_str='x./max(abs(x)).*common.Vpp./2'; %user defined function with argument x
user_function_para.interface=1;                    %Output interface

% eml

eml_para.mode=5;                                   %mode
eml_para.lambda=1550;                              %Emitted wavelength in [nm]
eml_para.f=0;                                      %Emitted frequency in [Hz]; Set either wavelength or frequency to zero
eml_para.power=0;                                  %Emitted power in [dBm]
eml_para.interface=2;                              %Output interface
eml_para.linewidth=0;                              %Laser linewidth in [Hz]
eml_para.ampl_imbal=0;                             %Amplitude Imbalance (-1:1) (only valid for IQ modulators)
eml_para.pha_imbal=0.07854;                        %Phase Imbalance (-pi/2:pi/2) (only valid for IQ modulators)
eml_para.fs=common. f_sim;                         %sampling frequency in [Hz]
                                                   %only required for nonzero linewidth
eml_para.bias=-1;                                  %bias
eml_para.U_pi=1;                                   %U_pi

% optatten

optatten_para.atten=-10;                           %Attenuation in [dB]
optatten_para.mode=2;                              %attenuation [dB] or output power [dBm]
optatten_para.optout=2;                            %type of optical output signal

% edfa

edfa_para.gain=common.plaunch;                     %small-signal gain in [dB]
edfa_para.mode=2;                                  %gain [dB] or output power [dBm]
edfa_para.sat=1;                                   %ideal gain or saturated gain
edfa_para.Psat=10;                                 %saturating output power in [dBm]
edfa_para.nf=5;                                    %noise figure of EDFA in [dB]
edfa_para.optout=2;                                %type of optical output signal
edfa_para.ase_enable=1;                            %1: enable ASE or 0: diable ASE
edfa_para.signal_off=0;                            %1: set signal to zero or 0: DONT
edfa_para.fs=common.f_sim;                         %sampling frequency in [Hz]

% spans_no_disp2

spans_no_disp2_para.number_of_spans=common.num_span; %Number of Spans
spans_no_disp2_para.pre_amp_enable=0;              %Enable pre-Amplifier
spans_no_disp2_para.pre_ase_conv_enable=0;         %Enable pre-ASE conv.
spans_no_disp2_para.pre_opt_filt_enable=0;         %Enable pre-optical filter.
spans_no_disp2_para.fiber1_enable=0;               %Enable 1st Fiber (eg. Main Fiber)
spans_no_disp2_para.fiber_adapt_enable=1;          %Adaptive Fiber
spans_no_disp2_para.post_amp_enable=1;             %Enable Post-Amplifier
spans_no_disp2_para.post_ase_conv_enable=0;        %Enable post-ASE conv.
spans_no_disp2_para.post_opt_filt_enable=0;        %Enable post-optical filter.
spans_no_disp2_para.pre_amp_gain=30;               %small-signal gain in [dB]
spans_no_disp2_para.pre_amp_mode=1;                %gain [dB] or output power [dBm]
spans_no_disp2_para.pre_amp_sat=1;                 %ideal gain or saturated gain
spans_no_disp2_para.pre_amp_Psat=10;               %saturating output power in [dBm]
spans_no_disp2_para.pre_amp_nf=5;                  %noise figure of EDFA in [dB]
spans_no_disp2_para.pre_amp_optout=2;              %type of optical output signal
spans_no_disp2_para.pre_amp_ase_enable=0;          %1: enable ASE or 0: diable ASE
spans_no_disp2_para.pre_amp_signal_off=0;          %1: set signal to zero or 0: DONT
spans_no_disp2_para.pre_amp_fs=320000000000;       %sampling frequency in [Hz]
spans_no_disp2_para.pre_ase_fs=320000000000;       %sampling frequency in [Hz]
spans_no_disp2_para.pre_ase_bw=1;                  %Determines, if generated noise covers full or half sampling bandwidth
spans_no_disp2_para.pre_ase_conv=0;                %linear (false) or cyclic (true) convolution for half bandwidth noise
spans_no_disp2_para.pre_ase_coeff=4;               %number of coefficients for interpolation filter
spans_no_disp2_para.pre_ase_sigdelay=1;            %delay data signal according to delay of interpolation filter
spans_no_disp2_para.pre_ase_modes=1;               %number of ASE modes
spans_no_disp2_para.pre_ase_optout=2;              %type of optical output signal
spans_no_disp2_para.pre_filt_fs=320000000000;      %sampling frequency in [Hz]
spans_no_disp2_para.pre_filt_conv=1;               %linear (false) or cyclic (true) convolution
spans_no_disp2_para.pre_filt_fast=2;               %in case of linear convolution, use fast, optimized fast or time-domain convolution
spans_no_disp2_para.pre_filt_symreal=0;            %transfer function is symmetric (i.e. impulse response is real) 
                                                   %results in faster computation for time domain convolution in C-simulation
spans_no_disp2_para.pre_filt_filtype=1;            %filter type
spans_no_disp2_para.pre_filt_shift=0;              %frequency shift of transfer function in [Hz]
spans_no_disp2_para.pre_filt_interface=2;          %Output interface
spans_no_disp2_para.pre_filt_B=10000000000;        %3dB-bandwidth for gaussian filter in [Hz]
spans_no_disp2_para.pre_filt_mgauss=1;             %order of Gaussian filter, must be a multiple of 1/2
spans_no_disp2_para.pre_filt_Brect=10000000000;    %two-sided bandwidth for rectangular filter in [Hz]
spans_no_disp2_para.pre_filt_fbglength=0.005;      %Length of Bragg Gratting in [m]
spans_no_disp2_para.pre_filt_fbgkappa=200;         %kappa in [1/m]
spans_no_disp2_para.pre_filt_fbgvg=200000000;      %Group velocity in [m/s]
spans_no_disp2_para.pre_filt_m=200;                %number of coefficients for impulse response
spans_no_disp2_para.fst_fiber_fa=320000000000;     %Sampling frequency in Hz
spans_no_disp2_para.fst_fiber_conv=1;              %linear (false) or cyclic (true) convolution
spans_no_disp2_para.fst_fiber_len=100;             %length of fiber in km
spans_no_disp2_para.fst_fiber_alpha=0.2;           %attenuation coefficient in dB/km
spans_no_disp2_para.fst_fiber_D=17;                %dispersion coefficient in ps/(nm*km)
spans_no_disp2_para.fst_fiber_Ds=0.06;             %dispersion slope in ps/(nm²*km)
spans_no_disp2_para.fst_fiber_lambda0=1550;        %reference wavelength for dispersion coefficient and dispersion slope
spans_no_disp2_para.fst_fiber_nlin=1;              %include non-linearities
spans_no_disp2_para.fst_fiber_optout=2;            %type of optical output signal
spans_no_disp2_para.fiber_b3=1;                    %include beta_3
spans_no_disp2_para.fiber_ssfmstep=100;            %step size for Split-Step-Fourier in m
spans_no_disp2_para.fiber_n2=3.2e-20;              %Nonlinear refractive Index n2 in m²/W
spans_no_disp2_para.fiber_Aeff=8e-11;              %Effective Core Area in m²
spans_no_disp2_para.fiber_gamma_1=0;               %nonlinear coefficient of gamma, [1/w/m]
spans_no_disp2_para.fiber_spm_dgl=1;               %true:  consider SPM
                                                   %false: neglect SPM
spans_no_disp2_para.fiber_xpm_dgl=3;               %determines, if XPM should be calculated for separated channels approach
spans_no_disp2_para.fiber_xpm_dgl_eff=0;           %determines the number of XPM neighbours taken into account for 'selected tones'
spans_no_disp2_para.fiber_linmethod=1;             %method for linear convolution
spans_no_disp2_para.fiber_firlength=1024;          %number of samples for impulse response
spans_no_disp2_para.fiber_iirgdslope=0.4;          %normalized group delay slope per IIR-filter
spans_no_disp2_para.fiber_iircascades=16;          %maximal number of cascaded IIR-Filters
spans_no_disp2_para.fiber_adapt_fa=common. f_sim;  %Sampling frequency in Hz
spans_no_disp2_para.fiber_adapt_L=100;             %Length of Fiber
spans_no_disp2_para.fiber_adapt_lambda=1550;       %reference wavelength
spans_no_disp2_para.fiber_adapt_gamma=0.0013;      %nonlinear coefficient [1/W/m]
spans_no_disp2_para.fiber_adapt_alpha=0.2;         %alpha
spans_no_disp2_para.fiber_adapt_D=17;              %dispersion coefficient [ps/(nm*km)]
spans_no_disp2_para.fiber_adapt_Ds=0.06;           %dispersion slope [ps/(nm²*km)]
spans_no_disp2_para.fiber_adapt_SS_dphimax=0.0005; %maximum nonlinear rotation angle in nonlinear step
spans_no_disp2_para.fiber_adapt_optout=2;          %type of optical output signal
spans_no_disp2_para.post_amp_gain=common. plaunch; %small-signal gain in [dB]
spans_no_disp2_para.post_amp_mode=2;               %gain [dB] or output power [dBm]
spans_no_disp2_para.post_amp_sat=1;                %ideal gain or saturated gain
spans_no_disp2_para.post_amp_Psat=10;              %saturating output power in [dBm]
spans_no_disp2_para.post_amp_nf=5;                 %noise figure of EDFA in [dB]
spans_no_disp2_para.post_amp_optout=2;             %type of optical output signal
spans_no_disp2_para.post_amp_ase_enable=1;         %1: enable ASE or 0: diable ASE
spans_no_disp2_para.post_amp_signal_off=0;         %1: set signal to zero or 0: DONT
spans_no_disp2_para.post_amp_fs=common. f_sim;     %sampling frequency in [Hz]
spans_no_disp2_para.post_ase_fs=320000000000;      %sampling frequency in [Hz]
spans_no_disp2_para.post_ase_bw=1;                 %Determines, if generated noise covers full or half sampling bandwidth
spans_no_disp2_para.post_ase_conv=0;               %linear (false) or cyclic (true) convolution for half bandwidth noise
spans_no_disp2_para.post_ase_coeff=4;              %number of coefficients for interpolation filter
spans_no_disp2_para.post_ase_sigdelay=1;           %delay data signal according to delay of interpolation filter
spans_no_disp2_para.post_ase_modes=1;              %number of ASE modes
spans_no_disp2_para.post_ase_optout=2;             %type of optical output signal
spans_no_disp2_para.post_filt_fs=common. f_sim;    %sampling frequency in [Hz]
spans_no_disp2_para.post_filt_conv=1;              %linear (false) or cyclic (true) convolution
spans_no_disp2_para.post_filt_fast=2;              %in case of linear convolution, use fast, optimized fast or time-domain convolution
spans_no_disp2_para.post_filt_symreal=0;           %transfer function is symmetric (i.e. impulse response is real) 
                                                   %results in faster computation for time domain convolution in C-simulation
spans_no_disp2_para.post_filt_filtype=2;           %filter type
spans_no_disp2_para.post_filt_shift=0;             %frequency shift of transfer function in [Hz]
spans_no_disp2_para.post_filt_interface=2;         %Output interface
spans_no_disp2_para.post_filt_B=10000000000;       %3dB-bandwidth for gaussian filter in [Hz]
spans_no_disp2_para.post_filt_mgauss=1;            %order of Gaussian filter, must be a multiple of 1/2
spans_no_disp2_para.post_filt_Brect=common. f_sym*1.4; %two-sided bandwidth for rectangular filter in [Hz]
spans_no_disp2_para.post_filt_fbglength=0.005;     %Length of Bragg Gratting in [m]
spans_no_disp2_para.post_filt_fbgkappa=200;        %kappa in [1/m]
spans_no_disp2_para.post_filt_fbgvg=200000000;     %Group velocity in [m/s]
spans_no_disp2_para.post_filt_m=200;               %number of coefficients for impulse response

% fil_opt

fil_opt_para.fs=common.f_sim;                      %sampling frequency in [Hz]
fil_opt_para.conv=1;                               %linear (false) or cyclic (true) convolution
fil_opt_para.fast=2;                               %in case of linear convolution, use fast, optimized fast or time-domain convolution
fil_opt_para.symreal=0;                            %transfer function is symmetric (i.e. impulse response is real) 
                                                   %results in faster computation for time domain convolution in C-simulation
fil_opt_para.filtype=2;                            %filter type
fil_opt_para.shift=0;                              %frequency shift of transfer function in [Hz]
fil_opt_para.interface=2;                          %Output interface
fil_opt_para.B=10000000000;                        %3dB-bandwidth for gaussian filter in [Hz]
fil_opt_para.mgauss=1;                             %order of Gaussian filter, must be a multiple of 1/2
fil_opt_para.Brect=common.f_sym*1.4;               %two-sided bandwidth for rectangular filter in [Hz]
fil_opt_para.fbglength=0.005;                      %Length of Bragg Gratting in [m]
fil_opt_para.fbgkappa=200;                         %kappa in [1/m]
fil_opt_para.fbgvg=200000000;                      %Group velocity in [m/s]
fil_opt_para.m=200;                                %number of coefficients for impulse response

% dp_co_rx

dp_co_rx_para.lo_power=0;                          %Emitted power of local laser in [dBm]
dp_co_rx_para.f_offset=0;                          %Offset of LO frequency in Hz [dBm]
dp_co_rx_para.linewidth=0;                         %Laser linewidth in [Hz]
dp_co_rx_para.fsim=common.f_sim;                   %sampling frequency in [Hz]
                                                   %only required for nonzero linewidth
dp_co_rx_para.output=2;                            %output method
dp_co_rx_para.responsivity=1;                      %responsivity of photo diodes in [A/W]
dp_co_rx_para.shot=0;                              %include shot noise
dp_co_rx_para.PSD_therm=0;                         %single-sided power spectral density of thermal noise in [A^2/Hz]
dp_co_rx_para.ampl_imbal=[0 0 ];                   %Amplitude Imbalance (-1:1)
dp_co_rx_para.pha_imbal=[0 0 ];                    %Phase Imbalance (-pi/2:pi/2)
dp_co_rx_para.active=0;                            %active
dp_co_rx_para.conv=1;                              %conv
dp_co_rx_para.fixed_delay=1;                       %fixed_delay
dp_co_rx_para.delay=128;                           %delay
dp_co_rx_para.filter=3;                            %filter
dp_co_rx_para.pass=1;                              %pass
dp_co_rx_para.fc=25000000000;                      %fc
dp_co_rx_para.n=3;                                 %n
dp_co_rx_para.Rp=0.5;                              %Rp
dp_co_rx_para.Rs=20;                               %stopband ripple in [dB]
dp_co_rx_para.df=0.01;                             %transition width
dp_co_rx_para.f=[0 0.0313 0.0313 1 ];              %frequency vector
dp_co_rx_para.m=[1 1 0 0 ];                        %magnitude vector
dp_co_rx_para.skew=0;                              %skew on/off
dp_co_rx_para.round=0;                             %1: rounds the delay time towards the next multiple of the sampling period; 
                                                   %0: exact realization of the given delay
dp_co_rx_para.out1_skew=0;                         %out1_skew
dp_co_rx_para.out2_skew=0;                         %out2_skew
dp_co_rx_para.out3_skew=0;                         %out3_skew
dp_co_rx_para.out4_skew=0;                         %out4_skew

% sampling

sampling_para.conv=1;                              %cyclic (1) or linear (0) convolution
sampling_para.randel=0;                            %random sampling delay (Activate/Deactivate)
sampling_para.freq_in=common.f_sim;                %frequency of the input signal [Hz]
sampling_para.freq_out=common.f_adc;               %desired frequency of the output signal [Hz]
sampling_para.freq_offset=0;                       %frequency offset of sampler [Hz]
sampling_para.samp_delay=0;                        %constant sampling delay [s]
sampling_para.samp_jitter=0;                       %random sampling jitter RMS value [s]
sampling_para.P_out=-1;                            %normalize output power to P_out

% matched

matched_para.enabled=1;                            %enable/disable
matched_para.f_sym=common.f_sym;                   %symbol frequency (input) 
matched_para.fs=common.f_adc;                      %sampling frequency
matched_para.conv=1;                               %0: linear convolution with overlap-add; 1: cyclic convolution
matched_para.filter=1;                             %Filter type
matched_para.K_MF=1;                               %Matched Filter Gain
matched_para.alpharacos=common. alpha;             %Roll-off factor (alpha) of raised-cosine (racos) filter
matched_para.truncracos=common. trunc;             %truncate the impulse response to a number of symbol periods according to this parameter
matched_para.recwidth=1;                           %normalized width of rectangle pulse
matched_para.risetrapez=0.7;                       %normalized rise time
matched_para.durtrapez=0.3;                        %normalized duration time
matched_para.falltrapez=0.7;                       %normalized fall time

% edc_module

edc_module_para.enable=1;                          %Turn NFT on/off
edc_module_para.fsim=common.f_adc;                 %Sampling Frequency in Hz
edc_module_para.spanlen=100*common.num_span;       %length of fiber in km
edc_module_para.lambda_T=1550;                     %mean wavelength in [nm]
edc_module_para.D=17;                              %dispersion coefficient in [ps/(nm*km)[
edc_module_para.Ds=0.06;                           %dispersion slope in [ps/(nm²*km)]
edc_module_para.lambda0=1550;                      %Reference wavelength in [nm]

% normalize__2

normalize__2_para.active=1;                        %Activate/Deactivate module
normalize__2_para.rms=1;                           %Normalize blockwise to RMS
normalize__2_para.indx=0;                          %Normalize to make a specific value in the vector with magnitude of one
normalize__2_para.loc=-1;                          %index of value upon the whole signal is to be normalized (if "index" was chosen), if -1 is chosen normalize the max
normalize__2_para.factor=1;                        %multiply signal with this factor

% eq_ffe_volterra_function_coeff_corr__2

eq_ffe_volterra_function_coeff_corr__2_para.active=1; %Activate/Deactivate Module
eq_ffe_volterra_function_coeff_corr__2_para.Ne=4;  %Number of delay taps - 1st order
eq_ffe_volterra_function_coeff_corr__2_para.Nc=48; %Estimated number of channel coefficinets
eq_ffe_volterra_function_coeff_corr__2_para.delay=-1; %Delay of feed forward equalizer
eq_ffe_volterra_function_coeff_corr__2_para.FSE=common.K; %K - Number of samples per symbol (sampling period T/K)
eq_ffe_volterra_function_coeff_corr__2_para.conv=0; %Cyclic or Linear Convolution - if 2nd or 3rd order of Volterra is active, only linear convolution works
eq_ffe_volterra_function_coeff_corr__2_para.conv_blocksize=128; %Blocksize for parallel computing of the linear convolition
eq_ffe_volterra_function_coeff_corr__2_para.ch_dly=0; %Delay of the received signal (in samples!!! not symbols!!!)
eq_ffe_volterra_function_coeff_corr__2_para.train=1; %frequentness of the training sequence
eq_ffe_volterra_function_coeff_corr__2_para.Nt=common. volt_train; %Number of training symbols
eq_ffe_volterra_function_coeff_corr__2_para.save=0; %saving coefficients x loops (x number of loops)
eq_ffe_volterra_function_coeff_corr__2_para.filen=''; %Name of the file to be save in the workspace, which will contain the feed forward and feedback filter coefficients
eq_ffe_volterra_function_coeff_corr__2_para.field_ffe='e'; %Name of the field containing the FFE coefficients
eq_ffe_volterra_function_coeff_corr__2_para.plot_coeff=0; %plot coefficients in figure 1001
eq_ffe_volterra_function_coeff_corr__2_para.store_coeff=0; %store coefficients
eq_ffe_volterra_function_coeff_corr__2_para.nl2=0; %Include nonlinear coefficients of 2nd order
eq_ffe_volterra_function_coeff_corr__2_para.Ne2=common. taps_NL2; %Number of delay taps - 2nd order
eq_ffe_volterra_function_coeff_corr__2_para.nl3=0; %Include nonlinear coefficients of 3rd order
eq_ffe_volterra_function_coeff_corr__2_para.Ne3=common. taps_NL3/2; %Number of delay taps - 3rd order
eq_ffe_volterra_function_coeff_corr__2_para.Ne3_pow=0.001; %Max. power of coefficients
eq_ffe_volterra_function_coeff_corr__2_para.twoReal=0; %Seperation of Real and Imaginary Part of the signal

% user_function__5

user_function__5_para.enable=1;                    %enable?
user_function__5_para.func_str='(x(common.buffer:end-common.buffer))'; %user defined function with argument x
user_function__5_para.interface=1;                 %Output interface

% eq_ffe_volterra_function_coeff_corr

eq_ffe_volterra_function_coeff_corr_para.active=1; %Activate/Deactivate Module
eq_ffe_volterra_function_coeff_corr_para.Ne=4;     %Number of delay taps - 1st order
eq_ffe_volterra_function_coeff_corr_para.Nc=48;    %Estimated number of channel coefficinets
eq_ffe_volterra_function_coeff_corr_para.delay=-1; %Delay of feed forward equalizer
eq_ffe_volterra_function_coeff_corr_para.FSE=common.K; %K - Number of samples per symbol (sampling period T/K)
eq_ffe_volterra_function_coeff_corr_para.conv=0;   %Cyclic or Linear Convolution - if 2nd or 3rd order of Volterra is active, only linear convolution works
eq_ffe_volterra_function_coeff_corr_para.conv_blocksize=128; %Blocksize for parallel computing of the linear convolition
eq_ffe_volterra_function_coeff_corr_para.ch_dly=0; %Delay of the received signal (in samples!!! not symbols!!!)
eq_ffe_volterra_function_coeff_corr_para.train=1;  %frequentness of the training sequence
eq_ffe_volterra_function_coeff_corr_para.Nt=common. volt_train; %Number of training symbols
eq_ffe_volterra_function_coeff_corr_para.save=0;   %saving coefficients x loops (x number of loops)
eq_ffe_volterra_function_coeff_corr_para.filen=''; %Name of the file to be save in the workspace, which will contain the feed forward and feedback filter coefficients
eq_ffe_volterra_function_coeff_corr_para.field_ffe='e'; %Name of the field containing the FFE coefficients
eq_ffe_volterra_function_coeff_corr_para.plot_coeff=0; %plot coefficients in figure 1001
eq_ffe_volterra_function_coeff_corr_para.store_coeff=0; %store coefficients
eq_ffe_volterra_function_coeff_corr_para.nl2=1;    %Include nonlinear coefficients of 2nd order
eq_ffe_volterra_function_coeff_corr_para.Ne2=common. taps_NL2; %Number of delay taps - 2nd order
eq_ffe_volterra_function_coeff_corr_para.nl3=1;    %Include nonlinear coefficients of 3rd order
eq_ffe_volterra_function_coeff_corr_para.Ne3=common. taps_NL3; %Number of delay taps - 3rd order
eq_ffe_volterra_function_coeff_corr_para.Ne3_pow=0.001; %Max. power of coefficients
eq_ffe_volterra_function_coeff_corr_para.twoReal=0; %Seperation of Real and Imaginary Part of the signal

% user_function__3

user_function__3_para.enable=1;                    %enable?
user_function__3_para.func_str='(x(common.buffer:end-common.buffer))'; %user defined function with argument x
user_function__3_para.interface=1;                 %Output interface

% user_function__6

user_function__6_para.enable=1;                    %enable?
user_function__6_para.func_str='(x(common.buffer:end-common.buffer))'; %user defined function with argument x
user_function__6_para.interface=1;                 %Output interface

% dp_anal__2

dp_anal__2_para.mode=1;                            %MODE:
                                                   %0:off
                                                   %1:on
dp_anal__2_para.delay=0;                           %delay in samples
dp_anal__2_para.skip=0;                            %loops to skip
dp_anal__2_para.norm=1;                            %Normalize inputs (data and ref)
dp_anal__2_para.dc_block=0;                        %Normalize inputs
dp_anal__2_para.lab_mode=1;                        %delay in every loop .. delayed symbols are skipped
dp_anal__2_para.carrier=1;                         %Number of used subcarriers (1 = without multicarrier)
dp_anal__2_para.fft=1024;                          %FFT length in samples (for OFDM)
dp_anal__2_para.order=common.dimension;            %modulation order
dp_anal__2_para.dec=3;                             %nearest neighbor
                                                   %kmeans
                                                   %threashold of digidemod
dp_anal__2_para.MOD_Pr=0;                          %MOD_Pr
dp_anal__2_para.MOD=2;                             %MOD
dp_anal__2_para.Polarity=1;                        %polarity of mapping from dp flexmod for correct demodulation of PAM. important for multidimensional optimisations
dp_anal__2_para.cluster=0;                         %Estimate IRR
dp_anal__2_para.fom_evm=1;                         %Estimate EVM
dp_anal__2_para.fom_BER=1;                         %Estimate BER
dp_anal__2_para.fom_Q=0;                           %Estimate Q
dp_anal__2_para.fom_IRR=0;                         %Estimate IRR
dp_anal__2_para.plot_mode=1;                       %correlation
dp_anal__2_para.figure_ind=0;                      %figure index
dp_anal__2_para.nblock=[3 4 ];                     %screen sections
dp_anal__2_para.tile=[0 0 ];                       %figure position (y,x)
dp_anal__2_para.monitor=1;                         %monitor index
dp_anal__2_para.choicew=1;                         %choicew
dp_anal__2_para.hide_axis=1;                       %hide_axis
dp_anal__2_para.plot=0;                            %activate plots
dp_anal__2_para.plot_evm=0;                        %Ploting evm in 2D
dp_anal__2_para.plot_errors=0;                     %Ploting errors in demod and mse over time
dp_anal__2_para.cluster_plot=0;                    %Plot cluster analyses
dp_anal__2_para.cluster_texts=0;                   %Plot cluster analyses
dp_anal__2_para.plot_delay=0;                      %correlation
dp_anal__2_para.plot_flexanal=0;                   %correlation
dp_anal__2_para.berWin=1;                          %Stat window
dp_anal__2_para.store_plots=0;                     %store plots
dp_anal__2_para.sto=1;                             %Point of storage in simulation
dp_anal__2_para.win_update=1;                      %when the display figure is updated
dp_anal__2_para.sto_mode=1;                        %stores cumulated evm values or evms of every loop!
dp_anal__2_para.disp_2D=0;                         %disp_2D
dp_anal__2_para.disp_sym=0;                        %disp_sym
dp_anal__2_para.disp_bits=0;                       %disp_bits
dp_anal__2_para.disp_errs=1;                       %disp_errs
dp_anal__2_para.disp_g=0;                          %disp_g
dp_anal__2_para.disp_phi=0;                        %disp_phi
dp_anal__2_para.FEC_active=0;                      %FEC_active
dp_anal__2_para.FEC_algo=1;                        %FEC Encoding (for each bitstream separately - only implemented for realmod with random data)
dp_anal__2_para.minloops=1e+30;                    %minimal number of simulation loops before module can break the simulation
dp_anal__2_para.minerrors=1e+30;                   %minimal number of bit errors before module breaks the simulatio
dp_anal__2_para.break_ber=1e-05;                   %ber break uloop


% constellation__3

constellation__3_para.plot=1;                      %true: plot eye diagram; false: do not plot eye diagram
constellation__3_para.skip='0';                    %samples to skip
constellation__3_para.plot_color=1;                %color of plot
constellation__3_para.figure_ind=0;                %0=private figure, >0=index of common figure handle
constellation__3_para.monitor=1;                   %monitor index
constellation__3_para.plot_style=0;                %style of plot
constellation__3_para.plotmode=1;                  %plots:
                                                   %-only the current block
                                                   %-all blocks since the last plot
                                                   %-the whole signal calculated by now
constellation__3_para.plothist=1;                  %plot histogram of data
constellation__3_para.color_period=1;              %color bar periodicity
constellation__3_para.cofmagy=4;                   %n
                                                   %µ
                                                   %m
                                                   %1
                                                   %k
                                                   %M
                                                   %G
constellation__3_para.min_y_val=0;                 %min. y-value in plot
constellation__3_para.max_y_val=0;                 %max. y-value in plot
constellation__3_para.fig_plo=0;                   %0 to plot on the same figure
                                                   % 1 for separeted figures
constellation__3_para.markersize=3;                %size of const marker
constellation__3_para.grd=1;                       %determines, if a grid should be plotted
constellation__3_para.hide_axis=0;                 %hide labels from axis and rescale axis
constellation__3_para.plot_ref=1;                  %Plots reference constellation
constellation__3_para.nblock=[3 4 ];               %screen sections
constellation__3_para.tile=[0 0 ];                 %figure position (y,x)
constellation__3_para.choicew=1;                   %window proportion
constellation__3_para.interval=1;                  %show plots every x loop (x=0: at end of simulation)
constellation__3_para.stp=0;                       %stop after each plot
constellation__3_para.plotsize=[0 0 ];             %plot size(x,y)
constellation__3_para.histpoints=301;              %histogram dimension
constellation__3_para.hist_max_val=0;              %histogram max value (0==autoscaling)
constellation__3_para.norm=0;                      %Normalize constellation to 1
constellation__3_para.hist_colorbar=0;             %plot colorbarof histogram
constellation__3_para.store=1;                     %store figures with dcStore
constellation__3_para.storeName='const';           %name of
constellation__3_para.gif='';                      %if folder given stroes gif
constellation__3_para.fs=320000000000;             %Sampling frequency 
constellation__3_para.sym_rate=10000000000;        %symbol rate in [1/s]
constellation__3_para.start_pos=1;                 %determines the sampling delay
                                                   %depending on parameter 'indexing' the delay is given relatively within the bit period
                                                   %or it is given as index into the input vector

% user_function__8

user_function__8_para.enable=1;                    %enable?
user_function__8_para.func_str='(x(common.buffer:end-common.buffer))'; %user defined function with argument x
user_function__8_para.interface=1;                 %Output interface








% dp_anal

dp_anal_para.mode=1;                               %MODE:
                                                   %0:off
                                                   %1:on
dp_anal_para.delay=0;                              %delay in samples
dp_anal_para.skip=0;                               %loops to skip
dp_anal_para.norm=1;                               %Normalize inputs (data and ref)
dp_anal_para.dc_block=0;                           %Normalize inputs
dp_anal_para.lab_mode=1;                           %delay in every loop .. delayed symbols are skipped
dp_anal_para.carrier=1;                            %Number of used subcarriers (1 = without multicarrier)
dp_anal_para.fft=1024;                             %FFT length in samples (for OFDM)
dp_anal_para.order=common.dimension;               %modulation order
dp_anal_para.dec=3;                                %nearest neighbor
                                                   %kmeans
                                                   %threashold of digidemod
dp_anal_para.MOD_Pr=0;                             %MOD_Pr
dp_anal_para.MOD=2;                                %MOD
dp_anal_para.Polarity=1;                           %polarity of mapping from dp flexmod for correct demodulation of PAM. important for multidimensional optimisations
dp_anal_para.cluster=0;                            %Estimate IRR
dp_anal_para.fom_evm=1;                            %Estimate EVM
dp_anal_para.fom_BER=1;                            %Estimate BER
dp_anal_para.fom_Q=0;                              %Estimate Q
dp_anal_para.fom_IRR=0;                            %Estimate IRR
dp_anal_para.plot_mode=1;                          %correlation
dp_anal_para.figure_ind=0;                         %figure index
dp_anal_para.nblock=[3 4 ];                        %screen sections
dp_anal_para.tile=[0 0 ];                          %figure position (y,x)
dp_anal_para.monitor=1;                            %monitor index
dp_anal_para.choicew=1;                            %choicew
dp_anal_para.hide_axis=1;                          %hide_axis
dp_anal_para.plot=0;                               %activate plots
dp_anal_para.plot_evm=0;                           %Ploting evm in 2D
dp_anal_para.plot_errors=0;                        %Ploting errors in demod and mse over time
dp_anal_para.cluster_plot=0;                       %Plot cluster analyses
dp_anal_para.cluster_texts=0;                      %Plot cluster analyses
dp_anal_para.plot_delay=0;                         %correlation
dp_anal_para.plot_flexanal=0;                      %correlation
dp_anal_para.berWin=1;                             %Stat window
dp_anal_para.store_plots=0;                        %store plots
dp_anal_para.sto=1;                                %Point of storage in simulation
dp_anal_para.win_update=1;                         %when the display figure is updated
dp_anal_para.sto_mode=1;                           %stores cumulated evm values or evms of every loop!
dp_anal_para.disp_2D=0;                            %disp_2D
dp_anal_para.disp_sym=0;                           %disp_sym
dp_anal_para.disp_bits=0;                          %disp_bits
dp_anal_para.disp_errs=0;                          %disp_errs
dp_anal_para.disp_g=0;                             %disp_g
dp_anal_para.disp_phi=0;                           %disp_phi
dp_anal_para.FEC_active=0;                         %FEC_active
dp_anal_para.FEC_algo=1;                           %FEC Encoding (for each bitstream separately - only implemented for realmod with random data)
dp_anal_para.minloops=1e+30;                       %minimal number of simulation loops before module can break the simulation
dp_anal_para.minerrors=1e+30;                      %minimal number of bit errors before module breaks the simulatio
dp_anal_para.break_ber=1e-05;                      %ber break uloop







% constellation__2

constellation__2_para.plot=1;                      %true: plot eye diagram; false: do not plot eye diagram
constellation__2_para.skip='0';                    %samples to skip
constellation__2_para.plot_color=1;                %color of plot
constellation__2_para.figure_ind=0;                %0=private figure, >0=index of common figure handle
constellation__2_para.monitor=1;                   %monitor index
constellation__2_para.plot_style=0;                %style of plot
constellation__2_para.plotmode=1;                  %plots:
                                                   %-only the current block
                                                   %-all blocks since the last plot
                                                   %-the whole signal calculated by now
constellation__2_para.plothist=1;                  %plot histogram of data
constellation__2_para.color_period=1;              %color bar periodicity
constellation__2_para.cofmagy=4;                   %n
                                                   %µ
                                                   %m
                                                   %1
                                                   %k
                                                   %M
                                                   %G
constellation__2_para.min_y_val=0;                 %min. y-value in plot
constellation__2_para.max_y_val=0;                 %max. y-value in plot
constellation__2_para.fig_plo=0;                   %0 to plot on the same figure
                                                   % 1 for separeted figures
constellation__2_para.markersize=3;                %size of const marker
constellation__2_para.grd=1;                       %determines, if a grid should be plotted
constellation__2_para.hide_axis=0;                 %hide labels from axis and rescale axis
constellation__2_para.plot_ref=1;                  %Plots reference constellation
constellation__2_para.nblock=[3 4 ];               %screen sections
constellation__2_para.tile=[0 0 ];                 %figure position (y,x)
constellation__2_para.choicew=1;                   %window proportion
constellation__2_para.interval=1;                  %show plots every x loop (x=0: at end of simulation)
constellation__2_para.stp=0;                       %stop after each plot
constellation__2_para.plotsize=[0 0 ];             %plot size(x,y)
constellation__2_para.histpoints=301;              %histogram dimension
constellation__2_para.hist_max_val=0;              %histogram max value (0==autoscaling)
constellation__2_para.norm=0;                      %Normalize constellation to 1
constellation__2_para.hist_colorbar=0;             %plot colorbarof histogram
constellation__2_para.store=1;                     %store figures with dcStore
constellation__2_para.storeName='const';           %name of
constellation__2_para.gif='';                      %if folder given stroes gif
constellation__2_para.fs=320000000000;             %Sampling frequency 
constellation__2_para.sym_rate=10000000000;        %symbol rate in [1/s]
constellation__2_para.start_pos=1;                 %determines the sampling delay
                                                   %depending on parameter 'indexing' the delay is given relatively within the bit period
                                                   %or it is given as index into the input vector

% dp_saser_monitor

dp_saser_monitor_para.active=1;                    %Activate/Deactivate Module
dp_saser_monitor_para.dc=0;                        %monitor dc value
dp_saser_monitor_para.pow_dBm=1;                   %monitor power value [dBm]
dp_saser_monitor_para.var=0;                       %monitor variance value
dp_saser_monitor_para.std=0;                       %monitor standard deviation value
dp_saser_monitor_para.osnr=0;                      %monitor OSNR like a boss!
dp_saser_monitor_para.scr=0;                       %monitor SCR value [dB]
dp_saser_monitor_para.accumulate=1;                %Accumulate monitored value over simulation blocks (true) or calculate independently per block (false)
dp_saser_monitor_para.delay=0;                     %Delay of the received signal (ignore meaningless zeros)
dp_saser_monitor_para.berWin=1;                    %Stat window
dp_saser_monitor_para.win_update=3;                %Point of storage in simulation
dp_saser_monitor_para.format_dc=2;                 %display format for dc value
dp_saser_monitor_para.format_pow_dBm=2;            %display format for pow_dBm value
dp_saser_monitor_para.format_var=2;                %display format for var value
dp_saser_monitor_para.format_std=2;                %display format for std value
dp_saser_monitor_para.format_osnr=2;               %display format for OSNR value
dp_saser_monitor_para.format_scr=2;                %display format for scr value
dp_saser_monitor_para.precision=2;                 %display precision (number of digits after decimal point)
dp_saser_monitor_para.store_val=1;                 %store value?
dp_saser_monitor_para.sto=2;                       %Point of storage in simulation
dp_saser_monitor_para.tot_pow=0;                   %Display total power of the dual-polarization signal
dp_saser_monitor_para.bandwidth=12500000000;       %Resolution bandwidth for OSNR monitoring [Hz]
dp_saser_monitor_para.points=501;                  %number of points for the power spectrum



% dc_sto_data

dc_sto_data_para.store=1;                          %Store data?
dc_sto_data_para.path='';                          %Storage path
dc_sto_data_para.name='DUL_Training_6dBm_512';     %Name of current simulation
dc_sto_data_para.subDirStr='';                     %Subdir string
dc_sto_data_para.pause=0;                          %Pause after user loop (1,2,3... - put index of uloops, which should be paused) and ask for input?
dc_sto_data_para.sendmail=0;                       %sendmail
dc_sto_data_para.recs='';                          %recipients
dc_sto_data_para.wvar='user_function__3_out;user_function__5_out;normalize__2_out;prms_out'; %Store variable from ws at loopmax. Use ';' to seperate variables and ',' behind a var to seperate it from a name
dc_sto_data_para.dispvar='common.taps_NL2,NL2;common.taps_NL3,NL3;common.volt_train,Volterra_Train'; %dis number from ws at loopmax. Use ';' to seperate variables and ',' behind a var to seperate it from a name
dc_sto_data_para.dzip=0;                           %Zip model and module data
dc_sto_data_para.exh=1;                            %Exhausting storage (every uloop iteration)=> if cancelation or errer is possible choose this
dc_sto_data_para.model=0;                          %Store model
dc_sto_data_para.figs=0;                           %Store figures
dc_sto_data_para.module=0;                         %Store modules
dc_sto_data_para.data=1;                           %Stores data
dc_sto_data_para.dirnames='';                      %Names for storage directorys
dc_sto_data_para.common=1;                         %store commons in loopmax
dc_sto_data_para.parameters=1;                     %store parameters
dc_sto_data_para.updateMode=1;                     %Defines when statistic window is updated!
dc_sto_data_para.long_meas=0;                      %long meas
