%% Shell Script for the 2D Facet-based Radar Altimeter Echo Model for Sea Ice

% Shell script controlling application of the facet-based model for
% pulse-limited or SAR altimeter echo over snow-covered sea ice

% (c) J.C. Landy, UiT The Arctic University of Norway, 2024

% Included Codes:
% Facet_Echo_Model_2D.m
% RelDielConst_Brine.m
% rsgene2D_anisotrop.m
% synthetic_topo_shell.m
% snow_backscatter.m
% ice_backscatter.m

% Uses the following codes from external sources:
% I2EM_Backscatter_model.m (Fawwaz Ulaby)
% MieExtinc_DrySnow.m (Fawwaz Ulaby)
% RelDielConst_PureIce.m (Fawwaz Ulaby)
% RelDielConst_DrySnow.m (Fawwaz Ulaby)
% RelDielConst_SalineWater.m (Fawwaz Ulaby)
% TVBmodel_HeterogeneousMix.m (Fawwaz Ulaby)

% clear; clc

%% Model Variables (MODIFIABLE)

% Up to 3 variables can be vectors


% Geophysical parameters
sigma_si = [0.0001:0.0001:0.0014 0.0015:0.0005:0.006]; % sea ice rms height (default = 0.002 m)
l_si = 0.02; % sea ice correlation length (default = 0.02 m)
sigma_surf_si = [0:0.001:0.009 0.01:0.005:0.09 0.1:0.025:1]; % sea ice large-scale rms roughness height (default = 0.1 m)
T_si = -15; % sea ice bulk temperature (default = -15 C)
S_si = 6; % sea ice bulk salinity (default = 6 ppt)

sigma_s = 0.001; % snow rms height (default = 0.001 m)
l_s = 0.04; % snow correlation length (default = 0.04 m)
sigma_surf_s = sigma_surf_si*1.5; % snow large-scale rms roughness height (default = x1.5 sea ice rms height)
T_s = -20; % snow bulk temperature (default = -20 C)
rho_s = 350; % snow bulk density (default = 350 kg/m^3)
r_s = 0.001; % snow grain size (normal range from 0.0001 to 0.004 m, default 1 mm)
h_s = 0; % snow depth, m

% Antenna parameters
c = 299792458; % speed of light, m/s
lambda = c/13.5e9;  % ku-band
GP = whos('*'); % all parameters controlling scattering signatures

op_mode = 2; % operational mode: 1 = pulse-limited, 2 = SAR
beam_weighting = 1; % weighting on the beam-wise azimuth FFT: 1 = rectangular, 2 = Hamming (default = Hamming)
P_T = 2.188e-5; % transmitted peak power (default = 2.188e-5 watts)

pitch = 0; % antenna bench pitch counterclockwise (up to ~0.01 rads)
roll = 0; % antenna bench roll counterclockwise (up to ~0.005 rads)

h = 698900; % satellite altitude (default = 720000 m)
v = 7508; % satellite velocity (default = 7500 m/s)
N_b = 64; % no. beams in synthetic aperture / pulses in a burst (default = 64)
if op_mode==1 % single beam for PL echo
    N_b = 1;
else
end
            
prf = 15119; % pulse-repetition frequency (default = 15119 Hz, CRISTAL)
bandwidth = 500*10^6; % antenna bandwidth
G_0 = 42.3; % peak antenna gain, dB, ku-band
D_0 = 20*log10(N_b); % synthetic beam gain

beam_div_1 = 0.94;
beam_div_2 = 1.18;
gamma1 = (beam_div_1*pi/180)/(2*sqrt(log(2))); % CRISTAL along-track term, ku-band
gamma2 = (beam_div_2*pi/180)/(2*sqrt(log(2))); % CRISTAL across-track term, ku-band

% Number of range bins
N_tb = 60; % (default = 70)

% Range bin at mean scattering surface, i.e. time = 0
t_0 = 15; % (default = 15)

% Time oversampling factor on original bandwidth
t_sub = 10;

% Parameters of synthetic topography
topo_type = 2; % type of surface: 1 = Gaussian, 2 = lognormal
l_surf = 15; % large-scale correlation length (default = 5 m)
dx = 5; % resolution of grid for deriving the distribution of surface slopes, m (default = 1 m)

% Antenna Geometry
epsilon_b = lambda/(2*N_b*v*(1/prf)); % angular resolution of beams from full look crescent (beam separation angle) 

%% Loop Echo Model

% Use parallel processing
% parpool

% Time domain
t = (1/bandwidth)*((1:(1/t_sub):N_tb) - t_0);

save('FEM2D_Simulations_CRISTAL_Ku_Lognormal_Rectangular');

% Loop model over vector variables
P_t_full_range = cell(length(sigma_si),length(sigma_surf_si));
P_t_ml_range = cell(length(sigma_si),length(sigma_surf_si));
P_t_full_comp_range = cell(length(sigma_si),length(sigma_surf_si));
P_t_ml_comp_range = cell(length(sigma_si),length(sigma_surf_si));
EM_bias_range = NaN(length(sigma_si),length(sigma_surf_si));
for i = 1:length(sigma_si)
    
    % Effective width of angular extent of coherent component (TUNING PARAMETER FOR LEADS)
    beta_c = epsilon_b;
    % beta_c = 0.0221/(2*64*7500*(1/18182)); % same as Cryosat-2 angular beam resolution, i.e.: 0.0221/(2*64*7500*(1/18182))
    % c = 299792458; % speed of light, m/s
    % beta_c = sqrt((c*1/bandwidth)/h); % from Carsey et al 1992 (pp 117)
    % beta_c = 0.5*lambda/1.15; % from Fung & Eom 1983 eq 10, far-field condition based on antenna diameter

    % Surface & volume backscattering properties
    [theta,sigma_0_snow_surf,sigma_0_snow_vol,kappa_e,tau_snow,c_s,epsr_ds] = snow_backscatter(lambda,sigma_s,l_s,T_s,rho_s,r_s,h_s,beta_c);
    [~,sigma_0_ice_surf,~] = ice_backscatter(lambda,sigma_si(i),l_si,T_si,S_si,h_s,beta_c,epsr_ds);
    
    for j = 1:length(sigma_surf_si)
        
        fprintf(['Simulation ' num2str((length(sigma_surf_si)*(i-1) + j)) '/' num2str(length(sigma_si)*length(sigma_surf_si)) '\n']);
        
        % if sigma_si(i) >= (sigma_surf_si(j)*5e-3) %% Ignore unrealistic combos

            [P_t_full,P_t_ml,P_t_full_comp,P_t_ml_comp,em_bias] = Facet_Echo_Model_2D(op_mode,lambda,bandwidth,P_T,h,v,pitch,roll,prf,beam_weighting,G_0,D_0,gamma1,gamma2,N_b,t,sigma_0_snow_surf,sigma_0_snow_vol,kappa_e,tau_snow,c_s,h_s,sigma_0_ice_surf,sigma_surf_s(j),sigma_surf_si(j),l_surf,topo_type,dx);
            
            P_t_full_range{i,j} = P_t_full;
            P_t_ml_range{i,j} = P_t_ml;
            P_t_full_comp_range{i,j} = P_t_full_comp;
            P_t_ml_comp_range{i,j} = P_t_ml_comp;
            EM_bias_range(i,j) = em_bias;

        % else
        %     P_t_full_range{i,j} = NaN(N_b,length(t));
        %     P_t_ml_range{i,j} = NaN(1,length(t));
        %     P_t_full_comp_range{i,j} = NaN(length(t),N_b,3);
        %     P_t_ml_comp_range{i,j} = NaN(length(t),3);
        % end
        
        % Add echo to results
        save('FEM2D_Simulations_CRISTAL_Ku_Lognormal_Rectangular','P_t_full_range','P_t_ml_range','P_t_full_comp_range','P_t_ml_comp_range','EM_bias_range','beta_c','-append');

    end
   
end

