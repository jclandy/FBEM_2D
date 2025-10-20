function [P_t_full,P_t_ml,P_t_full_comp,P_t_ml_comp,em_bias] = Facet_Echo_Model_2D(op_mode,lambda,bandwidth,P_T,h,v,pitch,roll,prf,beam_weighting,G_0,D_0,gamma1,gamma2,N_b,t,sigma_0_snow_surf,sigma_0_snow_vol,kappa_e,tau_snow,c_s,h_s,sigma_0_ice_surf,sigma_surf_s,sigma_surf_si,l_surf,topo_type,dx)

%% 2D Facet-based Radar Altimeter Echo Model for Sea Ice

% Simulates the backscattered echo response of a pulse-limited or synthetic
% aperture radar altimeter from snow-covered sea ice, over a slope
% distribution representing an idealized statistical surface

%% Input (reference values)
% op_mode = operational mode: 1 = pulse-limited, 2 = SAR (PL-mode only feasible on high memory machines)
% lambda = radar wavelength, 0.0221 m
% bandwidth = antenna bandwidth, Hz
% P_T = transmitted peak power, 2.188e-5 watts
% h = satellite altitude, 720000 m
% v = satellite velocity, 7500 m/s
% pitch = antenna bench pitch counterclockwise, rads (up to ~0.005 rads)
% roll = antenna bench roll counterclockwise, rads (up to ~0.005 rads)
% prf = pulse-repitition frequency, Hz
% beam_weighting = weighting function on beam pattern (1 = rectangular, 2 =
% Hamming)
% G_0 = peak antenna gain, dB
% gamma1 = along-track antenna parameter, 0.0116 rads
% gamma2 = across-track antenna parameter, 0.0129 rads
% N_b = no. beams in synthetic aperture, 64 (1 in PL mode)
% t = time, s
% theta = angular sampling of scattering signatures, rads
% sigma_0_snow_surf = backscattering coefficient of snow surface, dB
% sigma_0_snow_vol = backscattering coefficient of snow volume, dB 
% kappa_e = extinction coefficient of snow volume, Np/m
% tau_snow = transmission coefficient at air-snow interface
% c_s = speed of light in snowpack, m/s
% h_s = snow depth, m
% sigma_0_ice_surf = backscattering coefficient of ice surface, dB
% sigma_0_lead_surf = backscattering coefficient of lead surface, dB

%% Output
% P_t_full = delay-Doppler map (DDM) of single look echoes, watts
% P_t_ml = multi-looked power waveform, watts
% P_t_full_comp = DDM for individual snow surface, snow volume, and ice
% surface components, watts
% P_t_ml_comp = multi-looked power waveforms for individual snow surface,
% snow volume, and ice surface components, watts

% Based on model equations introduced in Landy et al, TGARS, 2019
% Builing on theory of Wingham et al 2006, Giles et al 2007,
% Makynen et al 2009, Ulaby et al 2014, Halimi et al 2014, Wingham et al
% 2021

% Uses the following codes from external sources:

% (c) Jack Landy, UiT The Arctic University of Norway, 2024


%% Antenna parameters

c = 299792458; % speed of light, m/s
Re = 6371*10^3; % earth's radius, m

f_c = c/lambda; % radar frequency, Hz
k0 = (2*pi)/lambda; % wavenumber

delta_x = v/prf; % distance between coherent pulses in synthetic aperture, m

delta_x_dopp = (h*prf*c)/(2*N_b*v*f_c); % along-track doppler-beam limited footprint diameter, m
delta_x_pl = 2*sqrt(c*(h/((Re+h)/Re))*(1/bandwidth)); % across-track pulse-limited footprint diameter, m
beam_div = 1.1992*pi/180; % across-track CS2 beam divergence
delta_x_bl = 2*h*tan(beam_div/2); % across-track beam-limited footprint diameter, m
% A_pl = pi*(delta_x_pl/2)^2; % area of each range ring (after waveform peak), m

epsilon_b = lambda/(2*N_b*v*(1/prf)); % angular resolution of beams from full look crescent (beam separation angle) 

% Antenna look geometry
m = -(N_b-1)/2:(N_b-1)/2;

%% Basic geometrical functions
% Antenna location (x_0,y_0) with respect to footprint centre

% Range
R = @(x,y,x_0,y_0) sqrt(h^2 + ((x-x_0).^2 + (y-y_0).^2)*(1 + h/Re));

% Antenna illumination angles
THETA_G = @(x,y,x_0,y_0) pi/2 + atan2(-R(x,y,x_0,y_0), sqrt((x-x_0+h*tan(pitch)).^2 + (y-y_0+h*tan(roll)).^2));
PHI_G = @(x,y,x_0,y_0) atan2((y-y_0+h*tan(roll)),(x-x_0+h*tan(pitch)));

% Look angle of synthetic Doppler beam
theta_l = @(x,x_0) atan(-(x-x_0+h*tan(pitch))./-h); 

%% Integration limits

x_lim = delta_x_pl*4;
y_lim = delta_x_pl*4;

%% Gain function
    
% Antenna gain pattern
G = @(x,y,x_0,y_0) G_0*exp(-THETA_G(x,y,x_0,y_0).^2.*(cos(PHI_G(x,y,x_0,y_0)).^2/gamma1^2 + sin(PHI_G(x,y,x_0,y_0)).^2/gamma2^2));

%% Surface height distribution functions

if topo_type==1
    % Gaussian height distribution function
    mu = @(sigma) 0;
    var = @(sigma) sigma^2;
    % rough0 = @(x,mu,var) (1/(sqrt(var)*sqrt(2*pi))).*exp(-0.5*((-(x + 1 - mu)/sqrt(var)).^2)); % gaussian
    rough0 = @(x,mu,var) normpdf(-x-1,mu,sqrt(var));
    
elseif topo_type==2
    % Lognormal height distribution function
    mu = @(sigma) log(1/sqrt(1 + sigma^2));
    var = @(sigma) log(1 + sigma^2/1^2);
    % rough0 = @(x,mu,var) IF(x<0,1./(-x*sqrt(var)*sqrt(2*pi)).*exp(-log(-x - mu).^2/(2*var)),0); % lognormal
    rough0 = @(x,mu,var) lognpdf(-x,mu,sqrt(var));
    
end

%% Sigma0 functions

% Slope distribution std dev from random statistical surface
[x,y,z] = synthetic_topo_shell(1,topo_type,sigma_surf_si,l_surf,0,dx);

[FX,FY] = gradient(z,x(1,:),y(:,1));
S = atan(sqrt(FX.^2 + FY.^2));
sigma_dxdyh_si = rms(S(:))/sqrt(2); % rms slope that matches the rician pdf scale parameter equilavent to histogram(S(:))
sigma_dxdyh_s = rms(S(:)*sigma_surf_s/sigma_surf_si)/sqrt(2); % scale slope distribution for snow topography

% sigma_dxdyh_si = sigma_surf_si/l_surf; % effective rms slope for exponentially-correlated surface (not applicable for non-Gaussian surface?)
% sigma_dxdyh_s = sigma_surf_s/l_surf;

% Slope PDF for angle theta from antenna location, follows the Rice distribution
S_PDF = @(slope,theta,sigma) pdf('Rician',slope,theta,sigma);

% Surface backscatter coefficients as a function of local slope wrt nadir
MU_S_si = @(slope) 10.^(ppval(sigma_0_ice_surf,slope)/10).*ppval(tau_snow,slope).^2*exp(-kappa_e*h_s/2);
MU_S_s = @(slope) 10.^(ppval(sigma_0_snow_surf,slope)/10)*(h_s~=0);

% Volume backscatter coefficient as a function of off-nadir angle
VU_THETA_t_s = @(theta,t) 10.^(ppval(sigma_0_snow_vol,theta)/10)*kappa_e.*exp(-c_s*kappa_e*(t + 2*h_s/c_s));

%% Apply electromagnetic bias
% In the case of surfaces with non-Gaussian height distribution, the
% surface height and slope can be correlated and this is accounted for
% here by weighting the height distribution for the backscatter
% contribution as a function of slope

% If sigma_dxdyh is explicitly calculated from a random surface slope distribution
covmat_si = cov(z(:),S(:)); % sea ice
cov_h_dhdx_si = covmat_si(1,2);
covmat_s = cov(z(:)*sigma_surf_s/sigma_surf_si,S(:)*sigma_surf_s/sigma_surf_si); % scale slope distribution for snow topography
cov_h_dhdx_s = covmat_s(1,2);

% Mean of slope distribution, when noncentrality = 0
% https://stats.libretexts.org/Bookshelves/Probability_Theory/Probability_Mathematical_Statistics_and_Stochastic_Processes_(Siegrist)/04%3A_Expected_Value/4.05%3A_Covariance_and_Correlation
% eq 4.5.34
S_mean = @(sigma_dxdyh) sigma_dxdyh*sqrt(pi/2);

% Mean slope of distribution as a function of surface height z,
% with z reversed in terms of range distance
S_z = @(z,sigma_surf,sigma_dxdyh,cov_h_dhdx) S_mean(sigma_dxdyh) + cov_h_dhdx/sigma_surf^2*-z;

% Rician scale parameter as a function of surface height z
sigma_dxdyh_z = @(z,sigma_surf,sigma_dxdyh,cov_h_dhdx) S_z(z,sigma_surf,sigma_dxdyh,cov_h_dhdx)/sqrt(pi/2);

% Surface backscattering distribution as a function of surface height z
MU_z_si = @(z) integral(@(slope) MU_S_si(slope).*S_PDF(slope,0,sigma_dxdyh_z(z,sigma_surf_si,sigma_dxdyh_si,cov_h_dhdx_si)), 0, sigma_dxdyh_si*10, 'arrayvalued', true);
MU_z_s = @(z) integral(@(slope) MU_S_s(slope).*S_PDF(slope,0,sigma_dxdyh_z(z,sigma_surf_s,sigma_dxdyh_s,cov_h_dhdx_s)), 0, sigma_dxdyh_s*10, 'arrayvalued', true);

% Discretized range bins
tsamp = round(1/bandwidth/(t(2)-t(1)));
L = 1049;
Z0 = 512;
Z = ((1/tsamp)*((1:L)-1-Z0)*1/bandwidth)*c/2;

% Ice and snow surface backscattering distributions in discretized range bins
warning('off','all')
MU_Z_si = MU_z_si(Z);
MU_Z_si(isnan(MU_Z_si)) = 0;
MU_Z_s = MU_z_s(Z);
MU_Z_s(isnan(MU_Z_s)) = 0;
warning('on','all')

%% Estimate radar electromagnetic bias for rough sea ice surface

% Generally EM bias for sea ice will be higher than for ocean with same
% height stddev, because the sea ice has a more skewed/lognormal height
% distribution than the skewed-gaussian ocean distribution
em_bias = sum(Z.*rough0(Z-1,mu(sigma_surf_si),var(sigma_surf_si)).*MU_Z_si)/sum(rough0(Z-1,mu(sigma_surf_si),var(sigma_surf_si)).*MU_Z_si);

%% Transmitted power envelope
    
% Time offset including slant-range time correction
T = @(x,y,t,x_0,y_0,tc) t + 2*h/c + tc - 2*(R(x,y,x_0,y_0)/c);

% Power envelope
P_t = @(x,y,t,x_0,y_0,tc) sinc(bandwidth*T(x,y,t,x_0,y_0,tc)).^2;
    
%% Beam weighting function
if op_mode==1
    beam_weighting = 1;
else
end

%% Radar Simulator Loop
% Formulated for parallel processing

P_t_full_si_surf = zeros(length(m),length(t));
P_t_full_s_surf = zeros(length(m),length(t));
P_t_full_s_vol = zeros(length(m),length(t));
parfevalOnAll(@warning,0,'off','all');
parfor ii = 1:length(m)
    
    %% Antenna location

    X_0 = h*m(ii)*epsilon_b + h*tan(pitch); % along track
    Y_0 = h*tan(roll); % across track
    
    %% Synthetic beam gain function
    
    P_m = [];
    if beam_weighting == 1
        % Rectangular
        P_m = @(x) D_0*sin(N_b*(k0*delta_x*sin(theta_l(x,X_0) + m(ii)*epsilon_b))).^2./(N_b*sin(k0*delta_x*sin(theta_l(x,X_0) + m(ii)*epsilon_b))).^2; 
    elseif beam_weighting == 2
        % Apply hamming window to azimuthal response function
        n = 0:(N_b-1);
        P_m = @(x) abs(sum((0.54 - 0.46*cos((2*pi*n)/(N_b - 1))) .* exp(2*1j*k0*v/prf*(theta_l(x(:),X_0) + m(ii)*epsilon_b).*(n - (N_b-1)/2)),2).^2);
    end

    %% Slant-range time correction

    TC = IF(op_mode==1,0,2*(sqrt(X_0.^2*(1 + h/Re)+h^2)-h)/c);
    
    %% Backscatter response for surfaces
        
    % Calculate mean sigma0 as a function of theta
    the = (pi/2 - atan2(h, sqrt(((0:1:delta_x_bl/2)-X_0).^2 + ((0:1:delta_x_bl/2)-Y_0).^2)))';
    
    % Sea ice surface
    MU_THETA_si = @(theta) integral(@(slope) MU_S_si(slope).*S_PDF(slope,theta,sigma_dxdyh_si), 0, max(the), 'arrayvalued', true, 'abstol',1e-5,'reltol',1e-5);
    mu_the = IF(sigma_surf_si*ones(size(the))==0,10.^(ppval(sigma_0_ice_surf,the)/10).*ppval(tau_snow,the).^2*exp(-kappa_e*h_s/2),MU_THETA_si(the));
    
    % Gaussian function approximating sea ice surface sigma0 distribution with theta
    % (error typically <0.01% and enables x10+ speedup)
    [mu_the_fit,gof] = fit(the(~isnan(mu_the)),mu_the(~isnan(mu_the)),'a*exp(-((x-b)/c)^2) + d','StartPoint',[nanmax(mu_the) 0 the(find(mu_the<(nanmax(mu_the)-range(mu_the)*0.5),1)) 0],'Upper',[Inf 0 Inf Inf], 'robust','LAR','display','off');
    MU_THETA_FIT_si = @(x,y) mu_the_fit.a*exp(-((THETA_G(x,y,X_0,Y_0) - mu_the_fit.b)/mu_the_fit.c).^2) + mu_the_fit.d;
    
    % Snow surface
    MU_THETA_s = @(theta) integral(@(slope) MU_S_s(slope).*S_PDF(slope,theta,sigma_dxdyh_s), 0, max(the), 'arrayvalued', true, 'abstol',1e-5,'reltol',1e-5);
    mu_the = IF(sigma_surf_s*ones(size(the))==0,10.^(ppval(sigma_0_snow_surf,the)/10)*(h_s~=0),MU_THETA_s(the));
    
    % Gaussian function approximating snow surface sigma0 distribution with theta
    if h_s > 0
        [mu_the_fit,gof] = fit(the(~isnan(mu_the)),mu_the(~isnan(mu_the)),'a*exp(-((x-b)/c)^2) + d','StartPoint',[nanmax(mu_the) 0 the(find(mu_the<(nanmax(mu_the)-range(mu_the)*0.5),1)) 0],'Upper',[Inf 0 Inf Inf], 'robust','LAR','display','off');
        MU_THETA_FIT_s = @(x,y) mu_the_fit.a*exp(-((THETA_G(x,y,X_0,Y_0) - mu_the_fit.b)/mu_the_fit.c).^2) + mu_the_fit.d;
    else
        MU_THETA_FIT_s = @(x,y) 0;
    end

    %% Backscatter response for snow volume
    
    VU_THETA_T_s = @(theta,t) IF(t>=-2*h_s/c_s & t<0,VU_THETA_t_s(theta,t),0);

    %% Integrate Flat-Surface Impulse Response, Transmitted Echo and Backscatter

    % Volume echo integration limits
    if X_0 < 0 && op_mode == 2
        x_low = @(t) min([-delta_x_dopp, -delta_x_dopp + X_0 + real((Re*c*h*(t + TC + h_s/c_s)/(Re + h))^(1/2))]);
        x_high = @(t) min([delta_x_dopp, X_0 + real((Re*c*h*(t + TC + 2*h_s/c_s)/(Re + h))^(1/2))]);
    elseif X_0 >= 0 && op_mode == 2
        x_low = @(t) max([-delta_x_dopp, X_0 - real((Re*c*h*(t + TC + 2*h_s/c_s)/(Re + h))^(1/2))]);
        x_high = @(t) max([delta_x_dopp, delta_x_dopp + X_0 - real((Re*c*h*(t + TC + h_s/c_s)/(Re + h))^(1/2))]);
    else
        x_low = @(t) X_0 - real((Re*c*h*(t + TC + 2*h_s/c_s)/(Re + h))^(1/2));
        x_high = @(t) X_0 + real((Re*c*h*(t + TC + 2*h_s/c_s)/(Re + h))^(1/2));
    end

    % Radar equation integral per range bin
    % warning('off','all')
    FSIR_S0_P_t_si_surf = zeros(size(t));
    FSIR_S0_P_t_s_surf = zeros(size(t));
    FSIR_S0_P_t_s_vol = zeros(size(t));
    for jj = 1:length(t)
        
        % sea ice surface
        % si_surf_tracer(jj) = quad2d(@(x,y) P_t(x,y,t(jj)).*G(x,y,X_0,Y_0).^2.*P_m(x).*MU_THETA_si(THETA_G(x,y,X_0,Y_0)), -x_lim + X_0, x_lim + X_0, -y_lim + Y_0, y_lim + Y_0, 'abstol',1e-8,'reltol',1e-8,'MaxFunEvals',1500)*1/h^4; % exact
        FSIR_S0_P_t_si_surf(jj) = quad2d(@(x,y) P_t(x,y,t(jj),X_0,Y_0,TC).*G(x,y,X_0,Y_0).^2.*P_m(x).*MU_THETA_FIT_si(x,y), -x_lim + X_0, x_lim + X_0, -y_lim + Y_0, y_lim + Y_0, 'abstol',1e-8,'reltol',1e-8,'MaxFunEvals',1500)*1/h^4; % approximated
        
        if h_s > 0
            % snow surface
            % s_surf_tracer(jj) = quad2d(@(x,y) P_t(x,y,t(jj)) + 2*h_s/c_s.*G(x,y,X_0,Y_0).^2.*P_m(x).*MU_THETA_s(THETA_G(x,y,X_0,Y_0)), -x_lim + X_0, x_lim + X_0, -y_lim + Y_0, y_lim + Y_0, 'abstol',1e-8,'reltol',1e-8,'MaxFunEvals',1500)*1/h^4; % exact
            FSIR_S0_P_t_s_surf(jj) = quad2d(@(x,y) P_t(x,y,t(jj) + 2*h_s/c_s,X_0,Y_0,TC).*G(x,y,X_0,Y_0).^2.*P_m(x).*MU_THETA_FIT_s(x,y), -x_lim + X_0, x_lim + X_0, -y_lim + Y_0, y_lim + Y_0, 'abstol',1e-8,'reltol',1e-8,'MaxFunEvals',1500)*1/h^4; % approximated
            
            % snow volume (symmetric along y-axis about origin [X_0,Y_0] except when roll ~= 0)
            y_low = @(x) real(sqrt(real(sqrt((h/((Re+h)/Re))*c*(t(jj) + TC))).^2 - (x - X_0).^2)) + Y_0;
            y_high = @(x) real(sqrt(real(sqrt((h/((Re+h)/Re))*c*(t(jj) + TC + 2*h_s/c_s))).^2 - (x - X_0).^2)) + Y_0;
        
            FSIR_S0_P_t_s_vol(jj) = 2*quad2d(@(x,y) P_t(x,y,t(jj) + 2*h_s/c_s,X_0,Y_0,TC).*G(x,y,X_0,Y_0).^2.*P_m(x).*VU_THETA_T_s(THETA_G(x,y,X_0,Y_0),T(x,y,t(jj),X_0,Y_0,TC)), x_low(t(jj)), x_high(t(jj)), y_low, y_high,'abstol',1e-8)*1/h^4; % exact
        else
            FSIR_S0_P_t_s_surf(jj) = 0;
            FSIR_S0_P_t_s_vol(jj) = 0;
        end

    end
    FSIR_S0_P_t_si_surf(isnan(FSIR_S0_P_t_si_surf)) = 0;
    FSIR_S0_P_t_s_surf(isnan(FSIR_S0_P_t_s_surf)) = 0;
    FSIR_S0_P_t_s_vol(isnan(FSIR_S0_P_t_s_vol)) = 0;
    % FSIR_S0_P_t = FSIR_S0_P_t_si_surf + FSIR_S0_P_t_s_surf + FSIR_S0_P_t_s_vol;
    % warning('on','all')
    
    %% Convolution of FSIR and P_t over the sigma0 weighted surface height distribution
    
    % Weight height distribution by covariance with sigma0
    if sigma_surf_si==0
        
        % rough1 = kroneckerDelta(Z,sym(0));
        P_r_extend_s_surf = FSIR_S0_P_t_s_surf;
        P_r_extend_s_vol = FSIR_S0_P_t_s_vol;
        P_r_extend_si_surf = FSIR_S0_P_t_si_surf;
        
    else
        
        %%% Snow surface
        rough1 = rough0(Z-1,mu(sigma_surf_s),var(sigma_surf_s)).*MU_Z_s; rough1 = rough1/trapz(Z,rough1);
    
        % Prepare vectors for FFT
        rough1_spl = [rough1(Z0:L) rough1(1:Z0-1)];
    
        zoos = zeros(1,length(length(t)+1:L));
        FSIR_S0_P_t_s_surf_extend = [FSIR_S0_P_t_s_surf zoos];
        
        % Convolution through FFT
        f_rough1_spl = sqrt(L)*10^(-10)*dft(rough1_spl,0,1);
        f_FSIR_S0_P_t_s_surf_extend = sqrt(L)*10^(-10)*dft(FSIR_S0_P_t_s_surf_extend,0,1);
        
        f_rough1_FSIR_S0_P_t_s_surf = f_rough1_spl.*f_FSIR_S0_P_t_s_surf_extend;
        
        P_r_extend_s_surf = real((1/sqrt(L))*(1/10^(-10))*dft(f_rough1_FSIR_S0_P_t_s_surf,0,-1));
        
        %%% Snow volume
        rough1 = rough0(Z-1,mu(sigma_surf_s),var(sigma_surf_s)); % no em bias
    
        % Prepare vectors for FFT
        rough1_spl = [rough1(Z0:L) rough1(1:Z0-1)];
    
        zoos = zeros(1,length(length(t)+1:L));
        FSIR_S0_P_t_s_vol_extend = [FSIR_S0_P_t_s_vol zoos];
    
        % Convolution through FFT
        f_rough1_spl = sqrt(L)*10^(-10)*dft(rough1_spl,0,1);
        f_FSIR_S0_P_t_s_vol_extend = sqrt(L)*10^(-10)*dft(FSIR_S0_P_t_s_vol_extend,0,1);
        
        f_rough1_FSIR_S0_P_t_s_vol = f_rough1_spl.*f_FSIR_S0_P_t_s_vol_extend;

        P_r_extend_s_vol = real((1/sqrt(L))*(1/10^(-10))*dft(f_rough1_FSIR_S0_P_t_s_vol,0,-1));
        
        %%% Sea ice
        rough1 = rough0(Z-1,mu(sigma_surf_si),var(sigma_surf_si)).*MU_Z_si; rough1 = rough1/trapz(Z,rough1);
    
        % Prepare vectors for FFT
        rough1_spl = [rough1(Z0:L) rough1(1:Z0-1)];
    
        zoos = zeros(1,length(length(t)+1:L));
        FSIR_S0_P_t_si_surf_extend = [FSIR_S0_P_t_si_surf zoos];

        % Convolution through FFT
        f_rough1_spl = sqrt(L)*10^(-10)*dft(rough1_spl,0,1);
        f_FSIR_S0_P_t_si_surf_extend = sqrt(L)*10^(-10)*dft(FSIR_S0_P_t_si_surf_extend,0,1);
        
        f_rough1_FSIR_S0_P_t_si_surf = f_rough1_spl.*f_FSIR_S0_P_t_si_surf_extend;

        P_r_extend_si_surf = real((1/sqrt(L))*(1/10^(-10))*dft(f_rough1_FSIR_S0_P_t_si_surf,0,-1));
        
    end

    %% Normalize with radar equation (absolute numbers not trustworthy)
    
    P_r_si_surf = ((lambda^2*P_T)/(4*pi)^3)*(0.5*c*h)*P_r_extend_si_surf((1:length(t)));
    P_r_si_surf(P_r_si_surf<0) = 0;
    P_t_full_si_surf(ii,:) = P_r_si_surf;
    
    P_r_s_surf = ((lambda^2*P_T)/(4*pi)^3)*(0.5*c*h)*P_r_extend_s_surf((1:length(t)));
    P_r_s_surf(P_r_s_surf<0) = 0;
    P_t_full_s_surf(ii,:) = P_r_s_surf;
    
    P_r_s_vol = ((lambda^2*P_T)/(4*pi)^3)*(0.5*c*h)*P_r_extend_s_vol((1:length(t)));
    P_r_s_vol(P_r_s_vol<0) = 0;
    P_t_full_s_vol(ii,:) = P_r_s_vol;
       
end
parfevalOnAll(@warning,0,'on','all');

% Full DDM
P_t_full_comp = permute(cat(3,P_t_full_s_surf,P_t_full_s_vol,P_t_full_si_surf),[2 1 3]);
P_t_full = (P_t_full_si_surf + P_t_full_s_surf + P_t_full_s_vol)';

% Multi-looked echoes
P_t_ml_comp = permute(nansum(P_t_full_comp,2),[1 3 2]);
P_t_ml = nansum(P_t_full,2);


end

