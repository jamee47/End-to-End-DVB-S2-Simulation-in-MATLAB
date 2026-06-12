function Ap_dB = calcRainAttenuationP618(chanParams, p)
    % CALCRAINATTENUATIONP618  ITU-R P.618-13 rain attenuation at exceedance p (%).
    f_GHz     = chanParams.Frequency_GHz;
    theta_deg = chanParams.ElevationAngle_deg;
    lat_deg   = chanParams.Latitude_deg;
    hs_km     = chanParams.Altitude_m / 1000;
    R001      = chanParams.R001_mmph;
    h0_km     = chanParams.IsothermHeight_km;
    
    hR_km = h0_km + 0.36;
    Re_km = 8500;
    
    if R001 <= 0 || hR_km <= hs_km
        Ap_dB = 0; return;
    end
    
    if theta_deg >= 5
        Ls = (hR_km - hs_km) / sind(theta_deg);
    else
        Ls = 2*(hR_km - hs_km) / ...
             (sqrt(sind(theta_deg)^2 + 2*(hR_km-hs_km)/Re_km) + sind(theta_deg));
    end
    LG = Ls * cosd(theta_deg);
    
    [k, alpha] = ituP838_k_alpha(f_GHz, theta_deg, chanParams.Polarization);
    gamma_R    = k * (R001^alpha);
    
    r001 = 1 / (1 + 0.78*sqrt((LG*gamma_R)/f_GHz) - 0.38*(1 - exp(-2*LG)));
    
    zeta = atand((hR_km - hs_km) / (LG*r001));
    if zeta > theta_deg
        LR = (LG*r001) / cosd(theta_deg);
    else
        LR = (hR_km - hs_km) / sind(theta_deg);
    end
    
    chi  = max(0, 36 - abs(lat_deg));
    v001 = 1 / (1 + sqrt(sind(theta_deg)) * ...
           (31*(1 - exp(-theta_deg/(1+chi))) * sqrt(LR*gamma_R)/(f_GHz^2) - 0.45));
    
    LE   = LR * v001;
    A001 = gamma_R * LE;
    
    if abs(p - 0.01) < 1e-12
        Ap_dB = A001; return;
    end
    
    if p >= 1 || abs(lat_deg) >= 36
        beta = 0;
    elseif theta_deg >= 25
        beta = -0.005*(abs(lat_deg) - 36);
    else
        beta = -0.005*(abs(lat_deg) - 36) + 1.8 - 4.25*sind(theta_deg);
    end
    
    if A001 <= 0
        Ap_dB = 0; return;
    end
    
    exponent = -(0.655 + 0.033*log(p) - 0.045*log(A001) ...
                 - beta*(1-p)*sind(theta_deg));
    Ap_dB    = A001 * (p/0.01)^exponent;
end