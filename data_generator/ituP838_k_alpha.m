function [k, alpha] = ituP838_k_alpha(f_GHz, theta_deg, polarization)
    % ITUP838_K_ALPHA  Specific attenuation model coefficients (ITU-R P.838-3).
    x = log10(f_GHz);
    
    akH = [-5.33980, -0.35351, -0.23789, -0.94158];
    bkH = [-0.10008,  1.26970,  0.86036,  0.64552];
    ckH = [ 1.13098,  0.45400,  0.15354,  0.16817];
    
    akV = [-3.80595, -3.44965, -0.39902,  0.50167];
    bkV = [ 0.56934, -0.22911,  0.73042,  1.07319];
    ckV = [ 0.81061,  0.51059,  0.11899,  0.27195];
    
    logkH = sum(akH .* exp(-((x-bkH)./ckH).^2)) - 0.18961*x + 0.71147;
    logkV = sum(akV .* exp(-((x-bkV)./ckV).^2)) - 0.16398*x + 0.63297;
    kH = 10^logkH;
    kV = 10^logkV;
    
    aaH = [-0.14318,  0.29591,  0.32177, -5.37610, 16.1721];
    baH = [ 1.82442,  0.77564,  0.63773, -0.96230, -3.29980];
    caH = [-0.55187,  0.19822,  0.13164,  1.47828,  3.43990];
    
    aaV = [-0.07771,  0.56727, -0.20238, -48.2991, 48.5833];
    baV = [ 2.33840,  0.95545,  1.14520,  0.791669, 0.791459];
    caV = [-0.76284,  0.54039,  0.26809,  0.116226, 0.116479];
    
    alphaH = sum(aaH .* exp(-((x-baH)./caH).^2)) + 0.67849*x - 1.95537;
    alphaV = sum(aaV .* exp(-((x-baV)./caV).^2)) - 0.053739*x + 0.83433;
    
    pol = upper(string(polarization));
    if     pol == "H", tau_deg = 0;
    elseif pol == "V", tau_deg = 90;
    else,              tau_deg = 45;
    end
    
    k     = (kH + kV + (kH-kV)*cosd(theta_deg)^2*cosd(2*tau_deg)) / 2;
    alpha = (kH*alphaH + kV*alphaV + ...
            (kH*alphaH - kV*alphaV)*cosd(theta_deg)^2*cosd(2*tau_deg)) / (2*k);
end