function f = test_enthalpy_derivative
f = figure;
p = logspace(log10(0.000611657),log10(22.064),10000)'; % [MPa] pressure (p) range
hL = IAPWS_IF97('hL_p',p);
pL = p(~isnan(hL));hL = hL(~isnan(hL));
dhLdp_avg = diff(hL)./diff(pL);
pL_avg = (pL(1:end-1)+pL(2:end))/2;
dhLdp = IAPWS_IF97('dhLdp_p',pL);
subplot(1,2,1),plot(pL,dhLdp,pL_avg,dhLdp_avg)

hV = IAPWS_IF97('hV_p',p);
pV = p(~isnan(hV));hV = hV(~isnan(hV));
dhVdp_avg = diff(hV)./diff(pV);
pV_avg = (pV(1:end-1)+pV(2:end))/2;
dhVdp = IAPWS_IF97('dhVdp_p',pV);
subplot(1,2,2),plot(pV,dhVdp,pV_avg,dhVdp_avg)
end