function f = test_density_derivative
f = figure;
p = logspace(log10(0.000611657),log10(22.064),10000)'; % [MPa] pressure (p) range
vL = IAPWS_IF97('vL_p',p);
pL = p(~isnan(vL));vL = vL(~isnan(vL));
dvLdp_avg = diff(vL)./diff(pL);
pL_avg = (pL(1:end-1)+pL(2:end))/2;
dvLdp = IAPWS_IF97('dvLdp_p',pL);
subplot(1,2,1),plot(pL,dvLdp,pL_avg,dvLdp_avg)

vV = IAPWS_IF97('vV_p',p);
pV = p(~isnan(vV));vV = vV(~isnan(vV));
dvVdp_avg = diff(vV)./diff(pV);
pV_avg = (pV(1:end-1)+pV(2:end))/2;
dvVdp = IAPWS_IF97('dvVdp_p',pV);
subplot(1,2,2),plot(pV,dvVdp,pV_avg,dvVdp_avg)
end