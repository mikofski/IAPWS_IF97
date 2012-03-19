function k = k_pT(p,T)
% k = k_pT(p,T)
%   thermal conductivity, k [W/m/K], as a function of pressure, p [MPa], and temperature, T [K]
% based on Revised Release on the IAPWS Formulation 1985 for the Thermal Conductivity of Ordinary Water Substance, 2008
% Reference: http://www.iapws.org/
% June 23, 2009
% Mark Mikofski
%% size of inputs
dim = size(p);
initnan = NaN(dim);
k = initnan;
%% constants and calculated
Tstar = 647.26; % [K]
rhostar = 317.7; % [K]
kstar = 1; % [W/m/K]
a = [0.0102811, 0.0299621, 0.0156146, -0.00422464];
b = [-0.397070, 0.400302, 1.060000];
B = [-0.171587, 2.392190];
d = [0.0701309, 0.0118520, 0.00169937, -1.0200];
C = [0.642857, -4.11717, -6.17937, 0.00308976, 0.0822994, 10.0932];
Tmin = 273.16; % [K] minimum temperature is triple point
TB13 = 623.15; % [K] temperature at boundary between region 1 and 3
TB23 = 863.15; % [K] temperature of boundary between region 2 and 3
Tmax = 1073.15; % [K] maximum valid temperature
pmin = psat_T(Tmin); % [MPa] minimum pressure is 611.657 Pa
pmax = 100; % [MPa] maximum valid pressure
psat = psat_T(T); % [MPa] saturation pressures
pB23 = initnan; valid = T>=TB13 & T<=TB23;
pB23(valid) = pB23_T(T(valid)); % [MPa] pressure on boundary between region 2 and region 3
%% valid ranges
valid1 = p>=psat & p<=pmax & T>=Tmin & T<=TB13; % valid range for region 1, include B13 in region 1
valid2 = p>=pmin & ((T>=Tmin & T<=TB13 & p<=psat) | (T>TB13 & T<=TB23 & p<=pB23) | (T>TB23 & T<=Tmax & p<=pmax)); % valid range for region 2, include B23 in region 2
valid3 = p>pB23 & p<=pmax & T>TB13 & T<TB23;
if any(any(valid1))
    T1 = T(valid1);rho1bar = 1./v1_pT(p(valid1),T1)/rhostar;T1bar = T1/Tstar;
    k0 = sqrt(T1bar).*(a(1) + (a(2) + (a(3) + a(4).*T1bar).*T1bar).*T1bar);
    k1 = b(1) + b(2)*rho1bar + b(3)*exp(B(1)*(rho1bar + B(2)).^2);
    deltaTbar = abs(T1bar-1)+C(4);Q = 2 + C(5)./deltaTbar.^0.6;S = (1./deltaTbar).*(T1bar>=1) + (C(6)./deltaTbar.^0.6).*(T1bar<1);
    k2 = (d(1)./T1bar.^10 + d(2)).*rho1bar.^1.8.*exp(C(1)*(1-rho1bar.^2.8)) + d(3)*S.*rho1bar.^Q.*exp(Q./(1+Q).*(1-rho1bar.^(1+Q))) + d(4)*exp(C(2)*T1bar.^1.5 + C(3)./rho1bar.^5);
    k(valid1) = (k0+k1+k2)*kstar;
end
if any(any(valid2))
    T2 = T(valid2);rho2bar = 1./v2_pT(p(valid2),T2)/rhostar;T2bar = T2/Tstar;
    k0 = sqrt(T2bar).*(a(1) + (a(2) + (a(3) + a(4).*T2bar).*T2bar).*T2bar);
    k1 = b(1) + b(2)*rho2bar + b(3)*exp(B(1)*(rho2bar + B(2)).^2);
    deltaTbar = abs(T2bar-1)+C(4);Q = 2 + C(5)./deltaTbar.^0.6;S = (1./deltaTbar).*(T2bar>=1) + (C(6)./deltaTbar.^0.6).*(T2bar<1);
    k2 = (d(1)./T2bar.^10 + d(2)).*rho2bar.^1.8.*exp(C(1)*(1-rho2bar.^2.8)) + d(3)*S.*rho2bar.^Q.*exp(Q./(1+Q).*(1-rho2bar.^(1+Q))) + d(4)*exp(C(2)*T2bar.^1.5 + C(3)./rho2bar.^5);
    k(valid2) = (k0+k1+k2)*kstar;
end
if any(any(valid3))
    T3 = T(valid3);rho3bar = 1./v_pT(p(valid3),T3)/rhostar;T3bar = T3/Tstar;
    k0 = sqrt(T3bar).*(a(1) + (a(2) + (a(3) + a(4).*T3bar).*T3bar).*T3bar);
    k1 = b(1) + b(2)*rho3bar + b(3)*exp(B(1)*(rho3bar + B(2)).^2);
    deltaTbar = abs(T3bar-1)+C(4);Q = 2 + C(5)./deltaTbar.^0.6;S = (1./deltaTbar).*(T3bar>=1) + (C(6)./deltaTbar.^0.6).*(T3bar<1);
    k2 = (d(1)./T3bar.^10 + d(2)).*rho3bar.^1.8.*exp(C(1)*(1-rho3bar.^2.8)) + d(3)*S.*rho3bar.^Q.*exp(Q./(1+Q).*(1-rho3bar.^(1+Q))) + d(4)*exp(C(2)*T3bar.^1.5 + C(3)./rho3bar.^5);
    k(valid3) = (k0+k1+k2)*kstar;
end
end