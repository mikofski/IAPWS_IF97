function k = k_ph(p,h)
% k = k_ph(p,h)
%   thermal conductivity, k [W/m/K], as a function of pressure, p [MPa], and enthalpy, h [kJ/kg]
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
T2bcsat = 554.485; % [K] saturation temperature at 5.85 kJ/kg/K isentropic line between region 2b and 2c
TB13 = 623.15; % [K] temperature at boundary between region 1 and 3
Tmax = 1073.15; % [K] maximum temperature
pmin = psat_T(Tmin); % [MPa] minimum pressure is 611.657 Pa
p2ab = 4; % [MPa] pressure along boundary between region 2a and 2b
p2bcsat = psat_T(T2bcsat); % [MPa] saturation pressure at 5.85 kJ/kg/K isentropic line between region 2b and 2c
pB13sat = psat_T(TB13); % [MPa] saturation pressure at boundary between region 1 and 3, 16.5291643 MPa
pmax = 100; % [MPa] maximum pressure
h1B13L = h1_pT(pB13sat,TB13); % [kJ/kg] saturated liquid enthalpy at boundary between region 1, region 3 and region 4
h2B13V = h2_pT(pB13sat,TB13); % [kJ/kg] saturated vapor enthalpy at boundary between region 2, region 3 and region 4
%% calculated matrices
Tsat = Tsat_p(p); % [K] saturation temperatures
h1min = initnan;h2max = initnan;valid = p>=pmin & p<=pmax;pvalid = p(valid); % initialize matricies with NaN and set valid range of parameters
if any(any(valid))
    Tmin = Tmin*ones(dim);Tmax = Tmax*ones(dim); % copy to matrix of size dim
    h1min(valid) = h1_pT(pvalid,Tmin(valid)); % [kJ/kg] minimum enthalpies
    h2max(valid) = h2_pT(pvalid,Tmax(valid)); % [kJ/kg] maximum enthalpies in region 2
end
h1L = initnan;h2V = initnan;valid = p>=pmin & p<=pB13sat;pvalid = p(valid); % initialize matricies with NaN and set valid range of parameters
if any(any(valid))
    h1L(valid) = h1_pT(pvalid,Tsat(valid)); % [kJ/kg] saturated liquid enthalpies in region 1
    h2V(valid) = h2_pT(pvalid,Tsat(valid)); % [kJ/kg] saturated vapor enthalpies in region 2
end
h1B13 = initnan;h3ab = initnan;h2B23 = initnan;valid = p>=pB13sat & p<=pmax;pvalid = p(valid); % initialize matricies with NaN and set valid range of parameters
if any(any(valid))
    TB13 = TB13*ones(dim); % copy to matrix of size dim
    h1B13(valid) = h1_pT(pvalid,TB13(valid)); % [kJ/kg] enthalpies on boundary between region 1 and region 3
    h3ab(valid) = h3ab_p(pvalid); % [kJ/kg] enthalpies on critical entropy isentropic line between regions 3a and region 3b
    h2B23(valid) = h2_pT(pvalid,TB23_p(pvalid)); % [kJ/kg] enthalpies on boundary between region 2 and region 3
end
h2bc = initnan;valid = p>=p2bcsat & p<=pmax; % initialize matricies with NaN and set valid range of parameters
h2bc(valid) = h2bc_p(p(valid)); % [kJ/kg] enthalpies on boundary between region 2b and 2c
p3sat = pB13sat*ones(dim);valid = h>=h1B13L & h<=h2B13V; % % do NOT use NaN to initialize p3sat, b/c for h<h1B13L or h>h2B13V p>NaN = 0, instead use pB13sat
if any(any(valid))
    p3sat(valid) = p3sat_h(h(valid)); % [MPa] saturation pressure on boundary between region 3 and 4
end
%% valid ranges
valid1 = (p>=pmin & p<=pB13sat & h>=h1min & h<=h1L) | (p>pB13sat & p<=pmax & h>=h1min & h<=h1B13); % valid range for region 1
valid2a = p>=pmin & p<=p2ab & h>h2V & h<=h2max; % valid range for region 2a
valid2b = (p>p2ab & p<=p2bcsat & h>h2V & h<=h2max) | (p>p2bcsat & p<=pmax & h>h2bc & h<=h2max); % valid range for region 2b
valid2c = (p>p2bcsat & p<=pB13sat & h>h2V & h<=h2bc) | (p>pB13sat & p<=pmax & h>h2B23 & h<=h2bc); % valid range for region 2c
valid3a = p>p3sat & p<=pmax & h>h1B13 & h<=h3ab; % valid range for region 3a
valid3b = p>p3sat & p<=pmax & h>h3ab & h<=h2B23; % valid range for region 3b
valid4a = p>=pmin & p<=pB13sat & h>h1L & h<=h2V; % valid range for region 4a
valid4b = p>pB13sat & p<=p3sat & h>h1B13L & h<=h2B13V; % valid range for region 4b
if any(any(valid1))
    p1 = p(valid1);T1 = T1_ph(p1,h(valid1));rho1bar = 1./v1_pT(p1,T1)/rhostar;T1bar = T1/Tstar;
    k0 = sqrt(T1bar).*(a(1) + (a(2) + (a(3) + a(4).*T1bar).*T1bar).*T1bar);
    k1 = b(1) + b(2)*rho1bar + b(3)*exp(B(1)*(rho1bar + B(2)).^2);
    deltaTbar = abs(T1bar-1)+C(4);Q = 2 + C(5)./deltaTbar.^0.6;S = (1./deltaTbar).*(T1bar>=1) + (C(6)./deltaTbar.^0.6).*(T1bar<1);
    k2 = (d(1)./T1bar.^10 + d(2)).*rho1bar.^1.8.*exp(C(1)*(1-rho1bar.^2.8)) + d(3)*S.*rho1bar.^Q.*exp(Q./(1+Q).*(1-rho1bar.^(1+Q))) + d(4)*exp(C(2)*T1bar.^1.5 + C(3)./rho1bar.^5);
    k(valid1) = (k0+k1+k2)*kstar;
end
if any(any(valid2a))
    p2a = p(valid2a);T2a = T2a_ph(p2a,h(valid2a));rho2abar = 1./v2_pT(p2a,T2a)/rhostar;T2abar = T2a/Tstar;
    k0 = sqrt(T2abar).*(a(1) + (a(2) + (a(3) + a(4).*T2abar).*T2abar).*T2abar);
    k1 = b(1) + b(2)*rho2abar + b(3)*exp(B(1)*(rho2abar + B(2)).^2);
    deltaTbar = abs(T2abar-1)+C(4);Q = 2 + C(5)./deltaTbar.^0.6;S = (1./deltaTbar).*(T2abar>=1) + (C(6)./deltaTbar.^0.6).*(T2abar<1);
    k2 = (d(1)./T2abar.^10 + d(2)).*rho2abar.^1.8.*exp(C(1)*(1-rho2abar.^2.8)) + d(3)*S.*rho2abar.^Q.*exp(Q./(1+Q).*(1-rho2abar.^(1+Q))) + d(4)*exp(C(2)*T2abar.^1.5 + C(3)./rho2abar.^5);
    k(valid2a) = (k0+k1+k2)*kstar;
end
if any(any(valid2b))
    p2b = p(valid2b);T2b = T2b_ph(p2b,h(valid2b));rho2bbar = 1./v2_pT(p2b,T2b)/rhostar;T2bbar = T2b/Tstar;
    k0 = sqrt(T2bbar).*(a(1) + (a(2) + (a(3) + a(4).*T2bbar).*T2bbar).*T2bbar);
    k1 = b(1) + b(2)*rho2bbar + b(3)*exp(B(1)*(rho2bbar + B(2)).^2);
    deltaTbar = abs(T2bbar-1)+C(4);Q = 2 + C(5)./deltaTbar.^0.6;S = (1./deltaTbar).*(T2bbar>=1) + (C(6)./deltaTbar.^0.6).*(T2bbar<1);
    k2 = (d(1)./T2bbar.^10 + d(2)).*rho2bbar.^1.8.*exp(C(1)*(1-rho2bbar.^2.8)) + d(3)*S.*rho2bbar.^Q.*exp(Q./(1+Q).*(1-rho2bbar.^(1+Q))) + d(4)*exp(C(2)*T2bbar.^1.5 + C(3)./rho2bbar.^5);
    k(valid2b) = (k0+k1+k2)*kstar;
end
if any(any(valid2c))
    p2c = p(valid2c);T2c = T2c_ph(p2c,h(valid2c));rho2cbar = 1./v2_pT(p2c,T2c)/rhostar;T2cbar = T2c/Tstar;
    k0 = sqrt(T2cbar).*(a(1) + (a(2) + (a(3) + a(4).*T2cbar).*T2cbar).*T2cbar);
    k1 = b(1) + b(2)*rho2cbar + b(3)*exp(B(1)*(rho2cbar + B(2)).^2);
    deltaTbar = abs(T2cbar-1)+C(4);Q = 2 + C(5)./deltaTbar.^0.6;S = (1./deltaTbar).*(T2cbar>=1) + (C(6)./deltaTbar.^0.6).*(T2cbar<1);
    k2 = (d(1)./T2cbar.^10 + d(2)).*rho2cbar.^1.8.*exp(C(1)*(1-rho2cbar.^2.8)) + d(3)*S.*rho2cbar.^Q.*exp(Q./(1+Q).*(1-rho2cbar.^(1+Q))) + d(4)*exp(C(2)*T2cbar.^1.5 + C(3)./rho2cbar.^5);
    k(valid2c) = (k0+k1+k2)*kstar;
end
if any(any(valid3a))
    p3a = p(valid3a);h3a = h(valid3a);T3abar = T3a_ph(p3a,h3a)/Tstar;rho3abar = 1./v3a_ph(p3a,h3a)/rhostar;
    k0 = sqrt(T3abar).*(a(1) + (a(2) + (a(3) + a(4).*T3abar).*T3abar).*T3abar);
    k1 = b(1) + b(2)*rho3abar + b(3)*exp(B(1)*(rho3abar + B(2)).^2);
    deltaTbar = abs(T3abar-1)+C(4);Q = 2 + C(5)./deltaTbar.^0.6;S = (1./deltaTbar).*(T3abar>=1) + (C(6)./deltaTbar.^0.6).*(T3abar<1);
    k2 = (d(1)./T3abar.^10 + d(2)).*rho3abar.^1.8.*exp(C(1)*(1-rho3abar.^2.8)) + d(3)*S.*rho3abar.^Q.*exp(Q./(1+Q).*(1-rho3abar.^(1+Q))) + d(4)*exp(C(2)*T3abar.^1.5 + C(3)./rho3abar.^5);
    k(valid3a) = (k0+k1+k2)*kstar;
end
if any(any(valid3b))
    p3b = p(valid3b);h3b = h(valid3b);T3bbar = T3b_ph(p3b,h3b)/Tstar;rho3bbar = 1./v3b_ph(p3b,h3b)/rhostar;
    k0 = sqrt(T3bbar).*(a(1) + (a(2) + (a(3) + a(4).*T3bbar).*T3bbar).*T3bbar);
    k1 = b(1) + b(2)*rho3bbar + b(3)*exp(B(1)*(rho3bbar + B(2)).^2);
    deltaTbar = abs(T3bbar-1)+C(4);Q = 2 + C(5)./deltaTbar.^0.6;S = (1./deltaTbar).*(T3bbar>=1) + (C(6)./deltaTbar.^0.6).*(T3bbar<1);
    k2 = (d(1)./T3bbar.^10 + d(2)).*rho3bbar.^1.8.*exp(C(1)*(1-rho3bbar.^2.8)) + d(3)*S.*rho3bbar.^Q.*exp(Q./(1+Q).*(1-rho3bbar.^(1+Q))) + d(4)*exp(C(2)*T3bbar.^1.5 + C(3)./rho3bbar.^5);
    k(valid3b) = (k0+k1+k2)*kstar;
end
if any(any(valid4a))
    p4a = p(valid4a);Tsat4a = Tsat(valid4a);h1L4a = h1L(valid4a);
    x = (h(valid4a)-h1L4a)./(h2V(valid4a)-h1L4a); % quality
    v1L = v1_pT(p4a,Tsat4a); % [m^3/kg] saturated liquid specific volumes
    v2V = v2_pT(p4a,Tsat4a); % [m^3/kg] saturated vapor specific volumes
    T4abar = Tsat4a/Tstar;rho4abar = 1./(v1L + x.*(v2V - v1L))/rhostar;
    k0 = sqrt(T4abar).*(a(1) + (a(2) + (a(3) + a(4).*T4abar).*T4abar).*T4abar);
    k1 = b(1) + b(2)*rho4abar + b(3)*exp(B(1)*(rho4abar + B(2)).^2);
    deltaTbar = abs(T4abar-1)+C(4);Q = 2 + C(5)./deltaTbar.^0.6;S = (1./deltaTbar).*(T4abar>=1) + (C(6)./deltaTbar.^0.6).*(T4abar<1);
    k2 = (d(1)./T4abar.^10 + d(2)).*rho4abar.^1.8.*exp(C(1)*(1-rho4abar.^2.8)) + d(3)*S.*rho4abar.^Q.*exp(Q./(1+Q).*(1-rho4abar.^(1+Q))) + d(4)*exp(C(2)*T4abar.^1.5 + C(3)./rho4abar.^5);
    k(valid4a) = (k0+k1+k2)*kstar;
end
if any(any(valid4b))
    p4b = p(valid4b); Tsat4b = Tsat(valid4b);
    v3L = vL_p(p4b); v3V = vV_p(p4b);
    h3L = h3_rhoT(1./v3L,Tsat4b);
    h3V = h3_rhoT(1./v3V,Tsat4b);
    x = (h(valid4b)-h3L)./(h3V-h3L); % quality
    T4bbar = Tsat4b/Tstar;rho4bbar = 1./(v3L + x.*(v3V - v3L))/rhostar;
    k0 = sqrt(T4bbar).*(a(1) + (a(2) + (a(3) + a(4).*T4bbar).*T4bbar).*T4bbar);
    k1 = b(1) + b(2)*rho4bbar + b(3)*exp(B(1)*(rho4bbar + B(2)).^2);
    deltaTbar = abs(T4bbar-1)+C(4);Q = 2 + C(5)./deltaTbar.^0.6;S = (1./deltaTbar).*(T4bbar>=1) + (C(6)./deltaTbar.^0.6).*(T4bbar<1);
    k2 = (d(1)./T4bbar.^10 + d(2)).*rho4bbar.^1.8.*exp(C(1)*(1-rho4bbar.^2.8)) + d(3)*S.*rho4bbar.^Q.*exp(Q./(1+Q).*(1-rho4bbar.^(1+Q))) + d(4)*exp(C(2)*T4bbar.^1.5 + C(3)./rho4bbar.^5);
    k(valid4b) = (k0+k1+k2)*kstar;
end
end
