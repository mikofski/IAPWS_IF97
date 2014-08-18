function [p,T,resnorm,F,exitflag,output,jacob] = pT_uv(u,v,p0,T0,options)
% PT_UV pressure and temperature from internal energy and specific volume
%
% NOTE: Pressure is very sensitive to enthalpy and internal energy in the
% saturated liquid side. IE: Very large changes in pressuer results from
% very small changes in internal energy. Therefore use as many significant
% digits as possible for internal energy and use the relatively high tolerance
% for optimisation.
%
% [P,T] = PT_UV(U,V) calculate pressure (p [MPa]) and temperature (T [K])
% at given internal energy (u [kJ/kg]) and specific volume (v [m^3/kg])
%
% [P,T] = PT_UV(...,P0,T0) provide initial guess to improve conversion
%
% [P,T] = PT_UV(...,OPTIONS) provide optimisation options
%
% [...,RESNORM,F,EXITFLAG,OUTPUT,JACOB] = PT_UV(...) return optimisation
% information.
%
% See also optimset, optimget, IAPWS_IF97, newtonraphson
%
% Example:
% >> [p,T,resnorm] = pT_uv(230.3137076633380,0.001014500255655)
% p =
%     0.1013
% T =
%   328.1500
% resnorm =
%    9.9476e-13
%
% Copyright (c) 2014 Mark Mifofski

if nargin<3
    p0 = 0.1013;
    T0 = 55+273.15;
end
if nargin<5
    options = optimset('TolFun',1e-12,'TolX',1e-12,'Display','off'); % options
end
x0 = [p0,T0];
h_pT = @(p,T)IAPWS_IF97('h_pT',p,T);
v_pT = @(p,T)IAPWS_IF97('v_pT',p,T);
u_pT = @(p,T)h_pT(p,T)-p*v_pT(p,T); % [kJ/kg] internal energy
Fu_pT = @(p,T,u)u-u_pT(p,T);
Fv_pT = @(p,T,v)v-v_pT(p,T);
fun = @(x)[Fu_pT(x(1),x(2),u),Fv_pT(x(1),x(2),v)];
[x,resnorm,F,exitflag,output,jacob] = NewtonRaphson.newtonraphson(fun,x0,options);
p = x(1);
T = x(2);
end
