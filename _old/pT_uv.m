function [p,T,resnorm,F,exitflag,output,jacob] = pT_uv(u,v,p0,T0,options)
% PT_UV pressure and temperature from internal energy and specific volume
%
% NOTE: Pressure is very sensitive to enthalpy and internal energy in the
% compressed liquid state. IE: Very large changes in pressure results from
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
% Example 1: Compressed liquid
% >> [p,T,resnorm] = pT_uv(230.2110313177497,0.001014500244461827)
% p =
%     0.1013249999991855
% T =
%   328.1500000000002
% resnorm =
%    7.105427357601332e-13
%
% Example 2: Superheated steam
% >> [p,T,resnorm] = pT_uv(2881.2, 0.4742)
% p =
%     0.5999
% T =
%   622.9358
% resnorm =
%    9.0949e-13
%
% Copyright (c) 2014 Mark Mifofski

% if no initial guesses for p and T provided use critical point properties
if nargin<3
    p0 = 22.09; % [MPa] critical pressure
    T0 = 647.3; % [K] critical temperature
end
% default optimisation options if none provided
% maximum tolerance of the norm of the residuals, TolFun = 1e-12
% minimum tolerance of the relative maximum stepsize, TolX = 1e-12
% turn off display of optimisation steps, Display = 'off'
% default maximum iterations, MaxIter = 100
if nargin<5
    options = optimset('TolFun',1e-12,'TolX',1e-12,'Display','off'); % options
end
% backwards function
x0 = [p0,T0]; % initial guess
h_pT = @(p,T)IAPWS_IF97('h_pT',p,T); % [kJ/kg] enthalpy function handle
v_pT = @(p,T)IAPWS_IF97('v_pT',p,T); % [m^3/kg] specific volume function handle
u_pT = @(p,T)h_pT(p,T)-p*v_pT(p,T)*1000; % [kJ/kg] internal energy func handle
% residuals used in optimisation, F=0 when solution is found
Fu_pT = @(p,T,u)u-u_pT(p,T); % internal energy residual
Fv_pT = @(p,T,v)v-v_pT(p,T); % specific volumne residual
% backwards function to solve
fun = @(x)[Fu_pT(x(1),x(2),u),Fv_pT(x(1),x(2),v)];
% solve backwards function for x
[x,resnorm,F,exitflag,output,jacob] = newtonraphson(fun,x0,options);
if exitflag~=1
    warning(output.message)
    fprintf('Solver failed to converge.\n\t* Please check that inputs are in range,')
    fprintf('\n\t* specify initial guesses or\n\t* increase number of significant digits.\n')
end
% parse x for return
p = x(1); % [MPa] pressure
T = x(2); % [K] temperature
end
