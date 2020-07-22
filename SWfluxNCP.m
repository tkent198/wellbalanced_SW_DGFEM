function [Flux, VNC] = SWfluxNCP(Ul,Ur,Fr)

% 
% Function gives the intercell HLL flux value for use in FV Godunov numerical
% scheme for the 1D non-dimensional SWEs:
%
% U_t + F(U)_x + G(U) U_x = 0 with U = (h,hu,b) and 
% F = (hu, hu^2 + 1/2 gh^2,0) and G_23 = -gh
% 
%
% INPUT: Ul = (hl,hul,bl) and Ur = (hr,hur,br) and Fr.

hl = Ul(1); 
hr = Ur(1);

ul = Ul(2)/hl;
ur = Ur(2)/hr;

bl = Ul(3);
br = Ur(3);

g = Fr^-2;

% compute left and right wave speeds 
Sl = min(ul - sqrt(g*hl), ur - sqrt(g*hr)); % Sl = ul - al
Sr = max(ul + sqrt(g*hl), ur + sqrt(g*hr)); % Sr = ur + ar

VNC = [0; -g*0.5*(hl + hr)*(bl - br); 0];

if (Sl > 0)
    FluxL = [hl*ul; hl*ul^2 + 0.5*g*hl^2; 0];
    Flux = FluxL - 0.5*VNC;
elseif (Sr < 0)
    FluxR = [hr*ur; hr*ur^2 + 0.5*g*hr^2; 0];
    Flux = FluxR + 0.5*VNC;
else
    FluxL = [hl*ul; hl*ul^2 + 0.5*g*hl^2; 0];
    FluxR = [hr*ur; hr*ur^2 + 0.5*g*hr^2; 0];
    FluxHLL = (FluxL*Sr - FluxR*Sl + Sl*Sr*(Ur - Ul))/(Sr-Sl);
    Flux = FluxHLL - 0.5*VNC*(Sr + Sl)/(Sr - Sl);
end


% % star state (for checking rest flow)
% FluxL = [hl*ul; hl*ul^2 + 0.5*g*hl^2; 0];
% FluxR = [hr*ur; hr*ur^2 + 0.5*g*hr^2; 0];
% FluxHLL = (FluxL*Sr - FluxR*Sl + Sl*Sr*(Ur - Ul))/(Sr-Sl);
% Flux = FluxHLL - 0.5*VNC*(Sr + Sl)/(Sr - Sl);


