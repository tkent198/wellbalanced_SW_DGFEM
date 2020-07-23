function [ RHSm , RHSs ] = space_DG1_RHS( Um , Us , Fr )
% SWE space discretisation for DG1 : RHS

% input:
% > Um,s - means and slopes
% > Fr = Froude number


ar_size = size(Um);
Nk = ar_size(2);
Kk = 1/Nk;
Neq = ar_size(1);

Flux = zeros(Neq,Nk+1);
VNC = zeros(Neq,Nk+1);
INT_F = zeros(Neq,Nk);
INT_Gm = zeros(Neq,Nk);
INT_Gs = zeros(Neq,Nk);

FF = Fr*Fr;
g = FF^-1;

cm = 1/sqrt(3);

% Riemann data
Uminus = Um + Us; % left
Uplus = Um - Us; % right

% periodic BCs
for j = 1:Nk-1
    [Flux(:,j+1), VNC(:,j+1)] = SWfluxNCP(Uminus(:,j),Uplus(:,j+1),Fr);
end
[Flux(:,Nk+1), VNC(:,Nk+1)] = SWfluxNCP(Uminus(:,Nk),Uplus(:,1),Fr);
Flux(:,1) = Flux(:,Nk+1);
VNC(:,1) = VNC(:,Nk+1);

Pp = 0.5*VNC + Flux;
Pm = -0.5*VNC + Flux;

% compute integral for slope in local coord
% via 2-pt Gauss quadrature when appropriate

% F integral: gauss quadrature
INT_F(1,:) = 2*Um(2,:);
INT_F(2,:) = ((Um(2,:) + Us(2,:)*-cm).^2)./(Um(1,:) + Us(1,:)*-cm) + ...
    (0.5/FF)*((Um(1,:) + Us(1,:)*-cm).^2) + ...
    ((Um(2,:) + Us(2,:)*cm).^2)./(Um(1,:) + Us(1,:)*cm) + ...
    (0.5/FF)*((Um(1,:) + Us(1,:)*cm).^2);
INT_F(3,:) = 0;

% G integral: exact
INT_Gm(1,:) = 0;
INT_Gm(2,:) = g*2*Um(1,:).*Us(3,:);
INT_Gm(3,:) = 0;

INT_Gs(1,:) = 0;
INT_Gs(2,:) = g*(2/3)*Us(1,:).*Us(3,:);
INT_Gs(3,:) = 0;

% RHS of space discretisation for means and slopes
RHSm = -(Pp(:,2:Nk+1) - Pm(:,1:Nk) + INT_Gm)./Kk;
RHSs = -3*(Pp(:,2:Nk+1) + Pm(:,1:Nk) - INT_F + INT_Gs)./Kk;


end

