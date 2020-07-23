function [ RHS ] = space_DG0_RHS( U , Fr )
% space discretisation for DG1 : RHS

% input:
% > Um,s - means and slopes
% > Fr = Froude number


ar_size = size(U);
Nk = ar_size(2);
Kk = 1/Nk;
Neq = ar_size(1);

Flux = zeros(Neq,Nk+1);
VNC = zeros(Neq,Nk+1);


FF = Fr*Fr;

% periodic BCs
for j = 1:Nk-1
    [Flux(:,j+1), VNC(:,j+1)] = SWfluxNCP(U(:,j),U(:,j+1),Fr);
end
[Flux(:,Nk+1), VNC(:,Nk+1)] = SWfluxNCP(U(:,Nk),U(:,1),Fr);
Flux(:,1) = Flux(:,Nk+1);
VNC(:,1) = VNC(:,Nk+1);

% % neumann BCs
% for j = 1:Nk-1
%     [Flux(:,j+1), VNC(:,j+1)] = SWfluxNCP(U(:,j),U(:,j+1),Fr);
% end
% Flux(:,Nk+1) = Flux(:,Nk);
% VNC(:,Nk+1) = VNC(:,Nk);
% 
% Flux(:,1) = Flux(:,2);
% VNC(:,1) = VNC(:,2);


Pp = 0.5*VNC + Flux;
Pm = -0.5*VNC + Flux;

RHS = -(Pp(:,2:Nk+1) - Pm(:,1:Nk))./Kk;


end

