%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Space DGFEM scheme for shallow water - well-balanced tests at DG0/1
%
% Kent, T., & Bokhove, O. (2020). Ensuring 'well-balanced'shallow water
% flows via a discontinuous Galerkin finite element method: issues at
% lowest order. arXiv preprint arXiv:2006.03370.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script integrates the 1D SWEs using the space-DGFEM theory of
% Rhebergen et al. (2008) to elucidate the well-balanced issue at DG0
% reported in detail in arXiv:2006.03370.

% 1D SWEs: 

%         h_t + (hu)_x = 0
%         (hu)_t + (hu^2 + 1/2 gh^2)_x = -gh b_x
%         b_t = 0

% Rewrite as a nonconservative system:  

%         U_t + F(U)_x + G(U) U_x = 0

% with U = (h,hu,b)^T
%
% See eqs (1), (7), (8) in arXiv:2006.03370

% DGFEM, after Rhebergen et al. 2008: U(x,t) is approximated in each element
% by its mean only (DG0; piecewise constant) or mean and slope (DG1;
% piecewise linear).

% NOTE: for predictive simulations a slope limiter is required where 
% solution is discontinuous.


%% 
clear; clf;

%% choose order: 0 or 1
% DG = 0; %DG0?
DG = 1; %DG1?

%% Set up grid etc.
Nk=200; %number of cells - try 50, 100, 200.
Nf=Nk+1; %number of edges
L=1; %length of domain
Kk=L/Nk; %length of cell
x=0:Kk:L; %edge coordinates
xc=Kk/2:Kk:L-Kk/2; %cell center coordinates
Fr = 1.9;
FF = Fr*Fr;
cfl = 0.1; % CFL number
Neq = 3; % # of equations in system

%% parabolic ridge topography
% Project initial condition onto elements
% ICs as functions of x

bc = 0.5; %height 
xp = 0.5; %position in domain
a = 0.1; %width
% NOte manual intervention here to ensure b is continous where the parabolic ridge
% meets the constant topography -- works fine for known values but should automate this.
if (Nk == 200) 
    bmin = (25e-04)/12; %Nk=200
elseif (Nk == 100) 
    bmin = (25e-04)/3; %Nk=100
elseif (Nk == 50) 
	bmin = 1/300; % Nk= 50
else
    error('Choose Nk = 50, 100, or 200');
end

fun_b = @(x) max(bmin, bc*(1 - ((x - L*xp).^2)./a^2)); % parabolic ridge
fun_bexp = @(x) 0.5*exp(-500*(x-0.5).^2);
fun_h = @(x) 1 - fun_b(x); % h: h+b=1 (rest flow)
fun_hu = @(x) 0*x; % hu(x0) = 0 (rest flow)

[ bm_ic , bs_ic ] = init_cond_DGFEM_fun( x , Kk , fun_b);
[ hm_ic , hs_ic ] = init_cond_DGFEM_fun( x , Kk , fun_h);
[ hum_ic , hus_ic ] = init_cond_DGFEM_fun( x , Kk , fun_hu);

Um0 = zeros(Neq,Nk);
Us0 = zeros(Neq,Nk);

% MEAN
Um0(1,:) = 1-bm_ic; % h
Um0(2,:) = hum_ic; % hu
Um0(3,:) = bm_ic; % b
% SLOPE
Us0(1,:) = -bs_ic;
Us0(2,:) = hus_ic;
Us0(3,:) = bs_ic;

figure(11);
for k = 1:(length(x)-1)
   plot([x(k), x(k+1)],[Um0(3,k)-Us0(3,k),Um0(3,k)+Us0(3,k)],'k'); hold on;
end

%%
tmax = 10.;
tn = 0;
Um = Um0;
Us = Us0;

Nmeas = 100;
dtmeasure = tmax/Nmeas;
tmeasure = dtmeasure;


%%
% set up arrays for saving data

h_m = zeros(Nmeas+1,Nk); 
hu_m = zeros(Nmeas+1,Nk);
b_m = zeros(Nmeas+1,Nk);

h_m(1,:) = Um0(1,:); 
hu_m(1,:) = Um0(2,:); 
b_m(1,:) = Um0(3,:);

h_s = zeros(Nmeas+1,Nk); 
hu_s = zeros(Nmeas+1,Nk);
b_s = zeros(Nmeas+1,Nk);

h_s(1,:) = Us0(1,:); 
hu_s(1,:) = Us0(2,:); 
b_s(1,:) = Us0(3,:);

%%
axislim = 0.1;
isave = 2;

while (tn <= tmax)
    
    Uminus = Um + Us; % left
    Uplus = Um - Us; % right
        
        dt = cfl*Kk./max([abs(Uminus(2,:)./Uminus(1,:) + sqrt((Fr^-2)*Uminus(1,:))),...
            abs(Uminus(2,:)./Uminus(1,:) - sqrt((Fr^-2)*Uminus(1,:))),...
            abs(Uplus(2,:)./Uplus(1,:) + sqrt((Fr^-2)*Uplus(1,:))),...
            abs(Uplus(2,:)./Uplus(1,:) - sqrt((Fr^-2)*Uplus(1,:)))]);    

    % update time
    tn = tn + dt;
    if (tn > tmeasure)
        dt = dt - (tn - tmeasure) + 1.e-10;
        tn = tmeasure+1.e-10;
    end
    
    % step forward in time dt - use RK3
    [ UUm , UUs ] = RK3_DG1( dt, Um, Us, Fr);
    
    Um = UUm;
    Us = DG*UUs; % if DG0, set slopes to zero

    if (tn > tmeasure)
        
        tmeasure = tmeasure + dtmeasure;     
        
        h_m(isave,:) = Um(1,:);
        hu_m(isave,:) = Um(2,:);
        b_m(isave,:) = Um(3,:);
        
        h_s(isave,:) = Us(1,:);
        hu_s(isave,:) = Us(2,:);
        b_s(isave,:) = Us(3,:);
        
        zm = Um(1,:) + Um(3,:);
        zs = Us(1,:) + Us(3,:);
              
        figure(2); % plot at equal intervals
        
        subplot(Neq-1,1,1);
        for k = 1:(length(x)-1)
            plot([x(k), x(k+1)],[zm(k)-zs(k),zm(k)+zs(k)],'k-'); hold on;
            plot([x(k), x(k+1)],[Um(3,k)-Us(3,k),Um(3,k)+Us(3,k)],'k-'); hold on;
            plot([x(k), x(k+1)],[Um0(3,k)-Us0(3,k),Um0(3,k)+Us0(3,k)],'r:'); hold on;
        end
        hold off;
        text(0.75*L,0.75,['t=',num2str(tn)],'Interpreter','tex','fontsize',18);
        axis([0 1 -0.1 1.1]);
        xlabel('x','fontsize',18); ylabel('h+b','fontsize',18);
        title('1D SWEs', 'Interpreter','tex','fontsize',18);
        
        subplot(Neq-1,1,2);
        for k = 1:(length(x)-1)
            plot([x(k), x(k+1)],[Um(2,k)-Us(2,k),Um(2,k)+Us(2,k)],'k-'); hold on;
        end
        hold off;
        axis([0 1 min(hum_ic)-axislim max(hum_ic)+axislim]);
        xlabel('x','fontsize',18); ylabel('hu','fontsize',18);
        
        drawnow; pause(0.1);
        isave = isave+1;
        
    end % if
    
    
end % while

%% optional plot: hu and h+b at tmax to check well-balancedness
% figure(3);
% subplot(2,1,1);
% plot(xc,Um(2,:)); hold on;
% plot([0,1],10*[eps,eps],'k');hold on;
% plot([0,1], -10*[eps,eps],'k');
% 
% subplot(2,1,2);
% plot(xc,Um(1,:)+Um(3,:)-1); hold on;
% plot([0,1],10*[eps,eps],'k');hold on;
% plot([0,1],-10*[eps,eps],'k');


%% SAVE DATA
PDE.hm = h_m;
PDE.hum = hu_m;
PDE.bm = b_m;

PDE.hs = h_s;
PDE.hus = hu_s;
PDE.bs = b_s;

MESH.L = L;
MESH.Nk = Nk;
MESH.x = x;
MESH.xc = xc;

% make directory path
data_path = strcat(pwd, {'/'}, 'data/'); 
data_path = strjoin(data_path);

DG = num2str(DG);
res = num2str(Nk);
Tmax = strrep(num2str(tmax), '.', '_');
FR = strrep(num2str(Fr), '.','_');

hov_fname = sprintf('topog_DG=%s_res=%s_tmax=%s_Fr=%s',DG,res,Tmax,FR);

save(fullfile(data_path, hov_fname),'PDE','MESH');