%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Space DGFEM scheme for shallow water - well-balanced tests at DG0/1
%
% Kent, T., & Bokhove, O. (2020). Ensuring 'well-balanced'shallow water
% flows via a discontinuous Galerkin finite element method: issues at
% lowest order. arXiv preprint arXiv:2006.03370.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plotting routines for data saved from run_DGFEM_SW.m

%% working directory:

data_path = strcat(pwd, {'/'}, 'data/'); 
data_path = strjoin(data_path);
fig_path = strcat(pwd, {'/'}, 'figs/'); 
fig_path = strjoin(fig_path);


%% load data: 
% this needs choosing manually according to the simulation details. Check
% data dir.
fname1 = 'topog_DG=1_res=200_tmax=10_Fr=1_9.mat';
fname2 = 'topog_DG=0_res=200_tmax=10_Fr=1_9.mat';
% fname3 = 'topog_DG=1_res=200_tmax=10_Fr=0_7.mat';
% fname4 = 'topog_DG=0_res=200_tmax=10_Fr=0_7.mat';



DG1Fr19 = load(fullfile(data_path, fname1)); 
DG0Fr19 = load(fullfile(data_path, fname2)); 
% DG1Fr07 = load(fullfile(data_path, fname3)); 
% DG0Fr07 = load(fullfile(data_path, fname4));

% DATAn is a 2-structure
%            PDE: numerical sols h hu and b -- means and sleaps
%           MESH: mesh for PDE solver

%% unpack data from structures:

mesh = DG1Fr19.MESH;
Nk = mesh.Nk;
L = mesh.L;
x = mesh.x;
Kk = L/Nk;

% DG1
hm = DG1Fr19.PDE.hm;
hs = DG1Fr19.PDE.hs;
hum = DG1Fr19.PDE.hum;
hus = DG1Fr19.PDE.hus;
bm = DG1Fr19.PDE.bm;
bs = DG1Fr19.PDE.bs;
zm = hm + bm;
zs = hs + bs;

% DG0
hm0 = DG0Fr19.PDE.hm;
hs0 = DG0Fr19.PDE.hs;
hum0 = DG0Fr19.PDE.hum;
hus0 = DG0Fr19.PDE.hus;
bm0 = DG0Fr19.PDE.bm;
bs0 = DG0Fr19.PDE.bs;
zm0 = hm0 + bm0;
zs0 = hs0 + bs0;

%% fig: h+b at various times

f1 = figure(101);
subplot(1,2,2); % DG1
for k = 1:(length(x)-1)
    plot([x(k), x(k+1)],[zm(end,k)-zs(end,k),zm(end,k)+zs(end,k)],'k-'); hold on;
    plot([x(k), x(k+1)],[zm(1,k)-zs(1,k),zm(1,k)+zs(1,k)],'r:'); hold on;
    plot([x(k), x(k+1)],[bm(end,k)-bs(end,k),bm(end,k)+bs(end,k)],'k-'); hold on;
    plot([x(k), x(k+1)],[bm(1,k)-bs(1,k),bm(1,k)+bs(1,k)],'r:'); hold on;
end
text(0.75*L,0.85,'DG1','Interpreter','tex','fontsize',18);
hold off;
axis([0 1 -0.1 1.1]);
xlabel('x','fontsize',18);
% title('DG1', 'Interpreter','tex','fontsize',18);

subplot(1,2,1); % DG0
for k = 1:(length(x)-1)
    plot([x(k), x(k+1)],[zm0(end,k)-zs0(end,k),zm0(end,k)+zs0(end,k)],'k-'); hold on;
    plot([x(k), x(k+1)],[zm0(1,k)-zs0(1,k),zm0(1,k)+zs0(1,k)],'r:'); hold on;
    plot([x(k), x(k+1)],[bm0(1:25:101,k)-bs0(1:25:101,k),bm0(1:25:101,k)+bs0(1:25:101,k)],'k-'); hold on;
    plot([x(k), x(k+1)],[bm0(1,k)-bs0(1,k),bm0(1,k)+bs0(1,k)],'r:'); hold on;
end
text(0.75*L,0.85,'DG0','Interpreter','tex','fontsize',18);
xlabel('x','fontsize',18); ylabel('h+b, b','fontsize',18);
hold off;
axis([0 1 -0.1 1.1]);
% ylabel('h+b','fontsize',18);
% title('DG0', 'Interpreter','tex','fontsize',18);


f1_fname = 'fig1_res=200_tmax=10_Fr=1_9';

set(f1, 'PaperUnits', 'centimeters');
x_width = 20;
y_width = 10;
set(f1, 'PaperPosition', [0 0 x_width y_width]);
saveas(f1,fullfile(fig_path, f1_fname),'jpg')

%% fig: h+b and hu at tmax

f2 = figure(102);
subplot(2,2,1); % hu for DG0
for k = 1:(length(x)-1)
    plot([x(k), x(k+1)],[hum0(end,k)-hus0(end,k),hum0(end,k)+hus0(end,k)],'k'); hold on;
end
plot([0,1],50*[eps,eps],'k');hold on;
plot([0,1], -50*[eps,eps],'k');
text(0.75*L,25*eps,'DG0','Interpreter','tex','fontsize',18);
ylabel('hu','fontsize',18);
axis([0 1 -50*eps 50*eps]);

subplot(2,2,2); % hu for DG1
for k = 1:(length(x)-1)
    plot([x(k), x(k+1)],[hum(end,k)-hus(end,k),hum(end,k)+hus(end,k)],'k'); hold on;
end
plot([0,1],50*[eps,eps],'k');hold on;
plot([0,1], -50*[eps,eps],'k');
text(0.75*L,25*eps,'DG1','Interpreter','tex','fontsize',18);
% ylabel('hu','fontsize',18);
axis([0 1 -50*eps 50*eps]);

subplot(2,2,3); % h+b for DG0
for k = 1:(length(x)-1)
    plot([x(k), x(k+1)],[zm0(end,k)-zs0(end,k),zm0(end,k)+zs0(end,k)],'k'); hold on;
end
% plot([0,1],1+1*[eps,eps],'k');hold on;
% plot([0,1],1-1*[eps,eps],'k');
% text(0.75*L,0.85,'DG0','Interpreter','tex','fontsize',18);
xlabel('x','fontsize',18); ylabel('h+b','fontsize',18);
axis([0 1 1-0.5e-12 1+0.5e-12]);
gca.YAxis.Exponent = -14;

subplot(2,2,4); % h+b for DG1
for k = 1:(length(x)-1)
    plot([x(k), x(k+1)],[zm(end,k)-zs(end,k),zm(end,k)+zs(end,k)],'k'); hold on;
end
% plot([0,1],1+0.1*[eps,eps],'k');hold on;
% plot([0,1],1-0.1*[eps,eps],'k');
xlabel('x','fontsize',18);
axis([0 1 1-0.5e-12 1+0.5e-12]);
YAxis.Exponent = -14;

f2_fname = 'fig2_res=200_tmax=10_Fr=1_9';
set(f2, 'PaperUnits', 'centimeters');
x_width = 25;
y_width = 15;
set(f2, 'PaperPosition', [0 0 x_width y_width]);
saveas(f2,fullfile(fig_path, f2_fname),'jpg')