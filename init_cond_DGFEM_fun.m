function [ Um0 , Us0 ] = init_cond_DGFEM_fun( x , Kk , fun)

%%% Function computes DG initial data for means and slopes via Guass
%%% quadrature, given:
% > mesh x
% > step size Kk = x(k+1) - x(k)
% > a function fun defining the initial state for u(x0)

%%% OUTPUT
% > [Um0,Us0]: projection of the IC fun onto the basis functions in each
% element, [mean,slope]. Integrals computed via 2-pt Gauss quadrature.
    
    x_m = zeros(size(x)-1);
    x_p = zeros(size(x)-1);

    zeta1 = -1/sqrt(3);
    zeta2 = -zeta1;
    for k = 1:(length(x)-1)
        x_m(k) = 0.5*(x(k) + x(k+1) + zeta1*Kk);
        x_p(k) = 0.5*(x(k) + x(k+1) + zeta2*Kk);
    end
    
    Um0 = 0.5*(fun(x_m) + fun(x_p));
    Us0 = 1.5*(zeta1*fun(x_m) + zeta2*fun(x_p));
    
    for k = 1:(length(x)-1)
        plot([x(k), x(k+1)],[Um0(k)-Us0(k),Um0(k)+Us0(k)],'k-'); hold on;
        title('Initial state check','fontsize',18);
    end
    
end

