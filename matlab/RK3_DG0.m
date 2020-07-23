function [ U ] = RK3_DG0( dt, U, Fr)
% time discretisation for DG1 : RK3 of Shu and Osher (1999)

% input:
% > dt - time step
% > U - means 

[ RHS ] = space_DG0_RHS( U , Fr );

U1 = U - dt*RHS;

[ RHS ] = space_DG0_RHS( U1 , Fr );

U2 = 0.25*(3*U + U1 + dt*RHS);

[ RHS ] = space_DG0_RHS( U2 , Fr );

UU = (1/3)*(U + 2*U2 + 2*dt*RHS);


U = UU;


end

