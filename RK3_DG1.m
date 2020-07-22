function [ Um , Us ] = RK3_DG1( dt, Um, Us, Fr)
% time discretisation for DG1 : RK3 of Shu and Osher (1989)

% input:
% > dt - time step
% > Um,s - means and slopes

[ RHSm , RHSs ] = space_DG1_RHS( Um , Us , Fr );

Um1 = Um + dt*RHSm;
Us1 = Us + dt*RHSs;

[ RHSm , RHSs ] = space_DG1_RHS( Um1 , Us1 , Fr );

Um2 = 0.25*(3*Um + Um1 + dt*RHSm);
Us2 = 0.25*(3*Us + Us1 + dt*RHSs);

[ RHSm , RHSs ] = space_DG1_RHS( Um2 , Us2 , Fr );

UUm = (1/3)*(Um + 2*Um2 + 2*dt*RHSm);
UUs = (1/3)*(Us + 2*Us2 + 2*dt*RHSs);


Um = UUm;
Us = UUs;

end

