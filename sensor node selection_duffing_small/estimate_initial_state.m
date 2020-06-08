% - function that estimates the initial state for a given matrix Cmatrix
% - input parameters: 
%                   - output_sequence    - output sequence for estimation
%                   - initial_guess_state            - guess of the initial
%                   state
%                   - Cmatrix             - output matrix for estimation
%                   - h                   - discretization constant                  
%                   - fcnHandle - function handle that describes the system
%                   dynamics
% - output parameters:
%                   - solution         - estimate of the initial state
% Author: Aleksandar Haber
% December 2019 - February 2020

function [solution]=estimate_initial_state(output_sequence,initial_guess_state,Cmatrix,h,fcnHandle)

% M is the number of sensor nodes
[M,time_horizon]=size(output_sequence);
[~,N]=size(Cmatrix);

%keyboard
% algorithm can be "trust-region-reflective" or "interior-point" or "sqp"
options = optimoptions(@fminunc,'Algorithm','quasi-newton','Display','iter','MaxFunctionEvaluations',30000,'UseParallel',true,'OptimalityTolerance',1.0000e-8)
%options = optimoptions(@fminunc,'Algorithm','quasi-newton','MaxFunctionEvaluations',30000,'UseParallel',true,'OptimalityTolerance',1.0000e-8)

%keyboard
[solution,fval,exitflag,output] = fminunc(@(z)cost_function(z),initial_guess_state,options);

% cost function
function cost_function_value = cost_function(z)
% this variable contains the computed final states    
STATE=zeros(N,time_horizon);
    for o=1:time_horizon-1
        if o==1
           STATE(:,o)=z(1:N,1);
           STATE(:,o+1)=z(1:N,1)+h*fcnHandle(z(1:N,1));
        else 
           STATE(:,o+1)=STATE(:,o)+h*fcnHandle(STATE(:,o)); 
        end
    end
    kroneckerCopt=kron(eye(time_horizon,time_horizon),Cmatrix);
    cost_function_value=norm(output_sequence(:)-kroneckerCopt*STATE(:),2)^2;
end
end