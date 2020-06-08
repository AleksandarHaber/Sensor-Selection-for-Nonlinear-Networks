% - function that estimates the initial state for a given matrix Cmatrix
% - input parameters: 
%                   - output_sequence    - output sequence for estimation
%                   - initial_guess_state            - guess of the initial
%                   state
%                   - Cmatrix             - output matrix for estimation
%                   - h                   - discretization constant                  
%                   - fcnHandle - function handle that describes the system
%                   dynamics
%                   - fcnHandleGradient - function hangle that describes
%                   the gradient of the system dynamics
% - output parameters:
%                   - solution         - estimate of the initial state
% Author: Aleksandar Haber
% December 2019 - February 2020

function [solution]=estimate_initial_state_gradient(output_sequence,initial_guess_state,Cmatrix,h,fcnHandle,fcnHandleGradient)

% M is the number of sensor nodes
[M,time_horizon]=size(output_sequence);
[~,N]=size(Cmatrix);

%keyboard
% algorithm can be 'trust-region' or 'quasi-newton'.
options = optimoptions(@fminunc,'Algorithm','quasi-newton','SpecifyObjectiveGradient',true,'Display','iter','MaxFunctionEvaluations',30000,'UseParallel',true,'OptimalityTolerance',1.0000e-6)
%options = optimoptions(@fminunc,'Algorithm','quasi-newton','MaxFunctionEvaluations',30000,'UseParallel',true,'OptimalityTolerance',1.0000e-8)

%keyboard
[solution,fval,exitflag,output] = fminunc(@(z)cost_function(z),initial_guess_state,options);

% cost function
function [cost_function_value,cost_function_gradient] = cost_function(z)
% this variable contains the computed final states    
STATE=zeros(N,time_horizon);
    for o=1:time_horizon-1
        if o==1
           STATE(:,o)=z(1:N,1);
           STATE(:,o+1)=z(1:N,1)+h*fcnHandle(z(1:N,1));
           
           GR_ITER(:,(o-1)*N+1:o*N)=eye(N,N);
           GR_ITER(:,(o)*N+1:(o+1)*N)=GR_ITER(:,(o-1)*N+1:o*N)+h*GR_ITER(:,(o-1)*N+1:o*N)*fcnHandleGradient(z(1:N,1));
           
        else 
           STATE(:,o+1)=STATE(:,o)+h*fcnHandle(STATE(:,o)); 
           GR_ITER(:,(o)*N+1:(o+1)*N)=GR_ITER(:,(o-1)*N+1:o*N)+h*GR_ITER(:,(o-1)*N+1:o*N)*fcnHandleGradient(STATE(:,o));
        end
    end
        
    
    
    kroneckerCopt=kron(eye(time_horizon,time_horizon),Cmatrix);
    E=output_sequence(:)-kroneckerCopt*STATE(:);
    cost_function_value=norm(E,2)^2;
    
    %now compute the gradient
   
    
    G22=-1*GR_ITER*kroneckerCopt';
    cost_function_gradient=2*[G22]*E;
    
    
    
    
    
    
end
end