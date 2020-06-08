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

function [solution]=estimate_initial_state_gradient(output_sequence,initial_guess_state,Cmatrix,h,fcnHandle,fcnHandleGradient)

% M is the number of sensor nodes
[M,time_horizon]=size(output_sequence);
[~,n]=size(Cmatrix);

N=n/2;

In=eye(n,n);

%keyboard
% algorithm can be 'trust-region' or 'quasi-newton'.
options = optimoptions(@fminunc,'Algorithm','quasi-newton','SpecifyObjectiveGradient',true,'Display','iter','MaxFunctionEvaluations',30000,'UseParallel',true,'OptimalityTolerance',1.0000e-10)
%options = optimoptions(@fminunc,'Algorithm','quasi-newton','MaxFunctionEvaluations',30000,'UseParallel',true,'OptimalityTolerance',1.0000e-8)

%keyboard
[solution,fval,exitflag,output] = fminunc(@(z)cost_function(z),initial_guess_state,options);

% cost function
function [cost_function_value,cost_function_gradient] = cost_function(z)
    
    
    %function [cost_function_value] = cost_function(z)
        
    % this variable contains the computed final states    
      STATE=zeros(n,time_horizon);
      GR_ITER=zeros(n,time_horizon*n);
      
    % optimization options for fsolve
    options_fsolve = optimoptions('fsolve','Algorithm', 'trust-region-dogleg','Display','off','SpecifyObjectiveGradient',true,'UseParallel',true,'FunctionTolerance',1.0000e-10,'MaxIter',10000,'StepTolerance', 1.0000e-12);
    
    problem.options = options_fsolve;
    % this objective is for the nonlinear system of equations
    problem.objective = @objective_fun;
    problem.solver = 'fsolve';
    
     
    
     for o=1:time_horizon-1
        if o==1
           STATE(:,o)=z;
           tmp0=fcnHandle(STATE(:,o)); 
           problem.x0 = STATE(:,o)+h*tmp0;  % use the Forward Euler method to generate the initial guess
           STATE(:,o+1)=fsolve(problem);
                   
            GR_ITER(:,(o-1)*n+1:o*n)=In;
            
            tmp0gr=fcnHandleGradient(STATE(:,o));
            tmp1gr=fcnHandleGradient(STATE(:,o+1));
            
            GR_ITER(:,(o)*n+1:(o+1)*n)=inv(In-0.5*h*tmp1gr)*GR_ITER(:,(o-1)*n+1:o*n)*(In+0.5*h*tmp0gr);
            
        else 
           
           tmp0=fcnHandle(STATE(:,o)); 
           problem.x0 = STATE(:,o)+h*tmp0;
           STATE(:,o+1)=fsolve(problem);
            
           tmp0gr=fcnHandleGradient(STATE(:,o));
           tmp1gr=fcnHandleGradient(STATE(:,o+1));
            
           GR_ITER(:,(o)*n+1:(o+1)*n)=inv(In-0.5*h*tmp1gr)*GR_ITER(:,(o-1)*n+1:o*n)*(In+0.5*h*tmp0gr);
            
                    
        end
     end
    
    % form the C matrix
    
         
    kroneckerCopt=kron(eye(time_horizon,time_horizon),Cmatrix);
    E=output_sequence(:)-kroneckerCopt*STATE(:);
    cost_function_value=norm(E,2)^2;
    G22=-1*GR_ITER*kroneckerCopt';
    cost_function_gradient=2*[G22]*E;
      
function [f, jacobianz]=objective_fun(xk)
    tmp1=0.5*h*(fcnHandle(xk)+tmp0);% tmp0=fcnHandle(xk-1) is computed outside this function to speed up the computations (in the function above)
    f=xk-STATE(:,o)-tmp1;
    jacobianz=In-0.5*h*fcnHandleGradient(xk)'; % here we need to transpose since the function "fcnHandleGradient(z)" returns the gradient, and the Jacobian is its transpose
end
        
    
end
end