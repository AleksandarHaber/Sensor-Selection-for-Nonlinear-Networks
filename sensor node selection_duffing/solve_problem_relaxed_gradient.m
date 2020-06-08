% - function that solves the relaxed mixed integer nonlinear optimization
% problem that is used for sensor node selection
% this function solves a nonlinear optimization problem
%
% - input parameters: 
%                   - state_sequence   - simulated state sequence that
%                   is used to generate the output data for solving the
%                   problem
%                   - initial_guess    - initial guess for solving the
%                   optimization problem, first 2*N entries correspond to the
%                   guess of the initial state and the last N entries
%                   correspond to the initial guess of variables
%                   corresponding to the sensor nodes
%                   - no_observed_nodes - number of desired sensor nodes
%                   - h - discretization constant
%                   - fcnHandle - function handle that describes the system
%                   dynamics
%                   - fcnHandleGradient - function handle that describes the gradient of the system
%                   dynamics
% - output parameters:
%                   - solution - solution of the optimization problem
% Author: Aleksandar Haber 
% December 2019 - May 2020

function [solution]=solve_problem_relaxed_gradient(state_sequence,initial_guess,no_observed_nodes,h,fcnHandle,fcnHandleGradient)
% we need to retrieve the time_horizon - it is equal to the
% observation_horizon+2
% n is the total state dimension
% N is the number of nodes N=n/2
[n,time_horizon]=size(state_sequence);
N=n/2; 
C=zeros(N,n); % this matrix is used to define the output sequence for estimation
C1=[1 0];

for p=1:N
    C(p,2*(p-1)+1:2*p)=C1;
end

% identity matrix used later in the code
In=eye(n,n);

% summary of the optimization problem-general form
% fmincon attempts to solve problems of the form:
%     min F(z)  subject to:  A*z  <= B, Aeq*z  = Beq (linear constraints)
%      z                     C(z) <= 0, Ceq(z) = 0   (nonlinear constraints)
%                               lb <= z <= ub        (bounds)



% there are three possibilities to define the constraints
% I realized that the third option produces the best results in later steps
% Namely, the solution of the relaxed problem is used in the later steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first option - computed number of nodes is smaller or equal than "no_observed_nodes" 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A=[zeros(1,n), ones(1,N)];
% B=no_observed_nodes;
% Aeq=[]
% Beq=[]
% lb=[0*ones(n,1); zeros(N,1)];
% ub=[1*ones(n,1);   ones(N,1)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% second option - % computed number of sensor
% nodes is EQUAL to "no_observed_nodes" 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% A=[]
% B=[]
% Aeq=[zeros(1,n), ones(1,N)];
% Beq=no_observed_nodes;
% lb=[0*ones(n,1); zeros(N,1)];
% ub=[1*ones(n,1);   ones(N,1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% third option - % there are no constraints on the number of computed
% control nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A=[];
B=[];
Aeq=[];
Beq=[];
lb=[0*ones(n,1); zeros(N,1)];
ub=[1*ones(n,1);   ones(N,1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% keyboard
% algorithm can be 'interior-point', 'SQP','active set', and 'trust region reflective'
% set the options
options = optimoptions(@fmincon,'Algorithm','interior-point','SpecifyObjectiveGradient',true,'Display','iter','MaxFunctionEvaluations',200000,'UseParallel',true,'OptimalityTolerance',1.0000e-12, 'MaxIterations',10000,'ConstraintTolerance',1.0000e-8,'StepTolerance',1.0000e-8)
[solution,fval,exitflag,output] = fmincon(@(z)cost_function(z),initial_guess,A,B,Aeq,Beq,lb,ub,[],options);

% this is the cost function
function [cost_function_value, cost_function_gradient] = cost_function(z)
%function [cost_function_value] = cost_function(z)
        
    % this variable contains the computed final states    
    STATE=zeros(n,time_horizon);
    % this matrix is used to store gradients
    GR_ITER=zeros(n,time_horizon*n);
      
    % optimization options for fsolve, 'trust-region-dogleg',
    % 'trust-region', 'levenberg-marquardt'
    options_fsolve = optimoptions('fsolve','Algorithm',  'trust-region-dogleg','Display','off','SpecifyObjectiveGradient',true,'UseParallel',true,'FunctionTolerance',1.0000e-8,'MaxIter',10000,'StepTolerance', 1.0000e-8);
    problem.options = options_fsolve;
    problem.objective = @objective_fun;
    problem.solver = 'fsolve';
    
    for o=1:time_horizon-1
        if o==1
           STATE(:,o)=z(1:n,1);
           tmp0=fcnHandle(STATE(:,o)); 
           problem.x0 = STATE(:,o)+h*tmp0;  % use the Forward Euler method to generate the initial guess
           STATE(:,o+1)=fsolve(problem);
           
           % compute the outputs
           OUTPUT(:,o)=C*state_sequence(:,o);
           OUTPUT(:,o+1)=C*state_sequence(:,o+1);
        
           % compute the gradients
           GR_ITER(:,(o-1)*n+1:o*n)=In;
           tmp0gr=fcnHandleGradient(STATE(:,o));
           tmp1gr=fcnHandleGradient(STATE(:,o+1));
           GR_ITER(:,(o)*n+1:(o+1)*n)=inv(In-0.5*h*tmp1gr)*GR_ITER(:,(o-1)*n+1:o*n)*(In+0.5*h*tmp0gr);
            
        else 
           
           tmp0=fcnHandle(STATE(:,o)); 
           problem.x0 = STATE(:,o)+h*tmp0;  % use the Forward Euler method to generate the initial guess
           STATE(:,o+1)=fsolve(problem);
           
           % compute the outputs
           OUTPUT(:,o+1)=C*state_sequence(:,o+1);
           
           % compute the gradients
           tmp0gr=fcnHandleGradient(STATE(:,o));
           tmp1gr=fcnHandleGradient(STATE(:,o+1));
           GR_ITER(:,(o)*n+1:(o+1)*n)=inv(In-0.5*h*tmp1gr)*GR_ITER(:,(o-1)*n+1:o*n)*(In+0.5*h*tmp0gr);
           
        end
     end
    
    % form the parametrized C matrix - this matrix is used to define the output sequence for estimation
    Cz=zeros(N,n); 
    
    for p=1:N
        Cz(p,2*(p-1)+1:2*p)=[z(n+p,1) 0];
    end
     
    kroneckerCopt=kron(eye(time_horizon,time_horizon),Cz);
    E=OUTPUT(:)-kroneckerCopt*STATE(:);
    cost_function_value=norm(E,2)^2;
    
    
    
    %compute the gradient
    
    G11=[];
    for o=1:time_horizon
        G11=[G11,-1*diag(C*STATE(:,o))];
    end
    
% another method for computing G11 - you can uncomment this to double
% check the computations
% optimization variables
%   G11_second=zeros(N,2*N);
%     
%   for l=1:time_horizon
%        for g=1:N
%             G11_second(g,(l-1)*2*N+(g-1)*2+1:(l-1)*2*N+g*2)=STATE(2*(g-1)+1:2*g,l)';
%        end
%     end
%     
%     keyboard
%     
%     kroneckerCopt_2=kron(eye(time_horizon,time_horizon),C);
%     G11_second=-G11_second*kroneckerCopt_2';
    
       
    G22=-1*GR_ITER*kroneckerCopt';
    % this is the final gradient
    cost_function_gradient=2*[G22; G11]*E;
   
         
% function for solving nonlinear system of equations    
function [f, jacobianz]=objective_fun(xk)
    tmp1=0.5*h*(fcnHandle(xk)+tmp0);% tmp0=fcnHandle(xk-1) is computed outside this function to speed up the computations (in the function above)
    f=xk-STATE(:,o)-tmp1;
    jacobianz=In-0.5*h*fcnHandleGradient(xk)'; % here we need to transpose since the function "fcnHandleGradient(z)" returns the gradient, and the Jacobian is its transpose
end

    
    
end
end