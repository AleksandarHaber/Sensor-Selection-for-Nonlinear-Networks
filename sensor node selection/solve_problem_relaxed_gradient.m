% - function that solves the relaxed mixed integer nonlinear optimization
% problem that is used for sensor node selection, that is, this function
% solves a nonlinear optimization problem
% - input parameters: 
%                   - state_sequence   - simulated state sequence that
%                   is used to generate the output data for solving the
%                   problem
%                   - initial_guess    - initial guess for solving the
%                   optimization problem, first N entries correspond to the
%                   guess of the initial state and the last N entries
%                   correspond to the initial guess of variables
%                   corresponding to the sensor nodes
%                   - no_observed_nodes - number of desired sensor nodes
%                   - h - discretization constant
%                   - fcnHandle - function handle that describes the system
%                   dynamics
%                   - fcnHandleGradient - function handle that described
%                   the gradient of the system dynamics
% - output parameters:
%                   - solution - solution of the optimization problem
% Author: Aleksandar Haber 
% December 2019 - May 2020

function [solution]=solve_problem_relaxed_gradient(state_sequence,initial_guess,no_observed_nodes,h,fcnHandle,fcnHandleGradient)
% we need to retrieve the time_horizon - it is equal to the
% observation_horizon+2
% N is the total number of nodes
[N,time_horizon]=size(state_sequence);
C=eye(N,N); % this matrix is used to define the output sequence for estimation


% fmincon attempts to solve problems of the form:
%     min F(z)  subject to:  A*z  <= B, Aeq*z  = Beq (linear constraints)
%      z                     C(z) <= 0, Ceq(z) = 0   (nonlinear constraints)
%                               lb <= z <= ub        (bounds)

% defining the constraint matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first option - we consider that the computed number of
% nodes is smaller or equal than "no_observed_nodes" 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A=[zeros(1,N), ones(1,N)];
% B=no_observed_nodes;
% Aeq=[]
% Beq=[]
% lb=[-10*ones(N,1); zeros(N,1)];
% ub=[10*ones(N,1);   ones(N,1)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% second option - % computed number of sensor
% nodes is EQUAL to "no_observed_nodes" 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% A=[]
% B=[]
% Aeq=[zeros(1,N), ones(1,N)];
% Beq=no_observed_nodes;
% lb=[-5*ones(N,1); zeros(N,1)];
% ub=[5*ones(N,1);   ones(N,1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% third option - % there are no constraints on the number of computed
% control nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A=[]
B=[]
Aeq=[];
Beq=[];
lb=[-5*ones(N,1); zeros(N,1)];
ub=[5*ones(N,1);   ones(N,1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%keyboard
% algorithm can be 'interior-point', 'SQP','active set', and 'trust region reflective'
% set the options
options = optimoptions(@fmincon,'Algorithm','interior-point','SpecifyObjectiveGradient',true,'Display','iter','MaxFunctionEvaluations',200000,'UseParallel',true,'OptimalityTolerance',1.0000e-16, 'StepTolerance',1.0000e-16)

[solution,fval,exitflag,output] = fmincon(@(z)cost_function(z),initial_guess,A,B,Aeq,Beq,lb,ub,[],options);

% this is the cost function
function [cost_function_value, cost_function_gradient] = cost_function(z)
% z - variable 
% this variable contains the computed states    
STATE=zeros(N,time_horizon);
GR_ITER=zeros(N,time_horizon*N);
    for o=1:time_horizon-1
        if o==1
           %X=z(1:N,1)+h*(fcnHandle(z(1:N,1)));
           STATE(:,o)=z(1:N,1);
           STATE(:,o+1)=z(1:N,1)+h*fcnHandle(z(1:N,1));
           OUTPUT(:,o)=C*state_sequence(:,o);
           OUTPUT(:,o+1)=C*state_sequence(:,o+1);
        
           GR_ITER(:,(o-1)*N+1:o*N)=eye(N,N);
           GR_ITER(:,(o)*N+1:(o+1)*N)=GR_ITER(:,(o-1)*N+1:o*N)+h*GR_ITER(:,(o-1)*N+1:o*N)*fcnHandleGradient(z(1:N,1));
        else 
           STATE(:,o+1)=STATE(:,o)+h*fcnHandle(STATE(:,o));
           GR_ITER(:,(o)*N+1:(o+1)*N)=GR_ITER(:,(o-1)*N+1:o*N)+h*GR_ITER(:,(o-1)*N+1:o*N)*fcnHandleGradient(STATE(:,o));
           OUTPUT(:,o+1)=C*state_sequence(:,o+1);
        end
    end

    kroneckerCopt=kron(eye(time_horizon,time_horizon),diag(z(N+1:end,1)));
    E=OUTPUT(:)-kroneckerCopt*STATE(:);
    cost_function_value=norm(E,2)^2;
    
    %now compute the gradient
    
    G11=[];
    for o=1:time_horizon
        G11=[G11,-1*diag(STATE(:,o))];
    end
    G22=-1*GR_ITER*kroneckerCopt';
    cost_function_gradient=2*[G22; G11]*E;
    
    
    
end
end