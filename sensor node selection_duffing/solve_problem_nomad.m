% Function that selects the sensor nodes using the NOMAD solver

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

function [solution,fval,exitflag,info]=solve_problem_nomad(state_sequence,initial_guess,no_observed_nodes,h,fcnHandle,fcnHandleGradient)

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

nC = 2*N; %Number of Continuous Variables
nB = N;   %Number of Binary Variables

xtype = [repmat('C',1,nC),repmat('B',1,nB)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first option with equality constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tol=5e-1; % this is a tolerance for relaxing the equality constraints
% % the equality constraints are posed as double inequality constraints
% % this is a standard practice for relaxing the equality constraints
% % without the tolerance the solver can return infeasible solutions...
% A=[ zeros(1,nC), ones(1,nB) ; zeros(1,nC), -ones(1,nB) ]; b=[no_observed_nodes+tol; - no_observed_nodes+tol];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% second option with inequality constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A=[ zeros(1,nC), ones(1,nB)]; b=[no_observed_nodes];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% lower and upper bounds
lb=[0*ones(2*N,1); zeros(N,1)];
ub=[1*ones(2*N,1);   ones(N,1)];

% Options for the mixed integer optimization problem
% to get the complete set of options type: optiset('solver','nomad')
opts = optiset('solver','nomad','display','iter','maxfeval',100000,'maxiter', 50000,'maxtime',10000,'maxnodes',100000)
Opt = opti('fun',@cost_function_dynamic,'ineq', A,b,'bounds',lb,ub,'xtype',xtype,'options',opts);

%Solve the mixed integer optimization problem
[solution,fval,exitflag,info] = solve(Opt,initial_guess);


function cost_function_value = cost_function_dynamic(z)
   % this variable contains the computed final states    
    STATE=zeros(n,time_horizon);
         
    % optimization options for fsolve, 'trust-region-dogleg', 'trust-region'
    options_fsolve = optimoptions('fsolve','Algorithm', 'trust-region-dogleg','Display','off','SpecifyObjectiveGradient',true,'UseParallel',true,'FunctionTolerance',1.0000e-8,'MaxIter',10000,'StepTolerance', 1.0000e-8);
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
                  
        else 
           
           tmp0=fcnHandle(STATE(:,o)); 
           problem.x0 = STATE(:,o)+h*tmp0;  % use the Forward Euler method to generate the initial guess
           STATE(:,o+1)=fsolve(problem);
           
           % compute the outputs
           OUTPUT(:,o+1)=C*state_sequence(:,o+1);
           
           
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
        
         
% function for solving nonlinear system of equations    
function [f, jacobianz]=objective_fun(xk)
    tmp1=0.5*h*(fcnHandle(xk)+tmp0);% tmp0=fcnHandle(xk-1) is computed outside this function to speed up the computations (in the function above)
    f=xk-STATE(:,o)-tmp1;
    jacobianz=In-0.5*h*fcnHandleGradient(xk)'; % here we need to transpose since the function "fcnHandleGradient(z)" returns the gradient, and the Jacobian is its transpose
end

end
end