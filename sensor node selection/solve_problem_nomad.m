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
% - output parameters:
%                   - solution - solution of the optimization problem
% Author: Aleksandar Haber 
% December 2019 - May 2020
% December 2019 - February 2020
function [solution,fval,exitflag,info]=solve_problem_nomad(state_sequence,initial_guess,no_observed_nodes,h,fcnHandle)

% we need to retrieve the time_horizon - it is equal to the
% observation_horizon+2
% N is the total number of nodes
[N,time_horizon]=size(state_sequence);
C=eye(N,N); % this matrix is used to define the output sequence for estimation


nB = N; %Number of Binary Variables
nC = N; %Number of Continuous Variables

xtype = [repmat('C',1,nC),repmat('B',1,nB)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first option with equality constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tol=1e-1; % this is a tolerance for relaxing the equality constraints
% the equality constraints are posed as double inequality constraints
% this is a standard practice for relaxing the equality constraints
% without the tolerance the solver returns infeasible solutions...
% A=[ zeros(1,nC), ones(1,nB) ; zeros(1,nC), -ones(1,nB) ]; b=[no_observed_nodes+tol; - no_observed_nodes+tol];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% second option with inequality constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=[ zeros(1,nC), ones(1,nB)]; b=[no_observed_nodes];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% lower and upper bounds
lb=[-5*ones(N,1); zeros(N,1)];
ub=[5*ones(N,1);   ones(N,1)];

% Options for the mixed integer optimization problem
% to get the complete set of options type: optiset('solver','nomad')
opts = optiset('solver','nomad','display','iter','maxfeval',100000,'maxiter', 50000,'maxtime',100000,'maxnodes',100000)
Opt = opti('fun',@cost_function_dynamic,'ineq', A,b,'bounds',lb,ub,'xtype',xtype,'options',opts);

%Solve the mixed integer optimization problem
[solution,fval,exitflag,info] = solve(Opt,initial_guess);


function cost_function_value = cost_function_dynamic(z)
 % z - variable 
% this variable contains the computed states    
STATE=zeros(N,time_horizon);

    for o=1:time_horizon-1
        if o==1
           %X=z(1:N,1)+h*(fcnHandle(z(1:N,1)));
           STATE(:,o)=z(1:N,1);
           STATE(:,o+1)=z(1:N,1)+h*fcnHandle(z(1:N,1));
           OUTPUT(:,o)=C*state_sequence(:,o);
           OUTPUT(:,o+1)=C*state_sequence(:,o+1);
        
        else 
           STATE(:,o+1)=STATE(:,o)+h*fcnHandle(STATE(:,o));
           OUTPUT(:,o+1)=C*state_sequence(:,o+1);
        end
    end

    kroneckerCopt=kron(eye(time_horizon,time_horizon),diag(z(N+1:end,1)));
    E=OUTPUT(:)-kroneckerCopt*STATE(:);
    cost_function_value=norm(E,2)^2;
end
end