% - the function that approximates the vectors with entries between 0 and 1 
% by a binary vector with a number of non-zero entries smaller or equal (or
% just equal) to prescribed number of entries
% this function is used to compute the sensor nodes
% - input parameters: 
%                   - computed_state_sequence    - state sequence simulated
%                   from the estimated initial state 
%                   - output_sequence - output sequence used in estimation
%                   - no_observed_nodes   - total number of observed nodes,
%                   that is, the total number of non-zero entries
% - output parameters:
%                   - solution         - binary vector that approximates
%                     "initial_solution"
% Author: Aleksandar Haber
% December 2019 - May 2020

function [solution]=solve_problem_binary_state(computed_state_sequence,output_sequence, no_observed_nodes)
[N,time_horizon]=size(computed_state_sequence);

output=output_sequence(:);
%state=computed_state_sequence;

% min f'*z    subject to:  A*z  <= b
%             z                        Aeq*z = beq
%                                      lb <= z <= ub
%                                      z(i) integer, where i is in the index
%                                      vector intcon (integer constraints)


f=[ones(1,N*time_horizon),zeros(1,N)];
% defining integer variables
intcon=[(N*time_horizon+1):(N*time_horizon+N)];

% lower and upper bounds
lb=zeros(N*time_horizon+N,1);
ub=[Inf*ones(N*time_horizon,1); ones(N,1)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first option 
% the number of binary variables is smaller or equal than no_observed_nodes
% comment this option if you are not using it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A1=zeros(2*(N*time_horizon),N*time_horizon);
% A2=[];
% B1=[-1;-1];
% b=[];
% 
% for i=1:(N*time_horizon)
%    A1(2*(i-1)+1: 2*i,   i)  = B1;
%    b(2*(i-1)+1,1) =   -output(i,1);
%    b(2*i,1)       =    output(i,1);
% end
% 
% for i=1:time_horizon 
%    A2=[A2;kron(diag(computed_state_sequence(:,i)),[-1;1])];
% end
% 
% A=[A1,A2;zeros(1,N*time_horizon),ones(1,N)];
% b=[b;no_observed_nodes];
% 
% Aeq=[]
% beq=[]


% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % second option 
% % % the number of binary variables is equal to the number of observed nodes
% % % comment this if you are not using it
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A1=zeros(2*(N*time_horizon),N*time_horizon);
A2=[];
B1=[-1;-1];
b=[];

for i=1:(N*time_horizon)
   A1(2*(i-1)+1: 2*i,   i)  = B1;
   b(2*(i-1)+1,1) =   -output(i,1);
   b(2*i,1)       =    output(i,1);
end

for i=1:time_horizon 
   A2=[A2;kron(diag(computed_state_sequence(:,i)),[-1;1])];
end

A=[A1,A2];
b=[b];

Aeq=[zeros(1,N*time_horizon),ones(1,N)];
beq=[no_observed_nodes];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                   end of the second option
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options = optimoptions(@intlinprog,'IntegerTolerance',1e-06,'RelativeGapTolerance',1.0000e-06)
% here we solve the optimization problem
solution = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub,options);
solution=solution(N*time_horizon+1:N*time_horizon+N);
end