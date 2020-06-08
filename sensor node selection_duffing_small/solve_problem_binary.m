% - the function that approximates the vectors with entries between 0 and 1 
% by a binary vector with a number of non-zero entries smaller or equal (or
% just equal) to prescribed number of entries
% this function is used to compute the sensor nodes
% - input parameters: 
%                   - initial_solution    - initial vector of selected sensor nodes obtained by solving the relaxed solution 
%                   - no_observed_nodes   - total number of observed nodes,
%                   that is, the total number of non-zero entries
%                   
%                   
% - output parameters:
%                   - solution         - binary vector that approximates
%                    
% Author: Aleksandar Haber
% December 2019 - May 2020

function [solution]=solve_problem_binary(initial_solution,no_observed_nodes)
[N,~]=size(initial_solution);


% min f'*z    subject to:  A*z  <= b
%             z                        Aeq*z = beq
%                                      lb <= z <= ub
%                                      z(i) integer, where i is in the index
%                                      vector intcon (integer constraints)


f=[ones(1,N),zeros(1,N)]
% defining integer variables
intcon=[N+1:2*N]

% lower and upper bounds
lb=zeros(2*N,1);
ub=[Inf*ones(N,1); ones(N,1)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first option 
% the number of binary variables is smaller or equal than no_observed_nodes
% comment this option if you are not using it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A=zeros(2*N,2*N);
% B1=[-1;-1];
% B2=[1;-1];
% 
% 
% for i=1:N
%   
%    A(2*(i-1)+1: 2*i,   i)  = B1;
%    A(2*(i-1)+1: 2*i, i+N)  = B2;
%    b(2*(i-1)+1,1) =   initial_solution(i);
%    b(2*i,1)       =  -initial_solution(i);
% 
% end
% 
% b=[b;no_observed_nodes];
% A=[A;zeros(1,N),ones(1,N)];
% 
% Aeq=[]
% beq=[]


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % second option 
% % the number of binary variables is equal to the number of observed nodes
% % comment this if you are not using it
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A=zeros(2*N,2*N);
B1=[-1;-1];
B2=[1;-1];

for i=1:N
   A(2*(i-1)+1: 2*i,   i)  = B1;
   A(2*(i-1)+1: 2*i, i+N)  = B2;
   b(2*(i-1)+1,1) =   initial_solution(i);
   b(2*i,1)       =  -initial_solution(i);
end

Aeq=[zeros(1,N),ones(1,N)];
beq=[no_observed_nodes];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                   end of the second option
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% here we solve the optimization problem
solution = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub) 

end