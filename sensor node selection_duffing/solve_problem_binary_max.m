% - the function that approximates the vectors with entries between 0 and 1 
% by a binary vector with a number of non-zero entries smaller or equal (or
% just equal) to prescribed number of entries
% - this formulation is based on the minimization on the maximal difference
% between two vectors
% this function is used to compute the sensor nodes
% - input parameters: 
%                   - initial_solution    - initial vector of selected sensor nodes obtained by solving the relaxed solution 
%                   - no_observed_nodes   - total number of observed nodes,
%                   that is, the total number of non-zero entries
%                   
%                   
% - output parameters:
%                   - solution         - binary vector that approximates
%                     "initial_solution"
% Author: Aleksandar Haber
% December 2019 - May 2020

function [solution]=solve_problem_binary_max(initial_solution,no_observed_nodes)
[N,~]=size(initial_solution);


% min f'*z    subject to:  A*z  <= b
%             z                        Aeq*z = beq
%                                      lb <= z <= ub
%                                      z(i) integer, where i is in the index
%                                      vector intcon (integer constraints)


f=[1,zeros(1,N)];
% defining integer variables
intcon=[2:N+1];

% lower and upper bounds
lb=zeros(N+1,1);
ub=[Inf; ones(N,1)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first option 
% the number of binary variables is smaller or equal than no_observed_nodes
% comment this option if you are not using it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A=zeros(2*N,N);
% B1=[1;-1];
% 
% 
% for i=1:N
%   
%    A(2*(i-1)+1: 2*i,   i)  = B1;
%    b(2*(i-1)+1,1) =   initial_solution(i);
%    b(2*i,1)       =  -initial_solution(i);
% 
% end
% 
% b=[b;no_observed_nodes];
% A=[-1*ones(2*N,1),A; 0, ones(1,N)];
% 
% Aeq=[];
% beq=[];


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % second option 
% % the number of binary variables is equal to the number of observed nodes
% % comment this if you are not using it
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A=zeros(2*N,N);
B1=[1;-1];


for i=1:N
  
   A(2*(i-1)+1: 2*i,   i)  = B1;
   b(2*(i-1)+1,1) =   initial_solution(i);
   b(2*i,1)       =  -initial_solution(i);

end

b=[b;];
A=[-1*ones(2*N,1),A];

 Aeq=[0,ones(1,N)];
 beq=[no_observed_nodes];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                   end of the second option
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% here we solve the optimization problem
solution = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub) 

end