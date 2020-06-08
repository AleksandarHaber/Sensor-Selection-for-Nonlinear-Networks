%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - File that performs sensor node selection (and initial state estimation) for the Duffing network whose
% continuous-time dynamics is:
% \dot{x}=f(x) 
%      y =Cx   
% 
% where x is the state, f(x) is the vector fuction describing an
% uncontrolled dynamics, C is the output matrix, and y is the measurement
% vector
%
% - The discretization of the dynamics is performed using the Trapezoidal
% Implicit (TI) method
%
% - Before running this file, run  "generate_dynamics_duffing.m" file to
% generate two files:
%
%    1.) "duffing_network_dynamics.m" - function that computes f(x) for a given x
%    2.) "duffing_network_dynamics_gradient" - function that computes
%    \nabla f(x) (gradient) for a given x
%
% - Author: Aleksandar Haber
% December 2019 - May 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, pack, clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   parameter selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N is the number of subsystems, every subsystem is of the second order
% the precise value of N should be adjusted such that it matches the number
% of subsystems in "duffing_network_dynamics.m" and "duffing_network_dynamics_gradient"
N=10
% discretization constant 
h=0.0001

% number of time steps for sensor node selection - observation
% horizon-total observation horizon is observation_horizon+2
observation_horizon=100; 
time=0:h:observation_horizon*h;

% number of sensor nodes
no_sensor_nodes=4;

% select an initial state that we want to estimate
initial_state_true=rand(2*N,1);


rng('shuffle') 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   end of parameter selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   simulate the state sequence for estimation and sensor node selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % simulate the uncontrolled dynamics using the MATLAB built-in solver, this
% % is used to test the accuracy of the Trapezoidal Implicit (TI) method 
[time_tmp,STATE_ode45] = simulate_uncontrolled_ode45(time,initial_state_true,@duffing_network_dynamics);

% simulate the uncontrolled dynamics using the TI method- this state sequence is used for
% estimation
STATE_TI=simulate_uncontrolled_ti_fsolve_3(length(time),initial_state_true,h,@duffing_network_dynamics,@duffing_network_dynamics_gradient);

figure(1)
plot(STATE_TI(:,:)')
 
% % compare the TI method with ode45 simulation
STATE_ode45=STATE_ode45';
for i=1:length(time)
   error_simulation(i)= norm(STATE_TI(:,i)-STATE_ode45(:,i),2)/norm(STATE_ode45(:,i),2);
end
 
% figure(2)
% hold on
% plot(error_simulation,'k');
% figure(3)
% plot(STATE_ode45(2,:),'k')
% hold on 
% plot(STATE_TI(2,:),'m')
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               end of computation of state sequences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve the relaxed problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solving the relaxed mixed integer nonlinear optimization problem

% generate an initial guess
% first 2*N entries correspond to the initial state guess and the last N entries
% correspond to the sensor node selections

initial_guess_initial_state=rand(2*N,1);
initial_guess=[initial_guess_initial_state; rand(N,1)];

tic
% compute the relaxed problem
solution_relaxed = solve_problem_relaxed_gradient(STATE_TI,initial_guess,no_sensor_nodes,h,@duffing_network_dynamics,@duffing_network_dynamics_gradient)
time_solution_relaxed=toc;

% this is the initial state estimate
solution_relaxed(1:2*N);
% these variables correspond to the selected sensor nodes- they have to be
% transferred to other function below.
solution_relaxed(2*N+1:end)

% compute the errors
solution_optimal_error_relaxed=norm(solution_relaxed(1:2*N)-initial_state_true)/norm(initial_state_true)
initial_guess_error=norm(initial_guess_initial_state-initial_state_true)/norm(initial_state_true)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Option 1: select the sensor nodes using the NOMAD solver 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve the problem using the NOMAD solver
% NOTE THAT FOR THE NOMAD SOLVER, THE BEST RESULTS ARE OBTAINED BY SETTING
% THE SENSOR NODES TO BE SMALLER OR EQUAL THAN THE DESIRED NUMBER OF NODES
% IF YOU SET AN INEQUALITY, YOU MIGHT GET AN INFEASIBLE SOLUTION
tic
[solution_nomad,fval,exitflag,info]=solve_problem_nomad(STATE_TI,solution_relaxed,no_sensor_nodes,h,@duffing_network_dynamics,@duffing_network_dynamics_gradient)
time_nomad=toc

% form the optimal C matrix
selected_nodes_1=solution_nomad(2*N+1:end);
real_number_sensor_nodes_1= nnz(selected_nodes_1>10^(-1)) ;
Coptimal_1=zeros(real_number_sensor_nodes_1,2*N);
indx=1;
for i=1:N
   if(selected_nodes_1(i)>0)
       Coptimal_1(indx,2*(i-1)+1:2*i)=[1 0];
       indx=indx+1;
   end
end

% estimate the initial state using the selected Coptimal matrix
% form the output sequence
output_sequence_1=Coptimal_1*STATE_TI;
[solution_final_1]=estimate_initial_state_gradient(output_sequence_1,initial_guess_initial_state,Coptimal_1,h,@duffing_network_dynamics,@duffing_network_dynamics_gradient);

figure(4)
plot(solution_final_1)
hold on 
plot(initial_state_true,'r')
% this is the final estimation error
solution_optimal_error_1=norm(solution_final_1-initial_state_true)/norm(initial_state_true)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of the NOMAD solver selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Option 2: Determine the optimal selection of sensor nodes by
% approximating the solution of the relaxed problem by a binary vector- in
% this case, the estimated state vector is NOT used. 
tic
[solution_binary_1]=solve_problem_binary(solution_relaxed(2*N+1:end),no_sensor_nodes)
time_option2=toc
% check the result 
plot(solution_relaxed(2*N+1:end),'r')
hold on 
plot(solution_binary_1(N+1:2*N))

% estimate the initial state on the basis of the selected sensor nodes
% form the optimal C matrix
selected_nodes_2=solution_binary_1(N+1:2*N);
real_number_sensor_nodes_2= nnz(selected_nodes_2>10^(-1)) ;
Coptimal_2=zeros(real_number_sensor_nodes_2,2*N);

indx=1;
for i=1:N
   if(selected_nodes_2(i)>0)
       Coptimal_2(indx,2*(i-1)+1:2*i)=[1 0];
       indx=indx+1;
   end
end

% estimate the initial state
% form the output sequence
output_sequence_2=Coptimal_2*STATE_TI;
% initial guess of the initial state - this is the same for both options

[solution_final_2]=estimate_initial_state_gradient(output_sequence_2,initial_guess_initial_state,Coptimal_2,h,@duffing_network_dynamics,@duffing_network_dynamics_gradient)

figure(5)
plot(solution_final_2)
hold on 
plot(initial_state_true,'r')
% this is the final estimation error

solution_optimal_error_2=norm(solution_final_2-initial_state_true)/norm(initial_state_true)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Option 3: Determine the optimal selection of sensor nodes by
% approximating the solution of the relaxed problem by a binary vector- in
% this case the estimated state vector IS used. 

% simulate the state sequence from the state estimate that is obtained as
% the solution of the relaxed problem
STATE_fE_simulated_from_estimated_x0=simulate_uncontrolled_ti_fsolve_3(length(time),solution_relaxed(1:2*N),h,@duffing_network_dynamics,@duffing_network_dynamics_gradient);

% this matrix is used to define the output sequence for estimation
Call=zeros(N,2*N); 
C1=[1 0];
for p=1:N
    Call(p,2*(p-1)+1:2*p)=C1;
end

% compute the optimal location of sensor nodes
tic
solution_binary_2=solve_problem_binary_state(Call*STATE_fE_simulated_from_estimated_x0,Call*STATE_TI,no_sensor_nodes);
time_option3=toc
% form the optimal C matrix
selected_nodes_3=solution_binary_2;
real_number_sensor_nodes_3= nnz(solution_binary_2) ;
Coptimal_3=zeros(real_number_sensor_nodes_3,2*N);
indx=1;
for i=1:N
   if(selected_nodes_3(i)>0)
       Coptimal_3(indx,2*(i-1)+1:2*i)=[1 0];
       indx=indx+1;
   end
end


% estimate the initial state
% form the output sequence
output_sequence_3=Coptimal_3*STATE_TI;
% initial guess of the initial state - this is the same for both options
[solution_final_3]=estimate_initial_state_gradient(output_sequence_3,initial_guess_initial_state,Coptimal_3,h,@duffing_network_dynamics,@duffing_network_dynamics_gradient)

figure(6)
plot(solution_final_3)
hold on 
plot(initial_state_true,'r')
% this is the final estimation error
solution_optimal_error_3=norm(solution_final_3-initial_state_true)/norm(initial_state_true)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Option 4: Determine the optimal selection of sensor nodes by
% approximating the solution of the relaxed problem by a binary vector- in
% this case, the estimated state vector is NOT used. 
% This solution is based on min-max formulation
tic
[solution_binary_3]=solve_problem_binary_max(solution_relaxed(2*N+1:end),no_sensor_nodes)
time_option4=toc

% check the result 
plot(solution_relaxed(2*N+1:end),'r')
hold on 
plot(solution_binary_3(2:N+1))

% estimate the initial state on the basis of the selected sensor nodes
% form the optimal C matrix
selected_nodes_4=solution_binary_3(2:N+1);
real_number_sensor_nodes_4= nnz(selected_nodes_4>10^(-1)) ;
Coptimal_4=zeros(real_number_sensor_nodes_4,2*N);

indx=1;
for i=1:N
   if(selected_nodes_4(i)>0)
       Coptimal_4(indx,2*(i-1)+1:2*i)=[1 0];
       indx=indx+1;
   end
end

% estimate the initial state
% form the output sequence
output_sequence_4=Coptimal_4*STATE_TI;
% initial guess of the initial state - this is the same for both options

[solution_final_4]=estimate_initial_state_gradient(output_sequence_4,initial_guess_initial_state,Coptimal_4,h,@duffing_network_dynamics,@duffing_network_dynamics_gradient)

figure(7)
plot(solution_final_4)
hold on 
plot(initial_state_true,'r')
% this is the final estimation error

solution_optimal_error_4=norm(solution_final_4-initial_state_true)/norm(initial_state_true)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              test the method against random selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% number_of_random_experiments=100;
% 
% for i=1:number_of_random_experiments
% ss = sort(randsample(N,real_number_sensor_nodes_2));
% Crandom=zeros(real_number_sensor_nodes_2,2*N);
% for j=1:real_number_sensor_nodes_2
%     Crandom(j,(ss(j)-1)*2+1:ss(j)*2)=[1, 0];
% end
% % here uncomment if you want different initial guess in every iteration
% %initial_guess_initial_state=randn(N,1);
% output_sequence=Crandom*STATE_TI;
% [solution_random]=estimate_initial_state_gradient(output_sequence,initial_guess_initial_state,Crandom,h,@duffing_network_dynamics,@duffing_network_dynamics_gradient)
% error_random(i)=norm(solution_random-initial_state_true)/norm(initial_state_true)
% end
% 
% figure(7)
% hist(error_random)
% histogram(error_random, 10, 'Normalization','probability' )
% ytix = get(gca, 'YTick')
% set(gca, 'YTick',ytix, 'YTickLabel',ytix*100)
% hold on 
% line([solution_optimal_error_3 solution_optimal_error_3], [0 0.25],'color','r','linewidth',2);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is for smaller networks, generate all the possible selections of
% sensor nodes - exhaustive search 

v = 1:1:N;
% generate all the possible selections of "no_control_nodes" control nodes
% out of possible "N" control nodes
nodes_random_selection = nchoosek(v,real_number_sensor_nodes_2)


for i=1:numel(nodes_random_selection(:,1))
    ss = sort(nodes_random_selection(i,:))
    Crandom=zeros(real_number_sensor_nodes_2,2*N);
for j=1:real_number_sensor_nodes_2
    Crandom(j,(ss(j)-1)*2+1:ss(j)*2)=[1, 0];
end
% here uncomment if you want different initial guess in every iteration
%initial_guess_initial_state=randn(N,1);
output_sequence=Crandom*STATE_TI;
[solution_random]=estimate_initial_state_gradient(output_sequence,initial_guess_initial_state,Crandom,h,@duffing_network_dynamics,@duffing_network_dynamics_gradient)
error_random(i)=norm(solution_random-initial_state_true)/norm(initial_state_true)
end


figure(8)
%hist(error_random)
histogram(error_random, 20, 'Normalization','probability','FaceColor',[0.9290 0.6940 0.1250],'FaceAlpha',0.3 )
ytix = get(gca, 'YTick')
set(gca, 'YTick',ytix, 'YTickLabel',ytix*100)
hold on 
line([solution_optimal_error_1 solution_optimal_error_1], [0 0.25],'color','r','linewidth',3);
line([solution_optimal_error_3 solution_optimal_error_3], [0 0.25],'color','b','linewidth',3);
line([solution_optimal_error_4 solution_optimal_error_4], [0 0.25],'color','k','linewidth',3);

























