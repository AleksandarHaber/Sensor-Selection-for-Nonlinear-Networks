%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The file that performs sensor node selection (and an initial state estimation) for the memory network whose
% continuous-time dynamics is defined in the file:
% "memory_network_dynamics.m"
% Author: Aleksandar Haber
% December 2019 - May 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, pack, clc, rng('shuffle') 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   letter definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the letters
L=[1 -1 -1 -1 -1;
   1 -1 -1 -1 -1; 
   1 -1 -1 -1 -1;
   1 -1 -1 -1 -1;
   1  1  1  1  1;];

T=[ 1  1  1   1  1;  
   -1 -1  1  -1  -1;
   -1 -1  1  -1  -1;
   -1 -1  1  -1  -1;
   -1 -1  1  -1  -1;];
    
H=[ 1 -1 -1 -1  1 ;
    1 -1 -1 -1  1 ;
    1  1  1  1  1 ;
    1 -1 -1 -1  1 ;
    1 -1 -1 -1  1 ;
    ]
E=[1   1   1   1   1;
   1  -1  -1   -1 -1;
   1   1   1   -1 -1;
   1  -1  -1   -1 -1;
   1   1   1    1  1;];

F=[1   1   1    1   1;
   1  -1  -1   -1  -1;
   1   1   1   -1  -1;
   1  -1  -1   -1  -1;
   1  -1  -1   -1  -1;];

S=[1  1  1  1  1;
   1 -1 -1 -1 -1;
   1  1  1  1  1;
  -1  -1 -1 -1 1;
   1    1  1  1 1;];

O=[1   1   1   1  1;
   1  -1  -1  -1  1;
   1  -1  -1  -1  1;
   1  -1  -1  -1  1;
   1   1   1   1  1;];

Y=[1 -1 -1 -1  1;
   -1 1 -1  1 -1;
   -1 -1 1 -1 -1;
   -1 -1 1 -1 -1;
   -1 -1 1 -1 -1;];

X=[1 -1 -1 -1 1;
   -1 1 -1  1 -1;
   -1 -1 1  -1 -1;
   -1 1 -1  1 -1;
   1 -1 -1 -1  1;];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                end of the letter definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               select the main parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The size of the network is N=N1xN1 
N1=5; N=N1*N1;
% discretization constant for the forward Euler method
h=0.001

% number of time steps for sensor node selection - observation horizon
observation_horizon=20; 

% number of sensor nodes
no_sensor_nodes=22;


% the state to be estimated is a perturbed letter T
initial_state_true=T(:)+randn(N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               end of parameter definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  computation of the state trajectory for state estimation and comparison
%  with ode45()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first we compare ode45 with the forward Euler method
% simulate the dynamics using the ode45 function
time=0:h:observation_horizon*h;
[time_tmp,STATE_ode45]=simulate_uncontrolled_ode45(time,initial_state_true,@memory_network_dynamics);

% compare the forward Euler method simulation with ode45 simulation
% this state sequence is used for estimation
STATE_FE=simulate_uncontrolled_forward_Euler(length(time),initial_state_true,h,@memory_network_dynamics);

figure(1)
%plot(STATE_FE(:,:)')
plot(STATE_FE(1,:))
% compare the TI method with ode45 simulation
STATE_ode45=STATE_ode45';
for i=1:length(time)
   error_simulation(i)= norm(STATE_FE(:,i)-STATE_ode45(:,i),2)/norm(STATE_ode45(:,i),2);
end
 
% figure(2)
% hold on
% plot(error_simulation,'k');
% figure(3)
% plot(STATE_ode45(2,:),'k')
% hold on 
% plot(STATE_FE(2,:),'m')

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
initial_guess_initial_state=randn(N,1);
initial_guess=[initial_guess_initial_state; rand(N,1)];

solution_relaxed=solve_problem_relaxed_gradient(STATE_FE,initial_guess,no_sensor_nodes,h,@memory_network_dynamics,@memory_network_dynamics_gradient)
% this is the initial state estimate
solution_relaxed(1:N)
% these variables correspond to the selected sensor nodes- they have to be
% transferred to another function below.
solution_relaxed(N+1:end)

solution_optimal_error_relaxed=norm(solution_relaxed(1:N)-initial_state_true)
initial_guess_error=norm(initial_guess_initial_state-initial_state_true)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Option 1: select the sensor nodes using the NOMAD solver 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve the problem using the NOMAD solver
% NOTE THAT FOR THE NOMAD SOLVER, THE BEST RESULTS ARE OBTAINED BY SETTING
% THE SENSOR NODES TO BE SMALLER OR EQUAL THAN THE DESIRED NUMBER OF NODES
% IF YOU SET AN INEQUALITY, YOU MIGHT GET AN INFEASIBLE SOLUTION
[solution_nomad]=solve_problem_nomad(STATE_FE,solution_relaxed,no_sensor_nodes,h,@memory_network_dynamics)
selected_nodes_1=solution_nomad(N+1:2*N);
real_number_sensor_nodes_1= nnz(selected_nodes_1>10^(-1)) ;
Coptimal_1=zeros(real_number_sensor_nodes_1,N);

indx=1;
for i=1:N
   if(selected_nodes_1(i)>0)
       Coptimal_1(indx,i)=1;
       indx=indx+1;
   end
end

% estimate the initial state
% form the output sequence
output_sequence_1=Coptimal_1*STATE_FE;
% initial guess of the initial state - this is the same for both options

[solution_final_1]=estimate_initial_state_gradient(output_sequence_1,initial_guess_initial_state,Coptimal_1,h,@memory_network_dynamics,@memory_network_dynamics_gradient)

figure(2)
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

[solution_binary_1]=solve_problem_binary(solution_relaxed(N+1:end),no_sensor_nodes)



% check the result 
figure(3)
plot(solution_relaxed(N+1:end),'r')
hold on 
plot(solution_binary_1(N+1:2*N))

% form the optimal C matrix
selected_nodes_2=solution_binary_1(N+1:2*N);
real_number_sensor_nodes_2= nnz(selected_nodes_2>10^(-1)) ;
Coptimal_2=zeros(real_number_sensor_nodes_2,N)

indx=1;
for i=1:N
   if(selected_nodes_2(i)>0)
       Coptimal_2(indx,i)=1;
       indx=indx+1;
   end
end

% form the output sequence
output_sequence_2=Coptimal_2*STATE_FE;

% estimate the initial state on the basis of the selected sensor nodes
[solution_final_2]=estimate_initial_state_gradient(output_sequence_2,initial_guess_initial_state,Coptimal_2,h,@memory_network_dynamics,@memory_network_dynamics_gradient)


figure(3)
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

STATE_fE_simulated_from_estimated_x0=simulate_uncontrolled_forward_Euler(length(time),solution_relaxed(1:N),h,@memory_network_dynamics);
solution_binary_2=solve_problem_binary_state(STATE_fE_simulated_from_estimated_x0,eye(N,N)*STATE_FE,no_sensor_nodes);

% form the optimal C matrix
selected_nodes_3=solution_binary_2;
real_number_sensor_nodes_3= nnz(selected_nodes_3>10^(-1)) ;
Coptimal_3=zeros(real_number_sensor_nodes_3,N);

indx=1;
for i=1:N
   if(selected_nodes_3(i)>0)
       Coptimal_3(indx,i)=1;
       indx=indx+1;
   end
     
end

% form the output sequence
output_sequence_3=Coptimal_3*STATE_FE;
% initial guess of the initial state
[solution_final_3]=estimate_initial_state_gradient(output_sequence_3,initial_guess_initial_state,Coptimal_3,h,@memory_network_dynamics,@memory_network_dynamics_gradient)


figure(4)
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
% Here we approximate the binary vector on the basis of the min-max
% formulation

[solution_binary_3]=solve_problem_binary_max(solution_relaxed(N+1:end),no_sensor_nodes)



% check the result 
figure(5)
plot(solution_relaxed(N+1:end),'r')
hold on 
plot(solution_binary_3(2:N+1))

% form the optimal C matrix
selected_nodes_4=solution_binary_3(2:N+1);
real_number_sensor_nodes_4= nnz(selected_nodes_4>10^(-1)) ;
Coptimal_4=zeros(real_number_sensor_nodes_4,N);

indx=1;
for i=1:N
   if(selected_nodes_4(i)>0)
       Coptimal_4(indx,i)=1;
       indx=indx+1;
   end
end

% form the output sequence
output_sequence_4=Coptimal_4*STATE_FE;

% estimate the initial state on the basis of the selected sensor nodes
[solution_final_4]=estimate_initial_state_gradient(output_sequence_4,initial_guess_initial_state,Coptimal_4,h,@memory_network_dynamics,@memory_network_dynamics_gradient)


figure(6)
plot(solution_final_4)
hold on 
plot(initial_state_true,'r')
% this is the final estimation error
solution_optimal_error_4=norm(solution_final_4-initial_state_true)/norm(initial_state_true)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              test the method against random selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

number_of_random_experiments=1000;

for i=1:number_of_random_experiments
ss = sort(randsample(N,real_number_sensor_nodes_2));
Crandom=zeros(real_number_sensor_nodes_2,N);
for j=1:real_number_sensor_nodes_2
    Crandom(j,ss(j))=1;
end
% here uncomment if you want different initial guess in every iteration
%initial_guess_initial_state=randn(N,1);
output_sequence=Crandom*STATE_FE;
[solution_random]=estimate_initial_state_gradient(output_sequence,initial_guess_initial_state,Crandom,h,@memory_network_dynamics,@memory_network_dynamics_gradient)
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



