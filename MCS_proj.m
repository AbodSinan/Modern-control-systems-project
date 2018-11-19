%initialise the parameters (derivations will be shown elsewhere)
% the states include euler angles and displacement
A = [-0.0675 0 0 0 2.3694 0;
    0 -0.0675 0 0 0 -2.3694;
    1 0 0 0 0 0;
    0 2.7285 0 0 0 0;
    1 0 -0.02326 0 -0.03527 0;
    0 -1 0 0.02326 0 -0.03527]
B = [7.0084 0; 0 -7.0084; 0 0; 0 0; 0.05217 0; 0 0.05217]
C = [1 0 0 0 0 0;
    0 1 0 0 0 0;
    329.3374 0 -9.1181 0 -11.8469 0;
    0 -392.3374 0 9.1181 0 -11.8469]
D = [0 0; 0 0; 20.4512 0; 0 20.4512]
%Turn into a state-space representation
sys = ss(A, B, C, D)
%check if the system is ctrb/obsv
rank(ctrb(sys))
rank(obsv(sys))
%make a Q matrix(weights are all unity)
Q = eye(6)
% make an R matrix (weights are high because fuel is important)
R = 100*eye(2)
%obtain the gain matrix using LQR
K = lqr(sys, Q, R)
%Initialise some disturbances (will be better when implemented in simulink)
Vd = .001* eye(6);
Vn = .001;

%Obtain the gain fuction for the estimator
[L, P, K] = lqe(sys, Vd, Vn);
%obtain the state-space model of the kalman filter
Kal = ss(A-L*C, [B L], eye(6), 0*[B L]);
%The next step is to simulate with simulink