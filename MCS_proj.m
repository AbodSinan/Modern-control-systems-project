%This mini-Project is to design an altitude control system for a Vanguard Ballistic Missle
%Because the modeling of the missle is very complex, including aerodynamics and other things outside the scope of 
%this course I got help from Tewari chapter 13. utilizing the symmetry of The missle about the X and Z axes to 
%eliminate many terms

%I will try to make this program as dynamic as possible

%The following function is an interesting way to determine the cost of the parameters of the missile.
%I was thinking about how to allocate the cost of the parameters. So I had to think in-depth about what is important
%in this application. The fuel of the missle is a major concern in most of the flight, we don't care about how much
%time the missle takes to get back to its trajectory at the start. However, it is the opposite at the end of its
%flight because it is crucial for it to get quickly back to its trajectory ASAP to make sure it hits the target.
%And at this point, the remaining fuel does not matter because the missile is almost at its target. Which is why the
%cost varies with the time. I was initially planning to use a sigmoid function when I stumbled upon a more potent function
% called Gompertz Function which fits perfectly to this application. The function grows exponentially at a predetermined
%point. In this case, I decided the point is around half of the flight time, it grows from 1 to 100. The R values will be
%the inverse of the fuel cost

travel_time = 100
t = 75 %the current point of travel

%The following parameters are for the Gompertz function (obtained by trial and error)
a = 99
b = travel_time / 2
c = 0.07
 
fuel_cost = 1 + (a * exp(-b * exp(-c *t))) %the cost of fuel as a Gompertz function
acc_cost = 99 - fuel_cost %the cost of the accuracy relevant to the trajectory

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
Q = fuel_cost * eye(6)
% make an R matrix (weights are high because fuel is important)
R = acc_cost *eye(2)
%obtain the gain matrix using LQR
K = lqr(sys, Q, R)
%Initialise some disturbances
Vd = .001* eye(6);
Vn = 1 * eye(4); %sensor noise

BF = [B Vd 0*B 0*B] %the B matrix with disturbance and noise
val = zeros(4, 6)
sysC = ss(A, BF, C, [D val Vn]) % the system including noise and disturbances as inputs
val = zeros(4, 10)
sysFullOutput = ss(A,BF,C,[D val]) %the system with noise + dist = 0
%Obtain the gain fuction for the estimator
[L, P, K] = lqe(sys, Vd, Vn);
%obtain the state-space model of the kalman filter
Kal = ss(A-L*C, [B L], eye(6), 0*[B L]);

%the simulation of the system was done with the help of Control Bootcamp series on youtube
dt = .01
sim_time = dt:dt:travel_time
u = 0*sim_time
u(100:120) = 100;
u(1500:1520) = -100;

uDIST = randn (6, size(sim_time, 2)) % the disturbance input has 6 rows (one for each state)
uNOISE = randn (4, size(sim_time, 2))% the noise has 4 rows (one for each measured output
uAUG = [u; u; Vd*Vd*uDIST; uNOISE]  % the augmented input has 12 rows (6 for dist + 4 for noise + 2 inputs)

[xtrue,t] = lsim(sysFullOutput,uAUG',t) %to obtain the real values of x
[y, t] = lsim(sysC, uAUG', t); %to simulate the system without kalman filter and with noise+dist

[x, t] = lsim(Kal, [u;u; y']', t); %to simulate the system with kalman filter and noise+dist

figure
plot(t,y)
hold on
plot(t,xtrue(:,1),'r')
plot(t,x(:,1),'k--')





