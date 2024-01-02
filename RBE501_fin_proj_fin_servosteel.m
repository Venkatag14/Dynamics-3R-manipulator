%generating 2 trajectories
%trajectory 1 from pose 1 to pose 2 with velocity 0rad/s at pose 2

% Clear all the repos 
clc; close all;

% Number of Samples 
numsamples = 250;
% Time to run the total iterations
numsecond = 2.5;

% Num of joints to consider - Here 3 Link Manipulator  
numJoints = 3;

% Sampling of time 
t = linspace(0,numsecond,numsamples);

% Link Lengths
l1 = 0.12; 
l2 = 0.09;
l3 = 0.06;

% Home Configuration
S = [0 0 1 0 0 0;
     0 0 1 0 -l1 0;
     0 0 1 0 -(l1+l2) 0]';

% Mass and Inertia's of Links - Obtained from CAD Modelling 
m1 = 0.427764;
I1 = [ 0.000242 0 0;
      0 0.000744 0;
      0 0 0.000616];
G(:,:,1) = [I1 zeros(3,3);
            zeros(3,3) m1*eye(3,3)];
m2 =  0.34238111;
I2 = [0.00016354 0 0;
      0.0000260  0.00058778 0;
      0 0  0.00051727];
G(:,:,2) = [I2 zeros(3,3);
            zeros(3,3) m2*eye(3,3)];
m3 = 0.8683981;
I3 = [0.0004148 0 0;
      0 0.0010697 0;
       0 0  0.0010796];
G(:,:,3) = [I3 zeros(3,3);
            zeros(3,3) m3*eye(3,3)];

% Accleration due to Gravity
g = 9.81;

% Home Transformation Matrices
M1 = [1 0 0 0.030469;
      0 1 0 0;
      0 0 1 0;
      0 0 0 1];
M2 = [1 0 0 0.12+0.00817053;
      0 1 0 0;
      0 0 1 0;
      0 0 0 1];
M3 = [1 0 0 0.12+0.09+0.0343233;
      0 1 0 0;
      0 0 1 0;
      0 0 0 1];
M4 = [1 0 0 0.12+0.09+0.06;
      0 1 0 0;
      0 0 1 0;
      0 0 0 1];
M = zeros(4,4,4);
M(:,:,1) = M1;
M(:,:,2) = inverse(M1)*M2;
M(:,:,3) = inverse(M2)*M3;
M(:,:,4) = inverse(M3)*M4;

A = zeros(6,3);
A(:,1) = Scrsp_to_lo(M1,S(:,1));
A(:,2) = Scrsp_to_lo(M2,S(:,2));
A(:,3) = Scrsp_to_lo(M3,S(:,3));

% Initialization of Joint angles, Velocity and acceleration 
% and other parameters

theta = zeros(250,3);
theta_dot = zeros(250,3);
theth_ddot = zeros(250,3);
Vi0 = zeros(6,1);
Vi = zeros(6,3,250);
Vi_dot0 = [0,0,0,0,g,0]';
Vi_dot = zeros(6,3,250);
Ftip = [0,0,0,0,0,0]';
Fi = zeros(6,3,250);
tau = zeros(250,3);

% Goal Position 1 - POSE1
goal(:,:,1) = [1 0 0 0.1;
               0 1 0 0.05;
               0 0 1 0;
               0 0 0 1];
% Goal Position 2 - POSE2
goal(:,:,2) = [1 0 0 0.1;
               0 1 0 0.1;
               0 0 1 0;
               0 0 0 1];
% Goal Position 3 - POSE3
goal(:,:,3) = [0 1 0 0.15;
               -1 0 0 0.075;
               0 0 1 0;
               0 0 0 1];

% Angle 
q = [deg2rad(90) deg2rad(-120) deg2rad(30);
     deg2rad(110) deg2rad(-120) deg2rad(0);
     deg2rad(50) deg2rad(-30) deg2rad(-110)];

% Inverse Kinematics
for j = 1:3
    for i = 1:6
        Tbs = inverse(fkine(S, M4, q(j,:)));
        Twbd = logm(Tbs*goal(:,:,j));
        Vbd = [Twbd(3,2), Twbd(1,3), Twbd(2,1), Twbd(1,4), Twbd(2,4), Twbd(3,4)]';
        Jb = Adjoint(Tbs)*jacob0(S,q(j,:));
        %Jinv = pinv(Jb);
        Jinv = inv(Jb'*Jb)*Jb';
        Theta = q(j,:)';
        Theta = Theta+Jinv*Vbd;
        q(j,:) = Theta'; 
    end
end


disp(q);
disp(fkine(S,M4,q(1,:)));
disp(fkine(S,M4,q(2,:)));
disp(fkine(S,M4,q(3,:)));
pose1 = [q(1,1), 0, 0;
        q(1,2),0, 0;
         q(1,3), 0, 0];
pose2 = [q(2,1),0, 0;
        q(2,2),  0, 0;
        q(2,3), 0, 0];
pose3 = [q(3,1), 0, 0;
        q(3,2), 0, 0;
        q(3,3),0, 0];

% Initialize quinticpoly for trajectories
[coeff11] = quinticpoly(0, 1, pose1(1,1),pose1(1,2),pose1(1,3), pose2(1,1),pose2(1,2),pose2(1,3));
[coeff12] = quinticpoly(0, 1.5, pose2(1,1),pose2(1,2),pose2(1,3), pose3(1,1),pose3(1,2),pose3(1,3));
[coeff21] = quinticpoly(0, 1, pose1(2,1),pose1(2,2),pose1(2,3), pose2(2,1),pose2(2,2),pose2(2,3));
[coeff22] = quinticpoly(0, 1.5, pose2(2,1),pose2(2,2),pose2(2,3), pose3(2,1),pose3(2,2),pose3(2,3));
[coeff31] = quinticpoly(0, 1, pose1(3,1),pose1(3,2),pose1(3,3), pose2(3,1),pose2(3,2),pose2(3,3));
[coeff32] = quinticpoly(0, 1.5, pose2(3,1),pose2(3,2),pose2(3,3), pose3(3,1),pose3(3,2),pose3(3,3));

% Generation of intermediate position from Pose1 ==> Pose2
for i = 1:100
    a = coeff11(1);
    b = coeff11(2);
    c = coeff11(3);
    d = coeff11(4);
    e = coeff11(5);
    f = coeff11(6);
    theta(i,1) = a+b*t(i)+c*t(i)^2+d*t(i)^3+e*t(i)^4+f*t(i)^5;
    theta_dot(i,1) = b+2*c*t(i)+3*d*t(i)^2+4*e*t(i)^3+5*f*t(i)^4;
    theta_ddot(i,1) = 2*c+6*d*t(i)+12*e*t(i)^2+20*f*t(i)^3;
    a = coeff21(1);
    b = coeff21(2);
    c = coeff21(3);
    d = coeff21(4);
    e = coeff21(5);
    f = coeff21(6);
    theta(i,2) = a+b*t(i)+c*t(i)^2+d*t(i)^3+e*t(i)^4+f*t(i)^5;
    theta_dot(i,2) = b+2*c*t(i)+3*d*t(i)^2+4*e*t(i)^3+5*f*t(i)^4;
    theta_ddot(i,2) = 2*c+6*d*t(i)+12*e*t(i)^2+20*f*t(i)^3;
    a = coeff31(1);
    b = coeff31(2);
    c = coeff31(3);
    d = coeff31(4);
    e = coeff31(5);
    f = coeff31(6);
    theta(i,3) = a+b*t(i)+c*t(i)^2+d*t(i)^3+e*t(i)^4+f*t(i)^5;
    theta_dot(i,3) = b+2*c*t(i)+3*d*t(i)^2+4*e*t(i)^3+5*f*t(i)^4;
    theta_ddot(i,3) = 2*c+6*d*t(i)+12*e*t(i)^2+20*f*t(i)^3;    
end

% Generation of intermediate position from Pose2 ==> Pose3
for i = 101:250
    a = coeff12(1);
    b = coeff12(2);
    c = coeff12(3);
    d = coeff12(4);
    e = coeff12(5);
    f = coeff12(6);
    theta(i,1) = a+b*(t(i)-1)+c*(t(i)-1)^2+d*(t(i)-1)^3+e*(t(i)-1)^4+f*(t(i)-1)^5;
    theta_dot(i,1) = b+2*c*(t(i)-1)+3*d*(t(i)-1)^2+4*e*(t(i)-1)^3+5*f*(t(i)-1)^4;
    theta_ddot(i,1) = 2*c+6*d*(t(i)-1)+12*e*(t(i)-1)^2+20*f*(t(i)-1)^3;
    a = coeff22(1);
    b = coeff22(2);
    c = coeff22(3);
    d = coeff22(4);
    e = coeff22(5);
    f = coeff22(6);
    theta(i,2) = a+b*(t(i)-1)+c*(t(i)-1)^2+d*(t(i)-1)^3+e*(t(i)-1)^4+f*(t(i)-1)^5;
    theta_dot(i,2) = b+2*c*(t(i)-1)+3*d*(t(i)-1)^2+4*e*(t(i)-1)^3+5*f*(t(i)-1)^4;
    theta_ddot(i,2) = 2*c+6*d*(t(i)-1)+12*e*(t(i)-1)^2+20*f*(t(i)-1)^3;
    a = coeff32(1);
    b = coeff32(2);
    c = coeff32(3);
    d = coeff32(4);
    e = coeff32(5);
    f = coeff32(6);
    theta(i,3) = a+b*(t(i)-1)+c*(t(i)-1)^2+d*(t(i)-1)^3+e*(t(i)-1)^4+f*(t(i)-1)^5;
    theta_dot(i,3) = b+2*c*(t(i)-1)+3*d*(t(i)-1)^2+4*e*(t(i)-1)^3+5*f*(t(i)-1)^4;
    theta_ddot(i,3) = 2*c+6*d*(t(i)-1)+12*e*(t(i)-1)^2+20*f*(t(i)-1)^3;
end

x = zeros(250,3);
y = zeros(250,3);
for i = 1:250
    x(i,1) = l1*cos(theta(i,1));
    x(i,2) = l1*cos(theta(i,1))+l2*cos(theta(i,1)+theta(i,2));
    x(i,3) = l1*cos(theta(i,1))+l2*cos(theta(i,1)+theta(i,2))+l3*cos(theta(i,1)+theta(i,2)+theta(i,3));
    y(i,1) = l1*sin(theta(i,1));
    y(i,2) = l1*sin(theta(i,1))+l2*sin(theta(i,1)+theta(i,2));
    y(i,3) = l1*sin(theta(i,1))+l2*sin(theta(i,1)+theta(i,2))+l3*sin(theta(i,1)+theta(i,2)+theta(i,3));

end
for i = 1:250
    T10 = inverse(trans(A(:,1),M(:,:,1),theta(i,1)));
    Vi(:,1,i) = A(:,1)*theta_dot(i,1)+Adjoint(T10)*Vi0;
    Vi_dot(:,1,i) = A(:,1)*theta_ddot(i,1)+Adjoint(T10)*Vi_dot0+ad(Vi(:,1,i))*A(:,1)*theta_dot(i,1);
    for j = 2:3
        Tjj_1 = inverse(trans(A(:,j),M(:,:,j),theta(i,j)));
        Vi(:,j,i) = A(:,j)*theta_dot(i,j)+Adjoint(Tjj_1)*Vi(:,j-1,i);
        Vi_dot(:,j,i) = A(:,j)*theta_ddot(i,j)+Adjoint(Tjj_1)*Vi_dot(:,j-1,i)+ad(Vi(:,j,i))*A(:,j)*theta_dot(i,j);
    end
    T43 = inverse(M(:,:,4));
    Fi(:,3,i) = Adjoint(T43)'*Ftip+G(:,:,3)*Vi_dot(:,3,i)-ad(Vi(:,3,i))'*G(:,:,3)*Vi(:,3,i);
    tau(i,3) = Fi(:,3,i)'*A(:,3);
    for j = 2:-1:1
        Tj1j = inverse(trans(A(:,j+1),M(:,:,j+1),theta(i,j+1)));
        Fi(:,j,i) = Adjoint(Tj1j)'*Fi(:,j+1,i)+G(:,:,j)*Vi_dot(:,j,i)-ad(Vi(:,j,i))'*G(:,:,j)*Vi(:,j,i);
        tau(i,j) = Fi(:,j,i)'*A(:,j);
    end
end
% Realtime Simulation

% Create subplots for different plots
figure;

subplot(4, 2, [1, 3, 5]);
hold on;
position_2d_plot = gobjects(4, 1);  % Initialize an array to store line objects

% Define colors for each link
link_colors = {'k', 'r', 'g', 'b'};

for j = 1:4
    position_2d_plot(j) = plot(0, 0, '-o', 'LineWidth', 2, 'MarkerSize', 4, 'MarkerFaceColor', 'k');
    set(position_2d_plot(j), 'Color', link_colors{j});
end

xlabel('x'), ylabel('y'), title('Position');
axis([-0.2, 0.3, -0.1, 0.3]);
grid on;

% Other subplots
subplot(4, 2, 2);
torque_plot = plot(t, tau);
xlabel('time [s]'), ylabel('torque N-M'), title('Torque');

subplot(4, 2, 4);
position_plot = plot(t, theta);
xlabel('time [s]'), ylabel('q [rad]'), title('Angle');

subplot(4, 2, 6);
velocity_plot = plot(t, theta_dot);
xlabel('time [s]'), ylabel('q [rad/s]'), title('Joint angular Velocity');

subplot(4, 2, 7);
acceleration_plot = plot(t, theta_ddot);
xlabel('time [s]'), ylabel('q [rad/s^2]'), title('Joint angular Acceleration');

subplot(4,2,8);
ind_pos_plot = plot(x, y);
xlabel('X'), ylabel('Y'), title('Individual Link Positions');

% Display the size of theta
disp('Size of theta:');
disp(size(theta));
disp(size(theta_dot))
disp(size(x))
disp(size(y))
disp(size(tau))

% Simulation loop
for i = 1:numsamples
    % Compute joint positions
    
    joint_positions = [0, 0; 
                       l1*cos(theta(i, 1)), l1*sin(theta(i, 1)); 
                       l1*cos(theta(i, 1)) + l2*cos(theta(i, 1) + theta(i, 2)), l1*sin(theta(i, 1)) + l2*sin(theta(i, 1) + theta(i, 2));
                       l1*cos(theta(i, 1)) + l2*cos(theta(i, 1) + theta(i, 2)) + l3*cos(theta(i, 1) + theta(i, 2) + theta(i, 3)), l1*sin(theta(i, 1)) + l2*sin(theta(i, 1) + theta(i, 2)) + l3*sin(theta(i, 1) + theta(i, 2) + theta(i, 3))];

    % Update the stick plot
    for j = 1:4
        set(position_2d_plot(j), 'XData', joint_positions(1:j, 1), 'YData', joint_positions(1:j, 2));
    end

    % Update other plots
    set(torque_plot(1), 'XData', t(1:i), 'YData', tau(1:i, 1));
    set(torque_plot(2), 'XData', t(1:i), 'YData', tau(1:i, 2));
    set(torque_plot(3), 'XData', t(1:i), 'YData', tau(1:i, 3));
    set(position_plot(1), 'XData', t(1:i), 'YData', theta(1:i, 1));
    set(position_plot(2), 'XData', t(1:i), 'YData', theta(1:i, 2));
    set(position_plot(3), 'XData', t(1:i), 'YData', theta(1:i, 3));
    
    set(velocity_plot(1), 'XData', t(1:i), 'YData', theta_dot(1:i, 1));
    set(velocity_plot(2), 'XData', t(1:i), 'YData', theta_dot(1:i, 2));
    set(velocity_plot(3), 'XData', t(1:i), 'YData', theta_dot(1:i, 3));

    set(acceleration_plot(1), 'XData', t(1:i), 'YData', theta_ddot(1:i, 1));
    set(acceleration_plot(2), 'XData', t(1:i), 'YData', theta_ddot(1:i, 2));
    set(acceleration_plot(3), 'XData', t(1:i), 'YData', theta_ddot(1:i, 3));

    set(ind_pos_plot(1), 'XData', x(1:i), 'YData', y(1:i, 1));
    set(ind_pos_plot(2), 'XData', x(1:i), 'YData', y(1:i, 2));
    set(ind_pos_plot(3), 'XData', x(1:i), 'YData', y(1:i, 3));

    % Display simulated joint angles
    disp(['Simulated Joint Angles at Time ', num2str(t(i)), ': ', num2str(theta(i, :))]);

    % Pause for visualization (optional)
    pause(0.05);

    % Update the plots
    drawnow;
end


%% HELPER FUNCTIONS 
% SOME ARE IMPORTED FROM MODERN ROBOTICS OPEN SOURCE IMPLEMENTATION

function coefficients = quinticpoly(t0,tf,q0,qd0,qdd0,qf,qdf,qddf)
    % Compute coefficients vector

    % A matrix for polynomial to represent in Ax=B form
    A = [1, t0, t0^2, t0^3, t0^4, t0^5;
         0, 1, 2*t0, 3*t0^2, 4*t0^3, 5*t0^4;
         0, 0, 2, 6*t0, 12*t0^2, 20*t0^3;
         1, tf, tf^2, tf^3, tf^4, tf^5;
         0, 1, 2*tf, 3*tf^2, 4*tf^3, 5*tf^4;
         0, 0, 2, 6*tf, 12*tf^2, 20*tf^3];

    % Specify the column vector B with initial and final conditions (Position,Velocity and Acceleration)
    B = [q0; qd0; qdd0; qf; qdf; qddf];

    % Solve for the coefficients 
    coefficients = A\B;
end
function T = Scrsp_to_lo(M,S)
    T = Adjoint(inverse(M))*S;
end
function T = trans(S, M, q)
    T = eye(4);
    for i = 1:size(q)
        twist = S(:,i);
        v = twist(4:6);
        omega = twist(1:3);
        H = POE(v,omega,q(i));
        T = T*H;
    end
    function T = POE(v,omega,q)
        I = eye(3);
        Om1_hat = [0, -omega(3), omega(2); omega(3), 0, -omega(1); -omega(2), omega(1), 0];
        exp1 = eye(3)+Om1_hat*sin(q)+(Om1_hat^2)*(1-cos(q));
        T = [exp1,(I*q+(1-cos(q))*Om1_hat+(q-sin(q))*Om1_hat^2)*v;0,0,0,1];
    end 
    T = M*T;
end
function T = inverse(H)
    R = H(1: 3, 1: 3);
    p = H(1: 3, 4);
    T = [R.', -R.'*p;
        0, 0, 0, 1];
end
function AdT = Adjoint(T)
    [R, p] = TransToRp(T);
    AdT = [R, zeros(3); VecToso3(p) * R, R];
end
function [R, p] = TransToRp(T)
    R = T(1: 3, 1: 3);
    p = T(1: 3, 4);
end

function se3mat = VecTose3(V)
    se3mat = [VecToso3(V(1: 3)), V(4: 6); 0, 0, 0, 0];
end

function so3mat = VecToso3(omg)
    so3mat = [0, -omg(3), omg(2); omg(3), 0, -omg(1); -omg(2), omg(1), 0];
end

function omg = so3ToVec(so3mat)
    omg = [so3mat(3, 2); so3mat(1, 3); so3mat(2, 1)];
end

function V = se3ToVec(se3mat)
    V = [se3mat(3, 2); se3mat(1, 3); se3mat(2, 1); se3mat(1: 3, 4)];
end
function V_adj = ad(V1)
    omega = V1(1:3);
    v = V1(4:6);
    
    % Construct the Adjoint matrix for V
    V_adj = [skew(omega), zeros(3);
             skew(v), skew(omega)];
    
end
function S_skew = skew(Si)
    S_skew = [    0  -Si(3)  Si(2);
                Si(3)     0  -Si(1);
                -Si(2)  Si(1)     0];
end
function T = fkine(S, M, q)
    T = eye(4);
    for i = 1:size(S,2)
        twist = S(:,i);
        v = twist(4:6);
        omega = twist(1:3);
        H = POE(v,omega,q(i));
        T = T*H;
    end
    function T = POE(v,omega,q)
        I = eye(3);
        Om1_hat = [0, -omega(3), omega(2); omega(3), 0, -omega(1); -omega(2), omega(1), 0];
        exp1 = eye(3)+Om1_hat*sin(q)+(Om1_hat^2)*(1-cos(q));
        T = [exp1,(I*q+(1-cos(q))*Om1_hat+(q-sin(q))*Om1_hat^2)*v;0,0,0,1];
    end 
    T = T*M;
end
function Js = jacob0(Slist, thetalist)

    Js = Slist;
    T = eye(4);
    for i = 2: length(thetalist)
        T = T * MatrixExp6(VecTose3(Slist(:, i - 1) * thetalist(i - 1)));
	    Js(:, i) = Adjoint(T) * Slist(:, i);
    end
end
function T = MatrixExp6(se3mat)
omgtheta = so3ToVec(se3mat(1: 3, 1: 3));
if NearZero(norm(omgtheta))
    T = [eye(3), se3mat(1: 3, 4); 0, 0, 0, 1];
else
    [omghat, theta] = AxisAng3(omgtheta);
    omgmat = se3mat(1: 3, 1: 3) / theta; 
    T = [MatrixExp3(se3mat(1: 3, 1: 3)), ...
         (eye(3) * theta + (1 - cos(theta)) * omgmat ...
          + (theta - sin(theta)) * omgmat * omgmat) ...
            * se3mat(1: 3, 4) / theta;
         0, 0, 0, 1];
end
end

function judge = NearZero(near)
judge = norm(near) < 1e-6;
end

function [omghat, theta] = AxisAng3(expc3)
theta = norm(expc3);
omghat = expc3 / theta;
end

function  R = MatrixExp3(so3mat)
omgtheta = so3ToVec(so3mat);
if NearZero(norm(omgtheta))
    R = eye(3);
else
    [omghat, theta] = AxisAng3(omgtheta);
    omgmat = so3mat / theta;
    R = eye(3) + sin(theta) * omgmat + (1 - cos(theta)) * omgmat * omgmat;
end
end

    