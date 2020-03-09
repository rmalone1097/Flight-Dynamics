function Euler321_mod(tmax,psi0,theta0,phi0,omega,order,wtype)
% wtype: specify based on form of omega(t) specified, here, if: 
% wtype=0 (constant omega components)
% wtype=1 (components vary as polynomial functions of time)
% wtype=2 (components vary as sin(Oi t))

% omega is a 3x1 vector containing the components of body angular velocity, or
% coefficients in the case that omega(t)
w1 = omega(1);
w2 = omega(2);
w3 = omega(3);

% Order is a 3x1 vector that contains (a) order of the temporal polynomial 
% for body angular velocity components
% OR (b) frequencies of temporal oscillations, depending on "wtype"
O1 = order(1);
O2 = order(2);
O3 = order(3);

% Specify tolerance (optional)
options = odeset('RelTol',1e-5);

% Call ode45 to solve ODEs specified by function below
[T,Y] = ode45(@Euler321, [0, tmax],[psi0,theta0,phi0],options);

% ODEs to be solved in state-space form
    function [ dydt ] = Euler321( t,y )
    % Initialize 
    dydt = zeros(3,1); 
    
    % State-space form (Euler 3-2-1)
        if wtype==0 % Constant w
        dydt(1) = 1/cos(y(2))*(sin(y(3))*w2+cos(y(3))*w3);
        dydt(2) = 1/cos(y(2))*(cos(y(3))*cos(y(2))*w2-sin(y(3))*cos(y(2))*w3);
        dydt(3) = 1/cos(y(2))*(cos(y(2))*w1+sin(y(3))*sin(y(2))*w2+cos(y(3))*sin(y(2))*w3);
        
        elseif wtype==1 % Polynomial w
        dydt(1) = 1/cos(y(2))*(sin(y(3))*w2*t^O2+cos(y(3))*w3*t^O3);
        dydt(2) = 1/cos(y(2))*(cos(y(3))*cos(y(2))*w2*t^O2-sin(y(3))*cos(y(2))*w3*t^O3);
        dydt(3) = 1/cos(y(2))*(cos(y(2))*w1*t^O1+sin(y(3))*sin(y(2))*w2*t^O2...
            +cos(y(3))*sin(y(2))*w3*t^O3);
        
        elseif wtype==2 % Sine w
        dydt(1) = 1/cos(y(2))*(sin(y(3))*w2*sin(O2*t)+cos(y(3))*w3*sin(O3*t));
        dydt(2) = 1/cos(y(2))*(cos(y(3))*cos(y(2))*w2*sin(O2*t)-sin(y(3))*cos(y(2))*w3*sin(O3*t));
        dydt(3) = 1/cos(y(2))*(cos(y(2))*w1*sin(O1*t)+...
            sin(y(3))*sin(y(2))*w2*sin(O2*t)+cos(y(3))*sin(y(2))*w3*sin(O3*t));
        end
    end

    %Part B: Final DCM @ t=5
    psiF = (Y(length(Y),1));
    thetaF = (Y(length(Y),2));
    phiF = (Y(length(Y),3));
    dcmF = [cos(thetaF)*cos(psiF), cos(thetaF)*sin(psiF), -sin(thetaF);...
            sin(phiF)*sin(thetaF)*cos(psiF)-cos(phiF)*sin(psiF), sin(phiF)*sin(thetaF)*sin(psiF)+cos(phiF)*cos(psiF), sin(phiF)*cos(thetaF);...
            cos(phiF)*sin(thetaF)*cos(psiF)+sin(phiF)*sin(psiF), cos(phiF)*sin(thetaF)*sin(psiF)-sin(phiF)*cos(psiF), cos(phiF)*cos(thetaF)];
    
    fprintf('The DCM at t=5 is: \n');
    disp(dcmF);
    
    %Part C: DCM as funciton of time

    psiT = zeros(1,length(T))';
    thetaT = zeros(1, length(T))';
    phiT = zeros(1, length(T))';
     
    for i = 1:length(T)    
        dcm = [cos(Y(i,2))*cos(Y(i,1)), cos(Y(i,2))*sin(Y(i,1)), -sin(Y(i,2));...
               sin(Y(i,3))*sin(Y(i,2))*cos(Y(i,1))-cos(Y(i,3))*sin(Y(i,1)), sin(Y(i,3))*sin(Y(i,2))*sin(Y(i,1))+cos(Y(i,3))*cos(Y(i,1)), sin(Y(i,3))*cos(Y(i,2));...
               cos(Y(i,1))*sin(Y(i,2))*cos(Y(i,1))+sin(Y(i,3))*sin(Y(i,1)), cos(Y(i,3))*sin(Y(i,2))*sin(Y(i,1))-sin(Y(i,3))*cos(Y(i,1)), cos(Y(i,3))*cos(Y(i,2))];
        psiT(i,1) = atan(dcm(1,2)/dcm(1,1));
        thetaT(i,1) = -asin(dcm(1,3));
        phiT(i,1) = atan(dcm(2,3)/dcm(3,3));
    end
    
% Plot results
figure (1)
plot(T,psiT(:,1),'-',T,thetaT(:,1),'-.',T,phiT(:,1),'--')
legend('Yaw','Pitch','Roll')
title('Yaw, Pitch and Roll as Functions of Time');
xlabel('Time (s)');
ylabel('Angle (rad)');

end


