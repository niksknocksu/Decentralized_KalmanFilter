%{ 
DEKF for the original sensor network in the paper,
total no of DOA sensors is 5
total no of TOA sensors is 5
total comm nodes is 100
%}

% To construct wieght matrix and the sensor network
run sensor_net_nik
%Define all constants
N = 200;
M = 10; %Monte Carlo Runs
err_xpos = zeros(M,N);
dt = 1;
t = dt.*(0:N-1);
x_act = zeros(4,N);
F = [1,dt,0,0;...
     0,1,0,0;...
     0,0,1,dt;...
     0,0,0,1];
B = [dt*dt/2, 0;...
     dt,0;...
     0, dt*dt/2;...
     0, dt];
a = [0;0];
w_x = 1E-3;
w_y = 1E-3;
Q = [(w_x)*(dt^3)/3,(w_x)*(dt^2)/2,0,0;...
         (w_x) * (dt^2)/2 ,(w_x) *dt,0,0;...
         0,0,(w_y) * (dt^3)/3,(w_y) * (dt^2)/2;...
         0,0,(w_y) * (dt^2)/2,(w_y) *dt];
     
MATRIX = cell(noOfNodes,6);
for m = 1:M
%Start Monte Carlo iterations from here
    % 1st column of MATRIX is location of each sensor
    for i = 1:sen_1_Nodes  %Sen_1_Nodes = TOA Sensors
        MATRIX{i,1}(1) = sen1_Xloc(i); %first col of comm matrix is location
        MATRIX{i,1}(2) = sen1_Yloc(i);
    end
    
    for i = 1:(sen_2_Nodes) %Sen_2_Nodes = DOA Sensors
        MATRIX{(i+sen_1_Nodes),1}(1) = sen2_Xloc(i);
        MATRIX{(i+sen_1_Nodes),1}(2) = sen2_Yloc(i);
    end
    
    for i = 1:comNodes
        MATRIX{(i+sen_1_Nodes+sen_2_Nodes),1}(1) = comXloc(i);
        MATRIX{(i+sen_1_Nodes+sen_2_Nodes),1}(2) = comYloc(i);
    end
    
    % Generate Measurement Sequence for both TOA and DOA Sensor
    x_act(:,1) = [300;5;200;3];

    % Simulating Ground Truth around which measurements would be generated
    for i = 2:N
           x_act(:,i) = F*x_act(:,(i-1)) + (B*a) + mvnrnd([0,0,0,0],Q)';
    end
    
    %Measurement Simulation
    R_TOA = 50; %variance for TOA sensors
    z_TOA = zeros(sen_1_Nodes,N);
    for i = 1:sen_1_Nodes
        for k = 1:N
            %z_TOA = sqrt(x-x_i^2 + y-y_i^2) + Noise
            z_TOA(i,k) = func_TOA(x_act(1,k),x_act(3,k),MATRIX{i,1}(1),MATRIX{i,1}(2)) + ((mvnrnd(0,R_TOA))');
        end
    end
    
    z_DOA = zeros(sen_2_Nodes,N);
    R_DOA = 4*pi/180; %Variance for DOA Sensors in radians
    for i = 1:sen_2_Nodes
        for k = 1:N
            z_DOA(i,k) = func_DOA(x_act(1,k),x_act(3,k),MATRIX{i+sen_1_Nodes,1}(1),MATRIX{i+sen_1_Nodes,1}(2)) + ((mvnrnd(0,R_DOA))');
        end
    end
    
    %Leaving for dinner, continue updating the MATRIX contents when back-
    %30/7 8:10 PM
    
    % 2nd Column of MATRIX is state vector
    for i = 1:sen_1_Nodes
        MATRIX{i,2} = [z_TOA(i,1)*cos(z_DOA(i,1));3;z_TOA(i,1)*sin(z_DOA(i,1));1]; 
    end
    
    for i = 1:sen_2_Nodes
        MATRIX{(i+sen_1_Nodes),2} = [z_TOA(i,1)*cos(z_DOA(i,1));3;z_TOA(i,1)*sin(z_DOA(i,1));1];
    end
    
    for i = 1:comNodes
        MATRIX{(i+sen_1_Nodes+sen_2_Nodes),2} = [z_TOA(1,1)*cos(z_DOA(1,1));3;z_TOA(1,1)*sin(z_DOA(1,1));1];
    end
    
    % 4th Column of MATRIX is corresponding prior information MATRIX
    % Initial Information matrix is set to inv(P)
    for i = 1:noOfNodes
        MATRIX{i,4} = inv(diag([400,25,500,9]));
    end
    
    for i = 1:noOfNodes
        MATRIX{i,3} = MATRIX{i,4}*MATRIX{i,2};
    end
    
    for i = 1:noOfNodes
        MATRIX{i,5} = [0;0;0;0];
    end
    
    for i = 1:noOfNodes
        MATRIX{i,6} = zeros(4);
    end
    
  % Main DEKF Algorithm as given on 3rd page of the Paper    
  
   W = inv(Q);
   V_TOA = inv(R_TOA);
   V_DOA = inv(R_DOA);
   err_xpos(m,1) = x_act(1,1) - MATRIX{1,2}(1,1);
   C_TOA = zeros(sen_1_Nodes,4);% Jacobian for TOA
   C_DOA = zeros(sen_2_Nodes,4); % Jacobian for DOA
   y_TOA = zeros(sen_1_Nodes,1);
   y_DOA = zeros(sen_2_Nodes,1);
   for k = 2:N

        % calculate q and omega for this k
        % find prior estimates for q and omega
        % Prediction update :
            for i = 1:noOfNodes
                MATRIX{i,2}(:,k) = F*MATRIX{i,2}(:,k-1);   % x(t|t-1)
            end
        %prior Information Matrix :
            for i = 1:noOfNodes
                % MATRIX{i,4} = W - (W*F*(inv(MATRIX{i,4} + (F')*W*F))*(F')*W);
                %Replace above with the one given in wikipedia
                M_k = inv(F)' * MATRIX{i,4} * inv(F);
                C_k = M_k*inv(M_k + W);
                L_k = eye(4) - C_k;
                MATRIX{i,4} = (L_k * M_k *(L_k')) + C_k * W* (C_k');
            end
            
            for i = 1:noOfNodes
                MATRIX{i,3} = MATRIX{i,4} * MATRIX{i,2}(:,k);
            end
            
            for i = 1:sen_1_Nodes
                C_TOA(i,:) = [(MATRIX{i,2}(1,k) - MATRIX{i,1}(1))/(sqrt((MATRIX{i,2}(1,k) - MATRIX{i,1}(1))^2 + (MATRIX{i,2}(3,k) - MATRIX{i,1}(2))^2)),0,...
                              (MATRIX{i,2}(3,k) - MATRIX{i,1}(2))/(sqrt(((MATRIX{i,2}(1,k) - MATRIX{i,1}(1))^2) + ((MATRIX{i,2}(3,k) - MATRIX{i,1}(2))^2))),0];
                y_TOA(i) = z_TOA(i,k) - func_TOA(MATRIX{i,2}(1,k), MATRIX{i,2}(3,k), MATRIX{i,1}(1), MATRIX{i,1}(2)) + (C_TOA(i,:) * MATRIX{i,2}(:,k));
                MATRIX{i,5} = (C_TOA(i,:))'*V_TOA*y_TOA(i) ;
                MATRIX{i,6} = (C_TOA(i,:))'*V_TOA*C_TOA(i,:);
            end
                
            for i = sen_2_Nodes
                C_DOA(i,:) = [-1*(MATRIX{i+sen_1_Nodes,2}(3,k) - MATRIX{i+sen_1_Nodes,1}(2))/(((MATRIX{i+sen_1_Nodes,2}(1,k) - MATRIX{i+sen_1_Nodes,1}(1))^2) + ((MATRIX{i+sen_1_Nodes,2}(3,k) - MATRIX{i+sen_1_Nodes,1}(2))^2)),0,...
                          (MATRIX{i+sen_1_Nodes,2}(1,k) - MATRIX{i+sen_1_Nodes,1}(1))/(((MATRIX{i+sen_1_Nodes,2}(1,k) - MATRIX{i+sen_1_Nodes,1}(1))^2) + ((MATRIX{i+sen_1_Nodes,2}(3,k) - MATRIX{i+sen_1_Nodes,1}(2))^2)),0];
                y_DOA(i) = z_DOA(i,k) - func_DOA(MATRIX{i+sen_1_Nodes,2}(1,k), MATRIX{i+sen_1_Nodes,2}(3,k), MATRIX{i+sen_1_Nodes,1}(1), MATRIX{i+sen_1_Nodes,1}(2)) + (C_DOA(i,:) * MATRIX{i+sen_1_Nodes,2}(:,k));
                MATRIX{i+sen_1_Nodes,5} = (C_DOA(i,:))'*V_DOA*y_DOA(i) ;
                MATRIX{i+sen_1_Nodes,6} = (C_DOA(i,:))'*V_DOA*C_DOA(i,:);
            end
            
            % Apply consensus on both prior and novel information
            
            %Consensus on Prior Information State vectors  
        CONS = cell(noOfNodes);
        L = 5;
        for l = 1:L
            for i = 1:noOfNodes
             temp = [0;0;0;0];
                     for j = 1:noOfNodes
                           temp = (w(i,j).* MATRIX{j,3})+ temp;  
                     end
             CONS{i} = temp;
            end
            for i= 1:noOfNodes
                MATRIX{i,3} = CONS{i};
            end
        end
        
        %Consensus on prior information Matrices
        for l = 1:L
            for i = 1:noOfNodes
             temp = zeros(4);
                     for j = 1:noOfNodes
                           temp = (w(i,j).* MATRIX{j,4})+ temp;  
                     end
             CONS{i} = temp;
            end
            for i= 1:noOfNodes
                MATRIX{i,4} = CONS{i};
            end
        end
    
       %Consensus on Novel information vector
        for l = 1:L
            for i = 1:noOfNodes
             temp = [0;0;0;0];
                     for j = 1:noOfNodes
                           temp = (w(i,j).* MATRIX{j,5})+ temp;  
                     end
             CONS{i} = temp;
            end
            for i= 1:noOfNodes
                MATRIX{i,5} = CONS{i};
            end
        end
        
        %Consensus on Novel Information Matrix
    
        for l = 1:L
            for i = 1:noOfNodes
             temp = zeros(4);
                     for j = 1:noOfNodes
                           temp = (w(i,j).* MATRIX{j,6})+ temp;  
                     end
             CONS{i} = temp;
            end
            for i= 1:noOfNodes
                MATRIX{i,6} = CONS{i};
            end
        end
        
     % Correction
     % Use gamma to fuse novel and prior information and calculate
     % resulting state vector
      
       gamma = 110;
         for i = 1:noOfNodes
             MATRIX{i,3} = MATRIX{i,3} + gamma.*MATRIX{i,5};
             MATRIX{i,4} = MATRIX{i,4} + gamma.*MATRIX{i,6};
         end
         for i = 1:noOfNodes
             MATRIX{i,2}(:,k) = (inv(MATRIX{i,4}))*MATRIX{i,3};
         end
         err_xpos(m,k) = x_act(1,k) - MATRIX{1,2}(1,k);
   end
end
rms_err_xpos = rms(err_xpos);
figure(2)
plot(t,rms_err_xpos);
grid on
xlabel('Time in seconds')
ylabel('RMSE in meters')
title('200 Monte Carlo RMSE')

  

function g = func_TOA(x,y,x_i,y_i)
    g = sqrt(((x-x_i)^2) + ((y-y_i)^2));
end

function h = func_DOA(x,y,x_i,y_i)
    h = atan2((y-y_i),(x-x_i));
end