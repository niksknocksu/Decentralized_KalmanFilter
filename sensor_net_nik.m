clear;
noOfNodes  = 110; %Total Number of nodes
sen_1_Nodes = 5;  % Number of TOA Sensor Nodes
sen_2_Nodes = 5;  % Number of DOA Sensor Nodes
comNodes = 100; % Number of communication nodes
figure(1);
clf;
grid on
xticks(0:500:5000);
yticks(0:500:5500);
hold on;
R = 829; % maximum range in meters;
sen1_Xloc = xlsread('Sensor_Net.xlsx',1,'D3:D7'); % list of TOA Sensor X coord.
sen1_Yloc = xlsread('Sensor_Net.xlsx',1,'E3:E7'); % List of TOA sensor Nodes Y coordinate.
sen2_Xloc =  xlsread('Sensor_Net.xlsx',1,'G3:G7'); % list of DOA Sensor X coord.
sen2_Yloc =  xlsread('Sensor_Net.xlsx',1,'H3:H7'); % list of DOA Sensor Y coord.
comXloc = xlsread('Sensor_Net.xlsx',1,'A3:A102'); % List of Comm Nodes X coord.
comYloc = xlsread('Sensor_Net.xlsx',1,'B3:B102'); % List of Comm Nodes Y coord.
netXloc = [sen1_Xloc;sen2_Xloc;comXloc];  % Augment all lists 
netYloc = [sen1_Yloc;sen2_Yloc;comYloc];
for i = 1:sen_1_Nodes
    figure(1)
    p1 = plot(sen1_Xloc(i), sen1_Yloc(i), 'ko', 'MarkerFaceColor','red');
    %text(netXloc(i), netYloc(i), num2str(i));
end
hold on
for i = 1:sen_2_Nodes
    figure(1)
    p2 = plot(sen2_Xloc(i), sen2_Yloc(i), 'b^','MarkerFaceColor', 'yellow');
    %text(netXloc(i), netYloc(i), num2str(i));
end
hold on
for i = 1:comNodes
    figure(1)
    p3 = plot(comXloc(i), comYloc(i), 'ks');
    %text(netXloc(i), netYloc(i), num2str(i));hold on
end
hold on
for i = 1:noOfNodes
    for j = 1:noOfNodes
            distance = sqrt((netXloc(i) - netXloc(j))^2 + (netYloc(i) - netYloc(j))^2);
            if distance <= R
                matrix(i, j) = 1;   % there is a link;
              p4 =   line([netXloc(i) netXloc(j)], [netYloc(i) netYloc(j)],'Color','blue', 'Linestyle','-');
            else
                matrix(i, j) = inf;
            end
    end
end
legend([p1,p2,p3,p4],{'TOA','DOA', 'COM','Link'});

w = matrix; 
for i = 1:numel(w)
    if w(i) == Inf
        w(i) = 0;
    end
end
for i = 1:noOfNodes
    w(i,i) = 0;
    for j = 1:noOfNodes
            if ((i ~= j) && (w(i,j) ~= 0))
                X = matrix(i,:);
                Y = matrix(j,:); 
                w(i,j) = 1/(1+(max(numel(X(X~=Inf)), numel(Y(Y~= Inf)))));
            end
    end
end

for i = 1:noOfNodes
    w(i,i) = 1 - (sum(w(i,:)));
end

% Metropolis wieght matrix has been constructed and working fine.

