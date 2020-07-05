%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                               %%
%%                                                                               %%
%%               "Power Flow Program Using Gauss-Seidel Method"                  %%
%%                                                                               %%
%%                 Coded By: Mohammad Solki                                      %%
%%                 Student No: 9411412054                                        %%
%%                 Supervisor: Dr Moslem Salehi                                  %%
%%                 Course: Power System Analysis 1                               %%
%%                                                                               %%
%%     Lorestan University - Summer 2020                                         %%
%%     Source Code Also Available on: https://github.com/d4rk92/PowerFlow        %%
%%                                                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





clc
clear all;
close all;
V_base = 230e3;     % Base Voltage
MVA_base = 100;     % Base MVA 
alpha = 1.6;        % Accelarating Factor

%% Calling function for obtaining line and bus data.
linedata = input_line_data();
busdata = input_bus_data();
Ybus = Ybus_matrix();

%% Compute reactive power 'Q' for the voltage controlled buses.
bus_no = busdata(:,1);     % Bus Number
bus_type = busdata(:,2);   % Bus Type
GenMW = busdata(:,3);      % Active Power Generated
GenMVAR = busdata(:,4);    % Reactive Power Generated
LoadMW = busdata(:,5);     % Active Power Demanded
LoadMVAR = busdata(:,6);   % Reactive Power Demanded
V = busdata(:,7);          % Initial Bus Voltages
del = busdata(:,8);        % Initial Bus Angles
Qmin = busdata(:,9);       % Minimum limit on the reactive power
Qmax = busdata(:,10);      % Maximum limit on the reactive power

V_orig = V.*cos(del) + 1i*V.*sin(del);
V_new = V_orig;                       % V_orig for storing previous iteration values.
P = (GenMW - LoadMW)/MVA_base;        % Pi = PGi - PLi, Active Power at i'th bus.
Q = (GenMVAR - LoadMVAR)/MVA_base;    % Qi = QGi - QLi, Reactive Power at i'th bus.
        
%% Starting the iterations.
% The iteration will continue until the difference between two consecutive
% voltage values becomes less than 0.00001.

tot_buses = length(bus_no);
iter = 1;
disp("Per unit voltage of buses after each iteration are:")
while(V_orig - V_new < 0.00001)
for idx = 1:tot_buses
    temp1 = 0;
    % Computing new voltages for all the load buses.
    if bus_type(idx) == 3
        for b = 1:tot_buses
            if b ~= idx
                temp1 = temp1 + Ybus(idx,b)*V_new(b);
            end
        end
        Vidx_new = ((P(idx)-1i*Q(idx))/V_new(idx) - temp1)/Ybus(idx,idx);
        Vidx_new_acc = (1-alpha)*V_new(idx) + alpha*Vidx_new;
        V_new(idx) = Vidx_new_acc;
    end
    
    % Computing Q values for all the voltage controlled buses.
    for a = 1:tot_buses
        temp = 0;
        if bus_type(a) == 2
            for b = 1:tot_buses
                temp = temp + V_new(a)*Ybus(a,b)*V_new(b);
            end
            Q(a) = -imag(temp);
            if Q(a) < Qmin(a)
                Q(a) = Qmin(a);
            end
            if Q(a) > Qmax(a)
               Q(a) = Qmax(a);
            end
        end
    end

    % Computing new voltages for the voltage controlled buses.
    if bus_type(idx) == 2
        for b = 1:tot_buses
            if b ~= idx
                temp1 = temp1 + Ybus(idx,b)*V_new(b);
            end
        end
        Vidx_new = ((P(idx)-1i*Q(idx))/V_new(idx) - temp1)/Ybus(idx,idx);
        V_temp = Vidx_new*V_new(idx)/abs(Vidx_new);
        V_new(idx) = V_temp;
    end
end

disp(['Iteration no: ',num2str(iter)])
disp(V_orig)

if max(V_orig - V_new) > 0.00001
    V_orig = V_new;
end
iter = iter + 1;
end













function busdata = input_bus_data()
%          |Bus|Type|P(MW)G|Q(MVAR)G|P(MW)L| Q(MVAR)L |V,pu|delta|Qmin|Qmax|

busdata = [  1   4     0      0      0   0          1      0      0   0  ;
             4   6     0      0      0   0          1      0      0   0  ;
             4   5     0      0      0   0          1      0      0   0  ;
             6   9     0      0      0   0          1      0      0   0  ;
             5   7     0      0      0   0          1      0      0   0  ;
             9   8     0      0      0   0          1      0      0   0  ;
             9   3     0      0      0   0          1      0      0   0  ;
             7   2     0      0      0   0          1      0      0   0  ;
             7   8     0      0      0   0          1      0      0   0 ;];
 
end

function linedata = input_line_data()
% R - Resistance
% X - Reactance
% G - Conductance
% B - Susceptance
% Y/2 - Shunt admittance
%          |From|To |   R    |   X    |   G    |   B    |  Total  |  Y/2  |
%          |Bus |Bus|   pu   |   pu   |   pu   |   pu   | charging|  pu   |
linedata = [ 1    4   0        0        0        0         0       0      ;  
             4    6   0        0        0        0         0       0      ;
             4    5   0        0        0        0         0       0      ;
             6    9   0        0        0        0         0       0      ;
             5    7   0        0        0        0         0       0      ;
             9    8   0        0        0        0         0       0      ;
             9    3   0        0        0        0         0       0      ;
             7    2   0        0        0        0         0       0      ;
             7    8   0        0        0        0         0       0      ];

end
















function Ybus = Ybus_matrix()
linedata = input_line_data();

init_bus = linedata(:,1);       % From bus number
final_bus = linedata(:,2);      % To bus number
r = linedata(:,3);              % Resistance, R
x = linedata(:,4);              % Reactance, X
g = linedata(:,5);              % Conductance, G
b = linedata(:,6);              % Susceptance, B
s = linedata(:,8);              % Shunt or Ground Admittance
z = r + 1i*x;                   % Impedance
y = g + 1i*b;                   % Admittance

tot_buses = max(max(init_bus),max(final_bus));    % total no. of buses
tot_branches = length(init_bus);                  % no. of branches
Ybus = zeros(tot_buses,tot_branches);             % Initialising YBus
sbus = zeros(tot_buses,tot_branches);             % Initialising Shunt bus
for a = 1:tot_buses
    Ybus(init_bus(a),final_bus(a)) = -y(a);
    Ybus(final_bus(a),init_bus(a)) = -y(a);
end

for b = 1:tot_buses
    sbus(init_bus(b),final_bus(b)) = s(b);
    sbus(final_bus(b),init_bus(b)) = s(b);
end

for b = 1:tot_buses
    for c = 1:tot_buses
        if c ~= b
            Ybus(b,b) = Ybus(b,b) - Ybus(b,c) + sbus(b,c);
        end
    end
end
end














