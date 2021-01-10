%%  * 7--> 2

%% Loading paras

T = readtable('/Users/chentianchi/Desktop/code/circuit design/DEmodels/param_db/para_s4.csv');
% strcmp('AmeR',T(1,:).repressor{:})
row_names = strcat(T.repressor, '_', T.RBS);
T.Properties.RowNames = row_names;

para_name = {'Y_min','Y_max','K','n'};
HlyIIR = table2struct(T("HlyIIR_H1",para_name));
PhIF   = table2struct(T("PhIF_P1",para_name));
PsrA   = table2struct(T("PsrA_R1",para_name));
LmrA   = table2struct(T("LmrA_N1",para_name));
SrpR   = table2struct(T("SrpR_S4",para_name));
BM3R1  = table2struct(T("BM3R1_B3",para_name));
BetI   = table2struct(T("BetI_E1",para_name));

% ODE
xi = 0.025
gma =0.025

syms t x p
plasmid = @(t,x,p)[xi*(HlyIIR.Y_min + (HlyIIR.Y_max - HlyIIR.Y_min)*HlyIIR.K^HlyIIR.n/(HlyIIR.K^HlyIIR.n + (x(7) + p)^HlyIIR.n)) - gma*(x(1));
                   xi*(LmrA.Y_min + (LmrA.Y_max - LmrA.Y_min)*LmrA.K^LmrA.n/(LmrA.K^LmrA.n + (x(1) + p)^LmrA.n)) - gma*(x(2));
                   xi*(SrpR.Y_min + (SrpR.Y_max - SrpR.Y_min)*SrpR.K^SrpR.n/(SrpR.K^SrpR.n + (x(1) + x(7))^SrpR.n)) - gma*(x(3))
                   xi*(BM3R1.Y_min + (BM3R1.Y_max - BM3R1.Y_min)*BM3R1.K^BM3R1.n/(BM3R1.K^BM3R1.n + (x(2) + x(3))^BM3R1.n)) - gma*(x(4))
                   xi*(PhIF.Y_min + (PhIF.Y_max - PhIF.Y_min)*PhIF.K^PhIF.n/(PhIF.K^PhIF.n + (x(4))^PhIF.n)) - gma*(x(5))
                   xi*(PsrA.Y_min + (PsrA.Y_max - PsrA.Y_min)*PsrA.K^PsrA.n/(PsrA.K^PsrA.n + (x(5) + x(7))^PsrA.n)) - gma*(x(6))
                   xi*(BetI.Y_min + (BetI.Y_max - BetI.Y_min)*BetI.K^BetI.n/(BetI.K^BetI.n + (x(4) + x(6))^BetI.n)) - gma*(x(7))];


figure             
for rd = 1:50
    p=0;
    u0= randi([1 20],1,7);
    [t,xp] = ode45(@(t,x) plasmid(t,x,p),[0 1e3],u0);
%     figure
%     plot(t,xp(:,6:7),'LineWidth',2)
    plot(t,xp(:,6),'-.','LineWidth',2,"Marker",".","Color",'#77AC30','DisplayName','Gate6')
    
    legend('Gate6 : PsrA')
    hold on 
    plot(t,xp(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7 : BetI')
%     plot(t,xp,'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')

%     legend('Gate 7')

    xlabel("Time Steps");
    ylabel("Concentration");
    title("Randomized initial conditions yield stable states")
    ax = gca; % current axes
    ax.FontSize = 20;
    ax.FontName = "Apple Symbols";

end


%% Two attractors

figure           
for rd = 1:100
     p=0;
%     u0= rand([1 20],1,7)
    u0 = [rand(1,7)*20];% rand(1)+4 rand(1)*0.25]
    [t,xp] = ode45(@(t,x) plasmid(t,x,p),[0 1000.5],u0);
%     figure
    plot(xp(:,6),xp(:,7),'-','LineWidth',1);
    hold on;
    plot(xp(end,6),xp(end,7),'Marker' ,'h','MarkerSize',18,'MarkerFaceColor','#D95319')
    
end
axis([0 6 0 6])
xlabel('Gate 6')
ylabel('Gate 7')
title('Gate6 vs Gate7 all phase plane trajecteries convergence')
ax = gca; % current axes
ax.FontSize = 18;
ax.FontName = "Apple Symbols";


%% random seed
s = rng;

%% Sequential simulation

% 1. Simulation the system untill it goes to steady state
p=0;
rng(s)
u0= randi([1 20],1,7);
[t,xp] = ode45(@(t,x) plasmid(t,x,p),[0 1e3],u0);
figure
plot(t,xp(:,6),'LineWidth',2,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
hold on
plot(t,xp(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
legend('Gate 6','Gate 7')

p = 20;
u0_t1 = xp(end,:);
[t,xp] = ode45(@(t,x) plasmid(t,x,p),[0 3e3],u0_t1);
figure
plot(t,xp(:,6),'LineWidth',2,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
hold on
plot(t,xp(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
legend('Gate 6','Gate 7')
xlabel('Time')
ylabel('Concentration')
title('Circuit oscillatory switching dynamics with constant input signal')
ax = gca; % current axes
ax.FontSize = 20;
ax.FontName = "Apple Symbols";


%% Find the lb and ub of the P inpluse and the min time to wait for the next signal

% Find the local minimum index
locs2 =  find(islocalmin(xp(:,6)));
ind_min1 = locs2(2,1);

% Test the time starting from ti =1 to ti = first local min
t_wait = [];
ul_range =[];
t_range = [];
for ti = 1: ind_min1
    p = 0
    u0_t11 = xp(ti,:)
    fprintf('The duration of P is %d:', t(ti));
    [tt,xp1] = ode45(@(t,x) plasmid(t,x,p),[0 3e4],u0_t11);
    % How much time to wait until the system become statble again
    [~,locs]  =findpeaks(xp1(:,6));
    stable_t_ind = locs(1)
%     if t(ti) > 306 & t(ti) < 703% lb ; ub
%         t_wait = [t_wait, t(stable_t_ind)];
%     end
    if xp1(end, 7) > xp(end, 6)%t(ti) > 306 & t(ti) < 703% lb ; ub
        ul_range = [ul_range, ti];
        t_range = [t_range, t(ti)];
        t_wait = [t_wait, tt(stable_t_ind)];
    end
    
%     figure
%     plot(tt,xp1(:,6),'LineWidth',2,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
%     hold on
%     plot(tt,xp1(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
%     legend('Gate 6','Gate 7')
end
max(t_wait)% give the min time to wait(lower bound) for the next signal within the switching range
lb = t_range(1)
ub = t_range(end)

%% single example from above
u0_t11 = xp1(7,:);
p=0
[tt,xp1] = ode45(@(t,x) plasmid(t,x,p),[0 1.5e3],u0_t11);

figure
plot(tt,xp1(:,6),'LineWidth',2,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
hold on
plot(tt,xp1(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')

legend('Gate 6','Gate 7')

%% Now test the square wave inpulse

% set P duration= 400, and t_wait = 570
% Initialize vector of the final plot
P_set = [];
sol_set =[];


% 1. make system goes to steady state
p=0
rng(s)
u0= randi([1 20],1,7);
[t1,sol1] = ode45(@(t,x) plasmid(t,x,p),[0 1000],u0);
figure
plot(t1,sol1(:,6),'LineWidth',2,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
hold on
plot(t1,sol1(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
legend('Gate 6','Gate 7')
title('System goes to steady states')
p1 = zeros(size(t1));

% 2. add inpulse p=20 for 400s
p = 20;
u0_t1 = sol1(end,:);
[t2,sol2] = ode45(@(t,x) plasmid(t,x,p),[0 400],u0_t1);
figure
plot(t2,sol2(:,6),'LineWidth',2,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
hold on
plot(t2,sol2(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
legend('Gate 6','Gate 7')
title('Add inpulse 20 for 400s')
p2 = zeros(size(t2))+20;

% 3. turn off p, p=0, for 570s
p = 0;
u0_t2 = sol2(end,:);
[t3,sol3] = ode45(@(t,x) plasmid(t,x,p),[0 570],u0_t2);
figure
plot(t3,sol3(:,6),'LineWidth',2,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
hold on
plot(t3,sol3(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
legend('Gate 6','Gate 7')
title('Relax system for 570s')
p3 = zeros(size(t3));

% % repeat 2 & 3

% 4. add inpulse p=20 for 400s
p = 20;
u0_t3 = sol3(end,:);
[t4,sol4] = ode45(@(t,x) plasmid(t,x,p),[0 400],u0_t3);
figure
plot(t4,sol4(:,6),'LineWidth',2,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
hold on
plot(t4,sol4(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
legend('Gate 6','Gate 7')
title('Add inpulse 20 for 400s')
p4 = zeros(size(t4))+20;

% 5. turn off p, p=0, for 570s
p = 0;
u0_t4 = sol4(end,:);
[t5,sol5] = ode45(@(t,x) plasmid(t,x,p),[0 570],u0_t4);
figure
plot(t5,sol5(:,6),'LineWidth',2,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
hold on
plot(t5,sol5(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
legend('Gate 6','Gate 7')
title('Relax system for 570s')
p5 = zeros(size(t5));


%% Final plot
P_set = [p1;p2;p3;p4;p5];

t2n = t2+t1(end);
t3n = t3 + t2(end) + t1(end);
t4n = t4 + t3(end) + t2(end) + t1(end);
t5n = t5 + t4(end) + t3(end) + t2(end) + t1(end);

t_set = [t1;t2n;t3n;t4n;t5n];
sol_set =[sol1;sol2;sol3;sol4;sol5];

figure
plot(t_set,sol_set(:,6),'LineWidth',1,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
hold on;
plot(t_set,sol_set(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
hold on;
plot(t_set, P_set,"-.",'LineWidth',2,"Color",'#7E2F8E')
legend('Gate 6','Gate 7')
ylim([0 22])
title('Circuit Simulation')



%% For loop simulating more switching cycles
% set P duration= 300, and t_wait = 900
% Initialize vector of the final plot
P_set = [];
sol_set =[];
t_set = [];

% 1. make system goes to steady state
p=0
rng(s)
u0= randi([1 20],1,7);
[t1,sol_steady] = ode45(@(t,x) plasmid(t,x,p),[0 1e3],u0);
% figure
% plot(t1,sol_steady(:,6),'LineWidth',2,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
% hold on
% plot(t1,sol_steady(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
% hold on
% plot(t1,sol_steady(:,8),'LineWidth',2,"Marker","^","Color",'#0072BD','DisplayName','Gate8')
% legend('Gate 6','Gate 7','Gate8')
% title('System goes to steady states')
p1 = zeros(size(t1));
t_set = [t_set, t1];
% Set the time parameters
Time.p = 400;Time.wait = 1e3;
mag = 20
for cycle = 1:5
    %1. Add inpulse p=20 for 300s
    p = mag;
    if cycle ==1
        u0_t1 = sol_steady(end,:);
    else
        u0_t1 = sol_relax(end,:);
    end
    [t_inpulse,sol_inpulse] = ode45(@(t,x) plasmid(t,x,p),[0 Time.p],u0_t1);
%     figure
%     plot(t_inpulse,sol_inpulse(:,6),'LineWidth',2,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
%     hold on
%     plot(t_inpulse,sol_inpulse(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
%     hold on
%     plot(t_inpulse,sol_inpulse(:,8),'LineWidth',2,"Marker","^","Color",'#0072BD','DisplayName','Gate8')
%     legend('Gate 6','Gate 7','Gate8')
%     title('Add inpulse 20 for 300s')
    p_on = zeros(size(t_inpulse))+mag;
    
    if cycle ==1
        t_inpulse_n1 = t_inpulse + t1(end);% new t_inpulse 
        t_set = [t_set; t_inpulse_n1];
    else
        t_inpulse_n1 = t_inpulse + t_relax_n2(end);
        t_set = [t_set; t_inpulse_n1];
    end
    
    
    % 2. turn off p for 900s
    p = 0;
    u0_t2 = sol_inpulse(end,:);
    [t_relax,sol_relax] = ode45(@(t,x) plasmid(t,x,p),[0 Time.wait],u0_t2);
%     figure
%     plot(t_relax,sol_relax(:,6),'LineWidth',2,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
%     hold on
%     plot(t_relax,sol_relax(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
%     hold on
%     plot(t_relax,sol_relax(:,8),'LineWidth',2,"Marker","^","Color",'#0072BD','DisplayName','Gate8')
%     legend('Gate 6','Gate 7','Gate8')
%     title('Relax system for 900s')
    p_off = zeros(size(t_relax));
    t_relax_n2 = t_relax + t_inpulse_n1(end);
    t_set = [t_set; t_relax_n2];
    
    
    % Collect data vector for final plots
    % 1. inpulse duration vector 
    P_set = [P_set; p_on];
    P_set = [P_set; p_off];
    
    % 2. solution vector
    sol_set = [sol_set; sol_inpulse];
    sol_set = [sol_set; sol_relax];
   
    
end

% Collect 

t_set_f = t_set;
p_set_f = [p1;P_set];
sol_set_f = [sol_steady; sol_set];


% final plot for all cycles
figure
plot(t_set_f,sol_set_f(:,6),'--','LineWidth',1,"Marker",".","Color",'#77AC30','DisplayName','Gate6')
hold on;
plot(t_set_f,sol_set_f(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
hold on;

plot(t_set_f, p_set_f,"-.",'LineWidth',2,"Color",'#7E2F8E')
legend('Gate 6','Gate 7','Input')
ylim([0 22])
xlabel('Time')
ylabel('Concentration')
title('Circuit Simulation')
ax = gca; % current axes
ax.FontSize = 20;
ax.FontName = "Apple Symbols";






%%  * 6--> 2 (stable)
syms t x p
plasmid62 = @(t,x,p)[xi*(HlyIIR.Y_min + (HlyIIR.Y_max - HlyIIR.Y_min)*HlyIIR.K^HlyIIR.n/(HlyIIR.K^HlyIIR.n + (x(7) + p)^HlyIIR.n)) - gma*(x(1));
                   xi*(LmrA.Y_min + (LmrA.Y_max - LmrA.Y_min)*LmrA.K^LmrA.n/(LmrA.K^LmrA.n + (x(1) + x(7))^LmrA.n)) - gma*(x(2))
                   xi*(SrpR.Y_min + (SrpR.Y_max - SrpR.Y_min)*SrpR.K^SrpR.n/(SrpR.K^SrpR.n + (x(1) + p)^SrpR.n)) - gma*(x(3))
                   xi*(BM3R1.Y_min + (BM3R1.Y_max - BM3R1.Y_min)*BM3R1.K^BM3R1.n/(BM3R1.K^BM3R1.n + (x(2) + x(3))^BM3R1.n)) - gma*(x(4))
                   xi*(PhIF.Y_min + (PhIF.Y_max - PhIF.Y_min)*PhIF.K^PhIF.n/(PhIF.K^PhIF.n + (x(4))^PhIF.n)) - gma*(x(5))
                   xi*(PsrA.Y_min + (PsrA.Y_max - PsrA.Y_min)*PsrA.K^PsrA.n/(PsrA.K^PsrA.n + (x(5) + x(7))^PsrA.n)) - gma*(x(6))
                   xi*(BetI.Y_min + (BetI.Y_max - BetI.Y_min)*BetI.K^BetI.n/(BetI.K^BetI.n + (x(4) + x(6))^BetI.n)) - gma*(x(7))];




figure             
for rd = 1:50
    p=0;
    u0= randi([1 20],1,7);
    [t,xp] = ode45(@(t,x) plasmid62(t,x,p),[0 5000.5],u0);
%     figure
%     plot(t,xp(:,6:7),'LineWidth',2)
    plot(t,xp(:,6),'LineWidth',2,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
%     legend('x6')
    hold on
    plot(t,xp(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')

%     legend('x7')
    title("stable state")
end


%% random seed
s = rng;

%% Sequential simulation

% 1. Simulation the system untill it goes to steady state
p=0;
rng(s)
u0= randi([1 20],1,7);
[t,xp] = ode45(@(t,x) plasmid62(t,x,p),[0 1e4],u0);
figure
plot(t,xp(:,6),'LineWidth',2,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
hold on
plot(t,xp(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
legend('Gate 6','Gate 7')

p = 20;
u0_t1 = xp(end,:);
[t,xp] = ode45(@(t,x) plasmid62(t,x,p),[0 3000.5],u0_t1);
figure
plot(t,xp(:,6),'LineWidth',2,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
hold on
plot(t,xp(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
legend('Gate 6','Gate 7')
xlabel('Time')
ylabel('Concentration')
title('System continues switching after input is turned on')


%% Find the lb and ub of the P inpluse and the min time to wait for the next signal

% Find the local minimum index
locs2 =  find(islocalmin(xp(:,6)));
ind_min1 = locs2(2,1);

% Test the time starting from ti =1 to ti = first local min
t_wait = [];
ul_range =[];
t_range = [];
for ti = 1: ind_min1
    p = 0
    u0_t11 = xp(ti,:)
    fprintf('The duration of P is %d:', t(ti));
    [tt,xp1] = ode45(@(t,x) plasmid62(t,x,p),[0 3000.5],u0_t11);
    % How much time to wait until the system become statble again
    [~,locs]  =findpeaks(xp1(:,6));
    stable_t_ind = locs(1)
    
    if xp1(end, 7) > xp(end, 6)%t(ti) > 306 & t(ti) < 703% lb ; ub
        ul_range = [ul_range, ti];
        t_range = [t_range, t(ti)];
        t_wait = [t_wait, tt(stable_t_ind)];
    end
    
%     figure
%     plot(tt,xp1(:,6),'LineWidth',2,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
%     hold on
%     plot(tt,xp1(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
%     legend('Gate 6','Gate 7')
end
max(t_wait)% give the min time to wait(lower bound) for the next signal within the switching range
ul_range
lb = t_range(1)
ub = t_range(end)

%% single example from above
u0_t11 = xp1(52,:);
p=0
[tt,xp1] = ode45(@(t,x) plasmid(t,x,p),[0 1.5e3],u0_t11);

figure
plot(tt,xp1(:,6),'LineWidth',2,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
hold on
plot(tt,xp1(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')

legend('Gate 6','Gate 7')



%% For loop simulating more switching cycles
% set P duration= 300, and t_wait = 900
% Initialize vector of the final plot
P_set = [];
sol_set =[];
t_set = [];

% 1. make system goes to steady state
p=0
rng(s)
u0= randi([1 20],1,7);
[t1,sol_steady] = ode45(@(t,x) plasmid62(t,x,p),[0 1.5e3],u0);
% figure
% plot(t1,sol_steady(:,6),'LineWidth',2,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
% hold on
% plot(t1,sol_steady(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
% hold on
% plot(t1,sol_steady(:,8),'LineWidth',2,"Marker","^","Color",'#0072BD','DisplayName','Gate8')
% legend('Gate 6','Gate 7','Gate8')
% title('System goes to steady states')
p1 = zeros(size(t1));
t_set = [t_set, t1];
% Set the time parameters
Time.p = 600;Time.wait = 1000;
mag = 20
for cycle = 1:5
    %1. Add inpulse p=20 for 300s
    p = mag;
    if cycle ==1
        u0_t1 = sol_steady(end,:);
    else
        u0_t1 = sol_relax(end,:);
    end
    [t_inpulse,sol_inpulse] = ode45(@(t,x) plasmid62(t,x,p),[0 Time.p],u0_t1);
%     figure
%     plot(t_inpulse,sol_inpulse(:,6),'LineWidth',2,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
%     hold on
%     plot(t_inpulse,sol_inpulse(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
%     hold on
%     plot(t_inpulse,sol_inpulse(:,8),'LineWidth',2,"Marker","^","Color",'#0072BD','DisplayName','Gate8')
%     legend('Gate 6','Gate 7','Gate8')
%     title('Add inpulse 20 for 300s')
    p_on = zeros(size(t_inpulse))+mag;
    
    if cycle ==1
        t_inpulse_n1 = t_inpulse + t1(end);% new t_inpulse 
        t_set = [t_set; t_inpulse_n1];
    else
        t_inpulse_n1 = t_inpulse + t_relax_n2(end);
        t_set = [t_set; t_inpulse_n1];
    end
    
    
    % 2. turn off p for 900s
    p = 0;
    u0_t2 = sol_inpulse(end,:);
    [t_relax,sol_relax] = ode45(@(t,x) plasmid62(t,x,p),[0 Time.wait],u0_t2);
%     figure
%     plot(t_relax,sol_relax(:,6),'LineWidth',2,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
%     hold on
%     plot(t_relax,sol_relax(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
%     hold on
%     plot(t_relax,sol_relax(:,8),'LineWidth',2,"Marker","^","Color",'#0072BD','DisplayName','Gate8')
%     legend('Gate 6','Gate 7','Gate8')
%     title('Relax system for 900s')
    p_off = zeros(size(t_relax));
    t_relax_n2 = t_relax + t_inpulse_n1(end);
    t_set = [t_set; t_relax_n2];
    
    
    % Collect data vector for final plots
    % 1. inpulse duration vector 
    P_set = [P_set; p_on];
    P_set = [P_set; p_off];
    
    % 2. solution vector
    sol_set = [sol_set; sol_inpulse];
    sol_set = [sol_set; sol_relax];
   
    
end

% Collect 

t_set_f = t_set;
p_set_f = [p1;P_set];
sol_set_f = [sol_steady; sol_set];


% final plot for all cycles
figure
plot(t_set_f,sol_set_f(:,6),'LineWidth',1,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
hold on;
plot(t_set_f,sol_set_f(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
hold on;

plot(t_set_f, p_set_f,"-.",'LineWidth',2,"Color",'#7E2F8E')
legend('Gate 6','Gate 7','Input')
ylim([0 22])
title('Circuit Simulation : 8-->2')
%% -----------------------------------------------
%% ***** || Delay Case || ******
%% -----------------------------------------------



%% **Case 1:  8(AmtR) --> NOT -->2(pTet) 
% T
AmtR = table2struct(T("AmtR_A1",para_name));
PsrA = table2struct(T("PsrA_R1",para_name));
% AmeR = table2struct(T("AmeR_F1",para_name));
HlyIIR = table2struct(T("HlyIIR_H1",para_name));
BetI   = table2struct(T("BetI_E1",para_name));
LmrA   = table2struct(T("LmrA_N1",para_name));
BM3R1  = table2struct(T("BM3R1_B3",para_name));
PhIF   = table2struct(T("PhIF_P3",para_name));
SrpR   = table2struct(T("SrpR_S4",para_name));


plasmid8n2 = @(t,x,p)[xi*(AmtR.Y_min + (AmtR.Y_max - AmtR.Y_min)*AmtR.K^AmtR.n/(AmtR.K^AmtR.n + (x(8) + p)^AmtR.n)) - gma*(x(1));
                      xi*(PsrA.Y_min + (PsrA.Y_max - PsrA.Y_min)*PsrA.K^PsrA.n/(PsrA.K^PsrA.n + (x(1) + p )^PsrA.n)) - gma*(x(2));
                      xi*(HlyIIR.Y_min + (HlyIIR.Y_max - HlyIIR.Y_min)*HlyIIR.K^HlyIIR.n/(HlyIIR.K^HlyIIR.n + (x(8) + x(1))^HlyIIR.n)) - gma*(x(3));
                      xi*(BetI.Y_min + (BetI.Y_max - BetI.Y_min)*BetI.K^BetI.n/(BetI.K^BetI.n + (x(2) + x(3))^BetI.n)) - gma*(x(4));
                      xi*(LmrA.Y_min + (LmrA.Y_max - LmrA.Y_min)*LmrA.K^LmrA.n/(LmrA.K^LmrA.n + (x(4))^LmrA.n)) - gma*(x(5));
                      xi*(BM3R1.Y_min + (BM3R1.Y_max - BM3R1.Y_min)*BM3R1.K^BM3R1.n/(BM3R1.K^BM3R1.n + (x(4) + x(7))^BM3R1.n)) - gma*(x(6));
                      xi*(PhIF.Y_min + (PhIF.Y_max - PhIF.Y_min)*PhIF.K^PhIF.n/(PhIF.K^PhIF.n + (x(5) + x(6) )^PhIF.n)) - gma*(x(7));
                      xi*(SrpR.Y_min + (SrpR.Y_max - SrpR.Y_min)*SrpR.K^SrpR.n/(SrpR.K^SrpR.n + (x(7))^SrpR.n)) - gma*(x(8));];
                     
figure             
for rd = 1:50
    p=0;
    u0= randi([1 20],1,8);
    [t,xp] = ode45(@(t,x) plasmid8n2(t,x,p),[0 1e3],u0);
%     figure
%     plot(t,xp(:,6:7),'LineWidth',2)
    plot(t,xp(:,6),'LineWidth',2,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
%     legend('x6')
    hold on
    plot(t,xp(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
    hold on
    plot(t,xp(:,8),'LineWidth',2,"Marker","^","Color",'#0072BD','DisplayName','Gate8')
    legend('Gate 6','Gate 7','Gate8')
    title("stable state")
end

%% random seed
s = rng;
%% Sequential simulation

% 1. Simulation the system untill it goes to steady state
p=0;
rng(s)
u0= randi([1 20],1,8);
[t,xp] = ode45(@(t,x) plasmid8n2(t,x,p),[0 1e4],u0);
figure
plot(t,xp(:,6),'LineWidth',2,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
hold on
plot(t,xp(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
hold on
plot(t,xp(:,8),'LineWidth',2,"Marker","^","Color",'#0072BD','DisplayName','Gate8')
legend('Gate 6','Gate 7','Gate8')

p = 20;
u0_t1 = xp(end,:);
[t,xp] = ode45(@(t,x) plasmid8n2(t,x,p),[0 1.5e3],u0_t1);
figure
plot(t,xp(:,6),'LineWidth',2,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
hold on
plot(t,xp(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
hold on
plot(t,xp(:,8),'LineWidth',2,"Marker","^","Color",'#0072BD','DisplayName','Gate8')
legend('Gate 6','Gate 7','Gate8')
xlabel('Time')
ylabel('Concentration')
title('System continues switching after input is turned on')

%% Find the lb and ub of the P inpluse and the min time to wait for the next signal

% Find the local minimum index
locs2 =  find(islocalmin(xp(:,7)));
ind_min1 = locs2(2,1);

% Test the time starting from ti =1 to ti = first local min
t_wait = [];
ul_range =[];
t_range = [];
for ti = 1: ind_min1
    p = 0
    u0_t11 = xp(ti,:)
    fprintf('The duration of P is %d:', t(ti));
    [tt,xp1] = ode45(@(t,x) plasmid8n2(t,x,p),[0 1.5e3],u0_t11);
    % How much time to wait until the system become statble again
    [~,locs]  =findpeaks(xp1(:,7));
    stable_t_ind = locs(1);
    if xp1(end, 7) > xp(end, 6)%t(ti) > 306 & t(ti) < 703% lb ; ub
        ul_range = [ul_range, ti];
        t_range = [t_range, t(ti)];
        t_wait = [t_wait, tt(stable_t_ind)];
    end
    
%     figure
%     plot(tt,xp1(:,6),'LineWidth',2,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
%     hold on
%     plot(tt,xp1(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
%     hold on
%     plot(tt,xp1(:,8),'LineWidth',2,"Marker","^","Color",'#0072BD','DisplayName','Gate8')
%     legend('Gate 6','Gate 7','Gate8')

end
max(t_wait)% give the min time to wait(lower bound) for the next signal within the switching range

% switching index range
ul_range
% switching lower bound time
t_lb = t_range(1)
% switching upper bound time
t_ub = t_range(end)

%% single example from above
u0_t11 = xp(94,:);
p=0
[tt,xp1] = ode45(@(t,x) plasmid8n2(t,x,p),[0 1.5e3],u0_t11);

figure
plot(tt,xp1(:,6),'LineWidth',2,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
hold on
plot(tt,xp1(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
hold on
plot(tt,xp1(:,8),'LineWidth',2,"Marker","^","Color",'#0072BD','DisplayName','Gate8')
legend('Gate 6','Gate 7','Gate8')


%% Now test the square wave inpulse (need to finish this section)

% set P duration= 300, and t_wait = 900
% Initialize vector of the final plot
P_set = [];
sol_set =[];


% 1. make system goes to steady state
p=0
rng(s)
u0= randi([1 20],1,8);
[t1,sol1] = ode45(@(t,x) plasmid8n2(t,x,p),[0 1500],u0);
figure
plot(t1,sol1(:,6),'LineWidth',2,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
hold on
plot(t1,sol1(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
hold on
plot(t1,sol1(:,8),'LineWidth',2,"Marker","^","Color",'#0072BD','DisplayName','Gate8')
legend('Gate 6','Gate 7','Gate8')
title('System goes to steady states')
p1 = zeros(size(t1));

% 2. add inpulse p=20 for 300s
p = 20;
u0_t1 = sol1(end,:);
[t2,sol2] = ode45(@(t,x) plasmid8n2(t,x,p),[0 300],u0_t1);
figure
plot(t2,sol2(:,6),'LineWidth',2,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
hold on
plot(t2,sol2(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
hold on
plot(t2,sol2(:,8),'LineWidth',2,"Marker","^","Color",'#0072BD','DisplayName','Gate8')
legend('Gate 6','Gate 7','Gate8')
title('Add inpulse 20 for 300s')
p2 = zeros(size(t2))+20;

% 3. turn off p, p=0, for 900s
p = 0;
u0_t2 = sol2(end,:);
[t3,sol3] = ode45(@(t,x) plasmid8n2(t,x,p),[0 900],u0_t2);
figure
plot(t3,sol3(:,6),'LineWidth',2,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
hold on
plot(t3,sol3(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
hold on
plot(t3,sol3(:,8),'LineWidth',2,"Marker","^","Color",'#0072BD','DisplayName','Gate8')
legend('Gate 6','Gate 7','Gate8')
title('Relax system for 900s')
p3 = zeros(size(t3));

% % repeat 2 & 3

% 4. add inpulse p=20 for 300s
p = 20;
u0_t3 = sol3(end,:);
[t4,sol4] = ode45(@(t,x) plasmid8n2(t,x,p),[0 300],u0_t3);
figure
plot(t4,sol4(:,6),'LineWidth',2,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
hold on
plot(t4,sol4(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
hold on
plot(t4,sol4(:,8),'LineWidth',2,"Marker","^","Color",'#0072BD','DisplayName','Gate8')
legend('Gate 6','Gate 7','Gate8')
title('Add inpulse 20 for 300s')
p4 = zeros(size(t4))+20;

% 5. turn off p, p=0, for 900s
p = 0;
u0_t4 = sol4(end,:);
[t5,sol5] = ode45(@(t,x) plasmid8n2(t,x,p),[0 900],u0_t4);
figure
plot(t5,sol5(:,6),'LineWidth',2,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
hold on
plot(t5,sol5(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
hold on
plot(t5,sol5(:,8),'LineWidth',2,"Marker","^","Color",'#0072BD','DisplayName','Gate8')
legend('Gate 6','Gate 7','Gate8')
title('Relax system for 900s')
p5 = zeros(size(t5));



%% Final plot
P_set = [p1;p2;p3;p4;p5];

t2n = t2+t1(end);
t3n = t3 + t2(end) + t1(end);
t4n = t4 + t3(end) + t2(end) + t1(end);
t5n = t5 + t4(end) + t3(end) + t2(end) + t1(end);

t_set = [t1;t2n;t3n;t4n;t5n];
sol_set =[sol1;sol2;sol3;sol4;sol5];

figure
plot(t_set,sol_set(:,6),'LineWidth',1,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
hold on;
plot(t_set,sol_set(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
hold on;
plot(t_set,sol_set(:,8),'LineWidth',2,"Marker","^","Color",'#0072BD','DisplayName','Gate8')
hold on;
plot(t_set, P_set,"-.",'LineWidth',2,"Color",'#7E2F8E')
legend('Gate 6','Gate 7','Gate 8', 'Input')
ylim([0 22])
title('Circuit Simulation')






%% For loop simulating more switching cycles
% set P duration= 300, and t_wait = 900
% Initialize vector of the final plot
P_set = [];
sol_set =[];
t_set = [];

% 1. make system goes to steady state
p=0
rng(s)
u0= randi([1 20],1,8);
[t1,sol_steady] = ode45(@(t,x) plasmid8n2(t,x,p),[0 1500],u0);
% figure
% plot(t1,sol_steady(:,6),'LineWidth',2,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
% hold on
% plot(t1,sol_steady(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
% hold on
% plot(t1,sol_steady(:,8),'LineWidth',2,"Marker","^","Color",'#0072BD','DisplayName','Gate8')
% legend('Gate 6','Gate 7','Gate8')
% title('System goes to steady states')
p1 = zeros(size(t1));
t_set = [t_set, t1];
% Set the time parameters
Time.p = 300;Time.wait = 1000;

for cycle = 1:5
    %1. Add inpulse p=20 for 300s
    p = 20;
    if cycle ==1
        u0_t1 = sol_steady(end,:);
    else
        u0_t1 = sol_relax(end,:);
    end
    [t_inpulse,sol_inpulse] = ode45(@(t,x) plasmid8n2(t,x,p),[0 Time.p],u0_t1);
%     figure
%     plot(t_inpulse,sol_inpulse(:,6),'LineWidth',2,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
%     hold on
%     plot(t_inpulse,sol_inpulse(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
%     hold on
%     plot(t_inpulse,sol_inpulse(:,8),'LineWidth',2,"Marker","^","Color",'#0072BD','DisplayName','Gate8')
%     legend('Gate 6','Gate 7','Gate8')
%     title('Add inpulse 20 for 300s')
    p_on = zeros(size(t_inpulse))+20;
    
    if cycle ==1
        t_inpulse_n1 = t_inpulse + t1(end);% new t_inpulse 
        t_set = [t_set; t_inpulse_n1];
    else
        t_inpulse_n1 = t_inpulse + t_relax_n2(end);
        t_set = [t_set; t_inpulse_n1];
    end
    
    
    % 2. turn off p for 900s
    p = 0;
    u0_t2 = sol_inpulse(end,:);
    [t_relax,sol_relax] = ode45(@(t,x) plasmid8n2(t,x,p),[0 Time.wait],u0_t2);
%     figure
%     plot(t_relax,sol_relax(:,6),'LineWidth',2,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
%     hold on
%     plot(t_relax,sol_relax(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
%     hold on
%     plot(t_relax,sol_relax(:,8),'LineWidth',2,"Marker","^","Color",'#0072BD','DisplayName','Gate8')
%     legend('Gate 6','Gate 7','Gate8')
%     title('Relax system for 900s')
    p_off = zeros(size(t_relax));
    t_relax_n2 = t_relax + t_inpulse_n1(end);
    t_set = [t_set; t_relax_n2];
    
    
    % Collect data vector for final plots
    % 1. inpulse duration vector 
    P_set = [P_set; p_on];
    P_set = [P_set; p_off];
    
    % 2. solution vector
    sol_set = [sol_set; sol_inpulse];
    sol_set = [sol_set; sol_relax];
   
    
end

% Collect 

t_set_f = t_set;
p_set_f = [p1;P_set];
sol_set_f = [sol_steady; sol_set];


% final plot for all cycles
figure
plot(t_set_f,sol_set_f(:,6),'LineWidth',1,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
hold on;
plot(t_set_f,sol_set_f(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
hold on;
plot(t_set_f,sol_set_f(:,8),'LineWidth',2,"Marker","^","Color",'#0072BD','DisplayName','Gate8')
hold on;
plot(t_set_f, p_set_f,"-.",'LineWidth',2,"Color",'#7E2F8E')
legend('Gate 6','Gate 7','Gate 8', 'Input')
ylim([0 22])
title('Circuit Simulation : 8-->2')
xlabel('Time')
ylabel('Concentration')
% title('Circuit Simulation with delay( With an extra NOT gate )')
ax = gca; % current axes
ax.FontSize = 20;
ax.FontName = "Apple Symbols";


%% **Case 2:  8(AmtR) --> NOT -->1(pBAD) 
% T
AmeR = table2struct(T("AmeR_F1",para_name));
PhIF   = table2struct(T("PhIF_P3",para_name));
BetI   = table2struct(T("BetI_E1",para_name));
HlyIIR = table2struct(T("HlyIIR_H1",para_name));
LmrA   = table2struct(T("LmrA_N1",para_name));
BM3R1  = table2struct(T("BM3R1_B3",para_name));
SrpR   = table2struct(T("SrpR_S4",para_name));
AmtR = table2struct(T("AmtR_A1",para_name));

plasmid8n1 = @(t,x,p)[xi*(AmtR.Y_min + (AmtR.Y_max - AmtR.Y_min)*AmtR.K^AmtR.n/(AmtR.K^AmtR.n + (x(8) + p)^AmtR.n)) - gma*(x(1));
                      xi*(PsrA.Y_min + (PsrA.Y_max - PsrA.Y_min)*PsrA.K^PsrA.n/(PsrA.K^PsrA.n + (x(1) + x(8) )^PsrA.n)) - gma*(x(2));
                      xi*(HlyIIR.Y_min + (HlyIIR.Y_max - HlyIIR.Y_min)*HlyIIR.K^HlyIIR.n/(HlyIIR.K^HlyIIR.n + (p + x(1))^HlyIIR.n)) - gma*(x(3));
                      xi*(BetI.Y_min + (BetI.Y_max - BetI.Y_min)*BetI.K^BetI.n/(BetI.K^BetI.n + (x(2) + x(3))^BetI.n)) - gma*(x(4));
                      xi*(LmrA.Y_min + (LmrA.Y_max - LmrA.Y_min)*LmrA.K^LmrA.n/(LmrA.K^LmrA.n + (x(4))^LmrA.n)) - gma*(x(5));
                      xi*(BM3R1.Y_min + (BM3R1.Y_max - BM3R1.Y_min)*BM3R1.K^BM3R1.n/(BM3R1.K^BM3R1.n + (x(4) + x(7))^BM3R1.n)) - gma*(x(6));
                      xi*(PhIF.Y_min + (PhIF.Y_max - PhIF.Y_min)*PhIF.K^PhIF.n/(PhIF.K^PhIF.n + (x(5) + x(6) )^PhIF.n)) - gma*(x(7));
                      xi*(SrpR.Y_min + (SrpR.Y_max - SrpR.Y_min)*SrpR.K^SrpR.n/(SrpR.K^SrpR.n + (x(7))^SrpR.n)) - gma*(x(8));];
                     

figure             
for rd = 1:50
    p=0;
    u0= randi([1 20],1,8);
    [t,xp] = ode45(@(t,x) plasmid8n1(t,x,p),[0 1e3],u0);
%     figure
%     plot(t,xp(:,6:7),'LineWidth',2)
    plot(t,xp(:,6),'LineWidth',2,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
%     legend('x6')
    hold on
    plot(t,xp(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
    hold on
    plot(t,xp(:,8),'LineWidth',2,"Marker","^","Color",'#0072BD','DisplayName','Gate8')
    legend('Gate 6','Gate 7','Gate8')
    title("stable state")
end


%% Sequential simulation

% 1. Simulation the system untill it goes to steady state
p=0;
rng(s)
u0= randi([1 20],1,8);
[t,xp] = ode45(@(t,x) plasmid8n1(t,x,p),[0 1e4],u0);
figure
plot(t,xp(:,6),'LineWidth',2,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
hold on
plot(t,xp(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
hold on
plot(t,xp(:,8),'LineWidth',2,"Marker","^","Color",'#0072BD','DisplayName','Gate8')
legend('Gate 6','Gate 7','Gate8')

p = 20;
u0_t1 = xp(end,:);
[t,xp] = ode45(@(t,x) plasmid8n1(t,x,p),[0 1.5e3],u0_t1);
figure
plot(t,xp(:,6),'LineWidth',2,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
hold on
plot(t,xp(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
hold on
plot(t,xp(:,8),'LineWidth',2,"Marker","^","Color",'#0072BD','DisplayName','Gate8')
legend('Gate 6','Gate 7','Gate8')
xlabel('Time')
ylabel('Concentration')
title('System continues switching after input is turned on')


%% Find the lb and ub of the P inpluse and the min time to wait for the next signal

% Find the local minimum index
locs2 =  find(islocalmin(xp(:,7)));
ind_min1 = locs2(2,1);

% Test the time starting from ti =1 to ti = first local min
t_wait = [];
ul_range =[];
t_range = [];
for ti = 1: ind_min1
    p = 0
    u0_t11 = xp(ti,:)
    fprintf('The duration of P is %d:', t(ti));
    [tt,xp1] = ode45(@(t,x) plasmid8n1(t,x,p),[0 1.5e3],u0_t11);
    % How much time to wait until the system become statble again
    [~,locs]  =findpeaks(xp1(:,7));
    stable_t_ind = locs(1);
    if xp1(end, 7) > xp(end, 6)
        ul_range = [ul_range, ti];
        t_range = [t_range, t(ti)];
        t_wait = [t_wait, tt(stable_t_ind)];
    end
    
%     figure
%     plot(tt,xp1(:,6),'LineWidth',2,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
%     hold on
%     plot(tt,xp1(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
%     hold on
%     plot(tt,xp1(:,8),'LineWidth',2,"Marker","^","Color",'#0072BD','DisplayName','Gate8')
%     legend('Gate 6','Gate 7','Gate8')

end
max(t_wait)% give the min time to wait(lower bound) for the next signal within the switching range

% switching index range
ul_range
% switching lower bound time
t_lb = t_range(1)
% switching upper bound time
t_ub = t_range(end)



%% single example from above
u0_t11 = xp(75,:);
p=0
[tt,xp1] = ode45(@(t,x) plasmid8n1(t,x,p),[0 1.5e3],u0_t11);

figure
plot(tt,xp1(:,6),'LineWidth',2,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
hold on
plot(tt,xp1(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
hold on
plot(tt,xp1(:,8),'LineWidth',2,"Marker","^","Color",'#0072BD','DisplayName','Gate8')
legend('Gate 6','Gate 7','Gate8')


%% For loop simulating more switching cycles
% set P duration= 300, and t_wait = 900
% Initialize vector of the final plot
P_set = [];
sol_set =[];
t_set = [];

% 1. make system goes to steady state
p=0
rng(s)
u0= randi([1 20],1,8);
[t1,sol_steady] = ode45(@(t,x) plasmid8n1(t,x,p),[0 1500],u0);
p1 = zeros(size(t1));
t_set = [t_set, t1];
% Set the time parameters
Time.p = 300;Time.wait = 1000;

for cycle = 1:5
    %1. Add inpulse p=20 for 220s
    p = 20;
    if cycle ==1
        u0_t1 = sol_steady(end,:);
    else
        u0_t1 = sol_relax(end,:);
    end
    [t_inpulse,sol_inpulse] = ode45(@(t,x) plasmid8n1(t,x,p),[0 Time.p],u0_t1);
    p_on = zeros(size(t_inpulse))+20;
    
    if cycle ==1
        t_inpulse_n1 = t_inpulse + t1(end);% new t_inpulse 
        t_set = [t_set; t_inpulse_n1];
    else
        t_inpulse_n1 = t_inpulse + t_relax_n2(end);
        t_set = [t_set; t_inpulse_n1];
    end
    
    
    % 2. turn off p for 850s
    p = 0;
    u0_t2 = sol_inpulse(end,:);
    [t_relax,sol_relax] = ode45(@(t,x) plasmid8n1(t,x,p),[0 Time.wait],u0_t2);
    p_off = zeros(size(t_relax));
    t_relax_n2 = t_relax + t_inpulse_n1(end);
    t_set = [t_set; t_relax_n2];
    
    
    % Collect data vector for final plots
    % 1. inpulse duration vector 
    P_set = [P_set; p_on];
    P_set = [P_set; p_off];
    
    % 2. solution vector
    sol_set = [sol_set; sol_inpulse];
    sol_set = [sol_set; sol_relax];
   
    
end

% Collect 

t_set_f = t_set;
p_set_f = [p1;P_set];
sol_set_f = [sol_steady; sol_set];


% final plot for all cycles
figure
plot(t_set_f,sol_set_f(:,6),'LineWidth',1,"Marker","*","Color",'#77AC30','DisplayName','Gate6')
hold on;
plot(t_set_f,sol_set_f(:,7),'LineWidth',2,"Marker",".","Color",'#D95319','DisplayName','Gate7')
hold on;
plot(t_set_f,sol_set_f(:,8),'LineWidth',2,"Marker","^","Color",'#0072BD','DisplayName','Gate8')
hold on;
plot(t_set_f, p_set_f,"-.",'LineWidth',2,"Color",'#7E2F8E')
legend('Gate 6','Gate 7','Gate 8', 'Input')
ylim([0 22])
title('Circuit Simulation 8 --> 1')


