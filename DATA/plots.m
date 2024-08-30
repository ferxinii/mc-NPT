
kb_j = 1.38e-23 ; %J/K
kb_ev = 8.617333262e-5 %eV/K
eps_u = 119.8 * kb_ev; %eV
eps_p = 119.8 * kb_j; %J
sig = 3.405 ;





%% Sampling info

steps_pressure = 10;
p0 = 0.15;
pf = 15;

pressure_vec = [];
for ii = 1:steps_pressure
    pressure_vec = [pressure_vec p0 + (ii-1)*(pf-p0)/(steps_pressure-1)]
end


%% Pressure

p_mean = [];
p_std = [];
for ii = 1:10
    data = load("pressure"+string(ii)+".dat");
    
    p_mean = [p_mean mean(data(:,5))];
    p_std = [p_std std(data(:,5))];
end

figure;
errorbar(pressure_vec,p_mean,p_std, 'LineWidth',1,"LineStyle","none"); hold on;
plot(pressure_vec,p_mean,'x','LineWidth',2,'MarkerSize',10); hold on;
plot(pressure_vec,pressure_vec,'o','LineWidth',2,'color','black'); hold on;
grid on; xlabel("Reference P*"); ylabel("Mean virial Pressure*");
legend("","Virial P*","Reference P*");
sgtitle("Reference Pressure VS Mean virial Pressure");


p_units = p_mean * eps_p / sig^3; %J/A3 = Nm/A3
p_units = p_units * 1e30 % N/m2= Pa
%p_units = p_units * 9.86923267*1e-6 %atm


%% Energy

u_mean = [];
u_std = [];
for ii = 1:10
    data = load("energy"+string(ii)+".dat");
    
    u_mean = [u_mean mean(data(:,2))];
    u_std = [u_std std(data(:,2))];
end

figure;
errorbar(pressure_vec,u_mean,u_std, 'LineWidth',1,"LineStyle","none"); hold on;
plot(pressure_vec,u_mean,'x','LineWidth',2,'MarkerSize',10); hold on;
grid on; sgtitle("Mean potential energy"); xlabel("Reference P*"); ylabel("<U*>");

u_units = u_mean * eps_p;

%% Boxlength and RHO

L_mean = [];
L_std = [];
rho_mean = [];
rho_std = [];

for ii = 1:10
    data = load("boxlength"+string(ii)+".dat");
    
    L_mean = [L_mean mean(data(:,1))];
    L_std = [L_std std(data(:,1))];

    rho_mean = [rho_mean mean(500./(data(:,1).^3))];
    rho_std = [rho_std std(500./(data(:,1).^3))];
end

figure;
errorbar(pressure_vec,L_mean,L_std, 'LineWidth',1,"LineStyle","none"); hold on;
plot(pressure_vec,L_mean,'x','LineWidth',2,'MarkerSize',10); hold on;
grid on; sgtitle("Mean boxlength"); xlabel("Reference P*"); ylabel("<L*>");

figure;
errorbar(pressure_vec,rho_mean,rho_std, 'LineWidth',1,"LineStyle","none"); hold on;
plot(pressure_vec,rho_mean,'x','LineWidth',2,'MarkerSize',10); hold on;
grid on; sgtitle("Mean density"); xlabel("Reference P*"); ylabel("<rho*>");


rho_units = rho_mean / sig^3;

%% P and U as a function of rho

figure
errorbar(rho_mean,p_mean,p_std, 'LineWidth',1,"LineStyle","none"); hold on;
plot(rho_mean,p_mean,'x','LineWidth',2,'MarkerSize',10); hold on;
grid on; sgtitle("Mean pressure as a function of density"); xlabel("<rho*>"); ylabel("<Pvir*>");

figure;
errorbar(rho_mean,u_mean,u_std, 'LineWidth',1,"LineStyle","none"); hold on;
plot(rho_mean,u_mean,'x','LineWidth',2,'MarkerSize',10); hold on;
grid on; sgtitle("Mean potential energy as a function of density"); xlabel("<rho*>"); ylabel("<U*>");

%% Load P3 to compare
load("../P3/workspace.mat");

%%
figure
errorbar(rho_mean,u_mean,u_std, 'LineWidth',1,"LineStyle","none"); hold on;
plot(rho_mean,u_mean,'x','LineWidth',2,'MarkerSize',10); hold on;
errorbar(rho,energy_mean,energy_dev, 'LineWidth',1,"LineStyle","none",'Color',"#4DBEEE"); hold on;
plot(rho,energy_mean,'x','LineWidth',2,'Color','black'); hold on;

legend("","obtained","","P3 (N,V,T)"); grid on; xlabel("<rho*> or rho*"); ylabel("<U*>");
sgtitle("Comparison of <U*> with P3");

%%
figure
errorbar(rho_mean,p_mean,p_std, 'LineWidth',1,"LineStyle","none"); hold on;
plot(rho_mean,p_mean,'x','LineWidth',2,'MarkerSize',10); hold on;
errorbar(rho,press_mean,press_dev, 'LineWidth',1,"LineStyle","none",'Color',"#4DBEEE"); hold on;
plot(rho,press_mean,'x','LineWidth',2,'Color','black'); hold on;

legend("","obtained","","P3 (N,V,T)"); grid on; xlabel("<rho*> or rho*"); ylabel("<Pvir*>");
sgtitle("Comparison of <P*> with P3");

%% gdr

gdr = [];
r = [];
for ii = 1:10
    data = load("gdr"+string(ii)+".dat");
    
    gdr = [gdr, data(:,2)];
    r = [r,  data(:,1)];
end

%plotting
light = [0.5,0.7,0.9];
dark = [0.8,0.1,0.2];
GRADIENT = @(n,nn) interp1([1/nn 1],[light;dark],n/nn);

figure
for ii = flip(1:10)
    plot(r(1:end-100,ii), gdr(1:end-100,ii),'LineWidth',1.5,'color',GRADIENT(ii,10)); hold on;
end

legend(flip(string(pressure_vec)));
sgtitle("GDR for different reference P*"); grid on;
xlabel("r*"); ylabel("g(r*)");


%% acceptance info and delr
fid = fopen('output_console.dat', 'r');
file_data = textscan(fid, '%f %f %s %f %s %f %s %f %s %f %f');
fclose(fid);

data = [ file_data{1}, file_data{2}, file_data{4}, file_data{6}, file_data{8}, file_data{10}, file_data{11} ];



%% delr

figure
subplot(1,2,1);
for ii = 1:10
mask = data(:,1) == ii;
aux_data = data(mask,:);

mask2 = aux_data(:,7) ~= 2;
    if length(mask2)~=1
        plot(aux_data(1,2), aux_data(mask2,3), 'o', 'color', "#0072BD", 'LineWidth',2); hold on
        plot(aux_data(1,2), (aux_data(~mask2,3)), 'x', 'color', 'black','LineWidth',4,'MarkerSize',10); hold on
    else
        plot(aux_data(1,2), (aux_data(:,3)), 'x', 'color', 'black','LineWidth',4,'MarkerSize',10); hold on
    end
end
xlabel("Reference P*"); ylabel("∆r");
title("∆r"); grid on;

subplot(1,2,2);
for ii = 1:10
mask = data(:,1) == ii;
aux_data = data(mask,:);

mask2 = aux_data(:,7) ~= 2;
    if length(mask2)~=1
        plot(aux_data(1,2), log(aux_data(mask2,3)), 'o', 'color', "#0072BD", 'LineWidth',2); hold on
        plot(aux_data(1,2), log(aux_data(~mask2,3)), 'x', 'color', 'black','LineWidth',4,'MarkerSize',10); hold on
    else
        plot(aux_data(1,2), log(aux_data(:,3)), 'x', 'color', 'black','LineWidth',4,'MarkerSize',10); hold on
    end
end
xlabel("Reference P*"); ylabel("ln(∆r)");
title("ln(∆r)"); grid on;

sgtitle("∆r used (linear & ln scale)");

%% acc r
%{
subplot(2,1,1);
for ii = 1:10
mask = data(:,1) == ii;
aux_data = data(mask,:);

mask2 = aux_data(:,7) ~= 2;
    if length(mask2)~=1
        plot(aux_data(1,2), aux_data(mask2,3), 'o', 'color', "#0072BD", 'LineWidth',2); hold on
        plot(aux_data(1,2), (aux_data(~mask2,3)), 'x', 'color', 'black','LineWidth',4,'MarkerSize',10); hold on
    else
        plot(aux_data(1,2), (aux_data(:,3)), 'x', 'color', 'black','LineWidth',4,'MarkerSize',10); hold on
    end
end
xlabel("Reference P*"); ylabel("∆r");
title("∆r")

subplot(2,1,2);
%}

figure
for ii = 1:10
mask = data(:,1) == ii;
aux_data = data(mask,:);

mask2 = aux_data(:,7) ~= 2;
    if length(mask2)~=1
        plot(aux_data(1,2), aux_data(mask2,4), 'o', 'color', "#0072BD", 'LineWidth',2); hold on
        plot(aux_data(1,2), aux_data(~mask2,4), 'x', 'color', 'black','LineWidth',4,'MarkerSize',10); hold on
    else
        plot(aux_data(1,2), aux_data(:,4), 'x', 'color', 'black','LineWidth',4,'MarkerSize',10); hold on
    end
end

yline(0.39,'Color','#D95319','LineWidth',2);
yline(0.61,'Color','#D95319','LineWidth',2);
ylim([0,1]); xlabel("Reference P*"); ylabel("acc rate ∆ri"); grid on;
sgtitle("Acceptance rate of trial moves in positions");


%% delv

figure;
for ii = 1:10
mask = data(:,1) == ii;
aux_data = data(mask,:);

mask2 = aux_data(:,7) ~= 2;
    if length(mask2)~=1
        plot(aux_data(1,2), aux_data(mask2,5), 'o', 'color', "#0072BD", 'LineWidth',2); hold on
        plot(aux_data(1,2), (aux_data(~mask2,5)), 'x', 'color', 'black','LineWidth',4,'MarkerSize',10); hold on
    else
        plot(aux_data(1,2), (aux_data(:,5)), 'x', 'color', 'black','LineWidth',4,'MarkerSize',10); hold on
    end
end

xlabel("Reference P*"); ylabel("∆V*"); grid on;
sgtitle("∆V used");

%% acc V
figure;
for ii = 1:10
mask = data(:,1) == ii;
aux_data = data(mask,:);

mask2 = aux_data(:,7) ~= 2;
    if length(mask2)~=1
        plot(aux_data(1,2), aux_data(mask2,6), 'o', 'color', "#0072BD", 'LineWidth',2); hold on
        plot(aux_data(1,2), aux_data(~mask2,6), 'x', 'color', 'black','LineWidth',4,'MarkerSize',10); hold on
    else
        plot(aux_data(1,2), aux_data(:,6), 'x', 'color', 'black','LineWidth',4,'MarkerSize',10); hold on
    end
end

yline(0.39,'Color','#D95319','LineWidth',2);
yline(0.61,'Color','#D95319','LineWidth',2);
ylim([0,1]); xlabel("Reference P*"); ylabel("acc rate ∆V"); grid on;
sgtitle("Acceptance rate of trial moves in volume");
