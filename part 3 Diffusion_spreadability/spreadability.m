clc;clear all;
t=0:0.001:100;
spr1=zeros(size(t)); spr2=zeros(size(t)); spr3=zeros(size(t)); spr4=zeros(size(t)); 
chi5=readmatrix("TChi_k_vein.txt"); val1=zeros(size(chi5));
chi10=readmatrix("TChi_k_Poisson.txt"); val2=zeros(size(chi10));
chi15=readmatrix("TChi_k_RSA.txt"); val3=zeros(size(chi15));
chi20=readmatrix("TChi_k_real-vein.txt"); val4=zeros(size(chi20));

for m = 1:size(t,2)
    for i=1:size(chi5,1)
        for j=1:size(chi5,2)
            val1(i,j)=chi5(i,j)*exp(-((i-151)^2+(j-151)^2)*t(m));
            val2(i,j)=chi10(i,j)*exp(-((i-151)^2+(j-151)^2)*t(m));
            val3(i,j)=chi15(i,j)*exp(-((i-151)^2+(j-151)^2)*t(m));
            val4(i,j)=chi20(i,j)*exp(-((i-151)^2+(j-151)^2)*t(m));
        end
    end

    spr1(m)=sum(val1,"all")*(pi/150)^2;
    spr2(m)=sum(val2,"all")*(pi/150)^2;
    spr3(m)=sum(val3,"all")*(pi/150)^2;
    spr4(m)=sum(val4,"all")*(pi/150)^2;

end

V=0.34;
spr=[spr1;spr2;spr3;spr4];
spr=spr/(4*pi^2*V);


% Create a new figure
figure('Position', [0, 0, 500, 500]);
t_interest =[0.001,10]

% Create the main plot (full view)
loglog(t, spr(1, :), 'g--', t, spr(2, :), 'r:', t, spr(3, :), 'b-.', t, spr(4, :), 'k-','LineWidth',2);
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$S(\infty)-S(t)$', 'Interpreter', 'latex', 'FontSize', 16);
legend('vein', 'Poisson', 'RSA', 'real-vein', 'Interpreter', 'latex', 'FontSize', 10, 'Location', [0.15, 0.15, 0.1, 0.1]);

% Define the coordinates for the local zoomed-in plot
local_x = t(t >= t_interest(1) & t <= t_interest(2));
local_spr1 = spr1(t >= t_interest(1) & t <= t_interest(2));
local_spr2 = spr2(t >= t_interest(1) & t <= t_interest(2));
local_spr3 = spr3(t >= t_interest(1) & t <= t_interest(2));
local_spr4 = spr4(t >= t_interest(1) & t <= t_interest(2));

% Create a zoomed-in plot within the same axes
axes('Position', [0.35, 0.35, 0.35, 0.35]); % Define the position and size of the zoomed-in plot
loglog(local_x, local_spr1, 'g--', local_x, local_spr2, 'r:', local_x, local_spr3, 'b-.', local_x, local_spr4, 'k-');




saveas(gcf,"spr.png");