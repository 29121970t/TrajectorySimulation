clear all;
close all; 


N = 1e7;


Bon = 1; 


saveresults = 1;


% Магнетон Бора [Дж/с]
mu_b = 9.27*10^(-24); 

%Масса Ag [кг]
m = 1.7912* 10^(-25);

%Радиус аппертуры печи [м]
oven_inlet_radius =  5.6419e-04;

oven_y = -.0765; 

rslit1_inlet =  3.09e-05; 

slit1y = -5.15e-2;

slt_w = 8e-4;

slt_h = 3.5e-5;

slit2y = -1.8e-2;

slt_z =-1.9e-4 - 30e-6;




r = zeros(N,3);     % Координаты [м]
v = zeros(N,3);     % Скорости [м/с]
a = zeros(N,3);     % Ускорения [м/с^2]

r(:,2) = oven_y;
r(:,1) = oven_inlet_radius/3.*randn(N,1);
r(:,3) = oven_inlet_radius/3.*randn(N,1) + slt_z;


% Инициализация скоростей
alpha = rslit1_inlet/abs(2.5e-2);
theta = (2*rand(N,1)-1)* alpha/2;
phi = rand(N,1)* 2*pi; 
u = [0 slit1y slt_z] - r;
u = u./vecnorm(u')';
u(:,[1 2 3]) = u(:,[1 3 2]); 
ind = u(:,3) > 1-eps;
u_prime(ind,1) = sin(theta(ind)).*cos(phi(ind));
u_prime(ind,2) = sin(theta(ind)).*sin(phi(ind));
u_prime(ind,3) = sign(u(ind,3)).*cos(theta(ind));
ind = ~ind;
u_prime(ind,1) = sin(theta(ind)).*(u(ind,1).*u(ind,3).*cos(phi(ind)) - u(ind,2).*sin(phi(ind)) )./sqrt(1-u(ind,3).^2) + u(ind,1).*cos(theta(ind));
u_prime(ind,2) = sin(theta(ind)).*(u(ind,2).*u(ind,3).*cos(phi(ind)) + u(ind,1).*sin(phi(ind)) )./sqrt(1-u(ind,3).^2) + u(ind,2).*cos(theta(ind));
u_prime(ind,3) = -sqrt(1-u(ind,3).^2).*sin(theta(ind)).*cos(phi(ind)) + u(ind,3).*cos(theta(ind));
u_prime(:,[1 2 3]) = u_prime(:,[1 3 2]);


vstart = 625;                               % [м/с]
vend = 750;                                 % [м/с]
vabs = (vend-vstart)* rand(N,1) + vstart;   
vinit = vabs .* u_prime;                    
v = vinit;                                  


spin = randi(2,N,1)*2-3; 
% spin = rand(N,1)*2-1; 

alive = ones(N,1);

finished = zeros(N,1);


if Bon == 1
comsoldatagridded = importfile('Field.txt');



fieldmeshx = [-0.5:0.01:0.5]'/1000;          % [м]  
fieldmeshy = [-5.1:0.01:5]'/100;            % [м] 
fieldmeshz = [-0.5:0.01:0.5]'/100;          % [м]
meshsize = [length(fieldmeshx) length(fieldmeshy) length(fieldmeshz)];




x_comsol = reshape(comsoldatagridded(:,1),meshsize)/1000;     % [T] 
y_comsol = reshape(comsoldatagridded(:,2),meshsize)/1000;     % [T]
z_comsol = reshape(comsoldatagridded(:,3),meshsize)/1000;     % [T] 

dBzdz = reshape(comsoldatagridded(:,15),meshsize); % [Т/м]
dBxdz = reshape(comsoldatagridded(:,9),meshsize);
dBydz = reshape(comsoldatagridded(:,12),meshsize);


dBzdz_interpolant = griddedInterpolant(x_comsol,y_comsol,z_comsol,dBzdz,'linear','none');
dBxdz_interpolant = griddedInterpolant(x_comsol,y_comsol,z_comsol,dBxdz,'linear','none');
dBydz_interpolant = griddedInterpolant(x_comsol,y_comsol,z_comsol,dBydz,'linear','none');
end



tf = 1e-3;

dt = 1.5133e-6/2;
ts = 0:dt:tf;

if saveresults
    strnow = string(datetime('now','TimeZone','local','Format','y-M-d_HH-mm-ss_z'));
    mkdir(strnow)
    copyfile([mfilename('fullpath'),'.m'],strnow+"/executedscript.m");
end

reverseStr = '';
pause on

rtemp = r;
for it = 1:length(ts)
    t = ts(it);

    a(:,1:3) = 0;

    
    if Bon == 1
        dBxdz_atom = dBxdz_interpolant(rtemp(alive==1,1),rtemp(alive==1,2),rtemp(alive==1,3));
        dBydz_atom = dBydz_interpolant(rtemp(alive==1,1),rtemp(alive==1,2),rtemp(alive==1,3));
        dBzdz_atom = dBzdz_interpolant(rtemp(alive==1,1),rtemp(alive==1,2),rtemp(alive==1,3));
        dBxdz_atom(isnan(dBxdz_atom)) = 0;
        dBydz_atom(isnan(dBydz_atom)) = 0;
        dBzdz_atom(isnan(dBzdz_atom)) = 0;
        
        

        a(alive==1,3) = mu_b*spin(alive==1).*dBzdz_atom/m;
        a(alive==1,2) = mu_b*spin(alive==1).*dBydz_atom/m;
        a(alive==1,1) = mu_b*spin(alive==1).*dBxdz_atom/m;


        a(~(rtemp(:,2)> -.0175 & rtemp(:,2)<.0175),:) = 0;

    end

    v = v + a*dt;
    

    rtemp = r + v*dt;


    
    alive(rtemp(:,2)>slit1y & r(:,2)<slit1y & ...
       (rtemp(:,1)).^2+(rtemp(:,3)-slt_z).^2>rslit1_inlet^2 ...
        ) = 0;

    alive(rtemp(:,2)>slit2y & r(:,2)<slit2y & ...
       (p1(rtemp(:,1),slt_h,slt_w,slt_z)<rtemp(:,3 ) | ...
       p2(rtemp(:,1),slt_h,slt_w,slt_z)>rtemp(:,3)) ...
        ) = 0;


    finished(rtemp(:,2)> .0175) = 1;
    

    r(alive==1 & ~finished,:) = rtemp(alive==1 & ~finished,:);
    

    if size(finished(~finished& alive==1))/size(finished(finished& alive==1))<0.05
        alive(~finished& alive==1) = 0;
        fprintf(['\n']);
        disp("Less than 5% left running");
        break;
    end
    
    msg = sprintf('t = %f us', t*1e6);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    
    if ploton

        if ~exist("fig",'var')
            fig = figure(101);
            figs1 = scatter3(r(alive==1,1),r(alive==1,2),r(alive==1,3),5,'k');
            hold on;
            figs2 = scatter3(r(alive==0,1),r(alive==0,2),r(alive==0,3),1,'r');
            xlabel('X')
            ylabel('Y')
            zlabel('Z')
            hold off;
            title("t=  "+ string(t));
            xlim([-6e-4, 6e-4]);
            ylim([-8e-2, 2e-2]);
            drawnow;
        else
            figs1.XData = r(alive==1,1);
            figs1.YData = r(alive==1,2);
            figs1.ZData = r(alive==1,3);
            figs2.XData = r(alive==0,1);
            figs2.YData = r(alive==0,2);
            figs2.ZData = r(alive==0,3);
            fig.CurrentAxes.Title.String = "t =  "+ num2str(t,3) + ' s';
            drawnow;
        end
    end
end


figure;
scatter(r(alive==1,1).*1000,r(alive==1,3).*1000,.5,'k');
xlabel('x (mm)');
ylabel('z (mm)');
grid on;
axis equal;
box on;
drawnow

if saveresults
    saveas(gcf,strnow+"\fig1.png")
end

figure;
scatter(r(alive==1,1).*1000,r(alive==1,3).*1000,.5,v(alive==1,2));
colormap jet;
xlabel('x (mm)');
ylabel('z (mm)');
grid on;
axis equal;
ylim([-0.45 0.1])
xlim([-0.9 0.9])
box on;
clbr = colorbar;
clbr.Label.String = 'v_y (m/s)';
clbr.Label.FontSize = 12;
drawnow

if saveresults
    saveas(gcf,strnow+"\fig2.png")
    
    aalive = a(alive==1,:);
    ralive = r(alive==1,:);
    valive = v(alive==1,:);
    spinalive = spin(alive==1,:);

    save(strnow+"\workspace" + ".mat",'aalive','ralive','valive','spinalive')
end



function p1 = p1(x,a,b,c)
    p1 = -a/2/(b/2)^2 .* x.^2 + c +a/2;
end

function p2 = p2(x,a,b,c)
    p2 = a/2/(b/2)^2 .* x.^2 + c -a/2;
end 



function ComsolFieldExport = importfile(filename, dataLines)

if nargin < 2
    dataLines = [10, Inf];
end


opts = delimitedTextImportOptions("NumVariables", 28);


opts.DataLines = dataLines;
opts.Delimiter = " ";


opts.VariableNames = ["VarName1", "Version", "COMSOL", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "VarName11", "VarName12", "VarName13", "VarName14", "VarName15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25", "Var26", "Var27", "Var28"];
opts.SelectedVariableNames = ["VarName1", "Version", "COMSOL", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "VarName11", "VarName12", "VarName13", "VarName14", "VarName15"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string"];

opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

opts = setvaropts(opts, ["Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25", "Var26", "Var27", "Var28"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25", "Var26", "Var27", "Var28"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, ["VarName1", "Version", "COMSOL", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "VarName11", "VarName12", "VarName13", "VarName14", "VarName15"], "ThousandsSeparator", ",");


ComsolFieldExport = readtable(filename, opts);


ComsolFieldExport = table2array(ComsolFieldExport);
end


