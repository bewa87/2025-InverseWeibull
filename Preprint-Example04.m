% Step 0: Preparation of Velocity Vector

pkg load statistics;

fid = fopen('ReadData-01803.txt','r');
C   = textscan(fid,'%f32 %f32 %f32 %f32 %f32 %s','delimiter',',');
fclose(fid);

vel_proc       = cell2mat(C(4));
final_vel_proc = sort(vel_proc);
ind_proc       = find(final_vel_proc>=0.09);
final_vel_proc = final_vel_proc(ind_proc);

% Step 1: Data Vector and Its Logarithmic Vector

x_vec    = final_vel_proc;

y_vec    = (log(x_vec))';

z_vec    = x_vec(1)./x_vec;

G_infty  = sum(log(z_vec));

% Step 2: Definition of G, H, dH and dG from our Manuscript

G     = @(beta) (1./beta)-sum(log(x_vec))+(length(x_vec)*(((x_vec).^(-beta))*y_vec))/(sum((x_vec).^(-beta)));
H     = @(beta) length(x_vec)*(((x_vec).^(-beta))*y_vec)/(sum((x_vec).^(-beta)));
dH    = @(beta) length(x_vec)/(sum((x_vec).^(-beta))*sum((x_vec).^(-beta)))*((((x_vec).^(-beta))*y_vec)*(((x_vec).^(-beta))*y_vec)-(((x_vec).^(-beta))*(y_vec.^2))*(sum((x_vec).^(-beta))));
dG    = @(beta) -length(x_vec)/(beta*beta)+dH(beta);

% Step 3: Plotting of G and H

xx    = 0.1:0.1:34.9;
yy    = zeros(1,length(xx));
zz    = zeros(1,length(xx));

for j = 1:1:length(xx)
  yy(1,j) = G(xx(1,j));
  zz(1,j) = H(xx(1,j));
endfor

% Step 4: Newton-Iteration for G(beta_fin) = 0

N_max     = 200;
beta_s    = zeros(1,(N_max+1));
beta_s(1) = 0.1;

for j = 1:1:N_max
  beta_s(j+1) = beta_s(j) - G(beta_s(j))/dG(beta_s(j));
endfor

beta_fin  = beta_s(end)
alpha_fin = ((1/length(x_vec))*sum((x_vec).^(-beta_fin)))^(1/beta_fin)

f         = @(t) (beta_fin/(alpha_fin^(beta_fin))).*(t.^(-(beta_fin+1))).*exp(-(alpha_fin*t).^(-beta_fin));

% Step 5: Plotting of G, H and Relative Frequencies with Inverse Weibull Distribution

% Step 5.1: Plots of G and H

figure(1)
plot(xx,yy)
hold on
plot(xx,G_infty*ones(1,length(xx)),'linestyle','--')
hold on
plot(beta_fin,0,'r*')
title('Plot of function G({\beta})','Interpreter','tex','fontsize',18)
xlabel('{\beta}','Interpreter','tex','fontsize',14)
ylabel('G({\beta})','Interpreter','tex','fontsize',14)
xlim([0 15])

figure(2)
plot(xx,zz)
hold on
plot(xx,length(x_vec)*log(x_vec(1))*ones(1,length(xx)),'linestyle','--')
title('Plot of function H({\beta})','Interpreter','tex','fontsize',18)
xlabel('{\beta}','Interpreter','tex','fontsize',14)
ylabel('H({\beta})','Interpreter','tex','fontsize',14)
xlim([0 35])

% Step 5.2: Relative Frequencies with Inverse Weibull Distribution

x_help    = 0:0.5:35;
x_count   = zeros((length(x_help)-1),1);

for j = 1:1:(length(x_help)-1)
  for k = 1:1:length(x_vec)
    if (x_help(j)<=x_vec(k) && x_vec(k)<=x_help(j+1))
      x_count(j) = x_count(j) + 1;
    endif
  endfor
endfor

x_count = x_count./length(x_vec);

figure(3)
bar(0.25:0.25:34.75,x_count)
hold on
plot(x_vec,f(x_vec),'linestyle','-','linewidth',1)
legend({'Relative frequencies','Inverse Weibull'},'location','northeast')
xlabel('t')
ylabel('Relative frequencies')
