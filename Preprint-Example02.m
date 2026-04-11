% Step 1: Data Vector and Its Logarithmic Vector

x_vec    = [4 6 7 11 13];

y_vec    = (log(x_vec))';

z_vec    = x_vec(1)./x_vec;

G_infty  = sum(log(z_vec));

% Step 2: Definition of G, H, dH and dG from our Manuscript

G     = @(beta) (length(x_vec)./beta)-sum(log(x_vec))+(length(x_vec)*(((x_vec).^(-beta))*y_vec))/(sum((x_vec).^(-beta)));
H     = @(beta) length(x_vec)*(((x_vec).^(-beta))*y_vec)/(sum((x_vec).^(-beta)));
dH    = @(beta) length(x_vec)/(sum((x_vec).^(-beta))*sum((x_vec).^(-beta)))*((((x_vec).^(-beta))*y_vec)*(((x_vec).^(-beta))*y_vec)-(((x_vec).^(-beta))*(y_vec.^2))*(sum((x_vec).^(-beta))));
dG    = @(beta) -length(x_vec)/(beta*beta)+dH(beta);

% Step 3: Plotting of G and H

xx    = 0.1:0.1:14.9;
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

% Step 5: Plotting of G, H

figure(1)
plot(xx,yy)
hold on
plot(xx,G_infty*ones(1,length(xx)),'linestyle','--')
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
xlim([0 15])

% Step 5: Bisection-Iteration for G(beta_final) = 0

##a_guess   = 0.1;
##b_guess   = 1.1;
##x_guess   = (a_guess+b_guess)/2;
##epsilon   = 0.01;
##
##if ((G(a_guess)>0 & G(b_guess)<0))
##  while(abs(G(x_guess))>epsilon)
##    if ((G(a_guess)>0 & G(x_guess)<0))
##       x_guess = (a_guess+x_guess)/2;
##    elseif ((G(x_guess)>0 & G(b_guess)<0))
##       x_guess = (x_guess+b_guess)/2;
##    endif
##  endwhile
##elseif
##   x_guess = (a_guess+b_guess)/2;
##endif
