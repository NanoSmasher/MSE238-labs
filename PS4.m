clear all; close all; clc; warning('off','all');
%%% MSE238_W2015_NM_HW2
%%% Henry Lu :: 1000588534

%%%%%%%%%%%%%%%%%%%%
%% pause % Problem 1
%%%%%%%%%%%%%%%%%%%%

%%% Variables
x = @(R,t) R.*(cos(t)+sqrt(2.5^2 - sin(t).^2));
h = pi/36;
r = 0:h:pi;

%%% First Principles

% two point central and four point forward
x_tpf = @(R,t,h) (x(R,t+h)-x(R,t-h))/(2*h);
x_fpf = @(R,t,h) (-1*x(R,t+3*h)+4*x(R,t+2*h)-5*x(R,t+h)+2*x(R,t))/(h^2);

% plot
figure(1);
set(gcf,'name','Question 1: First Principles');
subplot(3,1,1);
plot(r,x(1,r)); title('Position wrt theta');
subplot(3,1,2);
plot(r,x_tpf(1,r,h)); title('Velocity wrt theta - 2 point central');
subplot(3,1,3);
plot(r,x_fpf(1,r,h)); title('Acceleration wrt theta - 4 point forward');

%%% Matlab In-Built

% derivatives
fx = x(1,r);
df = diff(fx)/h;
d2f = diff(df)/h;

% plot
figure(2);
set(gcf,'name','Question 1: Matlab In-Built');
subplot(3,1,1);
plot(r,x(1,r)); title('Position wrt theta');
subplot(3,1,2);
plot(r(:,1:length(df)),df); title('Velocity wrt theta');
subplot(3,1,3);
plot(r(:,1:length(d2f)),d2f); title('Acceleration wrt theta');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pause % Problem 2 part a
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Variables
w = @(wo,t,R) wo.*(cos(t).^2./sqrt(R.^2 -sin(t).^2));
a = 0;
b = pi/2;
h = pi/180;
r = 1:0.1:10;
A = [0.652145 0.652145 0.347855 0.347855];
ei = [-0.339981 0.339981 -0.861136 0.861136];

%%% First Principles

% Trapezoidal, simposons, and quadrature integration methods
Trap = @(a,b,h,wo,R) h/2*(w(wo,a,R)+w(wo,b,R)) + h*sum(w(wo,a+2*h:h:b,R));

Simp = @(a,b,h,wo,R) h/3*(w(wo,a,R)+w(wo,b,R)+4*sum(w(wo,a+2*h:h*2:b,R)) + 2*sum(w(wo,a+3*h:h*2:b-h,R)));

fx = @(a,b,t,wo,R) 1/2*(b-a) * w(wo,(a+b+t.*(b-a))/2,R); % the change of variable function
Quad = @(a,b,wo,R) sum(A(1:length(A)).*fx(a,b,ei(1:length(ei)),wo,R)); % the actual function

% compute values
for i = 1:length(r)
    f_t(i) = Trap(a,b,h,1,r(i));
    f_s(i) = Simp(a,b,h,1,r(i));
    f_q(i) = Quad(a,b,1,r(i));
end

% plot
figure(3);
set(gcf,'name','Question 2: First Principles');
subplot(3,1,1); 
plot(r,f_t);title('Integration using Trapezoidal');
subplot(3,1,2); 
plot(r,f_s); title('Integration using Simpsons 1/3');
subplot(3,1,3); 
plot(r,f_q); title('Integration using Gauss quadrature');

%%% Matlab In-Built

% compute values
for i = 1:length(r)
    b_t(i) = trapz(a:h:b,w(1,a:h:b,r(i)));
    b_s(i) = integral(@(x)w(1,x,r(i)),a,b);
    b_q(i) = quad(@(x) w(1,x,r(i)),a,b); % QUAD will be removed in a future release
end

% plot
figure(4);
set(gcf,'name','Question 2: Matlab In-Built');
subplot(3,1,1); 
plot(r,b_t);title('Integration using Trapezoidal');
subplot(3,1,2); 
plot(r,b_s); title('Integration using Integrate');
subplot(3,1,3); 
plot(r,b_q); title('Integration using Gauss quadrature');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pause % Problem 2 part b
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Variables
erf_t = @(t) 2/sqrt(pi)*exp(-(t.^2));
A = [0.652145 0.652145 0.347855 0.347855];
ei = [-0.339981 0.339981 -0.861136 0.861136];

fx = @(x,t) x/2 * erf_t((x+t.*x)/2);
Quad = @(x) sum(A(1:length(A)).*fx(x,ei(1:length(ei))));

% print
fprintf('\nQuestion 2 part b:\n');
fprintf('4 point erf(-1): %.6f\n',Quad(-1));
fprintf('lobatto erf(-1): %.6f\n',quadl(erf_t,0,-1)); % QUADL will be removed in a future release
fprintf('gauss quad erf(-1): %.6f\n',quad(erf_t,0,-1)); % QUAD will be removed in a future release
fprintf('4 point erf(2): %.6f\n',Quad(2));
fprintf('lobatto erf(2): %.6f\n',quadl(erf_t,0,2)); % QUADL will be removed in a future release
fprintf('gauss quad erf(2): %.6f\n',quad(erf_t,0,2)); % QUAD will be removed in a future release

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pause % Problem 2 part c
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Variables
s = @(x,y) sin(pi.*x);
xa = 0;
xb = 1;
ya = 0;
yb = 1;

fx = @(t,y) (xb-xa)/2 * s((xa+xb+t.*(xb-xa))/2,y);
fq = @(y) sum(A(1:length(A)).*fx(ei(1:length(ei)),y));
fy = @(v)(yb-ya)/2 * fq((ya+yb+v.*(yb-ya))/2);
Quad = sum(A(1:length(A)).*fy(ei(1:length(ei))));

%%% Print
fprintf('\nQuestion 2 part c:\n');
fprintf('First principles integration: %.6f\n',Quad);
fprintf('In-built integration: %.6f\n',dblquad(s,0,1,0,1));  % DBLQUAD will be removed in a future release

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROBLEM 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n\nSymbolic Math starts now. All variables and figures will be cleared. Hit ENTER to continue\n\n')
pause; clear all; close all; clc; 
fprintf('Symbolic Math starts now. All variables and figures will be cleared. Hit ENTER to continue\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pause % Problem 3 Problem 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

syms R t;

x = R.*(cos(t)+sqrt(2.5^2 - sin(t).^2));
h = pi/36;
r = 0:h:pi;

fx = subs(x,'t',r);
dfx = (subs(x,'t',r+h)-subs(x,'t',r))/h;
d2fx = (-1*subs(x,'t',r+3*h) + 4*subs(x,'t',r+2*h) - 5*subs(x,'t',r+h) + 2*subs(x,'t',r))/(h^2);

figure(1);
set(gcf,'name','Question 1: Symbolic Principles');
subplot(3,1,1);
plot(r,subs(fx,'R',1)); title('Position wrt theta');
subplot(3,1,2);
plot(r,subs(dfx,'R',1)); title('Velocity wrt theta - 2 point central');
subplot(3,1,3);
plot(r,subs(d2fx,'R',1)); title('Acceleration wrt theta - 4 point forward');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pause % Problem 3 Problem 2 part a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

syms wo t v R;

w = wo.*(cos(t).^2./sqrt(R.^2 -sin(t).^2));
a = 0;
b = pi/2;
h = pi/180;
r = 1:0.5:10; % Changed to 0.5 to increase speed
A = [0.652145 0.652145 0.347855 0.347855];
ei = [-0.339981 0.339981 -0.861136 0.861136];

Trap = h/2*(subs(w,'t',a)+subs(w,'t',a)) + h*sum(subs(w,'t',a+2*h:h:b));
T_wo = subs(Trap,'wo',1);

Simp = h/3*(subs(w,'t',a)+subs(w,'t',b) +4*sum(subs(w,'t',a+2*h:h*2:b)) + 2*sum(subs(w,'t',a+3*h:h*2:b-h)));
S_wo = subs(Simp,'wo',1);

fx = 1/2*(b-a) * subs(w,'t',(a+b+v.*(b-a))/2);
Quad = sum(A(1:length(A)).*subs(fx,'v',ei(1:length(ei)))); % the actual function
Q_wo = subs(Quad,'wo',1);

figure(2);
set(gcf,'name','Question 2: Symbolic Principles');
subplot(3,1,1);
plot(r,subs(T_wo,'R',r)); title('Integration using Trapezoidal');
subplot(3,1,2);
plot(r,subs(S_wo,'R',r)); title('Integration using Simpsons 1/3');
subplot(3,1,3);
plot(r,subs(Q_wo,'R',r)); title('Integration using Gauss quadrature');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pause % Problem 3 Problem 2 part b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

syms x t;
erf_t = 2/sqrt(pi)*exp(-(t.^2));
fx = x/2 * subs(erf_t,'t',(x+t.*x)/2);
Quad = sum(A(1:length(A)).*subs(fx,'t',ei(1:length(ei))));

fprintf('\nQuestion 2 part b:\n');
fprintf('4 point erf(-1): %.6f\n',double(subs(Quad,'x',-1)));
fprintf('4 point erf(2): %.6f\n',double(subs(Quad,'x',2)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pause % Problem 3 Problem 2 part c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

syms x y t v;
s = sin(pi.*x);
xa = 0;
xb = 1;
ya = 0;
yb = 1;

ft = (xb-xa)/2 * subs(s,'x',(xa+xb+t.*(xb-xa))/2);
fy = sum(A(1:length(A)).*subs(ft,'t',ei(1:length(ei))));
fv = (yb-ya)/2 * subs(fy,'y',(ya+yb+v.*(yb-ya))/2);
Dquad = sum(A(1:length(A)).*subs(fv,'v',ei(1:length(ei))));

fprintf('\nQuestion 2 part c:\n');
fprintf('First principles integration: %.6f\n',double(Dquad));

fprintf('\n\nEnd of file - Henry Lu\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF FILE
