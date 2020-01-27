function fig=ex9_1()
% fig=ex9_1()
%
% Executes Example 9.1 and generates one figure
%
% INPUTS
%   none
%
% OUTPUTS
%   fig     figure handle to generated graphic
%
% Nicholas O'Donoughue
% 1 July 2019

% Set up error covariance matrix
C = [10 -3;-3 5];
[V,D] = eig(C);
[lam_sort,i_sort] = sort(diag(D),'descend');
v_max = V(:,i_sort(1));
v_min = V(:,i_sort(2));

gamma = 4.601; % 90% confidence interval
a = sqrt(gamma*lam_sort(1));
b = sqrt(gamma*lam_sort(2));

numPts = 101;
th = linspace(0,2*pi,numPts);
x1 = a*cos(th);
x2 = b*sin(th);

alpha = atan2(v_max(2),v_max(1));
x = x1*cos(alpha)-x2*sin(alpha);
y = x1*sin(alpha)+x2*cos(alpha);

fig=figure;
plot(0,0,'k+','DisplayName','Bias Point');hold on;
text(-2,.25,'Bias Point','FontSize',12);
text(3.5,3,'90\% Error Ellipse','FontSize',12);
plot(x,y,'k-','DisplayName','Error Ellipse');

% Draw the semi-minor and semi-major axes
plot([0 a*cos(alpha+pi)],[0 a*sin(alpha+pi)],'k:');
plot([0 -b*sin(alpha+pi)],[0 b*cos(alpha+pi)],'k:');
text(4.5,-2,'$r_1=7.24$','FontSize',12);
text(1,2,'$r_2=4.07$','FontSize',12);

plot([0 3],[0,0],'k:');
plot(2*cosd(0:-.1:-25),2*sind(0:-.1:-25),'k-','LineWidth',.5)
text(2.1,-.75,'$\alpha = -25^\circ$','FontSize',12);
%xlabel('x');
%ylabel('y');
%legend('Location','NorthEast');
utils.setPlotStyle(gca,{'equal','clean','tight'});
