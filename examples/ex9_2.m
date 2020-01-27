function fig = ex9_2()
% fig=ex9_2()
%
% Executes Example 9.2 and generates one figure
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

cep50 = utils.computeCEP50(C);
fprintf('CEP50: %.2f \n',cep50);

ell = utils.drawErrorEllipse([0,0]',C,101,50);
ell_cep = utils.drawCEP50([0,0]',C,101);

fig=figure;
plot(0,0,'k^','DisplayName','Bias Point');hold on;
plot(ell(1,:),ell(2,:),'k--','DisplayName','Error Ellipse');
plot(ell_cep(1,:),ell_cep(2,:),'k-','DisplayName','$CEP_{50}$');
text(-1.3,.1,'Bias Point','FontSize',12);
text(-.3,1.4,'50\% Error Ellipse','FontSize',12);
text(2.2,2.3,'$CEP_{50}$','FontSize',12);

%xlabel('x');
%ylabel('y');
%legend('Location','NorthWest');
utils.setPlotStyle(gca,{'equal','clean','tight'});
