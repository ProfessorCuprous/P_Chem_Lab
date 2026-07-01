function overlap()
% This subroutine is written for Expt3 of Chem310L in Duke University
% overlap.m plots the overlapped figure and provides the solution of R_e 
% via an animation plot. This subroutine needs parameters in matmorse.mat,
% which is generated in morse.m.

% Revision Log
% V1.0 - Original 'morse2.m'
% V2.0 - Update for R11.1 UNIX, 01-03-00
% 
% V3.0 - Renamed as 'overlap.m', Nov 2019
% Rearrange the frame: morse.m -> matmorse.mat -> overlay.m
% Remoce the annoying name requirement
% Remove the print function for electronic notebook is widely uesd now
% Add data verification
% Add the solution for R_e
% Add animation to visualization
% Update some information of this course
% Update some Matlab functions to recommended new version
% Written by Dongtian Xu, e-mail: dongtianxu00@gmail.com

% load matmorse.mat
try
    load('matmorse.mat');
catch
    warning('matmorse.mat is not found and the literature values are used automatically. If you want to use your own data, please run morse.m again');
    De = 4112; a = 2; w = 132.11; x = 0.008;
    De2 = 12244; a2 = 1.87; % w2 = 214.5; x2 = 0.0031;
end

% display some information
clc
disp('    ****************************************************')
disp('    **         Chemistry 310L : Experiment 3          **')
disp('    **           Calculation of R''_e - R''''_e          **')
disp('    **  Overlapped Morse Potential for Iodine Vapor   **')
disp('    ****************************************************')
disp('Parameters are loaded from matmorse.m:')
paraSheet = {'ParaSheet','D_e','a','w_e','x_e'; ...
            'Excited State',De,a,w,x; ...
            'Ground State',De2,a2,w2,x2};
disp(paraSheet)

r = -0.45:0.01:3;
V = [-0.5,3,0,3e4];  
morse = De*(1-exp(-a*r)).^2;
morse2 = De2*(1-exp(-a2*r)).^2';
v = 25; % v' = 25 is the most probable excited state
G = ((w*(v+0.5))-((w*x)*((v+0.5).^2)));
leftR = (-log(1+sqrt(G./De))./a); % leftR is negrtive
rightR = (-log(1-sqrt(G./De))./a);
dr = -leftR;

% calculate R'_e - R''_e
figure
axis(V)
hold on
plot(r,morse2,'m','linewidth',1); % ground state morse potential
h1 = plot(r,morse,'b','linewidth',1); % 1st excited state morse potential
h2 = plot([leftR rightR],[G G],'k'); % draw the level line of v' = 25
plot([0 0],[0 2e4],'k'); % vertical line at r = 0
h3 = text(rightR*1.1,G,'v''=25');
legend('Ground State','First Excited State');
title('Overlapped Ground and First Excited State Morse Potential for Iodine Vapor');
xlabel('Interatomic Separation R (Angstroms)');
ylabel('Potential Energy V(R) (wavenumbers cm^{-1})');
% move the excited state
h4 = text(0.2,V(4)/2,'press any key to move the excited state curve');
pause
delete(h4);
refresh
ddr = dr/30;
for i = 1:30
    h1.XData = h1.XData + ddr;
    h2.XData = h2.XData + ddr;
    h3.Position(1) = h3.Position(1) + ddr;
    pause(0.05)
end
pause(0.3)
stem(0,G,'k')
stem(dr,G,'k')
h4 = plot([0,ddr],[G G],'r','linewidth',1.5);
legend('Ground State','First Excited State');
for i = 1:30
    h4.XData(2) = h4.XData(2) + ddr;
    pause(0.05)
end
legend('Ground State','First Excited State');
hold off
drtxt = num2str(dr);
text(dr,1.5*G,['$R''_e-R''''_e=$' drtxt],'interpreter','latex')

disp('Press any key to continue...')
pause

% Plot the final overlapped figure
clc
disp('    ****************************************************');
disp('    **         Plot the final overlapped figure       **');
disp('    ****************************************************');
disp(' ');
% Only T''_e is needed to be input here
Te = input('Enter calculated Excited State Electronic Term Energy (Te'') (cm-1) : ');
if (isempty(Te))
Te = 15730;
end

Re2 = 2.66;
Re = Re2 + dr;
Te2 = 0;
r = 2:0.01:6;
r2 = 2:0.01:6;
morse = Te+(De*(1-exp(-a*(r-Re))).^2);
morse(morse>2.5e4) = [];
r(Te+(De*(1-exp(-a*(r-Re))).^2)>2.5e4) = [];
morse2 = Te2+(De2*(1-exp(-a2*(r2-Re2))).^2)';
morse2(morse2>1.8e4) = [];
r2(Te2+(De2*(1-exp(-a2*(r2-Re2))).^2)>1.8e4) = [];
v = 0:50; % the number of levels is reduced 
G = (((w*(v+0.5))-((w*x)*((v+0.5).^2))));
G2 = (((w2*(v+0.5))-((w2*x2)*((v+0.5).^2))));
G(G>De) = [];
G2(G2>De2) = [];
levelr(1,1:length(G)) = (-log(1+sqrt(G./De))./a)+Re;
levelr(2,1:length(G)) = (-log(1-sqrt(G./De))./a)+Re;
levelr2(1,1:length(G2)) = (-log(1+sqrt(G2./De2))./a2)+Re2;
levelr2(2,1:length(G2)) = (-log(1-sqrt(G2./De2))./a2)+Re2;
G = G + Te;
G2 = G2 + Te2;
energy = [G; G];
energy2 = [G2; G2];

figure
% axis([2 5 0 2e4])
plot(r,morse,'linewidth',1);
hold on;
plot(r2,morse2,'linewidth',1);
plot(levelr,energy,'k-');
plot(levelr2,energy2,'k-');
title('I-I Ground and First Excited State Morse Potentials');
xlabel('Interatomic Separation R (Angstroms)');
ylabel('Potential Energy V(R) (wavenumbers cm^{-1})');
legend('First Excited State','Ground State')

end