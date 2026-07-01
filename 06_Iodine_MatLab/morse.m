function morse()
% This subroutine is written for Expt3 of Chem310L in Duke University
% morse.m plots Morse potential figures of the ground and excited state of 
% I_2 molecule and the parameters in matmorse.mat for overlap.m

% Revision Log
% V1.0 - Original
% V2.0 - Update for 5.2 PC, 11-24-99
% V2.1 - Save data for Morse2.m, 01-03-00
%
% V3.0 - Complete rearrangement, Nov 2019
% Rearrange the frame: morse.m -> matmorse.mat -> overlay.m
% Remoce the annoying name requirement
% Remove the print function for electronic notebook is widely uesd now
% Add data verification
% Add the solution for R_e
% Add animation to visualization
% Update some information of this course
% Update some Matlab functions to recommended new version
% Written by Dongtian Xu, e-mail: dongtianxu00@gmail.com

format compact
format shortG

% Display information of the subroutine
clc
disp('    ****************************************************')
disp('    **         Chemistry 310L : Experiment 3          **')
disp('    ** Calculation of Ground and First Excited State  **')
disp('    **     Morse Potential Wells for Iodine Vapor     **')
disp('    ****************************************************')
disp('You will be asked to enter your values of D_e, a, w_e and x_e');
disp('If you want to get an instant view of the potential figure, press');
disp('ENTER to skip the input part and literature values will be used.');
disp('If you want to plot figures for youe ELN, please make sure you use')
disp('your own data.')
disp('Press any key to continue...')
pause();

% Get excited state parameters
clc
disp('    ****************************************************');
disp('    **  Calculation of Excited State Morse Potential  **');
disp('    ****************************************************');
disp(' ')
De = input('Enter your Excited State Dissociation Energy (De'') (cm-1) : ');
if (isempty(De))
  De = 4112;
end
disp(' ');
a = input('Enter your Excited State Morse Parameter (a'') (1/Angstroms) : ');
if (isempty(a))
  a = 2;
end
disp(' ');
w = input('Enter your Excited State Harmonic Frequency (we'') (cm-1) : ');
if (isempty(w))
  w = 132.11;
end
disp(' ');
x = input('Enter your Excited State Anharmonicity Parameter (xe'') (cm-1) : ');
if (isempty(x))
  x = 0.008;
end
disp(' ');

% Get ground state parameters
clc
disp('    ****************************************************');
disp('    **   Calculation of Ground State Morse Potential  **');
disp('    ****************************************************');
disp(' ');
De2 = input('Enter your Ground State Dissociation Energy (De") (cm-1) : ');
if (isempty(De2))
  De2 = 12244;
end
disp(' ');
a2 = input('Enter your Ground State Morse Parameter (a") (1/Angstroms) : ');
if (isempty(a2))
  a2 = 1.87;
end
disp(' ');
w2 = input('Enter your Ground State Harmonic Frequency (we") (cm-1) : ');
if (isempty(w2))
  w2 = 214.5;
end
disp(' ');
x2 = input('Enter your Ground State Anharmonicity Parameter (xe") (cm-1) : ');
if (isempty(x2))
  x2 = 0.0031;
end
disp(' ');

% Verify input data
clc
disp('Your input parameters are listed below:')
paraSheet = {'ParaSheet','D_e','a','w_e','x_e'; ...
            'Excited State',De,a,w,x; ...
            'Ground State',De2,a2,w2,x2};
disp(paraSheet)

% Morse potential
r = -0.5:0.01:6;
morse = De*(1-exp(-a*r)).^2;
morse2 = De2*(1-exp(-a2*r)).^2';
% V = [-0.5 2 0 6000]; % axis range
% Energy levles
v = 0:80; 
G = ((w*(v+0.5))-((w*x)*((v+0.5).^2)));
G2 = ((w2*(v+0.5))-((w2*x2)*((v+0.5).^2)));
G(G > De) = []; % In case the energy level surpass the dissociation energy
G2(G2 > De2) = [];
Re = 0;
Re2 = 0;
levelr(1,1:length(G)) = (-log(1+sqrt(G./De))./a)+Re;
levelr(2,1:length(G)) = (-log(1-sqrt(G./De))./a)+Re;
levelr2(1,1:length(G2)) = (-log(1+sqrt(G2./De2))./a2)+Re2;
levelr2(2,1:length(G2)) = (-log(1-sqrt(G2./De2))./a2)+Re2;
energy = [G; G];
energy2 = [G2; G2];

% Plot excited state Morse potential
figure
hold on
% axis(V)
plot(r,morse,'linewidth',1);
plot(levelr,energy,'-k','linewidth',0.5);
title('Excited State Morse Potential for Iodine Vapor');
xlabel('Interatomic Separation R (Angstroms)');
ylabel('Potential Energy V(R) (wavenumbers cm^{-1})');

% Plot ground state Morse potential
figure
hold on
% axis(V)
plot(r,morse2,'linewidth',1);
plot(levelr2,energy2,'-k','linewidth',0.5);
title('Ground State Morse Potential for Iodine Vapor');
xlabel('Interatomic Separation R (Angstroms)');
ylabel('Potential Energy V(R) (wavenumbers cm^{-1})');

save matmorse De De2 a a2 w w2 x x2
