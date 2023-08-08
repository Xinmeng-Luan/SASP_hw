clear all
close all
clc

%% Import Input Audio Signal
[Vin,~] = audioread('ExpSweep.wav');

%% LTSpice Files for Ground-Truth
[OutLowSpice,~]=audioread('outlowsweep.wav');
[OutMidSpice,~]=audioread('outmidsweep.wav');
[OutHighSpice,FsLTSpice]=audioread('outhighsweep.wav');
TsLTSpice=1/FsLTSpice;

%% Sampling frequency (to be varied: FsLTSpice/downSampFact, with downSampFact={4,3,2})
downSampFact=2;
Fs =FsLTSpice/downSampFact; 

%% Downsample Input Signal
Vin=Vin([1:downSampFact:end]);

%% Sampling Period
Ts=1/Fs;
%% Number of Samples
Nsamp=length(Vin);
%% Simulated time
tstop=Nsamp*Ts;
%% Parameters of Dynamic Element
L1=0.35*10^(-3);
L2=0.35*10^(-3);
L3=3.5*10^(-3);
L4=3.5*10^(-3);
C1= 2.8*10^(-6);
C2= 2.8*10^(-6);
C3= 28*10^(-6);
C4= 4.7*10^(-6);
C5=28*10^(-6);
C6=47*10^(-6);
%% Resistive Parameters
R1=10;
RspkLow=8;
R2=10;
RspkMid=8;
RspkHigh=8;
%% WDF setting of free parameters (adaptation conditions)
%   Z = R               resistor
%   Z = (2*L) / Ts      inductor
%   Z = Ts / (2*C)      capacitor
%   
%   Z = Z + Z           series
%   Z = (Z*Z) / (Z+Z)   parallel

%--------------------------------------------------------
%   high
Z10 = RspkHigh;
Z12 = (2*L1) / Ts;
Z9 = Ts / (2*C1);

Z11 = (Z10*Z12) / (Z10+Z12);
Z7 = Z11;
Z8 = Z7 + Z9;

%--------------------------------------------------------
%   mid
Z30 = Ts / (2*C4);
Z28 = R1;
Z26 = RspkMid;
Z24 = (2*L3) / Ts;
Z21 = Ts / (2*C3);
Z18 = Ts / (2*C2);
Z13 = (2*L2) / Ts;

Z29 = Z28 + Z30;
Z27 = Z29;
Z25 = (Z26*Z27) / (Z26+Z27);
Z23 = Z25;
Z22 = (Z23*Z24) / (Z23+Z24);
Z20 = Z22;
Z19 = Z20 + Z21;
Z17 = Z19;
Z16 = (Z17*Z18) / (Z17+Z18);
Z15 = Z16;
Z14 = Z15 + Z13;

%--------------------------------------------------------
%   low
Z41 = Ts / (2*C6);
Z42 = R2;
Z39 = RspkLow;
Z36 = Ts / (2*C5);
Z31 = (2*L4) / Ts;

Z40 = Z41 + Z42;
Z38 = Z40;
Z37 = (Z38*Z39) / (Z38+Z39);
Z35 = Z37;
Z34 = (Z35*Z36) / (Z35+Z36);
Z33 = Z34;
Z32 = Z31 + Z33;

%--------------------------------------------------------
%   general circuit
Z4 = Z32;
Z6 = Z14;
Z5 = (Z4*Z6) / (Z4+Z6);
Z1 = Z8;
Z2 = Z5;
Z3 = (Z1*Z2) / (Z1+Z2);

%% Computation of Scattering Matrices
% We defined two different functions to compute the matrices

% 3-port series adaptor
Ss1 = series(Z7,Z8,Z9);
Ss2 = series(Z13,Z14,Z15);
Ss3 = series(Z19,Z20,Z21);
Ss4 = series(Z28,Z29,Z30);
Ss5 = series(Z31,Z32,Z33);
Ss6 = series(Z40,Z41,Z42);

% 3-port parallel adaptor
Sp1 = parallel(Z1,Z2,Z3);
Sp2 = parallel(Z4,Z5,Z6);
Sp3 = parallel(Z10,Z11,Z12);
Sp4 = parallel(Z16,Z17,Z18);
Sp5 = parallel(Z22,Z23,Z24);
Sp6 = parallel(Z25,Z26,Z27);
Sp7 = parallel(Z34,Z35,Z36);
Sp8 = parallel(Z37,Z38,Z39);


%% Initialization of Waves
a = zeros(42, Nsamp+1);
b = zeros(42, Nsamp+1);

%% Initialize Output Signals
% Low
VoutLow=zeros(size(Vin));
% Mid
VoutMid=zeros(size(Vin));
% High
VoutHigh=zeros(size(Vin));

ii=1;       % because of the first sample (padded)
while (ii<Nsamp)
    ii=ii+1;

    %% Manage Dynamic Elements
    % high 
    a(9, ii) = b(9, ii-1);
    a(12, ii) = -b(12, ii-1);
    a(10, ii) = 0;

    % mid 
    a(13, ii) = -b(13, ii-1);
    a(18, ii) = b(18, ii-1);
    a(21, ii) = b(21, ii-1);
    a(24, ii) = -b(24, ii-1);
    a(26, ii) = 0;
    a(28, ii) = 0;
    a(30, ii) = b(30, ii-1);

    % low
    a(31, ii) = -b(31, ii-1);
    a(36, ii) = b(36, ii-1);
    a(39, ii) = 0;
    a(42, ii) = 0;
    a(41, ii) = b(41, ii-1);

    %% Forward Scan
    % high
    b(11, ii) = Sp3(2,:) * a(10:12, ii);
    a(7, ii) = b(11, ii);
    b(8, ii) = Ss1(2,:) * a(7:9, ii);
    a(1, ii) = b(8, ii);
    
    % mid 
    b(29, ii) = Ss4(2,:) * a(28:30, ii);
    a(27, ii) = b(29, ii);
    b(25, ii) = Sp6(1,:) * a(25:27, ii);
    a(23, ii) = b(25, ii);
    b(22, ii) = Sp5(1,:) * a(22:24, ii);
    a(20, ii) = b(22, ii);
    b(19, ii) = Ss3(1,:) * a(19:21, ii);
    a(17, ii) = b(19, ii);
    b(16, ii) = Sp4(1,:) * a(16:18, ii);
    a(15, ii) = b(16, ii);
    b(14, ii) = Ss2(2,:) * a(13:15, ii);
    a(6, ii) = b(14, ii);

    % low
    b(40, ii) = Ss6(1,:) * a(40:42, ii);
    a(38, ii) = b(40, ii);
    b(37, ii) = Sp8(1,:) * a(37:39, ii);
    a(35, ii) = b(37, ii);
    b(34, ii) = Sp7(1,:) * a(34:36, ii);
    a(33, ii) = b(34, ii);
    b(32, ii) = Ss5(2,:) * a(31:33, ii);
    a(4, ii) = b(32, ii);

    % general circuit
    b(5, ii) = Sp2(2,:) * a(4:6, ii);
    a(2, ii) = b(5, ii);
    b(3, ii) = Sp1(3,:) * a(1:3, ii);

    %% Local Root Scattering
    a(3, ii) = 2*Vin(ii-1) - b(3, ii);

    %% Backward Scan
    % general circuit
    b(1, ii) = Sp1(1,:) * a(1:3, ii);
    a(8, ii) = b(1, ii);
    b(2, ii) = Sp1(2,:) * a(1:3, ii);
    a(5, ii) = b(2, ii);
    b(4, ii) = Sp2(1,:) * a(4:6, ii);
    a(32, ii) = b(4, ii);
    b(6, ii) = Sp2(3,:) * a(4:6, ii);
    a(14, ii) = b(6, ii);

    % high
    b(7, ii) = Ss1(1,:) * a(7:9, ii);
    a(11, ii) = b(7, ii);
    b(9, ii) = Ss1(3,:) * a(7:9, ii);
    b(10, ii) = Sp3(1,:) * a(10:12, ii);
    b(12, ii) = Sp3(3,:) * a(10:12, ii);

    % mid 
    b(13, ii) = Ss2(1,:) * a(13:15, ii);
    b(15, ii) = Ss2(3,:) * a(13:15, ii);
    a(16, ii) = b(13, ii);
    b(18, ii) = Sp4(3,:) * a(16:18, ii);
    b(17, ii) = Sp4(2,:) * a(16:18, ii);
    a(19, ii) = b(17, ii);
    b(21, ii) = Ss3(3,:) * a(19:21, ii);
    b(20, ii) = Ss3(2,:) * a(19:21, ii);
    a(22, ii) = b(20, ii);
    b(24, ii) = Sp5(3,:) * a(22:24, ii);
    b(23, ii) = Sp5(2,:) * a(22:24, ii);
    a(25, ii) = b(23, ii);
    b(26, ii) = Sp6(2,:) * a(25:27, ii);
    b(27, ii) = Sp6(3,:) * a(25:27, ii);
    a(29, ii) = b(27, ii);
    b(28, ii) = Ss4(1,:) * a(28:30, ii);
    b(30, ii) = Ss4(3,:) * a(28:30, ii);

    % low 
    b(31, ii) = Ss5(1,:) * a(31:33, ii);
    b(33, ii) = Ss5(3,:) * a(31:33, ii);
    a(34, ii) = b(33, ii);
    b(36, ii) = Sp7(3,:) * a(34:36, ii);
    b(35, ii) = Sp7(2,:) * a(34:36, ii);
    a(37, ii) = b(35, ii);
    b(39, ii) = Sp8(3,:) * a(37:39, ii);
    b(38, ii) = Sp8(2,:) * a(37:39, ii);
    a(40, ii) = b(38, ii);
    b(42, ii) = Ss6(3,:) * a(40:42, ii);
    b(41, ii) = Ss6(2,:) * a(40:42, ii);

    %% Read Output
    
    
end


%% Output Plots
figure
set(gcf, 'Color', 'w');
subplot(311)
plot(TsLTSpice*[1:length(OutLowSpice)],OutLowSpice,'r','Linewidth',2); hold on;
plot(Ts*[1:Nsamp],VoutLow,'b--','Linewidth',1); grid on; xlim([0,tstop]); 
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$V_{\mathrm{outLow}}$ [V]','Fontsize',16,'interpreter','latex');
legend('LTspice','WDF','Fontsize',16,'interpreter','latex');
title('Output Signals','Fontsize',18,'interpreter','latex');
subplot(312)
plot(TsLTSpice*[1:length(OutMidSpice)],OutMidSpice,'r','Linewidth',2); hold on;
plot(Ts*[1:Nsamp],VoutMid,'b--','Linewidth',1); grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$V_{\mathrm{outMid}}$ [V]','Fontsize',16,'interpreter','latex');
subplot(313)
plot(TsLTSpice*[1:length(OutHighSpice)],OutHighSpice,'r','Linewidth',2); hold on;
plot(Ts*[1:Nsamp],VoutHigh,'b--','Linewidth',1); grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$V_{\mathrm{outHigh}}$ [V]','Fontsize',16,'interpreter','latex');

%% Error Plots
figure
set(gcf, 'Color', 'w');
subplot(311)
plot(Ts*[1:Nsamp],OutLowSpice([1:downSampFact:end])-VoutLow,'k','Linewidth',1);grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$E_{\mathrm{outLow}}$ [V]','Fontsize',16,'interpreter','latex');
title(['Error Signals. $F_{\mathrm{s}}=$ ',num2str(Fs),' Hz'],'Fontsize',18,'interpreter','latex');
subplot(312)
plot(Ts*[1:Nsamp],OutMidSpice([1:downSampFact:end])-VoutMid,'k','Linewidth',1);grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$E_{\mathrm{outMid}}$ [V]','Fontsize',16,'interpreter','latex');
subplot(313)
plot(Ts*[1:Nsamp],OutHighSpice([1:downSampFact:end])-VoutHigh,'k','Linewidth',1);grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$E_{\mathrm{outHigh}}$ [V]','Fontsize',16,'interpreter','latex');

