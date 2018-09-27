%% Rio
load([dropboxPath,'/ratsEphys/Rio/Rio_txtDataFiles/cmbObjsV2/CMBH_2017-08-20_CircleTrack.mat'])

% 
[thSlp_byVel,thAmp_byVel, thFrq_byVel, vel_dim] = TW_spdMod(root,[6 10],0,18:22);

figure; plot(vel_dim,thSlp_byVel)
xlabel('Running speed'); ylabel('Pairwise phase difference (radians)');legend('1-2','2-3','3-4','4-5','5-6')
figure; plot(vel_dim,thAmp_byVel)
xlabel('Running speed'); ylabel('Theta power (log(A.U.))'); legend('1','2','3','4','5','6')
figure; plot(vel_dim,thFrq_byVel)
xlabel('Running speed'); ylabel('Theta freq (Hz)'); title('mean'); legend('1','2','3','4','5','6')
text(10,7.6,'REGIO','FontSize',22)

figure; plot(vel_dim,mean(thSlp_byVel,1),'-k','LineWidth',2);
xlabel('Running speed'); ylabel('Pairwise phase difference (radians)'); title('mean per channel')

keyboard;
%% Regio
comclear all
load('/Users/ehren/Dropbox (NewmanLab)/ratsEphys/Regio/Regio_txtDataFiles/cmbObjsV2/CMBH_CircleTrack_2017-12-14.mat');

[thSlp_byVel,thAmp_byVel, thFrq_byVel, vel_dim] = TW_spdMod(root,[6 10],0,1:6);

figure; plot(vel_dim,thSlp_byVel)
xlabel('Running speed'); ylabel('Pairwise phase difference (radians)');legend('1-2','2-3','3-4','4-5','5-6')
figure; plot(vel_dim,thAmp_byVel)
xlabel('Running speed'); ylabel('Theta power (log(A.U.))'); legend('1','2','3','4','5','6')
figure; plot(vel_dim,thFrq_byVel)
xlabel('Running speed'); ylabel('Theta freq (Hz)'); title('mean'); legend('1','2','3','4','5','6')
text(10,7.6,'REGIO','FontSize',22)

figure; plot(vel_dim,mean(thSlp_byVel,1),'-k','LineWidth',2);
xlabel('Running speed'); ylabel('Pairwise phase difference (radians)'); title('mean per channel')

%% Romo
clear all
load(fullfile(dropboxPath,array2path(path2array('ratsEphys/Romo/Romo_txtDataFiles/cmbObjsV2/CMBH_CircleTrack_2017-11-22.mat'))));

[thSlp_byVel,thAmp_byVel, thFrq_byVel, vel_dim] = TW_spdMod(root,[6 10],0,1:6);

figure; plot(vel_dim,thSlp_byVel)
xlabel('Running speed'); ylabel('Pairwise phase difference (radians)');legend('1-2','2-3','3-4','4-5','5-6')
figure; plot(vel_dim,thAmp_byVel)
xlabel('Running speed'); ylabel('Theta power (log(A.U.))'); legend('1','2','3','4','5','6')
figure; plot(vel_dim,thFrq_byVel)
xlabel('Running speed'); ylabel('Theta freq (Hz)'); title('mean'); legend('1','2','3','4','5','6')
text(10,7.6,'ROMO','FontSize',22)

figure; plot(vel_dim,mean(thSlp_byVel,1),'-k','LineWidth',2);
xlabel('Running speed'); ylabel('Pairwise phase difference (radians)'); title('mean per channel')

%% Roble
clear all
load(fullfile(dropboxPath,array2path(path2array('ratsEphys/Roble/Roble_txtDataFiles/cmbObjsV2/CMBH_CircleTrack_2018-04-18.mat'))));

[thSlp_byVel,thAmp_byVel, thFrq_byVel, vel_dim] = TW_spdMod(root,[6 10],0,1:6);

figure; plot(vel_dim,thSlp_byVel)
xlabel('Running speed'); ylabel('Pairwise phase difference (radians)');legend('1-2','2-3','3-4','4-5','5-6')
figure; plot(vel_dim,thAmp_byVel)
xlabel('Running speed'); ylabel('Theta power (log(A.U.))'); legend('1','2','3','4','5','6')
figure; plot(vel_dim,thFrq_byVel)
xlabel('Running speed'); ylabel('Theta freq (Hz)'); title('mean'); legend('1','2','3','4','5','6')
text(10,7.6,'ROBLE','FontSize',22)
figure; plot(vel_dim,mean(thSlp_byVel,1),'-k','LineWxidth',2);
xlabel('Running speed'); ylabel('Pairwise phase difference (radians)'); title('mean per channel')


figure; plot(acc_dim,thSlp_byAcc)
xlabel('Running speed'); ylabel('Pairwise phase difference (radians)');legend('1-2','2-3','3-4','4-5','5-6')
figure; plot(acc_dim,thAmp_byAcc)
xlabel('Running speed'); ylabel('Theta power (log(A.U.))'); legend('1','2','3','4','5','6')
figure; plot(acc_dim,thFrq_byAcc)
xlabel('Running speed'); ylabel('Theta freq (Hz)'); title('mean'); legend('1','2','3','4','5','6')
text(10,7.6,'ROBLE','FontSize',22)
figure; plot(acc_dim,mean(thSlp_byAcc,1),'-k','LineWxidth',2);
xlabel('Running speed'); ylabel('Pairwise phase difference (radians)'); title('mean per channel')


