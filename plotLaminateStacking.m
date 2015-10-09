% Plots deflections based on laminate configurations

load('laminae.mat')

figure(1)
hold on
plot(radtodeg(result(1).stack(:,1)),result(1).deflection,'k')
plot(radtodeg(result(2).stack(:,1)),result(2).deflection,'k--')
plot(radtodeg(result(3).stack(:,1)),result(3).deflection,'k:')
plot(radtodeg(abs(result(4).stack(:,1))),result(4).deflection,'k*')

xlabel('Ply angle [deg]')
ylabel('Deflection [mm]')
legend('Case 1','Case 2','Case 3','Case 4')

figure(2)
hold on
plot(radtodeg(result(1).stack(:,1)),result(1).deflection,'k')
plot(radtodeg(result(5).stack(:,4)),result(5).deflection,'k--')
plot(radtodeg(result(6).stack(:,4)),result(6).deflection,'k:')

xlabel('Ply angle [deg]')
ylabel('Deflection [mm]')
legend('Case 1','Case 5','Case 6')

figure(3)
hold on
plot(radtodeg(result(1).stack(:,1)),result(1).deflection,'k')
plot(radtodeg(abs(result(7).stack(:,4))),result(7).deflection,'k--')
plot(radtodeg(abs(result(8).stack(:,1))),result(8).deflection,'k*')

xlabel('Ply angle [deg]')
ylabel('Deflection [mm]')
legend('Case 1','Case 7','Case 8')