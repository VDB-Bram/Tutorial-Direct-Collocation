% Compare optimization approaches

Shooting    = load(fullfile(pwd,'Shooting.mat'));
Collocation = load(fullfile(pwd,'Collocation_opti.mat'));


figure();

% plot joint angles
subplot(2,3,1:3)
plot(Shooting.tsim,Shooting.qsim,'b'); hold on;
plot(Collocation.tsim,Collocation.qsim,'r');
plot(Shooting.tsim,Shooting.q_exp,'--k');

legend('Shooting','Collocation','experimental');
xlabel('Time [s]'); ylabel('angle [rad]');

% plot the gains

subplot(2,3,4);
b = bar(1,Shooting.B);      b.FaceColor = [0 0 1]; hold on;
b = bar(2,Collocation.B);   b.FaceColor = [1 0 0];
set(gca,'XTickLabel',{'Shooting','Collocation'});
ylabel('B');


subplot(2,3,5);
b = bar(1,Shooting.Tb);     b.FaceColor = [0 0 1]; hold on;
b = bar(2,Collocation.Tb);  b.FaceColor = [1 0 0];
set(gca,'XTickLabel',{'Shooting','Collocation'});
ylabel('Tb');



% compare objective value
subplot(2,3,6)
J_shoot = sumsqr(Shooting.qsim-Shooting.q_exp);
J_Coll = sumsqr(Collocation.qsim-Collocation.q_exp);
b = bar(1,J_shoot); b.FaceColor = [0 0 1]; hold on;
b = bar(2,J_Coll);  b.FaceColor = [1 0 0];
ylabel('Objective value');
set(gca,'XTickLabel',{'Shooting','Collocation'});