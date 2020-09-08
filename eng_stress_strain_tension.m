% prgoram to compare engineering stress and engineering strain
clear
clc

% open stress-strain data from reference
fid1 = fopen('./Sim.dat','r');
A = fscanf(fid1,'%f');
A = transpose(reshape(A,[2, length(A)/2]));
A(:,2) = A(:,2)./1e6; % pa --> Mpa

% open stress-strain data from Yang's model
fid2 = fopen('~/damin_research/NASA/Damin/DDEHM/UMAT_Hex_phase/Sim.dat','r');
B = fscanf(fid2,'%f');
B = transpose(reshape(B,[2, length(B)/2]));
B(:,2) = B(:,2)./1e6; % pa --> Mpa


% Compare two stress strain curve obtained from different element type
figure(1)
hold on
grid on
box on
plot(A(:,1),A(:,2),'-.r','LineWidth',2);
plot(B(:,1),B(:,2),'-.k','LineWidth',2);
legend('reference model','DDEHM','location','best')
title('RVE9 tension')
xlabel("Strain");
ylabel("Stress (Mpa)")
ax = gca;
ax.FontSize = 15;
%pbaspect([1.65 1 1])
%xlim([0 0.5]);
hold off

% Save figures
saveflag = 0;
if(saveflag == 1) 
    saveas(figure(1),'unitension','fig');
    saveas(figure(1),'unitension','jpg');
    saveas(figure(2),'separation','fig');
    saveas(figure(2),'separation','jpg');
    saveas(figure(3),'traction','fig');
    saveas(figure(3),'traction','jpg');
    saveas(figure(4),'traction_separation','fig');
    saveas(figure(4),'traction_separation','jpg');
end