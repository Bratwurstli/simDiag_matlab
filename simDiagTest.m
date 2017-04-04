%% simDiag test
% test script for simDiag

%%
% generate test matrices

xCoordinates = [1 5 2 2 3 3];
yCoordinates = [-1 2 -1 3 -1 2];
matrix1 = diag(xCoordinates);
matrix2 = diag(yCoordinates);

%%
% solution are 2d points given by xCoordinates and yCoordinates

plot(xCoordinates,yCoordinates,'ob','Markersize',10,'Linewidth',2);
xlim([0 6]);
ylim([-2 4]);
hold on;
title('simDiag testplot');
legend('solution');

%%
% generate random matrices randMat1 and randMat2 that have a simultaneous diagonalization

B = rand(length(xCoordinates));
randMat1 = B\matrix1*B
randMat2 = B\matrix2*B

%%
% time measure and simultaneous diagonalization of randMat1 and randMat2

tic 
[basis values] = simDiag([randMat1 randMat2]);
toc

%%
% plot solution

plot(values(:,1),values(:,2),'+r','Linewidth', 1.5,'Markersize', 8);
legend('solution','simDiag');
hold off;
