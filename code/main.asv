clear; close all;

%% Generate Data
currents = [0, 20, 50,20, 100];

BCs = [10, -5; 
       1, 1;  
      -10, 20;
        5,15;
      -12,14;
       20, 20];

sensors = [0, 0, 0 ; 
           0, 5, 10];

[Data, tlist, FE_Matrices, SVD_Matrices] = genData(currents, BCs, sensors);

%% Find Singular Values
sigma = zeros(length(SVD_Matrices.Sc),1);
for i = 1:length(SVD_Matrices.Sc)
    sigma(i) = SVD_Matrices.Sc(i,i);
end

figure 
plot(sigma);
xlabel('Index')
ylabel('Explanatory Power of \sigma')
title('Number of Modes vs. Explanatory power of \sigma')

%% Project number of modes chosen to FE Matrices
Urc = SVD_Matrices.Uc(:,1:8); % projection matrix
Reduced_FE_Matrices.Mrc = Urc' * FE_Matrices.Mc * Urc;
Reduced_FE_Matrices.Krc = Urc' * FE_Matrices.Kc * Urc;
Reduced_FE_Matrices.Lr = FE_Matrices.L * Urc;

Reduced_FE_Matrices.Mrc_inv = inv(Reduced_FE_Matrices.Mrc);

%% Project modes to data
Reduced_Data = cell(size(Data));
for i = 1:size(Data, 1)
    Reduced_Data(i,1) = Data(i,1);
    Reduced_Data{i,2} = Urc' * Data{i,2};
    Reduced_Data{i,3} = Urc' * Data{i,3};    
end

%% Saving data
save("data.mat", "Reduced_Data")
save("R_FE_Mats.mat","Reduced_FE_Matrices")
save("FOM_Data", "Data")

%% 
figure
u = Data{27,2};
u = FE_Matrices.L * u;
plot(tlist, u(1,:))
hold on
plot(tlist, u(2,:))
plot(tlist, u(3,:))
