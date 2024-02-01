% Thomas C. Day
% Read results from overlap sweep and output the mean radius of each
% overlap threshold used. Then, pick an overlap threshold to use to match
% the data.

%% PA results:
% Visit each folder, load the overlap threshold and the mean radius of the
% clusters simulated.
folder = 'E:\My Drive\01_Data\ProjectPloidy\OverlapSweep_Hookian\PA_params\';
folderlist = dir(folder);
folderlist = folderlist([folderlist.isdir]);
folderlist = folderlist(3:end);

for ff = 1:length(folderlist)
    cd([folder,folderlist(ff).name]);
    X = load('parameters.mat');
    OT(ff) = X.overlap_thresh(ff);
    MR(ff) = mean(X.RAD);
    SR(ff) = std(X.RAD);
    cd(folder);
end

% Plot:
figure; hold on; box on; set(gca,'linewidth',1);
errorbar(OT, MR, SR, 'k.','linewidth',1,'capsize',8);
xlabel('OverlapThresh'); ylabel('Mean Cluster Radius');
title('PA t0 parameters (CellVol, CellAR)');
line([0,1e3],[23.3,23.3],'linestyle','--','color','k','linewidth',1);

% From these results, I use the average of the two points around the black
% line to get the correct overlap threshold:
OT_chosen = (OT(4) + OT(5))/2

%% PM results:
% Visit each folder, load the overlap threshold and the mean radius of the
% clusters simulated.
folder = 'E:\My Drive\01_Data\ProjectPloidy\OverlapSweep_Hookian\PM_params\';
folderlist = dir(folder);
folderlist = folderlist([folderlist.isdir]);
folderlist = folderlist(3:end);

for ff = 1:length(folderlist)
    cd([folder,folderlist(ff).name]);
    X = load('parameters.mat');
    OT(ff) = X.overlap_thresh(ff);
    MR(ff) = mean(X.RAD);
    SR(ff) = std(X.RAD);
    cd(folder);
end

% Plot:
figure; hold on; box on; set(gca,'linewidth',1);
errorbar(OT, MR, SR, 'k.','linewidth',1,'capsize',8);
xlabel('OverlapThresh'); ylabel('Mean Cluster Radius');
title('PM t0 parameters (CellVol, CellAR)');
line([500,2.5e3],[31.2,31.2],'linestyle','--','color','k','linewidth',1);

% From these results, I use the average of the three points around the black
% line to get the correct overlap threshold:
OT_chosen = (OT(7) + OT(7) + OT(8))/3


%% Next experiment analysis:
% Here I simulated 50 clusters with the (CellVol, CellAR) parameters
% measured in experiments, using the threshold overlap discovered above.
% Now, I compare the size at fracture from simulation to the top radius
% measured in experiment.

% PA lines: Visit each folder, retrieve the mean radius of clusters of each
% timepoint, then plot the residual.
trueRads = [23.2, 53.2, 85.1, 202.8, 117.8;
            23.2, 49.3, 54.5, 381.0, 408.3;
            23.2, 56.3, 72.3, 118.4, 467.9;
            23.2, 50.4, 69.5, 210.3, 217.7;
            23.2, 47.8, 80.0, 243.2, 209.4];
folder = 'E:\My Drive\01_Data\ProjectPloidy\OverlapSweep_Hookian\PA_exp\';
cd(folder);
t0_results = load('t000.mat');
folders = {'PA01','PA02','PA03','PA04','PA05'};
simRads_avg = zeros(5, 5); simRads_err = simRads_avg;
simRads_avg(:,1) = mean(t0_results.RAD);
simRads_err(:,1) = std(t0_results.RAD);
for ff = 1:length(folders)
    cd(folders{ff});
    folderlist = dir;
    folderlist = folderlist(3:end);
    for jj = 1:length(folderlist)
        cd(folderlist(jj).name);
        X = load('Results.mat');
        simRads_avg(ff,jj+1) = mean(X.RAD);
        simRads_err(ff,jj+1) = std(X.RAD);
        cd ..;
    end
    cd ..;
end

% Show true divergence:
T = [0,200,400,600,1000];
Colors = summer(5);
figure; hold on; box on; set(gca,'linewidth',1);
for ff = 1:size(trueRads,1)
    plot(T, trueRads(ff,:), '.-','linewidth',1,'color',Colors(ff,:));
    errorbar(T, simRads_avg(ff,:), simRads_err(ff,:), 'x','linewidth',1,'color',Colors(ff,:),'capsize',8);
end
xlabel('Time (days)');
ylabel('Mean Size [\mum]');
title('PA lines');

% Show percentage off:
T = [0,200,400,600,1000];
P = trueRads./simRads_avg - 1;
figure; hold on; box on; set(gca,'linewidth',1);
for ff = 1:size(trueRads,1)
    plot(T, P(ff,:), '.-','linewidth',1,'markersize',10,'color',Colors(ff,:));
end
xlabel('Time (days)');
ylabel('% Off');
title('PA lines');

% -------------------------------------------------------------------------
% PM Lines: Visit each folder, retrieve the mean radius of clusters of each
% timepoint, then plot the residual.
trueRads = [31.2, 49.9, 52.4, 52.8, 55.9;
            31.2, 50.0, 50.8, 52.0, 49.3;
            31.2, 50.3, 72.5, 49.7, 48.3;
            31.2, 55.6, 57.0, 56.7, 60.7;
            31.2, 51.9, 58.1, 80.9, 64.0];
folder = 'E:\My Drive\01_Data\ProjectPloidy\OverlapSweep_Hookian\PM_exp\';
cd(folder);
t0_results = load('t000.mat');
folders = {'PM01','PM02','PM03','PM04','PM05'};
simRads_avg = zeros(5, 5); simRads_err = simRads_avg;
simRads_avg(:,1) = mean(t0_results.RAD);
simRads_err(:,1) = std(t0_results.RAD);
for ff = 1:length(folders)
    cd(folders{ff});
    folderlist = dir;
    folderlist = folderlist(3:end);
    for jj = 1:length(folderlist)
        cd(folderlist(jj).name);
        X = load('Results.mat');
        simRads_avg(ff,jj+1) = mean(X.RAD);
        simRads_err(ff,jj+1) = std(X.RAD);
        cd ..;
    end
    cd ..;
end

% Show:
T = [0,200,400,600,1000];
Colors = winter(5);
figure; hold on; box on; set(gca,'linewidth',1);
for ff = 1:size(trueRads,1)
    plot(T, trueRads(ff,:), '.-','linewidth',1,'color',Colors(ff,:));
    errorbar(T, simRads_avg(ff,:), simRads_err(ff,:), 'x','linewidth',1,'color',Colors(ff,:),'capsize',8);
end
xlabel('Time (days)');
ylabel('Mean Size [\mum]');
title('PM lines');

% Show percentage off:
T = [0,200,400,600,1000];
P = trueRads./simRads_avg - 1;
figure; hold on; box on; set(gca,'linewidth',1);
for ff = 1:size(trueRads,1)
    plot(T, P(ff,:), '.-','linewidth',1,'markersize',10,'color',Colors(ff,:));
end
xlabel('Time (days)');
ylabel('% Off');
title('PM lines');


























