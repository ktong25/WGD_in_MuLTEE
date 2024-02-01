% Thomas C. Day
% Measure size of clusters:
folder = 'H:\My Drive\01_Data\ProjectPloidy\202312\PM_overlap\';
cd(folder);

folderlist = dir();
folderlist = folderlist([folderlist.isdir]);
folderlist = folderlist(3:end);

% Preallocate space:
Nbar = zeros(length(folderlist), 11); Nerr = Nbar;
Vbar = Nbar; Verr = Nbar; Rbar = Vbar; Rerr = Verr;
Vc   = Nbar; AR = Vc;
PackFracbar = Nbar; PackFracerr = Nbar;
KaiData = []; % this will eventually be Kai's data structure
for ff = 1:length(folderlist)
    cd(folderlist(ff).name);
    fprintf([folderlist(ff).name,'\n']);
    
    % Now get every file in folder:
    filelist = dir('thomas-sim*.mat');
    P = load('parameters.mat');
    
    % Measurements:
    for ii = 1:length(filelist)
        X = load(filelist(ii).name);
        
        % Number of cells per cluster:
        N = zeros(length(X.cell_list),1);
        for jj = 1:length(X.cell_list)
            N(jj) = length(X.cell_list{jj});
        end
        Nbar(ff,ii) = mean(N);
        Nerr(ff,ii) = std(N);
        
        % Volume of cluster:
        MaxRad = P.AR(ii)*P.diam(ii)/2;
        V = zeros(length(X.cell_list),1);
        for jj = 1:length(X.cell_list)
            Xc = [X.cell_list{jj}.Center];
            COM = mean(Xc,2);
            Xc = Xc - COM;
            Yc = Xc./vecnorm(Xc);
            Xc = Xc + MaxRad*Yc;
            [~,V(jj)] = convhull(Xc');
            R(jj) = ((3/4)*(1/pi)*V(jj)).^(1/3);
        end
        Vbar(ff,ii) = mean(V);
        Verr(ff,ii) = std(V);
        Rbar(ff,ii) = mean(R);
        Rerr(ff,ii) = std(R);
        
        % Volume of single cells:
        Vc(ff,ii) = P.V;
        AR(ff,ii) = P.AR(ii);
        
        % Packing fraction:
        PackFrac = N*P.V ./ V;
        PackFracbar(ff,ii) = mean(PackFrac);
        PackFracerr(ff,ii) = std(PackFrac);
        
        % Saving for Kai:
        CellVol = Vc(ff,ii)*ones(length(N),1);
        CellAR  = P.AR(ii)*ones(length(N),1);
        KaiData_toadd = [CellVol, CellAR, V, N, PackFrac];
        KaiData = [KaiData; KaiData_toadd]; % append to list
    end   
    
    cd(folder);
end

% Save data for Kai:
T = array2table(KaiData,'VariableNames',{'CellVol','CellAR','V','N','PackFrac'});
writetable(T,'ForKai_all_sim_measurements_PM_overlap.csv');

% Save all other data to matlab file:
save('ForKai_all_sim_measurements_PM_overlap.mat','Nbar','Nerr','Vbar','Verr','PackFracbar','PackFracerr','Rbar','Rerr','Vc','AR');


%% Make plots - PA values:

cd('H:\My Drive\01_Data\ProjectPloidy\202312\PA_overlap\');
load('ForKai_all_sim_measurements_PA_overlap.mat');

% Ncells:
figure('units','centimeters','position',[3,3,15,8]);
subplot(1,2,1); hold on; box on; set(gca,'linewidth',1);
colors = copper(size(Nbar,1));
for ff = 1:size(Nbar,1)
    errorbar(AR(ff,:), Nbar(ff,:), Nerr(ff,:), '.-','markersize',8,'linewidth',1,'color',colors(ff,:),'capsize',8);
end
xlabel('Aspect Ratio');
ylabel('Ncells');
subplot(1,2,2); hold on; box on; set(gca,'linewidth',1);
colors = summer(size(Nbar,2));
for ff = 1:size(Nbar,2)
    errorbar(Vc(:,ff), Nbar(:,ff), Nerr(:,ff), '.-','markersize',8,'linewidth',1,'color',colors(ff,:),'capsize',8);
end
xlabel('Cell Volume');
ylabel('Ncells');
print('PA_params_Ncells_vs_AR-Vc','-dpng','-r600');
close(gcf);

% Volume:
figure('units','centimeters','position',[3,3,15,8]);
subplot(1,2,1); hold on; box on; set(gca,'linewidth',1);
colors = copper(size(Vbar,1));
for ff = 1:size(Vbar,1)
    errorbar(AR(ff,:), Vbar(ff,:), Verr(ff,:), '.-','markersize',8,'linewidth',1,'color',colors(ff,:),'capsize',8);
end
xlabel('Aspect Ratio');
ylabel('Volume');
subplot(1,2,2); hold on; box on; set(gca,'linewidth',1);
colors = summer(size(Vbar,2));
for ff = 1:size(Vbar,2)
    errorbar(Vc(:,ff), Vbar(:,ff), Verr(:,ff), '.-','markersize',8,'linewidth',1,'color',colors(ff,:),'capsize',8);
end
xlabel('Cell Volume');
ylabel('Volume');
print('PA_params_V_vs_AR-Vc','-dpng','-r600');
close(gcf);

% Radius:
figure('units','centimeters','position',[3,3,15,8]);
subplot(1,2,1); hold on; box on; set(gca,'linewidth',1);
colors = copper(size(Rbar,1));
for ff = 1:size(Rbar,1)
    errorbar(AR(ff,:), Rbar(ff,:), Rerr(ff,:), '.-','markersize',8,'linewidth',1,'color',colors(ff,:),'capsize',8);
end
xlabel('Aspect Ratio');
ylabel('Radius [\mum]');
subplot(1,2,2); hold on; box on; set(gca,'linewidth',1);
colors = summer(size(Rbar,2));
for ff = 1:size(Rbar,2)
    errorbar(Vc(:,ff), Rbar(:,ff), Rerr(:,ff), '.-','markersize',8,'linewidth',1,'color',colors(ff,:),'capsize',8);
end
xlabel('Cell Volume');
ylabel('Radius [\mum]');
print('PA_params_R_vs_AR-Vc','-dpng','-r600');
close(gcf);

% Packing Fraction:
figure('units','centimeters','position',[3,3,15,8]);
subplot(1,2,1); hold on; box on; set(gca,'linewidth',1);
colors = copper(size(PackFracbar,1));
for ff = 1:size(PackFracbar,1)
    errorbar(AR(ff,:), PackFracbar(ff,:), PackFracerr(ff,:), '.-','markersize',8,'linewidth',1,'color',colors(ff,:),'capsize',8);
end
xlabel('Aspect Ratio');
ylabel('PackingFraction');
subplot(1,2,2); hold on; box on; set(gca,'linewidth',1);
colors = summer(size(PackFracbar,2));
for ff = 1:size(PackFracbar,2)
    errorbar(Vc(:,ff), PackFracbar(:,ff), PackFracerr(:,ff), '.-','markersize',8,'linewidth',1,'color',colors(ff,:),'capsize',8);
end
xlabel('Cell Volume');
ylabel('Packing Fraction');
print('PA_params_Phi_vs_AR-Vc','-dpng','-r600');
close(gcf);

% Heatmaps:
figure('units','centimeters','position',[3,3,15,8]);
subplot(1,2,1); hold on; box on; set(gca,'linewidth',1);
imagesc([1,2],[14,65],Nbar);
cb = colorbar;
set(cb,'Linewidth',1);
xlabel('AR'); ylabel('V_c');

subplot(1,2,2); hold on; box on; set(gca,'linewidth',1);
imagesc([1,2],[14,65],Vbar);
colormap('gray');
cb = colorbar;
set(cb,'linewidth',1);
xlabel('AR'); ylabel('V_c');
print('PA_params_N-V_vs_AR-Vc','-dpng','-r600');
close(gcf);

%% Make plots - PM values:

cd('H:\My Drive\01_Data\ProjectPloidy\202312\PM_overlap\');
load('ForKai_all_sim_measurements_PM_overlap.mat');

% Ncells:
figure('units','centimeters','position',[3,3,15,8]);
subplot(1,2,1); hold on; box on; set(gca,'linewidth',1);
colors = copper(size(Nbar,1));
for ff = 1:size(Nbar,1)
    errorbar(AR(ff,:), Nbar(ff,:), Nerr(ff,:), '.-','markersize',8,'linewidth',1,'color',colors(ff,:),'capsize',8);
end
xlabel('Aspect Ratio');
ylabel('Ncells');
subplot(1,2,2); hold on; box on; set(gca,'linewidth',1);
colors = summer(size(Nbar,2));
for ff = 1:size(Nbar,2)
    errorbar(Vc(:,ff), Nbar(:,ff), Nerr(:,ff), '.-','markersize',8,'linewidth',1,'color',colors(ff,:),'capsize',8);
end
xlabel('Cell Volume');
ylabel('Ncells');
print('PM_params_Ncells_vs_AR-Vc','-dpng','-r600');
close(gcf);

% Volume:
figure('units','centimeters','position',[3,3,15,8]);
subplot(1,2,1); hold on; box on; set(gca,'linewidth',1);
colors = copper(size(Vbar,1));
for ff = 1:size(Vbar,1)
    errorbar(AR(ff,:), Vbar(ff,:), Verr(ff,:), '.-','markersize',8,'linewidth',1,'color',colors(ff,:),'capsize',8);
end
xlabel('Aspect Ratio');
ylabel('Volume');
subplot(1,2,2); hold on; box on; set(gca,'linewidth',1);
colors = summer(size(Vbar,2));
for ff = 1:size(Vbar,2)
    errorbar(Vc(:,ff), Vbar(:,ff), Verr(:,ff), '.-','markersize',8,'linewidth',1,'color',colors(ff,:),'capsize',8);
end
xlabel('Cell Volume');
ylabel('Volume');
print('PM_params_V_vs_AR-Vc','-dpng','-r600');
close(gcf);

% Radius:
figure('units','centimeters','position',[3,3,15,8]);
subplot(1,2,1); hold on; box on; set(gca,'linewidth',1);
colors = copper(size(Rbar,1));
for ff = 1:size(Rbar,1)
    errorbar(AR(ff,:), Rbar(ff,:), Rerr(ff,:), '.-','markersize',8,'linewidth',1,'color',colors(ff,:),'capsize',8);
end
xlabel('Aspect Ratio');
ylabel('Radius [\mum]');
subplot(1,2,2); hold on; box on; set(gca,'linewidth',1);
colors = summer(size(Rbar,2));
for ff = 1:size(Rbar,2)
    errorbar(Vc(:,ff), Rbar(:,ff), Rerr(:,ff), '.-','markersize',8,'linewidth',1,'color',colors(ff,:),'capsize',8);
end
xlabel('Cell Volume');
ylabel('Radius [\mum]');
print('PM_params_R_vs_AR-Vc','-dpng','-r600');
close(gcf);

% Packing Fraction:
figure('units','centimeters','position',[3,3,15,8]);
subplot(1,2,1); hold on; box on; set(gca,'linewidth',1);
colors = copper(size(PackFracbar,1));
for ff = 1:size(PackFracbar,1)
    errorbar(AR(ff,:), PackFracbar(ff,:), PackFracerr(ff,:), '.-','markersize',8,'linewidth',1,'color',colors(ff,:),'capsize',8);
end
xlabel('Aspect Ratio');
ylabel('PackingFraction');
subplot(1,2,2); hold on; box on; set(gca,'linewidth',1);
colors = summer(size(PackFracbar,2));
for ff = 1:size(PackFracbar,2)
    errorbar(Vc(:,ff), PackFracbar(:,ff), PackFracerr(:,ff), '.-','markersize',8,'linewidth',1,'color',colors(ff,:),'capsize',8);
end
xlabel('Cell Volume');
ylabel('Packing Fraction');
print('PM_params_Phi_vs_AR-Vc','-dpng','-r600');
close(gcf);

% Heatmaps:
figure('units','centimeters','position',[3,3,15,8]);
subplot(1,2,1); hold on; box on; set(gca,'linewidth',1);
imagesc([1,2],[14,65],Nbar);
cb = colorbar;
set(cb,'Linewidth',1);
xlabel('AR'); ylabel('V_c');

subplot(1,2,2); hold on; box on; set(gca,'linewidth',1);
imagesc([1,2],[14,65],Vbar);
colormap('gray');
cb = colorbar;
set(cb,'linewidth',1);
xlabel('AR'); ylabel('V_c');
print('PM_params_N-V_vs_AR-Vc','-dpng','-r600');
close(gcf);