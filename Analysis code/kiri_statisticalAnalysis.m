% kiri_statisticalAnalysis calls the functions which ananlyse positions and 
% angles to give kinematic variables and statistically analyses them. Summary 
% statistics are calculated for each trial, for each individual and for 
% total data.
%
% Kiri Pullar, masters thesis 2009


%% load files and analyse tracking results
clear all;
clf;
%whatever folder tracking .mat files are in
savedir = cd('F:\Personal Data\My Folders\CurrentUniWork\From small\Results'); 
[filename, pathname] = uigetfile('*.mat', 'Select ProjectFiles','multiselect','on');
cd(savedir);

% If only one file was selected, make it a cell array.
if ~iscell(filename)
    filename = {filename};
end

% ...for all selected files...
for i = 1:size(filename,2)
    % ...load the project...
    load([pathname,filename{1,i}]);
    % ...and display the loaded project as this is batch processing and
    % we want to see how far we have got when we return later.
    %     disp(['project file ' filename{1,i} ' loaded ' datestr(now)]);
    [data(i).individual,data(i).weight,data(i).bodylength,data(i).carapacelength,data(i).speedcm,data(i).maxspeedcm,data(i).minspeedcm,data(i).speedbodylength,data(i).maxspeedbodylength,data(i).minspeedbodylength,data(i).maxjointangles, data(i).minjointangles, data(i).rangejointangles, data(i).legangle, data(i).legangluarvelocity, data(i).legangularacceleration, data(i).swingduration,data(i).stanceduration,data(i).strideperiod,data(i).legsstrideperiod,data(i).L1phase,data(i).R1phase,data(i).L2phase,data(i).R2phase,data(i).L3phase,data(i).R3phase,data(i).L4phase,data(i).R4phase, data(i).legdist, data(i).angexcur,data(i).swinglegvel,data(i).stancelegvel] = kiri_kinematicAnalysis(inputfile,model,video,filtering,tracking);

    %free up memory
    inputfile=[];
    video=[];
    model=[];
    optimization=[];
    tracking=[];
    filtering=[];

    % Display the analysis was successful.
    display(['File ' filename{1,i} ' analysed']);
end

%% summary stats for all animals in all trials
trial.averagespeedcm=zeros(length(data),1);
trial.maxspeedcm=zeros(length(data),1);
trial.minspeedcm=zeros(length(data),1);
trial.averagespeedbodylength=zeros(length(data),1);
trial.maxspeedbodylength=zeros(length(data),1);
trial.minspeedbodylength=zeros(length(data),1);
trial.speedcm=[];
trial.maxjointangles=[];
trial.minjointangles=[];
trial.rangejointangles=[];


for i=1:length(data)
    trial.individual(i)=data(i).individual;
    trial.speedcm{i}=data(i).speedcm;
    trial.averagespeedcm(i)=mean(data(i).speedcm);
    trial.maxspeedcm(i)=data(i).maxspeedcm;
    trial.minspeedcm(i)=data(i).minspeedcm;
    trial.averagespeedbodylength(i)=mean(data(i).speedbodylength);
    trial.speedbodylength{i}=data(i).speedbodylength;
    trial.maxspeedbodylength(i)=data(i).maxspeedbodylength;
    trial.minspeedbodylength(i)=data(i).minspeedbodylength;
    %     figure(1)
    %     hold on
    %     plot(data(i).weight,mean(data(i).speedcm),'ko','MarkerFaceColor','k')
    %     hold off
    %     xlabel('Mass (g)')
    %     ylabel('Speed (cm/s)')
    %
    
    trial.maxjointangles=[trial.maxjointangles data(i).maxjointangles];
    trial.minjointangles=[trial.minjointangles data(i).minjointangles];
    trial.rangejointangles=[trial.rangejointangles data(i).rangejointangles];

    trial.strideperiod{i}=data(i).strideperiod';
    trial.stridefrequency{i}=1./data(i).strideperiod';
    trial.swingduration{i}=data(i).swingduration;
    trial.stanceduration{i}=data(i).stanceduration;
    trial.stridelength{i}=mean(data(i).speedcm)./(1./data(i).strideperiod)';
    trial.protraction{i}=data(i).stanceduration{3};
    trial.retraction{i}=data(i).swingduration{3};
    trial.L1phase{i}=data(i).L1phase;
    trial.R1phase{i}=data(i).R1phase;
    trial.L2phase{i}=data(i).L2phase;
    trial.R2phase{i}=data(i).R2phase;
    trial.L3phase{i}=data(i).L3phase;
    trial.R3phase{i}=data(i).R3phase;
    trial.L4phase{i}=data(i).L4phase;
    trial.R4phase{i}=data(i).R4phase;
    
    trial.leg1dist{i}=data(i).legdist{1};
    trial.leg2dist{i}=data(i).legdist{2};
    trial.leg3dist{i}=data(i).legdist{3};
    trial.leg4dist{i}=data(i).legdist{4};
    trial.leg1angexcur{i}=data(i).angexcur{1};
    trial.leg2angexcur{i}=data(i).angexcur{2};
    trial.leg3angexcur{i}=data(i).angexcur{3};
    trial.leg4angexcur{i}=data(i).angexcur{4};
    %     trial.dutyfactor{i}=data(i).stanceduration{3}./data(i).strideperiod;
    %     some may have more stances than strides
    
    trial.leg1strideperiod{i}=data(i).legsstrideperiod{1};
    trial.leg2strideperiod{i}=data(i).legsstrideperiod{2};
    trial.leg3strideperiod{i}=data(i).legsstrideperiod{3};
    trial.leg4strideperiod{i}=data(i).legsstrideperiod{4};
    trial.leg1stridefrequency{i}=1/data(i).legsstrideperiod{1};
    trial.leg2stridefrequency{i}=1/data(i).legsstrideperiod{2};
    trial.leg3stridefrequency{i}=1/data(i).legsstrideperiod{3};
    trial.leg4stridefrequency{i}=1/data(i).legsstrideperiod{4};
    
    trial.leg1swingvel{i}=data(i).swinglegvel{1};
    trial.leg2swingvel{i}=data(i).swinglegvel{2};
    trial.leg3swingvel{i}=data(i).swinglegvel{3};
    trial.leg4swingvel{i}=data(i).swinglegvel{4};
    trial.leg1stancevel{i}=data(i).stancelegvel{1};
    trial.leg2stancevel{i}=data(i).stancelegvel{2};
    trial.leg3stancevel{i}=data(i).stancelegvel{3};
    trial.leg4stancevel{i}=data(i).stancelegvel{4};

    trial.averagestrideperiod(i)=mean(trial.strideperiod{i});
    trial.averagestridefrequency(i)=mean(trial.stridefrequency{i});
    trial.averagestridelength(i)=mean(trial.stridelength{i});
    trial.averageprotraction(i)=mean(trial.protraction{i});
    trial.averageretraction(i)=mean(trial.retraction{i});
    trial.averagedutyfactor(i)=mean(data(i).stanceduration{3})/mean(data(i).strideperiod);
    trial.averagedutyfactor1(i)=mean(data(i).stanceduration{1})/mean(data(i).strideperiod);
    trial.averagedutyfactor2(i)=mean(data(i).stanceduration{3})/mean(data(i).strideperiod);
    trial.averagedutyfactor3(i)=mean(data(i).stanceduration{5})/mean(data(i).strideperiod);
    trial.averagedutyfactor4(i)=mean(data(i).stanceduration{7})/mean(data(i).strideperiod);
    trial.averagestabilitymargin(i)=(8/4)*trial.averagedutyfactor(i)- 3/4;
    [trial.averageL1phase{i},trial.stdL1phase{i},trial.pL1phase{i}]=cellfun(@kiri_circStats,data(i).L1phase);
    [trial.averageR1phase{i},trial.stdR1phase{i},trial.pR1phase{i}]=cellfun(@kiri_circStats,data(i).R1phase);
    [trial.averageL2phase{i},trial.stdL2phase{i},trial.pL2phase{i}]=cellfun(@kiri_circStats,data(i).L2phase);
    [trial.averageR2phase{i},trial.stdR2phase{i},trial.pR2phase{i}]=cellfun(@kiri_circStats,data(i).R2phase);
    [trial.averageL3phase{i},trial.stdL3phase{i},trial.pL3phase{i}]=cellfun(@kiri_circStats,data(i).L3phase);
    [trial.averageR3phase{i},trial.stdR3phase{i},trial.pR3phase{i}]=cellfun(@kiri_circStats,data(i).R3phase);
    [trial.averageL4phase{i},trial.stdL4phase{i},trial.pL4phase{i}]=cellfun(@kiri_circStats,data(i).L4phase);
    [trial.averageR4phase{i},trial.stdR4phase{i},trial.pR4phase{i}]=cellfun(@kiri_circStats,data(i).R4phase);
    trial.averageleg1dist(i)=mean(trial.leg1dist{i});
    trial.averageleg2dist(i)=mean(trial.leg2dist{i});
    trial.averageleg3dist(i)=mean(trial.leg3dist{i});
    trial.averageleg4dist(i)=mean(trial.leg4dist{i});
    trial.averageleg1angexcur(i)=mean(trial.leg1angexcur{i});
    trial.averageleg2angexcur(i)=mean(trial.leg2angexcur{i});
    trial.averageleg3angexcur(i)=mean(trial.leg3angexcur{i});
    trial.averageleg4angexcur(i)=mean(trial.leg4angexcur{i});
    
    trial.averageleg1strideperiod(i)=mean(trial.leg1strideperiod{i});
    trial.averageleg2strideperiod(i)=mean(trial.leg2strideperiod{i});
    trial.averageleg3strideperiod(i)=mean(trial.leg3strideperiod{i});
    trial.averageleg4strideperiod(i)=mean(trial.leg4strideperiod{i});
    
    trial.averageleg1stridefrequency(i)=mean(trial.leg1stridefrequency{i});
    trial.averageleg2stridefrequency(i)=mean(trial.leg2stridefrequency{i});
    trial.averageleg3stridefrequency(i)=mean(trial.leg3stridefrequency{i});
    trial.averageleg4stridefrequency(i)=mean(trial.leg4stridefrequency{i});

    trial.averageleg1swingvel(i)=mean(trial.leg1swingvel{i});
    trial.averageleg2swingvel(i)=mean(trial.leg2swingvel{i});
    trial.averageleg3swingvel(i)=mean(trial.leg3swingvel{i});
    trial.averageleg4swingvel(i)=mean(trial.leg4swingvel{i});
    trial.averageleg1stancevel(i)=mean(trial.leg1stancevel{i});
    trial.averageleg2stancevel(i)=mean(trial.leg2stancevel{i});
    trial.averageleg3stancevel(i)=mean(trial.leg3stancevel{i});
    trial.averageleg4stancevel(i)=mean(trial.leg4stancevel{i});
    
    %% arrange data into individuals
    subject=[];
    subject(trial.individual(1)).weight=data(1).weight;
    subject(trial.individual(1)).bodylength=data(1).bodylength;
    subject(trial.individual(1)).carapacelength=data(1).carapacelength;
    subject(trial.individual(1)).speedcm(:,:)=data(1).speedcm;
    subject(trial.individual(1)).speedbodylength(:,:)=data(1).speedbodylength;

    subject(trial.individual(1)).strideperiod=data(1).strideperiod;
    subject(trial.individual(1)).stridefrequency=1./data(1).strideperiod;
    subject(trial.individual(1)).stridelength=mean(data(1).speedcm)./(1./data(1).strideperiod);
    subject(trial.individual(1)).swingduration(:,:)=data(1).swingduration;
    subject(trial.individual(1)).stanceduration(:,:)=data(1).stanceduration;


    subject(trial.individual(1)).maxjointangles(:,:)=data(1).maxjointangles;
    subject(trial.individual(1)).minjointangles(:,:)=data(1).minjointangles;
    subject(trial.individual(1)).rangejointangles(:,:)=data(1).rangejointangles;


    for i=2:length(trial.individual)
        subject(trial.individual(i)).weight=data(i).weight;
        if trial.individual(i)==trial.individual(i-1)
            subject(trial.individual(i)).bodylength=[subject(trial.individual(i)).bodylength; data(i).bodylength];
            subject(trial.individual(i)).carapacelength=[subject(trial.individual(i)).carapacelength; data(i).carapacelength];
            subject(trial.individual(i)).speedcm=[subject(trial.individual(i)).speedcm; data(i).speedcm];
            subject(trial.individual(i)).speedbodylength=[subject(trial.individual(i)).speedbodylength; data(i).speedbodylength];

            subject(trial.individual(i)).strideperiod=[subject(trial.individual(i)).strideperiod; data(i).strideperiod];
            subject(trial.individual(i)).stridefrequency=[subject(trial.individual(i)).stridefrequency; 1./data(i).strideperiod];
            subject(trial.individual(i)).stridelength=[subject(trial.individual(i)).stridelength; mean(data(i).speedcm)./(1./data(i).strideperiod)];
            subject(trial.individual(i)).swingduration=[subject(trial.individual(i)).swingduration; data(i).swingduration];
            subject(trial.individual(i)).stanceduration=[subject(trial.individual(i)).stanceduration; data(i).stanceduration];

            subject(trial.individual(i)).maxjointangles=[subject(trial.individual(i)).maxjointangles data(i).maxjointangles];
            subject(trial.individual(i)).minjointangles=[subject(trial.individual(i)).minjointangles data(i).minjointangles];
            subject(trial.individual(i)).rangejointangles=[subject(trial.individual(i)).rangejointangles data(i).rangejointangles];
        else
            subject(trial.individual(i)).bodylength=data(i).bodylength;
            subject(trial.individual(i)).carapacelength=data(i).carapacelength;
            subject(trial.individual(i)).speedcm=data(i).speedcm;
            subject(trial.individual(i)).speedbodylength=data(i).speedbodylength;

            subject(trial.individual(i)).strideperiod=data(i).strideperiod;
            subject(trial.individual(i)).stridefrequency=1./data(i).strideperiod;
            subject(trial.individual(i)).stridelength=mean(data(i).speedcm)./(1./data(i).strideperiod);
            subject(trial.individual(i)).swingduration=data(i).swingduration;
            subject(trial.individual(i)).stanceduration=data(i).stanceduration;

            subject(trial.individual(i)).maxjointangles(:,:)=data(i).maxjointangles;
            subject(trial.individual(i)).minjointangles(:,:)=data(i).minjointangles;
            subject(trial.individual(i)).rangejointangles(:,:)=data(i).rangejointangles;
        end
    end
end

total.meanspeedcm=mean(trial.averagespeedcm);
total.stdevspeedcm=std(trial.averagespeedcm);
total.rangespeedcm=[min(trial.minspeedcm) max(trial.maxspeedcm)];
total.meanspeedbodylength=mean(trial.averagespeedbodylength);
total.stdevspeedbodylength=std(trial.averagespeedbodylength);
total.rangespeedbodylength=[min(trial.minspeedbodylength) max(trial.maxspeedbodylength)];

total.meanmaxjointangles=round(mean(trial.maxjointangles,2));
total.stdmaxjointangles=round(std(trial.maxjointangles,0,2));
total.meanminjointangles=round(mean(trial.minjointangles,2));
total.stdminjointangles=round(std(trial.minjointangles,0,2));
total.meanrangejointangles=round(mean(trial.rangejointangles,2));
total.stdrangejointangles=round(std(trial.rangejointangles,0,2));

total.l1_r1phase=[];
for i=[1:length(trial.L1phase)]
total.l1_r1phase=[total.l1_r1phase; trial.L1phase{1,i}{1,2}];
end
% figure;
% hist(total.l1_r1phase)
% ylabel('f [\phi]')
% box off
% h = findobj(gca,'Type','patch');  
% set(h,'FaceColor','w','EdgeColor','k'); 
[total.meanl1_r1phase, total.stdl1_r1phase, total.pl1_r1phase]=kiri_circStats(total.l1_r1phase);
% hold on
% plot([total.meanl1_r1phase total.meanl1_r1phase],[0 25],'k:')
% hold off

total.l2_r2phase=[];
for i=[1:length(trial.L2phase)]
total.l2_r2phase=[total.l2_r2phase; trial.L2phase{1,i}{1,4}];
end
% figure;
% hist(total.l2_r2phase)
% ylabel('f [\phi]')
% box off
% h = findobj(gca,'Type','patch');  
% set(h,'FaceColor','w','EdgeColor','k');
[total.meanl2_r2phase, total.stdl2_r2phase, total.pl2_r2phase]=kiri_circStats(total.l2_r2phase);
% hold on
% plot([total.meanl2_r2phase total.meanl2_r2phase],[0 25],'k:')
% hold off

total.l3_r3phase=[];
for i=[1:length(trial.L3phase)]
total.l3_r3phase=[total.l3_r3phase; trial.L3phase{1,i}{1,6}];
end
% figure;
% hist(total.l3_r3phase)
% ylabel('f [\phi]')
% box off
% h = findobj(gca,'Type','patch');  
% set(h,'FaceColor','w','EdgeColor','k');
[total.meanl3_r3phase, total.stdl3_r3phase, total.pl3_r3phase]=kiri_circStats(total.l3_r3phase);
% hold on
% plot([total.meanl3_r3phase total.meanl3_r3phase],[0 25],'k:')
% hold off

total.l4_r4phase=[];
for i=[1:length(trial.L4phase)]
total.l4_r4phase=[total.l4_r4phase; trial.L4phase{1,i}{1,8}];
end
% figure;
% hist(total.l4_r4phase)
% ylabel('f [\phi]')
% box off
% h = findobj(gca,'Type','patch');  
% set(h,'FaceColor','w','EdgeColor','k');
[total.meanl4_r4phase, total.stdl4_r4phase, total.pl4_r4phase]=kiri_circStats(total.l4_r4phase);
% hold on
% plot([total.meanl4_r4phase total.meanl4_r4phase],[0 25],'k:')
% hold off

total.l1_l2phase=[];
for i=[1:length(trial.L1phase)]
total.l1_l2phase=[total.l1_l2phase; trial.L1phase{1,i}{1,3}];
end
% figure;
% hist(total.l1_l2phase)
% ylabel('f [\phi]')
% box off
% h = findobj(gca,'Type','patch');  
% set(h,'FaceColor','w','EdgeColor','k');
[total.meanl1_l2phase, total.stdl1_l2phase, total.pl1_l2phase]=kiri_circStats(total.l1_l2phase);
% hold on
% plot([total.meanl1_l2phase total.meanl1_l2phase],[0 35],'k:')
% hold off

total.r1_r2phase=[];
for i=[1:length(trial.R1phase)]
total.r1_r2phase=[total.r1_r2phase; trial.R1phase{1,i}{1,4}];
end
% figure;
% hist(total.r1_r2phase)
% ylabel('f [\phi]')
% box off
% h = findobj(gca,'Type','patch');  
% set(h,'FaceColor','w','EdgeColor','k');
[total.meanr1_r2phase, total.stdr1_r2phase, total.pr1_r2phase]=kiri_circStats(total.r1_r2phase);
% hold on
% plot([total.meanr1_r2phase total.meanr1_r2phase],[0 35],'k:')
% hold off

total.l2_l3phase=[];
for i=[1:length(trial.L2phase)]
total.l2_l3phase=[total.l2_l3phase; trial.L2phase{1,i}{1,5}];
end
% figure;
% hist(total.l2_l3phase)
% ylabel('f [\phi]')
% box off
% h = findobj(gca,'Type','patch');  
% set(h,'FaceColor','w','EdgeColor','k');
[total.meanl2_l3phase, total.stdl2_l3phase, total.pl2_l3phase]=kiri_circStats(total.l2_l3phase);
% hold on
% plot([total.meanl2_l3phase total.meanl2_l3phase],[0 35],'k:')
% hold off

total.r2_r3phase=[];
for i=[1:length(trial.R2phase)]
total.r2_r3phase=[total.r2_r3phase; trial.R2phase{1,i}{1,6}];
end
% figure;
% hist(total.r2_r3phase)
% ylabel('f [\phi]')
% box off
% h = findobj(gca,'Type','patch');  
% set(h,'FaceColor','w','EdgeColor','k');
[total.meanr2_r3phase, total.stdr2_r3phase, total.pr2_r3phase]=kiri_circStats(total.r2_r3phase);
% hold on
% plot([total.meanr2_r3phase total.meanr2_r3phase],[0 35],'k:')
% hold off

total.l3_l4phase=[];
for i=[1:length(trial.L3phase)]
total.l3_l4phase=[total.l3_l4phase; trial.L3phase{1,i}{1,7}];
end
% figure;
% hist(total.l3_l4phase)
% ylabel('f [\phi]')
% box off
% h = findobj(gca,'Type','patch');  
% set(h,'FaceColor','w','EdgeColor','k');
[total.meanl3_l4phase, total.stdl3_l4phase, total.pl3_l4phase]=kiri_circStats(total.l3_l4phase);
% hold on
% plot([total.meanl3_l4phase total.meanl3_l4phase],[0 35],'k:')
% hold off

total.r3_r4phase=[];
for i=[1:length(trial.R3phase)]
total.r3_r4phase=[total.r3_r4phase; trial.R3phase{1,i}{1,8}];
end
% figure;
% hist(total.r3_r4phase)
% ylabel('f [\phi]')
% box off
% h = findobj(gca,'Type','patch');  
% set(h,'FaceColor','w','EdgeColor','k');
[total.meanr3_r4phase, total.stdr3_r4phase, total.pr3_r4phase]=kiri_circStats(total.r3_r4phase);
% hold on
% plot([total.meanr3_r4phase total.meanr3_r4phase],[0 35],'k:')
% hold off

total.l1_l3phase=[];
for i=[1:length(trial.L1phase)]
total.l1_l3phase=[total.l1_l3phase; trial.L1phase{1,i}{1,5}];
end
% figure;
% hist(total.l1_l3phase)
% ylabel('f [\phi]')
% box off
% h = findobj(gca,'Type','patch');  
% set(h,'FaceColor','w','EdgeColor','k');
[total.meanl1_l3phase, total.stdl1_l3phase, total.pl1_l3phase]=kiri_circStats(total.l1_l3phase);
% hold on
% plot([total.meanl1_l3phase total.meanl1_l3phase],[0 35],'k:')
% hold off

total.r1_r3phase=[];
for i=[1:length(trial.R1phase)]
total.r1_r3phase=[total.r1_r3phase; trial.R1phase{1,i}{1,6}];
end
% figure;
% hist(total.r1_r3phase)
% ylabel('f [\phi]')
% box off
% h = findobj(gca,'Type','patch');  
% set(h,'FaceColor','w','EdgeColor','k');
[total.meanr1_r3phase, total.stdr1_r3phase, total.pr1_r3phase]=kiri_circStats(total.r1_r3phase);
% hold on
% plot([total.meanr1_r3phase total.meanr1_r3phase],[0 35],'k:')
% hold off

total.l2_l4phase=[];
for i=[1:length(trial.L2phase)]
total.l2_l4phase=[total.l2_l4phase; trial.L2phase{1,i}{1,7}];
end
% figure;
% hist(total.l2_l4phase)
% ylabel('f [\phi]')
% box off
% h = findobj(gca,'Type','patch');  
% set(h,'FaceColor','w','EdgeColor','k');
[total.meanl2_l4phase, total.stdl2_l4phase, total.pl2_l4phase]=kiri_circStats(total.l2_l4phase);
% hold on
% plot([total.meanl2_l4phase total.meanl2_l4phase],[0 35],'k:')
% hold off

total.r2_r4phase=[];
for i=[1:length(trial.R2phase)]
total.r2_r4phase=[total.r2_r4phase; trial.R2phase{1,i}{1,8}];
end
% figure;
% hist(total.r2_r4phase)
% ylabel('f [\phi]')
% box off
% h = findobj(gca,'Type','patch');  
% set(h,'FaceColor','w','EdgeColor','k');
[total.meanr2_r4phase, total.stdr2_r4phase, total.pr2_r4phase]=kiri_circStats(total.r2_r4phase);
% hold on
% plot([total.meanr2_r4phase total.meanr2_r4phase],[0 35],'k:')
% hold off

total.l1_l4phase=[];
for i=[1:length(trial.L1phase)]
total.l1_l4phase=[total.l1_l4phase; trial.L1phase{1,i}{1,7}];
end
% figure;
% hist(total.l1_l4phase)
% ylabel('f [\phi]')
% box off
% h = findobj(gca,'Type','patch');  
% set(h,'FaceColor','w','EdgeColor','k');
[total.meanl1_l4phase, total.stdl1_l4phase, total.pl1_l4phase]=kiri_circStats(total.l1_l4phase);
% hold on
% plot([total.meanl1_l4phase total.meanl1_l4phase],[0 35],'k:')
% hold off

total.r1_r4phase=[];
for i=[1:length(trial.R1phase)]
total.r1_r4phase=[total.r1_r4phase; trial.R1phase{1,i}{1,8}];
end
% figure;
% hist(total.r1_r4phase)
% ylabel('f [\phi]')
% box off
% h = findobj(gca,'Type','patch');  
% set(h,'FaceColor','w','EdgeColor','k');
[total.meanr1_r4phase, total.stdr1_r4phase, total.pr1_r4phase]=kiri_circStats(total.r1_r4phase);
% hold on
% plot([total.meanr1_r4phase total.meanr1_r4phase],[0 35],'k:')
% hold off

subjectstats.weight=[];
subjectstats.bodylength=[];
subjectstats.carapacelength=[];
subjectstats.averagespeedcm=[];
subjectstats.maxspeedcm=[];
subjectstats.minspeedcm=[];
subjectstats.averagespeedbodylength=[];
subjectstats.maxspeedbodylength=[];
subjectstats.minspeedbodylength=[];
subjectstats.strideperiod=[];
subjectstats.averagestridefrequency=[];
subjectstats.averagestridelength=[];
subjectstats.maxstridefrequency=[];
subjectstats.maxstridelength=[];
subjectstats.minstridefrequency=[];
subjectstats.minstridelength=[];
subjectstats.maxjointangles=[];
subjectstats.minjointangles=[];
subjectstats.rangejointangles=[];
subjectstats.stdmaxjointangles=[];
subjectstats.stdminjointangles=[];
subjectstats.stdrangejointangles=[];

%% summary stats for each animal
for i=1:length(subject)
    subject(i).nruns=length(subject(i).bodylength);
    subject(i).meanbodylength=mean(subject(i).bodylength);
    subject(i).meancarapacelength=mean(subject(i).carapacelength);
    subject(i).stdevbodylength=std(subject(i).bodylength);
    subject(i).meanspeedcm=mean(subject(i).speedcm);
    subject(i).stdevspeedcm=std(subject(i).speedcm);
    subject(i).rangespeedcm=[min(subject(i).speedcm) max(subject(i).speedcm)];
    subject(i).meanspeedbodylength=mean(subject(i).speedbodylength);
    subject(i).stdevspeedbodylength=std(subject(i).speedbodylength);
    subject(i).rangespeedbodylength=[min(subject(i).speedbodylength) max(subject(i).speedbodylength)];

    subject(i).meanstridefrequency=mean(subject(i).stridefrequency);
    subject(i).stdevstridefrequency=std(subject(i).stridefrequency);
    subject(i).meanstrideperiod=mean(subject(i).strideperiod);
    subject(i).stdevstrideperiod=std(subject(i).strideperiod);
    subject(i).meanstridelength=mean(subject(i).stridelength);
    subject(i).stdevstridelength=std(subject(i).stridelength);
    subject(i).meanmaxjointangles=mean(subject(i).maxjointangles,2);
    subject(i).stdmaxjointangles=std(subject(i).maxjointangles,0,2);
    subject(i).meanminjointangles=mean(subject(i).minjointangles,2);
    subject(i).stdminjointangles=std(subject(i).minjointangles,0,2);
    subject(i).meanrangejointangles=mean(subject(i).rangejointangles,2);
    subject(i).stdrangejointangles=std(subject(i).rangejointangles,0,2);


    if ~isempty(subject(i).bodylength)
        subjectstats.weight=[subjectstats.weight; subject(i).weight];
        subjectstats.bodylength=[subjectstats.bodylength; subject(i).meanbodylength];
        subjectstats.carapacelength=[subjectstats.carapacelength; subject(i).meancarapacelength];
        subjectstats.averagespeedcm=[subjectstats.averagespeedcm; subject(i).meanspeedcm];
        subjectstats.maxspeedcm=[subjectstats.maxspeedcm; max(subject(i).speedcm)];
        subjectstats.minspeedcm=[subjectstats.minspeedcm; min(subject(i).speedcm)];
        subjectstats.averagespeedbodylength=[subjectstats.averagespeedbodylength; subject(i).meanspeedbodylength];
        subjectstats.maxspeedbodylength=[subjectstats.maxspeedbodylength; max(subject(i).speedbodylength)];
        subjectstats.minspeedbodylength=[subjectstats.minspeedbodylength; min(subject(i).speedbodylength)];
        subjectstats.strideperiod=[subjectstats.strideperiod; mean(subject(i).strideperiod)];
        subjectstats.averagestridefrequency=[subjectstats.averagestridefrequency; mean(subject(i).stridefrequency)];
        subjectstats.averagestridelength=[subjectstats.averagestridelength; mean(subject(i).stridelength)];
        subjectstats.maxstridefrequency=[subjectstats.maxstridefrequency; max(subject(i).stridefrequency)];
        subjectstats.maxstridelength=[subjectstats.maxstridelength; max(subject(i).stridelength)];
        subjectstats.minstridefrequency=[subjectstats.minstridefrequency; min(subject(i).stridefrequency)];
        subjectstats.minstridelength=[subjectstats.minstridelength; min(subject(i).stridelength)];
        subjectstats.maxjointangles=[subjectstats.maxjointangles max(subject(i).maxjointangles')'];
        subjectstats.minjointangles=[subjectstats.minjointangles min(subject(i).minjointangles')'];
        subjectstats.rangejointangles=[subjectstats.rangejointangles max(subject(i).maxjointangles')'-min(subject(i).minjointangles')'];
    end
                subjectstats.meanmaxjointangles=round(mean(subjectstats.maxjointangles')');
        subjectstats.meanminjointangles=round(mean(subjectstats.minjointangles')');
        subjectstats.meanrangejointangles=round(mean(subjectstats.rangejointangles')');
            subjectstats.stdmaxjointangles=round(std(subjectstats.maxjointangles')');
        subjectstats.stdminjointangles=round(std(subjectstats.minjointangles')');
        subjectstats.stdrangejointangles=round(std(subjectstats.rangejointangles')');
end

%     W=[];
%     W=[W;data.weight];
%     M=[min(W)/1000 max(W)/1000];
%     s=(5.5*M.^0.24)*(100000/(60*60));
%     sl=(0.35*M.^0.38)*100;
%     sf=(269*M.^-.14)/60;

%     % speed stride parameters vs mass length parameters plots
%     for i=1:length(subject)
%             figure(2)
%             hold on
%             errorbar(subject(i).weight,subject(i).meanspeedcm,subject(i).stdevspeedcm,'ko','MarkerFaceColor','k')
%             xlabel('Mass (g)')
%             ylabel('Speed (cm/s)')
%             hold off
%             figure(3)
%             hold on
%             h1=Zorgiebel_ploterr(subject(i).meanbodylength, subject(i).meanspeedcm, subject(i).stdevbodylength, subject(i).stdevspeedcm, 'ko','abshhx',.8,'abshhy',.04);
%             hold off
%             set(h1,'MarkerFaceColor','k')
%             xlabel('Body length (cm)')
%             ylabel('Speed (cm/s)')
%             figure(4)
%             hold on
%             errorbar(subject(i).weight,subject(i).meanspeedbodylength,subject(i).stdevspeedbodylength,'ko','MarkerFaceColor','k')
%             xlabel('Mass (g)')
%             ylabel('Speed (body lengths/s)')
%             hold off
%             figure(5)
%             hold on
%             h2=Zorgiebel_ploterr(subject(i).meanbodylength, subject(i).meanspeedbodylength, subject(i).stdevbodylength, subject(i).stdevspeedbodylength,'ko','abshhx',.5,'abshhy',.02);
%             hold off
%             set(h2,'MarkerFaceColor','k')
%             xlabel('Body length (cm)')
%             ylabel('Speed (body lengths /s)')
%             figure(6)
%             hold on
%             errorbar(subject(i).weight,subject(i).meanstridefrequency,subject(i).stdevstridefrequency,'ko','MarkerFaceColor','k')
%             xlabel('Mass (g)')
%             ylabel('Stride frequency (Hz)')
%             figure(7)
%             hold on
%             h3=Zorgiebel_ploterr(subject(i).meanbodylength, subject(i).meanstridefrequency, subject(i).stdevbodylength, subject(i).stdevstridefrequency,'ko','abshhx',.5,'abshhy',.02);
%             xlabel('Body length (cm)')
%             ylabel('Stride frequency (Hz)')
%             set(h3,'MarkerFaceColor','k')
%             hold off
%             figure(8)
%             hold on
%             errorbar(subject(i).weight,subject(i).meanstridelength,subject(i).stdevstridelength,'ko','MarkerFaceColor','k')
%             xlabel('Mass (g)')
%             ylabel('Stride length (cm)')
%             hold off
%             figure(9)
%             hold on
%             h4=Zorgiebel_ploterr(subject(i).meanbodylength, subject(i).meanstridelength, subject(i).stdevbodylength, subject(i).stdevstridelength,'ko','abshhx',.5,'abshhy',.02);
%             xlabel('Body length (cm)')
%             ylabel('Stride length (cm)')
%             set(h4,'MarkerFaceColor','k')
%             hold off
%     end

%     figure(2)
%     hold on
%     plot([min(W) max(W)],s)
%     hold off
%     figure(6)
%     hold on
%     plot([min(W) max(W)],sf)
%     hold off
%     figure(8)
%     hold on
%     plot([min(W) max(W)],sl)
%     hold off



trialstats.averagespeedcm=trial.averagespeedcm;
trialstats.averagespeedbodylength=trial.averagespeedbodylength;
trialstats.strideperiod=trial.averagestrideperiod;
trialstats.stridefrequency=trial.averagestridefrequency;
trialstats.stridelength=trial.averagestridelength;
trialstats.protraction=trial.averageprotraction;
trialstats.retraction=trial.averageretraction;
trialstats.dutyfactor=trial.averagedutyfactor;
trialstats.stabilitymargin=trial.averagestabilitymargin;

trial.l1_r1phase=[];
for i=[1:length(trial.averageL1phase)]
    trial.l1_r1phase=[trial.l1_r1phase; trial.averageL1phase{i}(2)];
end
% figure;
% kiri_markerScatter(trial.averagespeedcm,trial.l1_r1phase,[1:length(trial.averageL1phase)],'trial.averagespeedcm','trial.l1_r1phase')
% hold on
% plot([0 45], [0.75 0.75], 'k:')
% plot([0 45], [0.25 0.25], 'k:')
% xlabel('Speed (cm/s)')
% ylabel('\phi')
% hold off
% ylim([0,1])
trial.l2_r2phase=[];
for i=[1:length(trial.averageL2phase)]
    trial.l2_r2phase=[trial.l2_r2phase; trial.averageL2phase{i}(4)];
end
% figure;
% kiri_markerScatter(trial.averagespeedcm,trial.l2_r2phase,[1:length(trial.averageL1phase)],'trial.averagespeedcm','trial.l2_r2phase')
% hold on
% plot([0 45], [0.75 0.75], 'k:')
% plot([0 45], [0.25 0.25], 'k:')
% xlabel('Speed (cm/s)')
% ylabel('\phi')
% hold off
% ylim([0,1])
trial.l3_r3phase=[];
for i=[1:length(trial.averageL3phase)]
    trial.l3_r3phase=[trial.l3_r3phase; trial.averageL3phase{i}(6)];
end
% figure;
% kiri_markerScatter(trial.averagespeedcm,trial.l3_r3phase,[1:length(trial.averageL3phase)],'trial.averagespeedcm','trial.l3_r3phase')
% hold on
% plot([0 45], [0.75 0.75], 'k:')
% plot([0 45], [0.25 0.25], 'k:')
% xlabel('Speed (cm/s)')
% ylabel('\phi')
% hold off
% ylim([0,1])
trial.l4_r4phase=[];
for i=[1:length(trial.averageL4phase)]
    trial.l4_r4phase=[trial.l4_r4phase; trial.averageL4phase{i}(8)];
end
% figure;
% kiri_markerScatter(trial.averagespeedcm,trial.l4_r4phase,[1:length(trial.averageL4phase)],'trial.averagespeedcm','trial.l4_r4phase')
% hold on
% plot([0 45], [0.75 0.75], 'k:')
% plot([0 45], [0.25 0.25], 'k:')
% xlabel('Speed (cm/s)')
% ylabel('\phi')
% hold off
% ylim([0,1])
trial.l1_l2phase=[];
for i=[1:length(trial.averageL1phase)]
    trial.l1_l2phase=[trial.l1_l2phase; trial.averageL1phase{i}(3)];
end
trial.r1_r2phase=[];
for i=[1:length(trial.averageR1phase)]
    trial.r1_r2phase=[trial.r1_r2phase; trial.averageR1phase{i}(4)];
end
% figure;
% kiri_markerScatter(trial.averagespeedcm,trial.r1_r2phase,[1:length(trial.averageR1phase)],'trial.averagespeedcm','trial.r1_r2phase')
% hold on
% plot([0 45], [0.75 0.75], 'k:')
% plot([0 45], [0.25 0.25], 'k:')
% xlabel('Speed (cm/s)')
% ylabel('\phi')
% hold off
% ylim([0,1])
trial.l2_l3phase=[];
for i=[1:length(trial.averageL2phase)]
    trial.l2_l3phase=[trial.l2_l3phase; trial.averageL2phase{i}(5)];
end
trial.r2_r3phase=[];
for i=[1:length(trial.averageR2phase)]
    trial.r2_r3phase=[trial.r2_r3phase; trial.averageR2phase{i}(6)];
end
% figure;
% kiri_markerScatter(trial.averagespeedcm,trial.averagedutyfactor1,[1:length(trial.averageR2phase)],'trial.averagespeedcm','trial.r2_r3phase')
% hold on
% plot([0 45], [0.75 0.75], 'k:')
% plot([0 45], [0.25 0.25], 'k:')
% xlabel('Speed (cm/s)')
% ylabel('\phi')
% hold off
% ylim([0,1])
trial.l3_l4phase=[];
for i=[1:length(trial.averageL3phase)]
    trial.l3_l4phase=[trial.l3_l4phase; trial.averageL3phase{i}(7)];
end
trial.r3_r4phase=[];
for i=[1:length(trial.averageR3phase)]
    trial.r3_r4phase=[trial.r3_r4phase; trial.averageR3phase{i}(8)];
end
% figure;
% kiri_markerScatter(trial.averagespeedcm,trial.averagedutyfactor2,[1:length(trial.averageR3phase)],'trial.averagespeedcm','trial.r3_r4phase')
% hold on
% plot([0 45], [0.75 0.75], 'k:')
% plot([0 45], [0.25 0.25], 'k:')
% xlabel('Speed (cm/s)')
% ylabel('\phi')
% hold off
% ylim([0,1])
trial.l1_l3phase=[];
for i=[1:length(trial.averageL1phase)]
    trial.l1_l3phase=[trial.l1_l3phase; trial.averageL1phase{i}(5)];
end
trial.r1_r3phase=[];
for i=[1:length(trial.averageR1phase)]
    trial.r1_r3phase=[trial.r1_r3phase; trial.averageR1phase{i}(6)];
end
% figure;
% kiri_markerScatter(trial.averagespeedcm,trial.averagedutyfactor3,[1:length(trial.averageR1phase)],'trial.averagespeedcm','trial.r1_r3phase')
% hold on
% plot([0 45], [0.75 0.75], 'k:')
% plot([0 45], [0.25 0.25], 'k:')
% xlabel('Speed (cm/s)')
% ylabel('\phi')
% hold off
% ylim([0,1])
trial.l1_l4phase=[];
for i=[1:length(trial.averageL1phase)]
    trial.l1_l4phase=[trial.l1_l4phase; trial.averageL1phase{i}(7)];
end
trial.r1_r4phase=[];
for i=[1:length(trial.averageR1phase)]
    trial.r1_r4phase=[trial.r1_r4phase; trial.averageR1phase{i}(8)];
end
% figure;
% kiri_markerScatter(trial.averagespeedcm,trial.averagedutyfactor4,[1:length(trial.averageR1phase)],'trial.averagespeedcm','trial.r1_r4phase')
% hold on
% plot([0 45], [0.75 0.75], 'k:')
% plot([0 45], [0.25 0.25], 'k:')
% xlabel('Speed (cm/s)')
% ylabel('\phi')
% hold off
% ylim([0,1])
trial.l2_l4phase=[];
for i=[1:length(trial.averageL2phase)]
    trial.l2_l4phase=[trial.l2_l4phase; trial.averageL2phase{i}(7)];
end
trial.r2_r4phase=[];
for i=[1:length(trial.averageR2phase)]
    trial.r2_r4phase=[trial.r2_r4phase; trial.averageR2phase{i}(8)];
end
% figure
% kiri_markerScatter(trial.averagespeedcm,trial.r2_r4phase,[1:length(trial.averageR2phase)],'trial.averagespeedcm','trial.r2_r4phase')
% hold on
% plot([0 45], [0.75 0.75], 'k:')
% plot([0 45], [0.25 0.25], 'k:')
% xlabel('Speed (cm/s)')
% ylabel('\phi')
% hold off
% ylim([0,1])
% 
% % stride parameters vs. speed plots
%         figure
%         kiri_markerScatter(trial.averagespeedcm,trial.averagestridefrequency,trial.individual,'trial.averagespeedcm','trial.averagestridefrequency')
%         xlabel('Speed (cm/s)')
%         ylabel('Stride frequency (Hz)')
%         W=[];
%         W=[W;data.weight];
%         M=[mean(W)/1000];
%         s=(5.5*M.^0.24)*(100000/(60*60));
%         sf=(269*M.^-.14)/60;
%         x=s*ones(1,50);
%         y=linspace(0,sf,50);
%         hold on
%         plot(x,y,'k:')
%         y=sf*ones(1,50);
%         x=linspace(0,s,50);
%         plot(x,y,'k:')
%         hold off
%         figure
%         hold on
%         kiri_markerScatter(trial.averagespeedbodylength,trial.averagestridefrequency, trial.individual,'trial.averagespeedbodylength','trial.averagestridefrequency')
%         xlabel('Speed (body lengths/s)')
%         ylabel('Stride frequency (Hz)')
%         hold off
%         figure
%         kiri_markerScatter(trial.averagespeedcm,trial.averagestridelength, trial.individual,'trial.averagespeedcm','trial.averagestridelength')
%         xlabel('Speed (cm/s)')
%         ylabel('Stride length (cm)')
%         sl=(0.35*M.^0.38)*100;
%         x=s*ones(1,50);
%         y=linspace(0,sl,50);
%         hold on
%         plot(x,y,'k:')
%         y=sl*ones(1,50);
%         x=linspace(0,s,50);
%         plot(x,y,'k:')
%         hold off
%         figure
%         hold on
%         kiri_markerScatter(trial.averagespeedbodylength,trial.averagestridelength, trial.individual,'trial.averagespeedbodylength','trial.averagestridelength')
%         xlabel('Speed (body lengths /s)')
%         ylabel('Stride length (cm)')
%         hold off
%         figure
%         hold on
%         kiri_markerScatter(trial.averagestridefrequency,trial.averageprotraction, trial.individual,'trial.averagestridefrequency','trial.averageprotraction')
%         xlabel('Stride frequency (Hz))')
%         ylabel('Protraction (s)')
%         hold off
%         figure
%         hold on
%         kiri_markerScatter(trial.averagestridefrequency,trial.averageretraction, trial.individual,'trial.averagestridefrequency','trial.averageretraction')
%         xlabel('Stride frequency (Hz)')
%         ylabel('Retraction (s)')
%         hold off
%         figure
%         hold on
%         kiri_markerScatter(trial.averagespeedcm,trial.averageprotraction, trial.individual,'trial.averagespeedcm','trial.averageprotraction')
%         xlabel('Speed (cm/s)')
%         ylabel('Protraction (s)')
%         hold off
%         figure
%         hold on
%         kiri_markerScatter(trial.averagespeedcm,trial.averageretraction, trial.individual,'trial.averagespeedcm','trial.averageretraction')
%         xlabel('Speed (cm/s)')
%         ylabel('Retraction (s)')
%         hold off
%         figure
%         hold on
%         kiri_markerScatter(trial.averagespeedbodylength,trial.averageprotraction, trial.individual,'trial.averagespeedbodylength','trial.averageprotraction')
%         xlabel('Speed (body lengths/s)')
%         ylabel('Protraction (s)')
%         hold off
%         figure
%         hold on
%         kiri_markerScatter(trial.averagespeedbodylength,trial.averageretraction, trial.individual,'trial.averagespeedbodylength','trial.averageretraction')
%         xlabel('Speed (body lengths/s)')
%         ylabel('Retraction (s)')
%         hold off
%         figure
%         kiri_markerScatter(trial.averagespeedcm,trial.averagedutyfactor, trial.individual,'trial.averagespeedcm','trial.averagedutyfactor')
%         hold on
%         plot([0 45], [0.375 0.375],'k:')
%         xlabel('Speed (cm/s)')
%         ylabel('Duty factor')
%         hold off
%         figure
%         kiri_markerScatter(trial.averagespeedbodylength,trial.averagedutyfactor, trial.individual,'trial.averagespeedbodylength','trial.averagedutyfactor')
%         hold on
%         plot([0 25], [0.375 0.375],'k:')
%         xlabel('Speed (body lengths /s)')
%         ylabel('Duty factor')
%         hold off
%         figure
%         hold on
%         kiri_markerScatter(trial.averageleg1strideperiod,trial.averageleg1stancevel,trial.individual,'trial.averageleg1strideperiod','trial.averageleg1stancevel')
%         xlabel('Stride period (s)')
%         ylabel('Angular velocity (degrees/s)')
%         hold off
%         figure
%         hold on
%         kiri_markerScatter(trial.averageleg2strideperiod,trial.averageleg2stancevel,trial.individual,'trial.averageleg2strideperiod','trial.averageleg2stancevel')
%         xlabel('Stride period (s)')
%         ylabel('Angular velocity (degrees/s)')
%         hold off
%         figure
%         hold on
%         kiri_markerScatter(trial.averageleg3strideperiod,trial.averageleg3stancevel,trial.individual,'trial.averageleg3strideperiod','trial.averageleg3stancevel')
%         xlabel('Stride period (s)')
%         ylabel('Angular velocity (degrees/s)')
%         hold off
%         figure
%         hold on
%         kiri_markerScatter(trial.averageleg4strideperiod,trial.averageleg4stancevel,trial.individual,'trial.averageleg4strideperiod','trial.averageleg4stancevel')
%         xlabel('Stride period (s)')
%         ylabel('Angular velocity (degrees/s)')
%         hold off
%         
%         figure
%         hold on
%         kiri_markerScatter(trial.averageleg1strideperiod,trial.averageleg1swingvel,trial.individual,'trial.averageleg1strideperiod','trial.averageleg1swingvel')
%         xlabel('Stride period (s)')
%         ylabel('Angular velocity (degrees/s)')
%         hold off
%         figure
%         hold on
%         kiri_markerScatter(trial.averageleg2strideperiod,trial.averageleg2swingvel,trial.individual,'trial.averageleg2strideperiod','trial.averageleg2swingvel')
%         xlabel('Stride period (s)')
%         ylabel('Angular velocity (degrees/s)')
%         hold off
%         figure
%         hold on
%         kiri_markerScatter(trial.averageleg3strideperiod,trial.averageleg3swingvel,trial.individual,'trial.averageleg3strideperiod','trial.averageleg3swingvel')
%         xlabel('Stride period (s)')
%         ylabel('Angular velocity (degrees/s)')
%         hold off
%         figure
%         hold on
%         kiri_markerScatter(trial.averageleg4strideperiod,trial.averageleg4swingvel,trial.individual,'trial.averageleg4strideperiod','trial.averageleg4swingvel')
%         xlabel('Stride period (s)')
%         ylabel('Angular velocity (degrees/s)')
%         hold off
%         figure
%         hold on
%         kiri_markerScatter(trial.averageleg1strideperiod,trial.averageleg1angexcur,trial.individual,'trial.averageleg1strideperiod','trial.averageleg1angexcur')
%         xlabel('Stride period (s)')
%         ylabel('Leg excursion (degrees)')
%         hold off
%         figure
%         hold on
%         kiri_markerScatter(trial.averageleg2strideperiod,trial.averageleg2angexcur,trial.individual,'trial.averageleg2strideperiod','trial.averageleg2angexcur')
%         xlabel('Stride period (s)')
%         ylabel('Leg excursion (degrees)')
%         hold off
%         figure
%         hold on
%         kiri_markerScatter(trial.averageleg3strideperiod,trial.averageleg3angexcur,trial.individual,'trial.averageleg3strideperiod','trial.averageleg3angexcur')
%         xlabel('Stride period (s)')
%         ylabel('Leg excursion (degrees)')
%         hold off
%         figure
%         hold on
%         kiri_markerScatter(trial.averageleg4strideperiod,trial.averageleg4angexcur,trial.individual,'trial.averageleg4strideperiod','trial.averageleg4angexcur')
%         xlabel('Stride period (s)')
%         ylabel('Leg excursion (degrees)')
%         hold off

%set up comparrisons for regression on log data
subjectvariables={'subjectstats.weight' 'subjectstats.averagespeedcm'; 'subjectstats.weight' 'subjectstats.minspeedcm'; 'subjectstats.weight' 'subjectstats.maxspeedcm';...
    'subjectstats.weight' 'subjectstats.averagespeedbodylength'; 'subjectstats.weight' 'subjectstats.minspeedbodylength'; 'subjectstats.weight' 'subjectstats.maxspeedbodylength';...
    'subjectstats.bodylength' 'subjectstats.averagespeedcm'; 'subjectstats.bodylength' 'subjectstats.minspeedcm'; 'subjectstats.bodylength' 'subjectstats.maxspeedcm';...
    'subjectstats.bodylength' 'subjectstats.averagespeedbodylength'; 'subjectstats.bodylength' 'subjectstats.minspeedbodylength'; 'subjectstats.bodylength' 'subjectstats.maxspeedbodylength';...
    'subjectstats.weight' 'subjectstats.averagestridefrequency';'subjectstats.bodylength' 'subjectstats.averagestridefrequency'; 'subjectstats.weight' 'subjectstats.averagestridelength';'subjectstats.bodylength' 'subjectstats.averagestridelength';...
    'subjectstats.weight' 'subjectstats.maxstridefrequency';'subjectstats.bodylength' 'subjectstats.maxstridefrequency'; 'subjectstats.weight' 'subjectstats.maxstridelength';'subjectstats.bodylength' 'subjectstats.maxstridelength';...
    'subjectstats.weight' 'subjectstats.minstridefrequency';'subjectstats.bodylength' 'subjectstats.minstridefrequency'; 'subjectstats.weight' 'subjectstats.minstridelength';'subjectstats.bodylength' 'subjectstats.minstridelength';};
trialvariables={'trialstats.averagespeedcm' 'trialstats.stridefrequency';'trialstats.averagespeedbodylength' 'trialstats.stridefrequency';...
    'trialstats.averagespeedcm' 'trialstats.stridelength';'trialstats.averagespeedbodylength' 'trialstats.stridelength';...
    'trialstats.averagespeedcm' 'trialstats.protraction';'trialstats.averagespeedbodylength' 'trialstats.protraction'; 'trial.averagestridefrequency' 'trialstats.protraction';...
    'trialstats.averagespeedcm' 'trialstats.retraction';'trialstats.averagespeedbodylength' 'trialstats.retraction'; 'trial.averagestridefrequency' 'trialstats.retraction';...
    'trialstats.averagespeedcm' 'trialstats.dutyfactor';'trialstats.averagespeedbodylength' 'trialstats.dutyfactor';...
    'trialstats.averagespeedcm' 'trialstats.stabilitymargin';'trialstats.averagespeedbodylength' 'trialstats.stabilitymargin';
    'trial.averageleg1strideperiod' 'trial.averageleg1angexcur';'trial.averageleg2strideperiod' 'trial.averageleg2angexcur';...
    'trial.averageleg3strideperiod' 'trial.averageleg3angexcur';'trial.averageleg4strideperiod' 'trial.averageleg4angexcur';...
    'trial.averageleg1stridefrequency' 'trial.averageleg1angexcur';'trial.averageleg2stridefrequency' 'trial.averageleg2angexcur';...
    'trial.averageleg3stridefrequency' 'trial.averageleg3angexcur';'trial.averageleg4stridefrequency' 'trial.averageleg4angexcur';...
    'trial.averageleg1strideperiod' 'trial.averageleg1dist';'trial.averageleg2strideperiod' 'trial.averageleg2dist';...
    'trial.averageleg3strideperiod' 'trial.averageleg3dist';'trial.averageleg4strideperiod' 'trial.averageleg4dist';...
    'trial.averageleg1stridefrequency' 'trial.averageleg1dist';'trial.averageleg2stridefrequency' 'trial.averageleg2dist';...
    'trial.averageleg3stridefrequency' 'trial.averageleg3dist';'trial.averageleg4stridefrequency' 'trial.averageleg4dist';...
    'trial.averageleg1strideperiod' 'trial.averageleg1swingvel';'trial.averageleg2strideperiod' 'trial.averageleg2swingvel';...
    'trial.averageleg3strideperiod' 'trial.averageleg3swingvel';'trial.averageleg4strideperiod' 'trial.averageleg4swingvel';...
    'trial.averageleg1stridefrequency' 'trial.averageleg1swingvel';'trial.averageleg2stridefrequency' 'trial.averageleg2swingvel';...
    'trial.averageleg3stridefrequency' 'trial.averageleg3swingvel';'trial.averageleg4stridefrequency' 'trial.averageleg4swingvel';...
     'trial.averageleg1strideperiod' 'trial.averageleg1stancevel';'trial.averageleg2strideperiod' 'trial.averageleg2stancevel';...
    'trial.averageleg3strideperiod' 'trial.averageleg3stancevel';'trial.averageleg4strideperiod' 'trial.averageleg4stancevel';...
    'trial.averageleg1stridefrequency' 'trial.averageleg1stancevel';'trial.averageleg2stridefrequency' 'trial.averageleg2stancevel';...
    'trial.averageleg3stridefrequency' 'trial.averageleg3stancevel';'trial.averageleg4stridefrequency' 'trial.averageleg4stancevel';...
    };

circvariables={'trialstats.averagespeedcm' 'trial.l1_r1phase';'trialstats.averagespeedbodylength' 'trial.l1_r1phase';...
    'trialstats.averagespeedcm' 'trial.l2_r2phase';'trialstats.averagespeedbodylength' 'trial.l2_r2phase';...
    'trialstats.averagespeedcm' 'trial.l3_r3phase';'trialstats.averagespeedbodylength' 'trial.l3_r3phase';...
    'trialstats.averagespeedcm' 'trial.l4_r4phase';'trialstats.averagespeedbodylength' 'trial.l4_r4phase';...
    'trialstats.averagespeedcm' 'trial.l1_l2phase';'trialstats.averagespeedbodylength' 'trial.l1_l2phase';...
    'trialstats.averagespeedcm' 'trial.l2_l3phase';'trialstats.averagespeedbodylength' 'trial.l2_l3phase';...
    'trialstats.averagespeedcm' 'trial.l3_l4phase';'trialstats.averagespeedbodylength' 'trial.l3_l4phase';...
    'trialstats.averagespeedcm' 'trial.l1_l3phase';'trialstats.averagespeedbodylength' 'trial.l1_l3phase';...
    'trialstats.averagespeedcm' 'trial.l1_l4phase';'trialstats.averagespeedbodylength' 'trial.l1_l4phase';...
    'trialstats.averagespeedcm' 'trial.l2_l4phase';'trialstats.averagespeedbodylength' 'trial.l2_l4phase';...
    'trialstats.averagespeedcm' 'trial.r1_r2phase';'trialstats.averagespeedbodylength' 'trial.r1_r2phase';...
    'trialstats.averagespeedcm' 'trial.r2_r3phase';'trialstats.averagespeedbodylength' 'trial.r2_r3phase';...
    'trialstats.averagespeedcm' 'trial.r3_r4phase';'trialstats.averagespeedbodylength' 'trial.r3_r4phase';...
    'trialstats.averagespeedcm' 'trial.r1_r3phase';'trialstats.averagespeedbodylength' 'trial.r1_r3phase';...
    'trialstats.averagespeedcm' 'trial.r1_r4phase';'trialstats.averagespeedbodylength' 'trial.r1_r4phase';...
    'trialstats.averagespeedcm' 'trial.r2_r4phase';'trialstats.averagespeedbodylength' 'trial.r2_r4phase';...
    'trial.averagestridefrequency' 'trial.l1_r1phase';'trial.averagestridefrequency' 'trial.l2_r2phase';...
    'trial.averagestridefrequency' 'trial.l3_r3phase';'trial.averagestridefrequency' 'trial.l4_r4phase';...
    'trial.averagestridefrequency' 'trial.l1_l2phase';'trial.averagestridefrequency' 'trial.r1_r2phase';...
    'trial.averagestridefrequency' 'trial.l2_l3phase';'trial.averagestridefrequency' 'trial.r2_r3phase';...
    'trial.averagestridefrequency' 'trial.l3_l4phase';'trial.averagestridefrequency' 'trial.r3_r4phase';...
    'trial.averagestridefrequency' 'trial.l1_l3phase';'trial.averagestridefrequency' 'trial.r1_r3phase';...
    'trial.averagestridefrequency' 'trial.l1_l4phase';'trial.averagestridefrequency' 'trial.r1_r4phase';};


for i=1:length(subjectvariables)
    subjectregression.subjectvariables(i,:)=subjectvariables(i,:);
    subjectregression.x(:,i)=eval(subjectvariables{i,1});
    subjectregression.y(:,i)=eval(subjectvariables{i,2});
    Y=subjectregression.y(:,i);
    X=[ones(size(subjectregression.x(:,i)),1) subjectregression.x(:,i)];
    [subjectregression.log.loga(i), subjectregression.log.a(i),subjectregression.log.b(i), subjectregression.log.bCI(:,i), subjectregression.log.SE(i), subjectregression.log.rsquared(i),subjectregression.log.p(i)] = kiri_correctedLogRegression(subjectregression.x(:,i),subjectregression.y(:,i));
    [b, bint,r,rint,stats] = regress(Y,X);
    subjectregression.linear.a(i)=b(1);
    subjectregression.linear.b(i)=b(2);
    subjectregression.linear.bCI(:,i)=bint(2,:);
    subjectregression.linear.rsquared(i)=stats(1);
    subjectregression.linear.p(i)=stats(3);

%     X=linspace(min(subjectregression.x(:,i)),max(subjectregression.x(:,i)));
%     Ylin=subjectregression.linear.a(i)+subjectregression.linear.b(i)*X;
%     Ylog=subjectregression.log.a(i)*(X.^subjectregression.log.b(i));

    % figure;
    % scatter(subjectregression.x(:,i),subjectregression.y(:,i),'ko','MarkerFaceColor','k')
    % xlabel(subjectvariables{i,1}); ylabel(subjectvariables{i,2});
    % axis tight
    % hold on
    % plot(X,Ylin,'r--')
    % plot(X,Ylog,'b--')
    % hold off

    % Y=(subjectregression.y(:,i));
    % X=[ones(size(subjectregression.x(:,i)),1) (subjectregression.x(:,i))];
    % [dataout{i} lowerLimit{i} upperLimit{i}] = lowess([(subjectregression.x(:,i)) (subjectregression.y(:,i))],1,0);
    % [B(i),BINT(i),R(i),RINT(i),STATS(i)] = REGRESS(Y,X);
    % s = regstats(Y,X,'linear','all');
    % S{i}=s;
    % figure;
    % scatter(s.yhat,s.r)
    % xlabel('Fitted Values'); ylabel('Residuals');
    % figure;
    % probplot('normal',abs(s.r))
    % figure;
    %
    % Y=log(subjectregression.y(:,i));
    % X=[ones(size(subjectregression.x(:,i)),1) log(subjectregression.x(:,i))];
    % % [B(i),BINT(i),R(i),RINT(i),STATS(i)] = REGRESS(Y,X);
    % s = regstats(Y,X,'linear','all');
    % S{i}=s;
    % figure;
    % scatter(s.yhat,s.r)
    % xlabel('Fitted Values'); ylabel('Residuals');
    % figure;
    % probplot('normal',abs(s.r))
    % figure;
end

for i=1:length(trialvariables)
    trialregression.trialvariables(i,:)=trialvariables(i,:);
    trialregression.x(:,i)=eval(trialvariables{i,1});
    trialregression.y(:,i)=eval(trialvariables{i,2});
    Y=(trialregression.y(:,i));
    X=[ones(size(trialregression.x(:,i)),1) (trialregression.x(:,i))];
    [trialregression.log.loga(i), trialregression.log.a(i),trialregression.log.b(i), trialregression.log.bCI(:,i), trialregression.log.SEE(i), trialregression.log.rsquared(i),trialregression.log.p(i)] = kiri_correctedLogRegression(trialregression.x(:,i),trialregression.y(:,i));
    [b, bint,r,rint,stats] = regress(Y,X);
    trialregression.linear.a(i)=b(1);
    trialregression.linear.b(i)=b(2);
    trialregression.linear.bCI(:,i)=bint(2,:);
    trialregression.linear.rsquared(i)=stats(1);
    trialregression.linear.p(i)=stats(3);
    % X=linspace(min(trialregression.x(:,i)),max(trialregression.x(:,i)));
    % Ylin=trialregression.linear.a(i)+trialregression.linear.b(i)*X;
    % Ylog=trialregression.log.a(i)*(X.^trialregression.log.b(i));
%  [dataout{i} lowerLimit{i} upperLimit{i}] = lowess([(trialregression.x(:,i)) (trialregression.y(:,i))],1,0);
    % figure;
    % scatter(trialregression.x(:,i),trialregression.y(:,i),'ko','MarkerFaceColor','k')
    % xlabel(trialvariables{i,1}); ylabel(trialvariables{i,2});
    % axis tight
    % hold on
    % plot(X,Ylin,'r--')
    % plot(X,Ylog,'b--')
    % hold off
end

for i=1:length(circvariables)
    circregression.variables(i,:)=circvariables(i,:);
    circregression.x(:,i)=eval(circvariables{i,1});
    circregression.y(:,i)=eval(circvariables{i,2});
    [circregression.rsquared(i),circregression.p(i)] = kiri_anglinRegression(circregression.y(:,i),circregression.x(:,i));
end

subjectregression.log.significantslope=find(subjectregression.log.p<0.05);
subjectregression.linear.significantslope=find(subjectregression.linear.p<0.05);
trialregression.log.significantslope=find(trialregression.log.p<0.05);
trialregression.linear.significantslope=find(trialregression.linear.p<0.05);
circregression.significantslope=find(circregression.p<0.05);