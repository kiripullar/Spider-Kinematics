% kiri_compareAnimals loads data files (which are .txt documents) and
% firstly compares Dolomedes aquaticus to other spider species using a
% scatter plot and regression on log transformed data and then repeats this
% process to compare results to terrestrial species, requires that either
% kiri_statisticalAnalysis has just been run or results have been opened
% from a .mat file to access the results from D. aquaticus. 
%
% Kiri Pullar, masters thesis 2009

%% Compare Dolomedes aquaticus to other spider species
fid=fopen('spiderfamilydata.txt');
c=textscan(fid,'%f%f%f', 'Delimiter', '\t'); %family, carapace length cm, speed cm/s
fclose(fid)

c{1}=[c{1}; 5];
c{2}=[c{2}; mean(subjectstats.carapacelength)];
c{3}=[c{3}; total.meanspeedcm];

for i=1:length(c{1})
    switch c{1}(i)
        case 1
            hold on
            scatter(c{2}(i),c{3}(i),'k','*')
            hold off
        case 2
            hold on
            scatter(c{2}(i),c{3}(i),'k','o')
            hold off
        case 3
            hold on
            scatter(c{2}(i),c{3}(i),'k','^')
            hold off
        case 4
            hold on
            scatter(c{2}(i),c{3}(i),'k','s')
            hold off
        case 5
            hold on
            scatter(c{2}(i),c{3}(i),'k','o','filled')
            hold off
    end
end

[loga, a, b, bCI, SE, rsquared, p] = regresslogdata(c{2},c{3});

%% Compare Dolomedes aquaticus to other terrestrial animals
figure
fid=fopen('mammaldata.txt');
c=textscan(fid,'%f%f%f', 'Delimiter', ' '); %mass kg, body length m, relative speed body length/s
fclose(fid)

mammalmasskg=c{1};
mammalspeedms=c{2}.*c{3};
hold on
loglog(mammalmasskg,mammalspeedms,'ks')
hold off
fid=fopen('reptiledata.txt');
c=textscan(fid,'%f%f%f', 'Delimiter', ' '); %mass g, speed ms
fclose(fid)

lizardmasskg=c{2}/1000;
lizardspeedms=c{1};
hold on
loglog(lizardmasskg,lizardspeedms,'k^')
hold off
fid=fopen('birddata.txt');
c=textscan(fid,'%f%f%f', 'Delimiter', ' '); %mass g, speed ms
fclose(fid)

birdmasskg=c{1}/1000;
birdspeedms=c{2};
hold on
loglog(birdmasskg,birdspeedms,'kd')
hold off
fid=fopen('amphibiandata.txt');
c=textscan(fid,'%f%f%f', 'Delimiter', ' '); %mass g, speed ms
fclose(fid)

amphibianmasskg=c{1}/1000;
amphibianspeedms=c{2};
hold on
loglog(amphibianmasskg,amphibianspeedms,'kp')
hold off
fid=fopen('arthropoddata.txt');
c=textscan(fid,'%f%f%f', 'Delimiter', ' '); %mass g, speed ms
fclose(fid)

arthropodmasskg=c{1}/1000;
arthropodspeedms=c{2};

mass=[];
for i=1:length(subject)
mass=[mass subject(i).weight];
end

hold on
loglog(arthropodmasskg,arthropodspeedms,'ko')
loglog(mean(mass)/1000, max(trial.averagespeedcm)/100,'ko','MarkerFaceColor','k')
hold off

totalmasskg=[mammalmasskg;lizardmasskg;birdmasskg;amphibianmasskg;arthropodmasskg;mean(mass)/1000];
totalspeedcm=[mammalspeedms;lizardspeedms;birdspeedms;amphibianspeedms;arthropodspeedms;max(trial.averagespeedcm)/100];
[loga,a,b,bCI,SE,rsquared,p] = regresslogdata(totalmasskg,totalspeedcm);

