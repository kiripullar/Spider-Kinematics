function [individual,weight,bodylength,carapacelength,speedcm,maxspeedcm,minspeedcm,speedbodylength,maxspeedbodylength,minspeedbodylength,maxjointangles, minjointangles, rangejointangles, legangle, legangularvelocity, legangularacceleration, swingduration,stanceduration,strideperiod, legstep, L1phase,R1phase,L2phase,R2phase,L3phase,R3phase,L4phase, R4phase, legexcur, angexcur,swinglegvel,stancelegvel] = kiri_kinematicAnalysis(inputfile,model,video,filtering,tracking)
% kiri_kinematicAnalysis(inputfile,model,video,filtering,tracking)
% requires inputfile, model, video, filtering, tracking fields from a .mat
% file containing results and outputs information relating to the size and
% kinematic parameters of the individual
%
% Kiri Pullar, masters thesis 2009

%% get individual, weight and length
str=inputfile.filename;
pat='r(\d*)';
t = regexp(str, pat, 'match');
individual=str2num(t{1}(2:length(t{1})));
weight=inputfile.spiderweight;
bodylength=sqrt((model.BodyPoint(1,3)-model.BodyPoint(1,1))^2+(model.BodyPoint(2,3)-model.BodyPoint(2,1))^2)*video.pix2cm;
carapacelength=sqrt((model.BodyPoint(1,2)-model.BodyPoint(1,1))^2+(model.BodyPoint(2,2)-model.BodyPoint(2,1))^2)*video.pix2cm;
%% overall speed calculations
position=filtering.positions*video.pix2cm;
t=(2:length(position)-1)*1/125;
velocityx=[];
for i=2:length(position)-1
    velocityx=[velocityx; (position(2,i+1)-position(2,i-1))/(2/125)];%velocity in x direction cm/s
end
% fig=figure;
% plot(t,velocityx, 'ko-')

velocityy=[];
for i=2:length(position)-1
    velocityy=[velocityy; (position(1,i+1)-position(1,i-1))/(2/125)];%velocity in y direction cm/s
end
% fig=figure;
% plot(t,velocityy, 'ko-')

speedcm=sqrt(velocityx.^2+velocityy.^2); %overall speed cm/s
[first,last]=kiri_findLocomotion(5,speedcm); %stop analysis if spider has slowed down too much i.e. is stopping
first=first+1;
last=last-1;
speedcm=sqrt(velocityx(first:last).^2+velocityy(first:last).^2);
maxspeedcm=max(speedcm);
minspeedcm=min(speedcm);
meanspeedcm=mean(speedcm);
% fig=figure;
% plot(t,speedcm, 'ko-')

speedbodylength=speedcm/bodylength; %overall speed bodylength/s
maxspeedbodylength=max(speedbodylength);
minspeedbodylength=min(speedbodylength);
meanspeedbodylength=mean(speedbodylength);
% fig=figure;
% plot(t,speedbodylength, 'ko-')

%% max and min angle of each joint
neutral=[35; 14; 18; 27];
filtering.angles(2,:)=filtering.angles(2,:)+neutral(1)*ones(1,length(filtering.angles(2,:)));
filtering.angles(9,:)=filtering.angles(9,:)+neutral(2)*ones(1,length(filtering.angles(9,:)));
filtering.angles(16,:)=filtering.angles(16,:)+neutral(3)*ones(1,length(filtering.angles(16,:)));
filtering.angles(23,:)=filtering.angles(23,:)+neutral(4)*ones(1,length(filtering.angles(23,:)));
filtering.angles(30,:)=filtering.angles(30,:)-neutral(1)*ones(1,length(filtering.angles(30,:)));
filtering.angles(37,:)=filtering.angles(37,:)-neutral(2)*ones(1,length(filtering.angles(37,:)));
filtering.angles(44,:)=filtering.angles(44,:)-neutral(3)*ones(1,length(filtering.angles(44,:)));
filtering.angles(51,:)=filtering.angles(51,:)-neutral(4)*ones(1,length(filtering.angles(51,:)));
filtering.angles=filtering.angles(2:57,first:last);

maxjoint=max(filtering.angles');
minjoint=min(filtering.angles');
rangejoint=maxjoint-minjoint;

maxjointangles=[];
for j=1:14:49
    for i=j:1:j+6
        maxjointangles=[maxjointangles; maxjoint(i) maxjoint(i+7)];
    end
end

minjointangles=[];
for j=1:14:49
    for i=j:1:j+6
        minjointangles=[minjointangles; minjoint(i) minjoint(i+7)];
    end
end

rangejointangles=[];
for j=1:14:49
    for i=j:1:j+6
        rangejointangles=[rangejointangles; rangejoint(i) rangejoint(i+7)];
    end
end

%% calculate total leg angles
filtering.bodypoints=filtering.bodypoints(:,:,first-1:last+1);
legangle=zeros(8,size(filtering.bodypoints,3));
for j=1:size(filtering.bodypoints,3)
    count=1;
    for i=[4:16:60]
        P1=filtering.bodypoints(:,i,j);
        P3=filtering.bodypoints(:,i+7,j);
        P2=filtering.bodypoints(:,i+8,j);
        legangle(count,j) = -kiri_rad2deg(atan2(det([P1-P2,P3-P2]),dot(P1-P2,P3-P2))); %-so that forwards is +ve
        count=count+2;
    end
    count=2;
    for i=[12:16:67]
        P1=filtering.bodypoints(:,i,j);
        P3=filtering.bodypoints(:,i+7,j);
        P2=filtering.bodypoints(:,i-8,j);
        legangle(count,j) = kiri_rad2deg(atan2(det([P1-P2,P3-P2]),dot(P1-P2,P3-P2)));
        count=count+2;
    end
end

try
    AEP=[];
    PEP=[];
    tibia=[11:8:67];
    for i=1:8
        %     figure;
        %     plot(legangle(i,:))
        %     hold on
        [maxang, minang] = Billauer_peakdet(legangle(i,:), 1);
        for j=1:size(minang,1)
            PEP(i,j)=minang(j,1);
            angex(i,j)=maxang(j,2)-minang(j,2);
            p1=filtering.bodypoints(:,tibia(i),maxang(j,1));
            p2=filtering.bodypoints(:,tibia(i),minang(j,1));
            legex(i,j)=sqrt((p2(2) - p1(2))^2+(p2(1) - p1(1))^2)*video.pix2cm;
        end
        for j=1:size(maxang,1)
            AEP(i,j)=maxang(j,1);
        end
        %
        %     plot(minang(:,1), minang(:,2), 'g*');
        %     plot(maxang(:,1), maxang(:,2), 'r*');
    end
    count=1;
    angexcur=[];
    legexcur=[];
    for i=1:2:7
        ang=[angex(i,:) angex(i+1,:)];
        leg=[legex(i,:) legex(i+1,:)];
        [r,c]=find(10<ang<50);
        angexcur{count}=ang(:,c);
        legexcur{count}=leg(:,c);
        count=count+1;
    end
    
    [r,c]=find(PEP==1|PEP==length(legangle));
    PEP(r,c)=0;
    [r,c]=find(AEP==1|AEP==length(legangle));
    AEP(r,c)=0;
    sortedPEP=Stillfried_arraysort(PEP,'sorted','nozeros');
    sortedPEP=sortedPEP(:,[1,3]);
    sortedAEP=Stillfried_arraysort(AEP,'sorted','nozeros');
    sortedAEP=sortedAEP(:,[1,3]);
    
    aeppep=[sortedAEP, ones(length(sortedAEP),1); sortedPEP, -ones(length(sortedPEP),1)];
    ind=Stillfried_arraysort(aeppep(:,2),'sorted','nozeros');
    ind=ind(:,1);%get indices for correct time ordering which are stored in first row
    timeordered=[];
    for j=1:length(ind)
        timeordered=[timeordered; aeppep(ind(j),:)];
    end
    
    for j=1:8
        ind=find(timeordered(:,1)==j);
        swing=[];
        stance=[];
        step=[];
        if timeordered(ind(1),3)==1
            step=[timeordered(ind(1),2)];
        end
        
        for i=2:length(ind)
            if timeordered(ind(i),3)==1 %a step is definded as two successive AEPs
                step=[step; timeordered(ind(i),2)];
                swing=[swing; timeordered(ind(i-1),2) timeordered(ind(i),2)];
            else if timeordered(ind(i),3)==-1
                    stance=[stance; timeordered(ind(i-1),2) timeordered(ind(i),2)];
                end
            end
        end
        steplegs{j}=step;
        swinglegs{j}=swing;
        stancelegs{j}=stance;
    end
    strides=max(cellfun(@length,steplegs))-1;
    strideduration=nan(strides,8);
    
    % kiri_plotStepPattern(stancelegs); %needs at least two steps
    
    for i=1:8
        swingduration{i}=(swinglegs{1,i}(:,2)-swinglegs{1,i}(:,1))*1/125;
        stanceduration{i}=(stancelegs{1,i}(:,2)-stancelegs{1,i}(:,1))*1/125;
        stepduration{i}=(diff(steplegs{1,i}))*1/125;
        for j=1:length(diff(steplegs{1,i}))
            strideduration(j,i)=stepduration{i}(j);
        end
    end
    
    count=1;
    legstep=[];
    for i=1:2:7
        leg=[strideduration(:,i) strideduration(:,i+1)];
        ind=~isnan(leg);
        legstep{count}=leg(ind);
        count=count+1;
    end
    
    j=~isnan(strideduration(:,3));
    strideperiod=strideduration(j,3);
    % p=anova1(strideduration,[],'off');
catch
    string=['** Error in ' str ' finding AEP and PEP **'];
    disp(string)
    swingduration=[];
    stanceduration=[];
    stepduration=[];
end

try
    aep=find(sortedAEP(:,1)==1);
    phaseL1=[];
    for i=1:length(aep)-1
        for j=aep(i)+1:aep(i+1)-1
            phaseL1=[phaseL1; sortedAEP(j,1) (sortedAEP(j,2)-sortedAEP(aep(i),2))/(sortedAEP(aep(i+1),2)-sortedAEP(aep(i),2))];
        end
    end
    for i=[2:8]
        ind=find(phaseL1(:,1)==i);
        L1phase{i}=phaseL1(ind,2);
    end
    aep=find(sortedAEP(:,1)==2);
    phaseR1=[];
    for i=1:length(aep)-1
        for j=aep(i)+1:aep(i+1)-1
            phaseR1=[phaseR1; sortedAEP(j,1) (sortedAEP(j,2)-sortedAEP(aep(i),2))/(sortedAEP(aep(i+1),2)-sortedAEP(aep(i),2))];
        end
    end
    for i=[1,3:8]
        ind=find(phaseR1(:,1)==i);
        R1phase{i}=phaseR1(ind,2);
    end
    aep=find(sortedAEP(:,1)==3);
    phaseL2=[];
    for i=1:length(aep)-1
        for j=aep(i)+1:aep(i+1)-1
            phaseL2=[phaseL2; sortedAEP(j,1) (sortedAEP(j,2)-sortedAEP(aep(i),2))/(sortedAEP(aep(i+1),2)-sortedAEP(aep(i),2))];
        end
    end
    for i=[1:2,4:8]
        ind=find(phaseL2(:,1)==i);
        L2phase{i}=phaseL2(ind,2);
    end
    aep=find(sortedAEP(:,1)==4);
    phaseR2=[];
    for i=1:length(aep)-1
        for j=aep(i)+1:aep(i+1)-1
            phaseR2=[phaseR2; sortedAEP(j,1) (sortedAEP(j,2)-sortedAEP(aep(i),2))/(sortedAEP(aep(i+1),2)-sortedAEP(aep(i),2))];
        end
    end
    for i=[1:3,5:8]
        ind=find(phaseR2(:,1)==i);
        R2phase{i}=phaseR2(ind,2);
    end
    aep=find(sortedAEP(:,1)==5);
    phaseL3=[];
    for i=1:length(aep)-1
        for j=aep(i)+1:aep(i+1)-1
            phaseL3=[phaseL3; sortedAEP(j,1) (sortedAEP(j,2)-sortedAEP(aep(i),2))/(sortedAEP(aep(i+1),2)-sortedAEP(aep(i),2))];
        end
    end
    for i=[1:4,6:8]
        ind=find(phaseL3(:,1)==i);
        L3phase{i}=phaseL3(ind,2);
    end
    aep=find(sortedAEP(:,1)==6);
    phaseR3=[];
    for i=1:length(aep)-1
        for j=aep(i)+1:aep(i+1)-1
            phaseR3=[phaseR3; sortedAEP(j,1) (sortedAEP(j,2)-sortedAEP(aep(i),2))/(sortedAEP(aep(i+1),2)-sortedAEP(aep(i),2))];
        end
    end
    for i=[1:5,7:8]
        ind=find(phaseR3(:,1)==i);
        R3phase{i}=phaseR3(ind,2);
    end
    aep=find(sortedAEP(:,1)==7);
    phaseL4=[];
    for i=1:length(aep)-1
        for j=aep(i)+1:aep(i+1)-1
            phaseL4=[phaseL4; sortedAEP(j,1) (sortedAEP(j,2)-sortedAEP(aep(i),2))/(sortedAEP(aep(i+1),2)-sortedAEP(aep(i),2))];
        end
    end
    for i=[1:6,8]
        ind=find(phaseL4(:,1)==i);
        L4phase{i}=phaseL4(ind,2);
    end
    aep=find(sortedAEP(:,1)==8);
    phaseR4=[];
    for i=1:length(aep)-1
        for j=aep(i)+1:aep(i+1)-1
            phaseR4=[phaseR4; sortedAEP(j,1) (sortedAEP(j,2)-sortedAEP(aep(i),2))/(sortedAEP(aep(i+1),2)-sortedAEP(aep(i),2))];
        end
    end
    for i=1:7
        ind=find(phaseR4(:,1)==i);
        R4phase{i}=phaseR4(ind,2);
    end
    
catch
    string=['** Error in ' str ' finding phase relationship **'];
    disp(string)
    R1phase={[]};
    R2phase={[]};
    R3phase={[]};
    R4phase={[]};
    L1phase={[]};
    L2phase={[]};
    L3phase={[]};
    L4phase={[]};
end

% figure
% plot([1:length(legangle)],legangle)

legangularvelocity=[];
for i=2:length(legangle)-1
    legangularvelocity=[legangularvelocity (legangle(:,i+1)-legangle(:,i-1))/(2/125)];%velocity in degrees/s
end
legangularvelocity=[nan(8,1) legangularvelocity nan(8,1)];
% figure
% plot([1:length(legangularvelocity)],legangularvelocity)
legangularacceleration=[];
for i=2:length(legangle)-1
    legangularacceleration=[legangularacceleration (legangle(:,i+1)- 2*legangle(:,i)+legangle(:,i-1))/((1/125)^2)];%acceleration in degrees/s^2
end
legangularacceleration=[nan(8,1) legangularacceleration nan(8,1)];

% figure
% plot([1:length(legangularacceleration)],legangularacceleration)
% for i=1:2:8
%     figure(i)
%     hold on
%     plot(legangle(i,5:length(legangle)-10),legangularvelocity(i,5:length(legangle)-10),'k')
%     plot(legangle(i+1,5:length(legangle)-10),legangularvelocity(i+1,5:length(legangle)-10),'k')
%     ylim([-1500 1500])
%     hold off
% end

legvel=[];
for i=1:8
    leg=[];
    for j=1:size(swinglegs{i},1)
        leg=[leg legangularvelocity(i,[swinglegs{i}(j,1):swinglegs{i}(j,2)])];
    end
    ind=~isnan(leg);
    legvel{i}=leg(ind);
end

swinglegvel=[];
count=1;
for i=1:2:7
    swinglegvel{count}=[legvel{i} legvel{i+1}];
    count=count+1;
end

legvel=[];
for i=1:8
    leg=[];
    for j=1:size(stancelegs{i},1)
        leg=[leg legangularvelocity(i,[stancelegs{i}(j,1):stancelegs{i}(j,2)])];
    end
    ind=~isnan(leg);
    legvel{i}=leg(ind);
end

stancelegvel=[];
count=1;
for i=1:2:7
    stancelegvel{count}=[legvel{i} legvel{i+1}];
    count=count+1;
end

% for i=1:2:8
%     figure(i)
%     for j=1:size(swinglegs{i},1)
%         [x,ind]=max(abs(legangularvelocity(i,[swinglegs{i}(j,1):swinglegs{i}(j,2)])));
%         hold on
%         plot(([-ind+1:length(legangularvelocity(i,[swinglegs{i}(j,1):swinglegs{i}(j,2)]))-ind]/length(legangularvelocity(i,[swinglegs{i}(j,1):swinglegs{i}(j,2)])))*100,legangularacceleration(i,[swinglegs{i}(j,1):swinglegs{i}(j,2)]),'k')
%         plot([0,0],[-800000,80000],'k')
%         plot([-80,80],[0,0],'k')
%         xlim([-80,80])
%         ylim([-80000,80000])
%         hold off
%     end
% end
%
% for i=1:2:8
%     figure(i+8)
%     for j=1:size(stancelegs{i},1)
%         [x,ind]=max(abs(legangularvelocity(i,[stancelegs{i}(j,1):stancelegs{i}(j,2)])));
%         hold on
%         plot(([-ind+1:length(legangularvelocity(i,[stancelegs{i}(j,1):stancelegs{i}(j,2)]))-ind]/length(legangularvelocity(i,[stancelegs{i}(j,1):stancelegs{i}(j,2)])))*100,legangularacceleration(i,[stancelegs{i}(j,1):stancelegs{i}(j,2)]),'k')
%         plot([0,0],[-80000,80000],'k')
%         plot([-80,80],[0,0],'k')
%         xlim([-80,80])
%         ylim([-80000,80000])
%         hold off
%     end
% end