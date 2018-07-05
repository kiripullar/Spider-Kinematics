function kiri_plotStepPattern(stancelegs)
% kiri_plotStepPattern takes touchdown and liftoff times of stance legs
% and creates a stepping pattern diagram with thick black lines
% representing stance phase (retraction), and white space representing a 
% swing phase (protraction)
%
% Kiri Pullar, masters thesis 2009

x=[];
for i=1:8
    for j=1:length(stancelegs{i})
        for k=stancelegs{i}(j,1):stancelegs{i}(j,2)
            x(k,i)=1;
        end
    end
end
[r,c]=find(x==0);
for i=1:length(r)
    x(r(i),c(i))=nan;
end

figure

count=8;
for i=2:2:8
    hold on
    plot(count*x(:,i),'k', 'LineWidth', 8)
    count=count-1;
    hold off
end
for i=1:2:8
    hold on
    plot(count*x(:,i),'k', 'LineWidth', 8)
    count=count-1;
    hold off
end

axis ij tight
ylim([0 9])
set(gca, 'ytick',[1 2 3 4 5 6 7 8], 'yticklabel',['R1'; 'R2'; 'R3'; 'R4'; 'L1'; 'L2'; 'L3'; 'L4'])
xlabel('Time (s)')
ylabel('Leg')

ax=get(gca,'xtick');
new=12.5:12.5:ax(length(ax));
set(gca,'xtick',new)
new=new/125;
set(gca,'xticklabel',new)
end
