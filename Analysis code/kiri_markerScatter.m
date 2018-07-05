function kiri_markerScatter(x,y,individual,xstring,ystring)
% kiri_markerScatter takes the same parameters as the scatter function in
% matlab but instead of representing individuals as different colours on
% graphs cycles through different marker shapes to make identification of
% individuals easier particularly in printed documents
%
% Kiri Pullar, masters thesis 2009

symbol=['v';'s';'o';'^';'p';'d';'*';'x';'+';'.'];

for i=1:length(individual)
switch individual(i)
case 1
hold on
scatter(x(i),y(i),'k',symbol(1),'filled','XDataSource',[xstring '(' num2str(i) ')'],'YDataSource',[ystring '(' num2str(i) ')'])
hold off
case 2
hold on
scatter(x(i),y(i),'k',symbol(2),'filled','XDataSource',[xstring '(' num2str(i) ')'],'YDataSource',[ystring '(' num2str(i) ')'])
hold off
case 3
hold on
scatter(x(i),y(i),'k',symbol(3),'filled','XDataSource',[xstring '(' num2str(i) ')'],'YDataSource',[ystring '(' num2str(i) ')'])
hold off
case 4
hold on
scatter(x(i),y(i),'k',symbol(4),'filled','XDataSource',[xstring '(' num2str(i) ')'],'YDataSource',[ystring '(' num2str(i) ')'])
hold off
case 5
hold on
scatter(x(i),y(i),'k',symbol(5),'filled','XDataSource',[xstring '(' num2str(i) ')'],'YDataSource',[ystring '(' num2str(i) ')'])
hold off
case 6
hold on
scatter(x(i),y(i),'k',symbol(6),'filled','XDataSource',[xstring '(' num2str(i) ')'],'YDataSource',[ystring '(' num2str(i) ')'])
hold off
case 7
hold on
scatter(x(i),y(i),'k',symbol(7),'XDataSource',[xstring '(' num2str(i) ')'],'YDataSource',[ystring '(' num2str(i) ')'])
hold off
case 8
hold on
scatter(x(i),y(i),'k',symbol(8),'XDataSource',[xstring '(' num2str(i) ')'],'YDataSource',[ystring '(' num2str(i) ')'])
hold off
case 9
hold on
scatter(x(i),y(i),'k',symbol(9),'XDataSource',[xstring '(' num2str(i) ')'],'YDataSource',[ystring '(' num2str(i) ')'])
hold off
case 10
hold on
scatter(x(i),y(i),'k',symbol(10),'XDataSource',[xstring '(' num2str(i) ')'],'YDataSource',[ystring '(' num2str(i) ')'])
hold off
case 11
hold on
scatter(x(i),y(i),'k',symbol(1),'XDataSource',[xstring '(' num2str(i) ')'],'YDataSource',[ystring '(' num2str(i) ')'])
hold off
case 12
hold on
scatter(x(i),y(i),'k',symbol(2),'XDataSource',[xstring '(' num2str(i) ')'],'YDataSource',[ystring '(' num2str(i) ')'])
hold off
case 13
hold on
scatter(x(i),y(i),'k',symbol(3),'XDataSource',[xstring '(' num2str(i) ')'],'YDataSource',[ystring '(' num2str(i) ')'])
hold off
case 14
hold on
scatter(x(i),y(i),'k',symbol(4),'XDataSource',[xstring '(' num2str(i) ')'],'YDataSource',[ystring '(' num2str(i) ')'])
hold off
case 15
hold on
scatter(x(i),y(i),'k',symbol(5),'XDataSource',[xstring '(' num2str(i) ')'],'YDataSource',[ystring '(' num2str(i) ')'])
hold off
case 16
hold on
scatter(x(i),y(i),'k',symbol(6),'XDataSource',[xstring '(' num2str(i) ')'],'YDataSource',[ystring '(' num2str(i) ')'])
hold off
case 17
hold on
scatter(x(i),y(i),'k',symbol(7),'XDataSource',[xstring '(' num2str(i) ')'],'YDataSource',[ystring '(' num2str(i) ')'])
hold off
case 18
hold on
scatter(x(i),y(i),'k',symbol(8),'XDataSource',[xstring '(' num2str(i) ')'],'YDataSource',[ystring '(' num2str(i) ')'])
hold off
case 19
hold on
scatter(x(i),y(i),'k',symbol(9),'XDataSource',[xstring '(' num2str(i) ')'],'YDataSource',[ystring '(' num2str(i) ')'])
hold off
case 20
hold on
scatter(x(i),y(i),'k',symbol(10),'XDataSource',[xstring '(' num2str(i) ')'],'YDataSource',[ystring '(' num2str(i) ')'])
hold off
end
end