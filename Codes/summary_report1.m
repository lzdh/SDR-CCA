 clear
 clc
% load('SliptHalf_pca_Cognition.mat')
% headers=header_Cog;
% refer to 'SubAnalysis.m' and 'SplitHalf_CV.m'

PaperDim=get(gcf,'PaperSize');
set(gcf,'PaperPosition',[0 0 PaperDim])
% Arbitrary Text, on background invisible axis
ht=axes('position',[0 0	1 1],'Visible','off');
text(ht,.3,.98,'Psychiatry','FontSize',20)
text(ht,.6,.98, '44 variables','FontSize',14)

text(ht,.1,.95,'Factor Definition','FontSize',14)
% text(ht,.1,.3,'Missingness Report','FontSize',14)

% Plotting axis
hf=axes('position',[.14 .71 .7 .2]);
dd_n=mean(d_n);
dom_n=sum(dd_n);
dom=sum(mean(d));
pca_n=0;
pca=0;
p=0;
for i=1:size(d,2)
    pca_n=dd_n(i)/dom_n*100;
    plot(i,pca_n,'go')
    hold on
    pca=dd(i)/dom*100;
    plot(i,pca,'bo')
    hold on
    p=p+pca;
    plot(i,p,'ro')
    hold on
end
xlabel('Number of Principal Components (PCs)');
xlim([1 size(dd,2)])
ylabel('% of variance explained by each PC');%# Add a label to the left y axis
ylim([0,100])
legend('Null eigen-spectrum', 'Covariance eigen-spectrum','Cumulative variance explained','Location','best')
yyaxis right
ylim([0,100])
ylabel('% of cumulative variance explained')
title('Scree plot & Null eigen-spectrum')

%------------------------------------
hf=axes('position',[.14 .44 .7 .2]);
dom=sum(mean(d));
dd=mean(d);
p=0;
for i=1:size(d,2)
    pca=dd(i)/dom*100;
    p=p+pca;
    plot(i,p,'ro')
    hold on
    plot(rbar*100,'bo')
end
xlim([1 size(dd,2)])
ylim([40 100])
legend({'In sample', 'Out sample'}, 'Location','SouthEast')
xlabel('Number of Principal Components (PCs)')
ylabel(gca,'% of cumulative variance explained');%# Add a label to the left y axis
title('PCs in split-half cross validation analysis')

%------------------------------------
text(ht,.1,.38,'PC summary table','FontSize',14)
% Table
uitable('Data', [1 dd(1)/sum(dd)*100 rbar(1)*100;2 dd(2)/sum(dd)*100 (rbar(2)-rbar(1))*100; 3 dd(3)/sum(dd)*100 (rbar(3)-rbar(2))*100],...
        'RowName',[],...
        'ColumnName', {'PC', 'Train %var', 'Test %var'},...
        'Units', 'Normalized',...
        'Position',[0.09 0.27 .39 .09]);

uitable('Data', [50 2 2;70 5 5; 90 13 14],...
        'RowName',[],...
        'ColumnName', {'>%var', '#PCs Train', '#PCs Test'},...
        'Units', 'Normalized',...
        'Position',[0.48 0.27 .39 .09]);    
      
% [bla,ind]=sort(abs([header_Alc{:,2}]), 'descend');
% weight=header_Alc(ind,:);

% uitable('Data', [weight(1,1) weight(2,1) weight(3,1);weight(1,2), weight(2,2) weight(3,2);...
%                  weight(4,1) weight(5,1) weight(6,1);weight(4,2), weight(5,2) weight(6,2);],...
%         'RowName',{'Ordered var name', 'Weights','Ordered var name', 'Weights'},...
%         'ColumnName', [],...
%         'Units', 'Normalized',...
%         'ColumnWidth', {107  107  107  107},...
%         'Position',[0.11 0.38 .8 .095]);   
       
         
print -dpng Report-Psy.png
