 clear
 clc
%  load('SliptHalf_pca_tobacco.mat')
%headers=header_Cog;

% for the final plot, run two_wayCV.m
load('Error_2wayCV.mat')

%sub_area={Psychiatric, Alcohol, Tobacco, Drug, Cognition, Emotion, FamHis, FemHeal, Demo, PhyHeal, Motor, Personality, Sensory};
Error1 = Error1{1};
Error2 = Error2{1};
%sub_header = erase(sub_header,"FamHist_");


PaperDim=get(gcf,'PaperSize');
set(gcf,'PaperPosition',[0 0 PaperDim])
% Arbitrary Text, on background invisible axis
ht=axes('position',[0 0	1 1],'Visible','off');
text(ht,.4,.99,'Psychiatry','FontSize',20)

text(ht,.1,.97,'Variable weights','FontSize',14)

hf=axes('position',[.11 .75 .82 .2]);
plot(weight(:,1:3),'linewidth',2);
hold on
stem(weight(:,1:3),':', 'linewidth',1)
set(gca,'XTickLabel',sub_header,'XTickLabelRotation',45,'XGrid','on',...
             'XTick',1:size(sub_header,1),'TickLabelInterpreter','none','FontName', 'Arial Narrow','Fontsize',9);
%ax = gca;
%ax.TickLabelInterpreter = 'latex';
%ax.XTickLabel{1} = ['\color{red}' ax.XTickLabel{1}];
xlim([0 size(dd,2)+1])
legend(strcat({'PC '},num2str([1:3]')),'Location','Northeast','Fontsize',7) %[0.45 0.15 0.05 0.03]
ylabel('PC weight')
title('First 3 Pricipal Loadings')  



hf=axes('position',[.11 .46 .82 .19]);
plot(weight_r,'linewidth',2);
hold on
stem(weight_r(:,1:3),':', 'linewidth',1)
set(gca,'XTickLabel',sub_header(:,1),'XTickLabelRotation',45,'XGrid','on',...
             'XTick',1:size(sub_header,1),'TickLabelInterpreter','none','FontName', 'Arial Narrow','Fontsize',9);
legend(strcat({'RC '},num2str([1:3]')),'Location','Northeast','Fontsize',7)
ylabel('RC weight')
xlim([0 size(dd,2)+1])
title('First 3 Rotated Component (RC) Loadings','Interpreter','none')


text(ht,.1,.37,'Optimal number of PCs','FontSize',14)

hf=axes('position',[.14 .15 .7 .19]);
plot(Error1(1:12), 'k.--')
hold on
plot(Error2(1:12), 'r.-')
legend({ 'Simple method', 'Two-way CV method'}, ...
    'Location', 'Northwest')
legend boxoff
set(gca, 'XTick', 1:12)
%set(gca, 'YTick', [])
xlabel('Number of PCs')
ylabel('Cross-validation error')
title('Error analysis using different number of PCs') 


print -dpng Report-Psy2.png  
