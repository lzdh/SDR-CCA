clear
clc

% Load data
vars0=importdata('../Data/MF_psy_vars.mat'); %original SM dataset with FS and all confounders in
%vars0=importdata('../Data/S815_unflipped.mat'); %original SM dataset with FS and all confounders in

%vars=importdata('Data/SDR_vars_PC.mat'); %SDR SM dataset of dim.55 with Demo
%vars=importdata('Data/SDR_SMnoDem.mat');
%vars=importdata('Data/deconf_SDR_SM.mat'); %SDR SM dataset of dim.75 with Demo
%vars=importdata('../Data/SDR_SM_predeconf.mat'); %SDR SM dataset of dim.70
vars1=importdata('Data/SDR_SM_NOpredeconf.mat'); %SDR SM dataset of dim.47
%vars=importdata('../Data/SDR_SM_dimRematch.mat'); %SDR SM dataset of dim.70
vars1=importdata('../Data/SDR_SM_dim60.mat'); %SDR SM dataset of dim.58
%vars=importdata('../Data/SDR_SM_dim60_unfilpped.mat'); %SDR SM dataset of dim.58


%BMvars=importdata('Data/SDR_NET1_PCnew.mat');%SDR BM of dim.1571
BMvars=importdata('../Data/SDR_NET1_newconf.mat');%SDR BM of dim.1571

NET=importdata('../Data/NET.mat');%original BM of dim.40000
%load('Data/SMvars.mat')%SM variables that were fed into SDR
%SMvars=importdata('Data/SMvars_noDemo.mat');%SM variables that were fed into SDR
%SMvars=importdata('../Data/SDR_CCAvars_3unflipped.mat');%SM variables that were fed into SDR
SMvars=importdata('../Data/SDR_CCAvars_3.mat');%SM variables that were fed into SDR

%SM_header=importdata('../SM_header_noDemo.mat');
%SM_header=importdata('Data/SDR_CCAheader_preDeconf.mat');
%SDR_CCAheader = importdata('Data/SDR_CCAheader_preDeconf.mat');
%SDR_CCAheader = importdata('../Header/.mat');

%SM_header=importdata('SM_header.txt');
header=importdata('../Header/MFvarnames.txt');
%header=importdata('../Header/headers.txt');

header_sdrvars=importdata('../Header/SDRHeader_SM213withsigns.mat');
%header_sdrvars=importdata('../Header/SDRHeader_SM213nosigns.mat');

%header_sdr=importdata('../Header/SDR_header58.mat');
header_sdr=importdata('../Header/SDR_header60_unflipped.mat');

%Prepare the permutation set
Nperm=10000;                                                                       % in the paper we used 100000 but 10000 should be enough
EB=hcp2blocks('../RESTRICTED_lzdh_3_20_2016_11_24_35.csv', [ ], false, vars0(:,1)); % change the filename to your version of the restricted file
PAPset=palm_quickperms([ ], EB, Nperm);                                            % the final matrix of permuations

Nkeep=100;

%%% setup confounds matrix
conf=palm_inormal([vars0(:,[2 3 4 7 14 15 22 23 25 479 480 481 486]) vars0(:,[265 266]).^(1/3) ]); % Gaussianise.  479 480 481 race variables 486 Ethnicity
conf(isnan(conf))=0;  % impute missing data as zeros
conf=nets_normalise([conf conf(:,5:end).^2]);  % add on squared terms and renormalise

% Remove bad variables which have more than half NaNs or more than 99% of
% their values are the same
%  badvars=[];
% for i=1:size(SMvars,2)
%   Y=SMvars(:,i); grotKEEP=~isnan(Y);  
%   %change the creterion when # of subjects changes
%   if (sum(grotKEEP)>400) & (max(sum(nets_class_vectomat(Y(grotKEEP))))/length(Y(grotKEEP))<0.95)
%       % the 3rd thing above is:  is the size of the largest equal-values-group too large?
%     i=i; % do nothing
%   else
%     [i sum(grotKEEP)  max(sum(nets_class_vectomat(Y(grotKEEP))))/length(Y(grotKEEP))]
%     badvars=[badvars i];
%   end 
% end

% varskeep=setdiff([1:size(SMvars,2)], badvars);  
% SMvars=SMvars(:,varskeep);
% SDR_CCAheader=SM_header(varskeep');

% Normalise and deconfound the SM variables
varsgrot=palm_inormal(SMvars);
for i=1:size(varsgrot,2)
  grot=(isnan(varsgrot(:,i))==0); grotconf=nets_demean(conf(grot,:)); varsgrot(grot,i)=nets_normalise(varsgrot(grot,i)-grotconf*(pinv(grotconf)*varsgrot(grot,i)));
end
% varsgrot(isnan(varsgrot))=0;

%%% prepare main BM matrix
NET1=nets_demean(NET);  NET1=NET1/std(NET1(:)); % no norm
grot=NET1;
NETd=nets_demean(grot-conf*(pinv(conf)*grot));   % deconfound and demean

%%% Do this if we want to feed PCA BM to CCA 
% [uu,ss,vv]=nets_svds(NETd,Nkeep); % SVD reduction
% BMpc=uu*ss;


%Apply PCA to further reduce the dimension of SDR BM
[uu,ss,vv]=nets_svds(BMvars,100); % SVD reduction
BMpc=uu*ss;
BMpc=nets_demean(BMpc-conf*(pinv(conf)*BMpc)); % Deconfound and demean SDR BM

% Normalise and deconfound SDR SM to feed into CCA
varsd=palm_inormal(vars);
for i=1:size(varsd,2)
  grot=(isnan(varsd(:,i))==0); grotconf=nets_demean(conf(grot,:)); varsd(grot,i)=nets_normalise(varsd(grot,i)-grotconf*(pinv(grotconf)*varsd(grot,i)));
end
varsd(isnan(varsd))=0;

%%% CCA/
[grotA,grotB,grotR,grotU,grotV,grotstats]=canoncorr(BMpc,varsd);


%%% CCA permutation testing
grotRp=[]; clear grotRpval;
for j=1:Nperm
  j
  [grotAr,grotBr,grotRp(j,:),grotUr,grotVr,grotstatsr]=canoncorr(BMpc,varsd(PAPset(:,j),:));
end
for i=1:size(grotRp,2);  % get FWE-corrected pvalues
  grotRpval(i)=(1+sum(grotRp(2:end,1)>=grotR(i)))/Nperm;
end
grotRpval
Ncca=sum(grotRpval<0.05)  % number of FWE-significant CCA components


%%%%%%%%%%%%%%% amount of variance between explained by canonical variates in original variables
for j=1:5
    TVE_SM(j,:)=corr(grotV(:,j),varsgrot,'rows','pairwise').^2;%varsgrot
    TVE_Net(j,:)=corr(grotU(:,j),NETd).^2; %NETd(:,1:size(NET,2))
end
AVE_SM=mean(TVE_SM,2,'omitnan');
AVE_Net=mean(TVE_Net,2,'omitnan');

for i=1:length(AVE_SM)
    cum(i)=sum(AVE_SM(1:i)');  
end

% Alternatively calculating the R-sqaured
% Note: if any variables in the original dataset (varsgrot) has too many NaNs, this method will lose precision. 
varsgrot(isnan(varsgrot))=0;
b1=grotV(:,1:3)\varsgrot;
yhat=grotV(:,1:3)*b1;
Rsq = 1 - sum(sum((varsgrot - yhat).^2))/sum(sum((varsgrot).^2))

b2=grotU(:,1:3)\NETd;
yhat2=grotU(:,1:3)*b2;
Rsq2 = 1 - sum((NETd - yhat2).^2)/sum((NETd).^2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% BM and SM weights for sig. CCA modes
varsgrot1=palm_inormal(vars0);
for i=1:size(varsgrot1,2)
  grot=(isnan(varsgrot1(:,i))==0); grotconf=nets_demean(conf(grot,:)); varsgrot1(grot,i)=nets_normalise(varsgrot1(grot,i)-grotconf*(pinv(grotconf)*varsgrot1(grot,i)));
end

for i=1:3

    grotAA(:,i) = corr(grotU(:,i),NET)';
     % or
    grotAAd(:,i) = corr(grotU(:,i),NETd)'; % weights after deconfounding    
    
    %SM
    grotBBd_sdr(:,i) = corr(grotV(:,i),varsgrot,'rows','pairwise')'; % weights after deconfounding 
    grotBBd(:,i) = corr(grotV(:,i),varsgrot1,'rows','pairwise')'; % weights after deconfounding
    grotBBd1(:,i) = corr(grotV(:,i),varsd,'rows','pairwise')'; % weights after deconfounding

end

grotAAd(isnan(grotAAd))=0;

%%% remove variables with all weights equal to 0 or too small so the
%%% rotation can be applied
grotBBd(isnan(grotBBd))=0;
rmvars=[];
for i=1:487
    for j=1:3
        if grotBBd(i,j)<10^(-3)
            grotBBd(i,j)=0; end
    end
    
    if grotBBd(i,1)~=0 || grotBBd(i,2)~=0 || grotBBd(i,3)~=0
       i=i;
    else
        rmvars=[rmvars i];
    end
    
end  


grotBBd(rmvars,:)=[];
header(rmvars)=[];    

[Bsdr_rot,T] = rotatefactors(grotBBd_sdr);
B_rot = rotatefactors(grotBBd);
Bt_rot = grotBBd*T;
B1_rot=rotatefactors(grotBBd1);

[~,I]=sort(abs(B_rot(:,1)),'descend');
[~,I0]=sort(abs(Bt_rot(:,1)),'descend');
[~,I1]=sort(abs(Bsdr_rot(:,1)),'descend');
[~,I2]=sort(abs(B1_rot(:,1)),'descend');
SM_rotw81=[header(I), num2cell(B_rot(I,1))];
SMt_rotw81=[header(I0), num2cell(Bt_rot(I0,1))];
SMsdr_rotw81=[header_sdrvars(I1), num2cell(Bsdr_rot(I1,1))];
SM1_rotw81=[header_sdr(I2), num2cell(B1_rot(I2,1))];

[~,I]=sort(abs(B_rot(:,2)),'descend');
[~,I0]=sort(abs(Bt_rot(:,2)),'descend');
[~,I1]=sort(abs(Bsdr_rot(:,2)),'descend');
[~,I2]=sort(abs(B1_rot(:,2)),'descend');
SM_rotw82=[header(I), num2cell(B_rot(I,2))];
SMt_rotw82=[header(I0), num2cell(Bt_rot(I0,2))];
SMsdr_rotw82=[header_sdrvars(I1), num2cell(Bsdr_rot(I1,2))];
SM1_rotw82=[header_sdr(I2), num2cell(B1_rot(I2,2))];

[~,I]=sort(abs(B_rot(:,3)),'descend');
[~,I0]=sort(abs(Bt_rot(:,3)),'descend');
[~,I1]=sort(abs(Bsdr_rot(:,3)),'descend');
[~,I2]=sort(abs(B1_rot(:,3)),'descend');
SM_rotw83=[header(I), num2cell(B_rot(I,3))];
SMt_rotw83=[header(I0), num2cell(Bt_rot(I0,3))];
SMsdr_rotw83=[header_sdrvars(I1), num2cell(Bsdr_rot(I1,3))];
SM1_rotw83=[header_sdr(I2), num2cell(B1_rot(I2,3))];

%%% Plot CCA weights for both methods: Note here the headers are different
%%% for each method in terms of the variable orders!

% % correlation method
% [~,idx1] = sort(abs(L1_r),'descend');
% figure(1)
% plot(grotBBd1(idx1(1:30)),'bo')
% set(gca,'XTickLabel',SDR_CCAheader(idx1(1:30)),'XTick',1:30,'XTickLabelRotation',45,'TickLabelInterpreter','none');
% 
% % plot(grotBBd1,'bo')
% % set(gca,'XTickLabel',SDR_CCAheader(),'XTick',1:70,'XTickLabelRotation',45,'TickLabelInterpreter','none');
% 
% [~,idx2] = sort(abs(L2_r),'descend');
% figure(2)
% plot(grotBBd2(idx2(1:30)),'bo')
% set(gca,'XTickLabel',SDR_CCAheader(idx2(1:30)),'XTick',1:30,'XTickLabelRotation',45,'TickLabelInterpreter','none');
% 
% % [~,idx3] = sort(abs(grotBBd3),'descend');
% % figure(3)
% % plot(grotBBd3(idx3(1:30)),'bo')
% % set(gca,'XTickLabel',SDR_CCAheader(idx1(1:30)),'XTick',1:30,'XTickLabelRotation',45,'TickLabelInterpreter','none');
% 
% %load('Data/SDR_pcaloadings.mat');
% load('../Data/PL_SM_predeconf.mat');
% %load('../Data/PL_SM_dimRematch.mat');
% 
% %SDR_CCAheader1= importdata('Data/SDR_CCAheader.mat');
%  SDR_CCAheader1=importdata('../Header_preDeconf_SDR_SM.mat');
% %SDR_CCAheader1=[SDR_CCAheader1; header([3:5 7 479:481 486])];
% 
% 
% CCAweight1=[];
% CCAweight2=[];
% % CCAweight3=[];
% k=0;
% for d=1:length(pcaloadings)
%     
%     w1=pcaloadings{d}*grotB((k+1):(k+size(pcaloadings{d},2)),1);
%     CCAweight1=[CCAweight1;w1]; 
%     
%     w2=pcaloadings{d}*grotB((k+1):(k+size(pcaloadings{d},2)),2);
%     CCAweight2=[CCAweight2;w2]; 
%     
% %     w3=pcaloadings{d}*grotB((k+1):(k+size(pcaloadings{d},2)),3);
% %     CCAweight3=[CCAweight3;w3]; 
%     
%     k=size(pcaloadings{d},2);
% end
% % CCAweight1=[CCAweight1; grotB(48:55,1)];
% % CCAweight2=[CCAweight2; grotB(48:55,2)];
% % CCAweight3=[CCAweight3; grotB(48:55,3)];
% 
% w_r=rotatefactors([CCAweight1, CCAweight2]);
% 
% [~,idx4] = sort(abs(w_r(:,1)),'descend');
% figure(4)
% plot(w_r(idx4(1:30),1),'bo')
% set(gca,'XTickLabel',SDR_CCAheader1(idx4(1:30)),'XTick',1:30,'XTickLabelRotation',45,'TickLabelInterpreter','none');
% 
% [~,idx5] = sort(abs(w_r(:,1)),'descend');
% figure(5)
% plot(w_r(idx5(1:30),2),'bo')
% set(gca,'XTickLabel',SDR_CCAheader1(idx5(1:30)),'XTick',1:30,'XTickLabelRotation',45,'TickLabelInterpreter','none');
% 
% % [~,idx6] = sort(abs(CCAweight3),'descend');
% % figure(6)
% % plot(CCAweight3(idx6(1:30)),'bo')
% % set(gca,'XTickLabel',SDR_CCAheader1(idx6(1:30)),'XTick',1:30,'XTickLabelRotation',45,'TickLabelInterpreter','none');

%% CCA weights for 2 CCA modes

n=30; %top n variables to look into for each CCA mode
grotBBd(isnan(grotBBd))=0;
rmvars=[];
for i=1:415
    for j=1:2
        if grotBBd(i,j)<10^(-3)
            grotBBd(i,j)=0; end
    end
    
    if grotBBd(i,1)~=0 || grotBBd(i,2)~=0
       i=i;
    else
        rmvars=[rmvars i];
    end
    
end  


grotBBd(rmvars,:)=[];
header(rmvars)=[];    

[Bsdr_rot,T] = rotatefactors(grotBBd_sdr(:,1:2));
B_rot = rotatefactors(grotBBd(:,1:2));
Bt_rot = grotBBd(:,1:2)*T;
B1_rot=rotatefactors(grotBBd1(:,1:2));

%%%%% Sphere area
BB_rot=sqrt(Bt_rot(:,1).^2+Bt_rot(:,2).^2);
[~,I0]=sort(BB_rot,'descend');
EDtopw=Bt_rot(I0(1:n),:);

fid = fopen('../Header/MF_VarNmsLongShort_plot.csv');
%fid = fopen('../Header/VarNmsLongShort_plot.csv');
tmp = textscan(fid, '%s%s%d%s%s%s%f%f%f','Delimiter',',','Headerlines',1);fclose(fid);
[~,Ikeep]=ismember(header(I0(1:n)),tmp{1});
for i=1:size(tmp,2)
    temp=tmp{i};
    tmp{i}=temp(Ikeep,:);
end
Nms=tmp{2};
ShowLab=double(tmp{3});
Categ=[tmp{4},tmp{6}];
CategGrp=tmp{5};
Color=double([tmp{7}, tmp{8}, tmp{9}]);
Color=num2cell(Color);
Categ=[Categ, Color];


figure(1)
mdsplot(EDtopw(:,1:2),{},[],[],Categ(),CategGrp(),'SDR CCA rotated loadings - Symbols & color by group', ...
         'Top 30 sphere SDR-CCA rotated loadings with symbols', {},{'CCA 1', 'CCA 2'})

figure(2)
mdsplot(EDtopw(:,1:2),Nms(),[],[],Categ(),CategGrp(),'SDR CCA rotated loadings - Names & color by group', ...
         'Top 30 sphere SDR-CCA rotated loadings with names', {},{'CCA 1', 'CCA 2'})


%% Analysis on rotated CCA weights (all variables)

n=30; %top n variables to look into for each CCA mode

%%%%% Sphere area
BB_rot=sqrt(Bt_rot(:,1).^2+Bt_rot(:,2).^2+Bt_rot(:,3).^2);
[~,I0]=sort(BB_rot,'descend');
EDtopw=Bt_rot(I0(1:n),:);

fid = fopen('../Header/MF_VarNmsLongShort_plot.csv');
%fid = fopen('../Header/VarNmsLongShort_plot.csv');
tmp = textscan(fid, '%s%s%d%s%s%s%f%f%f','Delimiter',',','Headerlines',1);fclose(fid);
[~,Ikeep]=ismember(header(I0(1:n)),tmp{1});
for i=1:size(tmp,2)
    temp=tmp{i};
    tmp{i}=temp(Ikeep,:);
end
Nms=tmp{2};
ShowLab=double(tmp{3});
Categ=[tmp{4},tmp{6}];
CategGrp=tmp{5};
Color=double([tmp{7}, tmp{8}, tmp{9}]);
Color=num2cell(Color);
Categ=[Categ, Color];


figure(1)
mdsplot3(EDtopw(:,1:3),{},[],[],Categ(),CategGrp(),'SDR CCA rotated loadings - Symbols & color by group', ...
         'Top 30 sphere SDR-CCA 3dim. rotated loadings with symbols', {},{'CCA 1', 'CCA 2','CCA 3'})
rotate3d on 

figure(2)
mdsplot3(EDtopw(:,1:3),Nms(),[],[],Categ(),CategGrp(),'SDR CCA rotated loadings - Naemes & color by group', ...
         'Top 30 sphere SDR-CCA 3dim. rotated loadings with names', {},{'CCA 1', 'CCA 2','CCA 3'})
rotate3d on    








%%%%% Cubic area

% name=union(SMt_rotw81(1:n,1), SMt_rotw82(1:n,1));
% name1=union(name, SMt_rotw83(1:n,1));
% 
% [~,loc1]=ismember(name1,SM_rotw81(:,1));
% [~,loc2]=ismember(name1,SM_rotw82(:,1));
% [~,loc3]=ismember(name1,SM_rotw83(:,1));
% 
% T_rotw8=[SM_rotw81(loc1,2) SM_rotw82(loc2,2) SM_rotw83(loc3,2)];
% 
% T_rotw8=cell2mat(T_rotw8);
% 
% fid = fopen('../Header/MF_VarNmsLongShort_plot.csv');
% tmp = textscan(fid, '%s%s%d%s%s%s%f%f%f','Delimiter',',','Headerlines',1);fclose(fid);
% [~,Ikeep]=ismember(name1',tmp{1});
% for i=1:size(tmp,2)
%     temp=tmp{i};
%     tmp{i}=temp(Ikeep,:);
% end
% Nms=tmp{2};
% ShowLab=double(tmp{3});
% Categ=[tmp{4},tmp{6}];
% CategGrp=tmp{5};
% Color=double([tmp{7}, tmp{8}, tmp{9}]);
% Color=num2cell(Color);
% Categ=[Categ, Color];
% 
% figure(1)
% mdsplot3(T_rotw8(:,1:3),{},[],[],Categ(),CategGrp(),'CCA rotated loadings - Symbols & color by group', ...
%          'SDR-CCA 3dim. rotated loadings with symbols', {},{'CCA 1', 'CCA 2','CCA 3'})
% rotate3d on 
% 
% figure(2)
% mdsplot3(T_rotw8(:,1:3),Nms(),[],[],Categ(),CategGrp(),'CCA rotated loadings - Naemes & color by group', ...
%          'SDR-CCA 3dim. rotated loadings with names', {},{'CCA 1', 'CCA 2','CCA 3'})
% rotate3d on    
%      
% figure
% pcaplot(T_rotw8(:,1:3),Nms(),[],[],Categ(),CategGrp(),'CCA rotated loadings - Names, color by group',... 
%     'SDR-CCA 2dim. rotated loadings with names','CCA3_1', 'CCA3_2')
% 
% pcaplot(T_rotw8(:,2:3),{},[],[],Categ(),CategGrp(),'CCA rotated loadings -  Symbols, color by group',... 
%     'SDR-CCA 2dim. rotated loadings with symbols','CCA3_2', 'CCA3_3')
% 

%% Analysis on rotated CCA weights (SDR varibles)

%%%%% Sphere area
BBsdr_rot=sqrt(Bsdr_rot(:,1).^2+Bsdr_rot(:,2).^2+Bsdr_rot(:,3).^2);
[~,I0]=sort(BBsdr_rot,'descend');
EDtopw=Bsdr_rot(I0(1:n),:);

%fid = fopen('../Header/VarNmsLongShort_plot.csv');
fid = fopen('../Header/MF_VarNmsLongShort_plot.csv');

tmp = textscan(fid, '%s%s%d%s%s%s%f%f%f','Delimiter',',','Headerlines',1);fclose(fid);
[~,Ikeep]=ismember(header_sdrvars(I0(1:n)),tmp{1});
for i=1:size(tmp,2)
    temp=tmp{i};
    tmp{i}=temp(Ikeep,:);
end
Nms=tmp{2};
ShowLab=double(tmp{3});
Categ=[tmp{4},tmp{6}];
CategGrp=tmp{5};
Color=double([tmp{7}, tmp{8}, tmp{9}]);
Color=num2cell(Color);
Categ=[Categ, Color];


figure(1)
mdsplot3(EDtopw(:,1:3),{},[],[],Categ(),CategGrp(),'SDRvars CCA rotated loadings - Symbols & color by group', ...
         'Top 30 sphere SDRvars-CCA 3dim. rotated loadings with symbols', {},{'CCA 1', 'CCA 2','CCA 3'})
rotate3d on 

figure(2)
mdsplot3(EDtopw(:,1:3),Nms(),[],[],Categ(),CategGrp(),'SDRvars CCA rotated loadings - Naemes & color by group', ...
         'Top 30 sphere SDRvars-CCA 3dim. rotated loadings with names', {},{'CCA 1', 'CCA 2','CCA 3'})
rotate3d on    


%%%%Cubic area
name=union(SMsdr_rotw81(1:n,1), SMsdr_rotw82(1:n,1));
name1=union(name, SMsdr_rotw83(1:n,1));

[~,loc1]=ismember(name1,SMsdr_rotw81(:,1));
[~,loc2]=ismember(name1,SMsdr_rotw82(:,1));
[~,loc3]=ismember(name1,SMsdr_rotw83(:,1));

T_rotw8=[SMsdr_rotw81(loc1,2) SMsdr_rotw82(loc2,2) SMsdr_rotw83(loc3,2)];

T_rotw8=cell2mat(T_rotw8);

fid = fopen('../Header/MF_VarNmsLongShort_plot.csv');
tmp = textscan(fid, '%s%s%d%s%s%s%f%f%f','Delimiter',',','Headerlines',1);fclose(fid);
[~,Ikeep]=ismember(name1',tmp{1});
for i=1:size(tmp,2)
    temp=tmp{i};
    tmp{i}=temp(Ikeep,:);
end
Nms=tmp{2};
ShowLab=double(tmp{3});
Categ=[tmp{4},tmp{6}];
CategGrp=tmp{5};
Color=double([tmp{7}, tmp{8}, tmp{9}]);
Color=num2cell(Color);
Categ=[Categ, Color];

figure(1)
mdsplot3(T_rotw8(:,1:3),{},[],[],Categ(),CategGrp(),'CCA rotated SDRvars loadings - Symbols & color by group', ...
         'SDR-CCA 3dim. rotated SDRvars loadings', {},{'CCA 1', 'CCA 2','CCA 3'})
rotate3d on 

figure(2)
mdsplot3(T_rotw8(:,1:3),Nms(),[],[],Categ(),CategGrp(),'CCA rotated SDRvars loadings - Names & color by group', ...
         'SDR-CCA 3dim. rotated SDRvars loadings', {},{'CCA 1', 'CCA 2','CCA 3'})
rotate3d on    
     
figure
pcaplot(T_rotw8(:,1:3),Nms(),[],[],Categ(),CategGrp(),'CCA rotated SDRvsrs loadings - Names, color by group',... 
    'SDR-CCA 2dim. rotated SDRvars loadings with names','CCA3_1', 'CCA3_2')




%%
n=10; %top n variables to look into for each CCA mode



%%%%% Sphere area
BB1_rot=sqrt(B1_rot(:,1).^2+B1_rot(:,2).^2+B1_rot(:,3).^2);
[~,I0]=sort(BB1_rot,'descend');
EDtopw=B1_rot(I0(1:n),:);

%fid = fopen('../Header/SDRheader_dim58_plot.csv');
fid = fopen('../Header/SDRheader_dim60_plot.csv');
tmp = textscan(fid, '%s%d%s%s%s%f%f%f','Delimiter',',','Headerlines',1);fclose(fid);
[~,Ikeep]=ismember(header_sdr(I0(1:n)),tmp{1});
for i=1:size(tmp,2)
    temp=tmp{i};
    tmp{i}=temp(Ikeep,:);
end
Nms=tmp{1};
ShowLab=double(tmp{2});
Categ=[tmp{3},tmp{5}];
CategGrp=tmp{4};
Color=double([tmp{6}, tmp{7}, tmp{8}]);
Color=num2cell(Color);
Categ=[Categ, Color];

figure(1)
mdsplot3(EDtopw(:,1:3),{},[],[],Categ(),CategGrp(),'SDR CCA rotated loadings - Symbols & color by group', ...
         'Top 30 sphere SDRSDR-CCA 3dim. rotated loadings with symbols', {},{'CCA 1', 'CCA 2','CCA 3'})
rotate3d on 

figure(2)
mdsplot3(EDtopw(:,1:3),Nms(),[],[],Categ(),CategGrp(),'SDRvars CCA rotated loadings - Naemes & color by group', ...
         'Top 30 sphere SDRSDRvars-CCA 3dim. rotated loadings with names', {},{'CCA 1', 'CCA 2','CCA 3'})
rotate3d on    







% 
% 
% name=union(SM1_rotw81(1:n,1), SM1_rotw82(1:n,1));
% name1=union(name, SM1_rotw83(1:n,1));
% 
% [~,loc1]=ismember(name1,SM1_rotw81(:,1));
% [~,loc2]=ismember(name1,SM1_rotw82(:,1));
% [~,loc3]=ismember(name1,SM1_rotw83(:,1));
% 
% T_rotw8=[SMsdr_rotw81(loc1,2) SMsdr_rotw82(loc2,2) SMsdr_rotw83(loc3,2)];
% 
% T_rotw8=cell2mat(T_rotw8);
% 
% fid = fopen('../Header/SDRheader_dim58_plot.csv');
% tmp = textscan(fid, '%s%d%s%s%s%f%f%f','Delimiter',',','Headerlines',1);fclose(fid);
% [~,Ikeep]=ismember(name1',tmp{1});
% for i=1:size(tmp,2)
%     temp=tmp{i};
%     tmp{i}=temp(Ikeep,:);
% end
% Nms=tmp{1};
% ShowLab=double(tmp{2});
% Categ=[tmp{3},tmp{5}];
% CategGrp=tmp{4};
% Color=double([tmp{6}, tmp{7}, tmp{8}]);
% Color=num2cell(Color);
% Categ=[Categ, Color];
% 
% 
% figure(1)
% mdsplot3(T_rotw8(:,1:3),{},[],[],Categ(),CategGrp(),'CCA rotated SDR loadings - Symbols & color by group', ...
%          'SDR-CCA 3dim. rotated SDR loadings', {},{'CCA 1', 'CCA 2','CCA 3'})
% rotate3d on 
% 
% figure(2)
% mdsplot3(T_rotw8(:,1:3),Nms(),[],[],Categ(),CategGrp(),'CCA rotated SDR loadings - Names & color by group', ...
%          'SDR-CCA 3dim. rotated SDR loadings', {},{'CCA 1', 'CCA 2','CCA 3'})
% rotate3d on    
% 
% 
% figure
% pcaplot(T_rotw8(:,1:3),Nms(),[],[],Categ(),CategGrp(),'CCA rotated SDR loadings - Names, color by group',... 
%     'SDR-CCA 2dim. rotated SDR loadings with names')
% 


     
     
     
%% Plotting CCA loadings
header=readtable('Header_SDR.csv', 'ReadVariableNames',false);
header=table2cell(header);
figure(1)
subplot(3,1,1);
plot(grotB(:,1),'LineWidth',1.5)
line([8.5,8.5],ylim,'LineStyle',':') 
line([16.5,16.5],ylim,'LineStyle',':')
line([17.5,17.5],ylim,'LineStyle',':')
line([18.5,18.5],ylim,'LineStyle',':')
line([21.5,21.5],ylim,'LineStyle',':')
line([32.5,32.5],ylim,'LineStyle',':')
line([35.5,35.5],ylim,'LineStyle',':')
line([36.5,36.5],ylim,'LineStyle',':')
line([39.5,39.5],ylim,'LineStyle',':')
line([42.5,42.5],ylim,'LineStyle',':')
line([43.5,43.5],ylim,'LineStyle',':')
line([46.5,46.5],ylim,'LineStyle',':')
line([47.5,47.5],ylim,'LineStyle',':')
refline(0,0)
set(gca,'xtick',[1:55],'XTickLabel',[],'XLim',[0,56])

%set(gca,'xtick',[1,7,15,17,18,21,32,34,37,40,43,45,48,50,58])
h2=subplot(3,1,2);
plot(grotB(:,2),'LineWidth',1.5)
line([8.5,8.5],ylim,'LineStyle',':') 
line([16.5,16.5],ylim,'LineStyle',':')
line([17.5,17.5],ylim,'LineStyle',':')
line([18.5,18.5],ylim,'LineStyle',':')
line([21.5,21.5],ylim,'LineStyle',':')
line([32.5,32.5],ylim,'LineStyle',':')
line([35.5,35.5],ylim,'LineStyle',':')
line([36.5,36.5],ylim,'LineStyle',':')
line([39.5,39.5],ylim,'LineStyle',':')
line([42.5,42.5],ylim,'LineStyle',':')
line([43.5,43.5],ylim,'LineStyle',':')
line([46.5,46.5],ylim,'LineStyle',':')
line([47.5,47.5],ylim,'LineStyle',':')
refline(0,0)
set(gca,'xtick',[1:55],'XTickLabel',[],'XLim',[0,56])

subplot(3,1,3);
plot(grotB(:,3),'LineWidth',1.5)
line([8.5,8.5],ylim,'LineStyle',':') 
line([16.5,16.5],ylim,'LineStyle',':')
line([17.5,17.5],ylim,'LineStyle',':')
line([18.5,18.5],ylim,'LineStyle',':')
line([21.5,21.5],ylim,'LineStyle',':')
line([32.5,32.5],ylim,'LineStyle',':')
line([35.5,35.5],ylim,'LineStyle',':')
line([36.5,36.5],ylim,'LineStyle',':')
line([39.5,39.5],ylim,'LineStyle',':')
line([42.5,42.5],ylim,'LineStyle',':')
line([43.5,43.5],ylim,'LineStyle',':')
line([46.5,46.5],ylim,'LineStyle',':')
line([47.5,47.5],ylim,'LineStyle',':')
refline(0,0)
set(gca,'xtick',[1:55],'XTickLabel',header(:,2),'XTickLabelRotation',90,'XLim',[0,56])


% Factor rotation
figure(2)
RW=rotatefactors(grotB(:,1:3),'Method','promax');

header=readtable('Header_SDR.csv', 'ReadVariableNames',false);
header=table2cell(header);
subplot(3,1,1);
plot(RW(:,1),'LineWidth',1.5)
line([8.5,8.5],ylim,'LineStyle',':') 
line([16.5,16.5],ylim,'LineStyle',':')
line([17.5,17.5],ylim,'LineStyle',':')
line([18.5,18.5],ylim,'LineStyle',':')
line([21.5,21.5],ylim,'LineStyle',':')
line([32.5,32.5],ylim,'LineStyle',':')
line([35.5,35.5],ylim,'LineStyle',':')
line([36.5,36.5],ylim,'LineStyle',':')
line([39.5,39.5],ylim,'LineStyle',':')
line([42.5,42.5],ylim,'LineStyle',':')
line([43.5,43.5],ylim,'LineStyle',':')
line([46.5,46.5],ylim,'LineStyle',':')
line([47.5,47.5],ylim,'LineStyle',':')
refline(0,0)
set(gca,'xtick',[1:54],'XTickLabel',[],'XLim',[0,55])

%set(gca,'xtick',[1,7,15,17,18,21,32,34,37,40,43,45,48,50,58])
h2=subplot(3,1,2);
plot(RW(:,2),'LineWidth',1.5)
line([8.5,8.5],ylim,'LineStyle',':') 
line([16.5,16.5],ylim,'LineStyle',':')
line([17.5,17.5],ylim,'LineStyle',':')
line([18.5,18.5],ylim,'LineStyle',':')
line([21.5,21.5],ylim,'LineStyle',':')
line([32.5,32.5],ylim,'LineStyle',':')
line([35.5,35.5],ylim,'LineStyle',':')
line([36.5,36.5],ylim,'LineStyle',':')
line([39.5,39.5],ylim,'LineStyle',':')
line([42.5,42.5],ylim,'LineStyle',':')
line([43.5,43.5],ylim,'LineStyle',':')
line([46.5,46.5],ylim,'LineStyle',':')
line([47.5,47.5],ylim,'LineStyle',':')
refline(0,0)
set(gca,'xtick',[1:54],'XTickLabel',[],'XLim',[0,55])

subplot(3,1,3);
plot(RW(:,3),'LineWidth',1.5)
line([8.5,8.5],ylim,'LineStyle',':') 
line([16.5,16.5],ylim,'LineStyle',':')
line([17.5,17.5],ylim,'LineStyle',':')
line([18.5,18.5],ylim,'LineStyle',':')
line([21.5,21.5],ylim,'LineStyle',':')
line([32.5,32.5],ylim,'LineStyle',':')
line([35.5,35.5],ylim,'LineStyle',':')
line([36.5,36.5],ylim,'LineStyle',':')
line([39.5,39.5],ylim,'LineStyle',':')
line([42.5,42.5],ylim,'LineStyle',':')
line([43.5,43.5],ylim,'LineStyle',':')
line([46.5,46.5],ylim,'LineStyle',':')
line([47.5,47.5],ylim,'LineStyle',':')
refline(0,0)
set(gca,'xtick',[1:54],'XTickLabel',header(:,2),'XTickLabelRotation',90,'XLim',[0,55])

