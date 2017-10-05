clear
clc

vars0=importdata('MF_psy_vars.mat');
NET=importdata('NET.mat');
%BMvars=importdata('SDR_NET1_PCnew.mat');
SDRvars=importdata('/Users/jessie_liu/Desktop/HCP900/Matlab/Data/SDR_CCAvars_preDeconf.mat');

%Prepare the permutation set
Nperm=10000;                                                                       % in the paper we used 100000 but 10000 should be enough
EB=hcp2blocks('../RESTRICTED_lzdh_3_20_2016_11_24_35.csv', [ ], false, vars0(:,1)); % change the filename to your version of the restricted file
PAPset=palm_quickperms([ ], EB, Nperm);                                            % the final matrix of permuations

Nkeep=100;

%%% setup confounds matrix
%conf=palm_inormal([vars0(:,[2 3 4 7 14 15 22 23 25]) vars0(:,[265 266]).^(1/3) ]);    % Gaussianise
conf=palm_inormal([vars0(:,[2 3 4 7 14 15 22 23 25 479 480 481 486]) vars0(:,[265 266]).^(1/3) ]);    % Gaussianise
conf(isnan(conf))=0;  % impute missing data as zeros
conf=nets_normalise([conf conf(:,5:end).^2]);  % add on squared terms and renormalise

% Remove bad variables
badvars=[];
for i=1:size(SDRvars,2)
  Y=SDRvars(:,i); grotKEEP=~isnan(Y);  
  %change the creterion when # of subjects changes
  if (sum(grotKEEP)>400) & (max(sum(nets_class_vectomat(Y(grotKEEP))))/length(Y(grotKEEP))<0.95)
      % the 3rd thing above is:  is the size of the largest equal-values-group too large?
    i=i; % do nothing
  else
    [i sum(grotKEEP)  max(sum(nets_class_vectomat(Y(grotKEEP))))/length(Y(grotKEEP))]
    badvars=[badvars i];
  end 
end

varskeep=setdiff([1:size(SDRvars,2)], badvars);  
SMvars=SDRvars(:,varskeep);

% Normalise and deconfound the original SM
varsgrot=palm_inormal(SMvars);
for i=1:size(varsgrot,2)
  grot=(isnan(varsgrot(:,i))==0); grotconf=nets_demean(conf(grot,:)); varsgrot(grot,i)=nets_normalise(varsgrot(grot,i)-grotconf*(pinv(grotconf)*varsgrot(grot,i)));
end
grot=varsgrot; grotI=double(~isnan(grot)); grot(isnan(grot))=0;
varsdCOV = grot*grot';
varsdCOV=nearestSPD(varsdCOV); % minor adjustment: project onto the nearest valid covariance matrix
[uu2,dd]=eigs(varsdCOV,100);  % SVD (eigs actually)
SMpc=uu2*sqrt(dd)-conf*(pinv(conf)*(uu2*sqrt(dd)));    % deconfound again just to be safe

%%% prepare main netmat matrix - we have a range of normalisation possibilities
NET1=nets_demean(NET);  NET1=NET1/std(NET1(:)); % no norm
grot=NET1; % concat horizontally
NETd=nets_demean(grot-conf*(pinv(conf)*grot));   % deconfound and demean
[uu,ss,vv]=nets_svds(NETd,100); % SVD reduction
BMpc=uu*ss;

% BMvars=nets_demean(BMvars-conf*(pinv(conf)*BMvars));   

% %Apply PCA to further reduce the dimension of SDR BM
% [uu,ss,vv]=nets_svds(BMvars,Nkeep); % SVD reduction
% BMpc=uu*ss;
 
%%% CCA/
[grotA,grotB,grotR,grotU,grotV,grotstats]=canoncorr(BMpc,SMpc);


%%% CCA permutation testing
grotRp=zeros(Nperm,Nkeep); clear grotRpval;
for j=1:Nperm
  j
  [grotAr,grotBr,grotRp(j,:),grotUr,grotVr,grotstatsr]=canoncorr(uu1,uu2(PAPset(:,j),:));
end
for i=1:Nkeep;  % get FWE-corrected pvalues
  grotRpval(i)=(1+sum(grotRp(2:end,1)>=grotR(i)))/Nperm;
end
grotRpval
Ncca=sum(grotRpval<0.05)  % number of FWE-significant CCA components

%%% netmat weights for CCA mode 1
grotAA = corr(grotU(:,1),NET)';
 % or
grotAAd = corr(grotU(:,1),NETd(:,1:size(NET,2)))'; % weights after deconfounding

%%% SM weights for CCA mode 1
grotBB = corr(grotV(:,1),palm_inormal(vars),'rows','pairwise');
 % or 
varsgrot=palm_inormal(vars);
for i=1:size(varsgrot,2)
  grot=(isnan(varsgrot(:,i))==0); grotconf=nets_demean(conf(grot,:)); varsgrot(grot,i)=nets_normalise(varsgrot(grot,i)-grotconf*(pinv(grotconf)*varsgrot(grot,i)));
end
grotBBd = corr(grotV(:,1),varsgrot,'rows','pairwise')'; % weights after deconfounding


%%% amount of variance between explained by canonical variates in original variables
for j=1:5
    TVE_SM(j,:)=corr(grotV(:,j),varsgrot,'rows','pairwise').^2;%varsgrot
    TVE_Net(j,:)=corr(grotU(:,j),NETd).^2; %NETd(:,1:size(NET,2))
end
AVE_SM=mean(TVE_SM,2,'omitnan');
AVE_Net=mean(TVE_Net,2,'omitnan');

for i=1:length(AVE_SM)
    cum(i)=sum(AVE_SM(1:i)');  
end


varsgrot(isnan(varsgrot))=0;
% Alternatively
b1=grotV(:,1:3)\varsgrot;
yhat=grotV(:,1:3)*b1;
Rsq = 1 - sum(sum((varsgrot - yhat).^2))/sum(sum((varsgrot).^2))

b2=grotU(:,1:3)\NETd;
yhat2=grotU(:,1:3)*b2;
Rsq2 = 1 - sum((NETd - yhat2).^2)/sum((NETd).^2)

%% Apply PCA method to 213 SDR variables

clear
clc

vars0=importdata('../Data/MF_psy_vars.mat');
NET=importdata('../Data/NET.mat');
%BMvars=importdata('SDR_NET1_PCnew.mat');
SDRvars=importdata('../Data/SDR_CCAvars_3.mat');
headers=importdata('../Header/MFvarnames.txt');
SDR_headers=importdata('../Header/SDRCCA_header_withsigns.mat');

%Prepare the permutation set
Nperm=10000;                                                                       % in the paper we used 100000 but 10000 should be enough
EB=hcp2blocks('../RESTRICTED_lzdh_3_20_2016_11_24_35.csv', [ ], false, vars0(:,1)); % change the filename to your version of the restricted file
PAPset=palm_quickperms([ ], EB, Nperm);                                            % the final matrix of permuations

Nkeep=100;

%%% setup confounds matrix
%conf=palm_inormal([vars0(:,[2 3 4 7 14 15 22 23 25]) vars0(:,[265 266]).^(1/3) ]);    % Gaussianise
conf=palm_inormal([vars0(:,[2 3 4 7 14 15 22 23 25 479 480 481 486]) vars0(:,[265 266]).^(1/3) ]);    % Gaussianise
conf(isnan(conf))=0;  % impute missing data as zeros
conf=nets_normalise([conf conf(:,5:end).^2]);  % add on squared terms and renormalise

% Normalise and deconfound the original SM
varsgrot=palm_inormal(SDRvars);
for i=1:size(varsgrot,2)
  grot=(isnan(varsgrot(:,i))==0); grotconf=nets_demean(conf(grot,:)); varsgrot(grot,i)=nets_normalise(varsgrot(grot,i)-grotconf*(pinv(grotconf)*varsgrot(grot,i)));
end
grot=varsgrot; grotI=double(~isnan(grot)); grot(isnan(grot))=0;
varsdCOV = grot*grot';
varsdCOV=nearestSPD(varsdCOV); % minor adjustment: project onto the nearest valid covariance matrix
[uu2,dd]=eigs(varsdCOV,Nkeep);  % SVD (eigs actually)
SMpc=uu2*sqrt(dd)-conf*(pinv(conf)*(uu2*sqrt(dd)));    % deconfound again just to be safe

%%% prepare main netmat matrix - we have a range of normalisation possibilities
NET1=nets_demean(NET);  NET1=NET1/std(NET1(:)); % no norm
grot=NET1; % concat horizontally
NETd=nets_demean(grot-conf*(pinv(conf)*grot));   % deconfound and demean
[uu,ss,vv]=nets_svds(NETd,Nkeep); % SVD reduction
BMpc=uu*ss;

%%% CCA/
[grotA,grotB,grotR,grotU,grotV,grotstats]=canoncorr(BMpc,SMpc);

%%% CCA permutation testing
grotRp=zeros(Nperm,Nkeep); clear grotRpval;
for j=1:Nperm
  j
  [grotAr,grotBr,grotRp(j,:),grotUr,grotVr,grotstatsr]=canoncorr(BMpc,SMpc(PAPset(:,j),:));
end
for i=1:Nkeep;  % get FWE-corrected pvalues
  grotRpval(i)=(1+sum(grotRp(2:end,1)>=grotR(i)))/Nperm;
end
grotRpval
Ncca=sum(grotRpval<0.05)  % number of FWE-significant CCA components

%%% BM and SM weights for sig. CCA modes
varsgrot1=palm_inormal(vars0);
for i=1:size(varsgrot1,2)
  grot=(isnan(varsgrot1(:,i))==0); grotconf=nets_demean(conf(grot,:)); varsgrot1(grot,i)=nets_normalise(varsgrot1(grot,i)-grotconf*(pinv(grotconf)*varsgrot1(grot,i)));
end

for i=1:Ncca

    grotAA(:,i) = corr(grotU(:,i),NET)';
     % or
    grotAAd(:,i) = corr(grotU(:,i),NETd)'; % weights after deconfounding    
    
    %SM
    grotBB_sdr(:,i) = corr(grotV(:,i),palm_inormal(SDRvars),'rows','pairwise')';
    grotBB(:,i) = corr(grotV(:,i),palm_inormal(vars0),'rows','pairwise')';

    grotBBd_sdr(:,i) = corr(grotV(:,i),varsgrot,'rows','pairwise')'; % weights after deconfounding 
    grotBBd(:,i) = corr(grotV(:,i),varsgrot1,'rows','pairwise')'; % weights after deconfounding
end

[~,I]=sort(abs(grotBBd_sdr(:,1)),'descend');
SMw8_sdr1=[SDR_headers(I), num2cell(grotBBd_sdr(I,1))];
[~,I]=sort(abs(grotBBd_sdr(:,2)),'descend');
SMw8_sdr2=[SDR_headers(I), num2cell(grotBBd_sdr(I,2))];
[~,I]=sort(abs(grotBBd_sdr(:,3)),'descend');
SMw8_sdr3=[SDR_headers(I), num2cell(grotBBd_sdr(I,3))];

[~,I1]=sort(abs(grotBBd(:,i)),'descend');
SMw81=[headers(I1), num2cell(grotBBd(I1,1))];
[~,I1]=sort(abs(grotBBd(:,i)),'descend');
SMw82=[headers(I1), num2cell(grotBBd(I1,2))];
[~,I1]=sort(abs(grotBBd(:,i)),'descend');
SMw83=[headers(I1), num2cell(grotBBd(I1,3))];



%%%Rotate sig. canonical variables
% grotAAd(isnan(grotAAd))=0;
% A_rot = rotatefactors(grotAAd);
grotBBd(isnan(grotBBd))=0;

for i=1:487
    if grotBBd 
        
B_rot = rotatefactors(grotBBd);
Bsdr_rot = rotatefactors(grotBBd_sdr);

[~,I]=sort(abs(B_rot(:,1),1),'descend');
[~,I1]=sort(abs(Bsdr_rot(:,1),1),'descend');
SM_rotw81=[headers(I), num2cell(B_rot(I,1))];
SMsdr_rotw81=[headers(I1), num2cell(Bsdr_rot(I1,1))];

[~,I]=sort(abs(B_rot(:,2),2),'descend');
[~,I1]=sort(abs(Bsdr_rot(:,2),2),'descend');
SM_rotw82=[headers(I), num2cell(B_rot(I,2))];
SMsdr_rotw82=[headers(I1), num2cell(Bsdr_rot(I1,2))];

[~,I]=sort(abs(B_rot(:,3),3),'descend');
[~,I1]=sort(abs(Bsdr_rot(:,3),3),'descend');
SM_rotw83=[headers(I), num2cell(B_rot(I,3))];
SMsdr_rotw83=[headers(I1), num2cell(Bsdr_rot(I1,3))];


%%%Variance Explained
for j=1:Ncca
    TVE_SM(j,:)=corr(grotV(:,j),varsgrot,'rows','pairwise').^2; %in & in
    TVE_Net(j,:)=corr(grotU(:,j),NETd).^2; %NETd(index1,1:size(NET,2))
  
    %%%Redundency index
%         RI_SM(j,:)=corr(grotU(:,j),varsgrot(index1,:),'rows','pairwise').^2;
%         RI_Net(j,:)=corr(grotV(:,j),NETd(index1,1:size(NET,2))).^2;
    
     CanVar=(P_in(:,j)+Q_in(:,j))./2;
     RI_SM(j,:)=corr(CanVar,varsgrot,'rows','pairwise').^2;
     RI_Net(j,:)=corr(CanVar,NETd).^2;

end

 AVE_SM(k,:)=mean(TVE_SM,2,'omitnan');
 AVE_Net(k,:)=mean(TVE_Net,2,'omitnan');

