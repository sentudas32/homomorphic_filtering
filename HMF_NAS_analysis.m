% This program is to extract time course from PCC, MPFC, rLPC and lLPC of
% DMN. Use Homomorphic filtering on Healthy rest

clc
clear all

filter=0;  % for BPF
emd1=1;
Band=[0.01 0.1]; % Frequency band(0.01- 0.08)
TR=2; %2.5 BH, 2 for rest, bh 2.5
isub=13;%5,6 
roi=3; %3,2  % For healthy rest 3 types o data 4,4,3,2,4,2,4 For DMN 4 roi% sub 2 roi 1,3%isub 13, 4,10,13,17,17,18,19
roival=[1 10 20 30];  %DMN= PCC-1,MPC-10,rLPC-20 lLPC-30
% roival=[1 2 3 4]; % For BB Breeath hold (BH) task, right frontal pole (RFP) - 1, LFP - 2, Right precentral gyrus (RPG) - 3, LPG-4

st=4;% How many scan will b discarded 
% tp=notp-st;
% nor=116; % no of ROI
TR=2;
hrf1=spm_hrf(TR);

fftv=[148 152 158 162 168 172 182 192 256];%[400 432 452 482 492 512];%;
% fftv=[250 256 280 286 290 296 300];%[400 432 452 482 492 512];
% crng=[3 4 5 6 7 8 9 11 13 15];
% fftv=[186 188 190 194 196]; %BB_BH 186, ft 236,bh 250 , rest 148;
crng=[3 4 5];
% crng=3;
% fftl=115;
% Ql=4;
% Qh=115;

maskp='/home/sukesh/fMRI/Data/Mask/BB_HealthyControl_DMNnAunDAMnCEN/BB_HealthyControl_DMN_1_10_20_30.nii'; % DMN with 4 nodes
datap='/home/sukesh/fMRI/Data/BB_healthyNDis/HC_Pre_4D';

% maskp='/home/sukesh/fMRI/Data/Mask/sub-A00057809_task-BH_mask_RFP1_LFP2_RPG3_LPG4.nii';
% datap='/home/sukesh/fMRI/Data/BB_BH/temp_10subj_BH/test';

% maskp='/home/sukesh/fMRI/Data/dataverse_files/test/CON3084_bh_mask_1_10_20_30.nii';
% datap='/home/sukesh/fMRI/Data/dataverse_files/test/4d';

% maskp='/home/sukesh/fMRI/Data/Mask/Vis_motor_sub_BB/pre/BB_1sub_per_rVis4.nii';
% datap='/home/sukesh/fMRI/Data/BB_event_Vis_Motor/5Sub_nifti/preprocessed/4D/per';


v=spm_vol(maskp);
brain=spm_read_vols(v);  % reading atlas


voxel_ind = find(brain==roival(roi)); %% if I am using greater than equal(zero) to means whole volume(row vector)
num_voxel = length(voxel_ind);

subs = dir(datap); 
subs(1:2)=[]; %this removes the "." and ".." 

sub=subs(isub).name

file_path=fullfile(datap,sub,strcat(sub,'.nii'));
%%% Data read %%%%%
data=spm_read_vols(spm_vol(file_path));

sd=size(data);
tp=sd(4)-st;
rsig = zeros(tp,num_voxel); % % Every column of it is a BOLD time series
for iscan=1:1:sd(4)-st   % Upto no of scan
    %%% Reading every scan %%%%
    data_scan=data(:,:,:,iscan+st);
    rsig(iscan,:)=data_scan(voxel_ind);
%     rsig(iscan,:)=data_scan(19,19,39);
    clear data_scan;
end
clear data;
sig1=mean(rsig,2);  % Taking avg of the voxels
if(filter==1)
    sig1=rest_IdealFilter(sig1,TR,Band);
end
if(emd1==1)
    imf=emd(sig1,'Display',0);
    if(size(imf,2)<3)
       imf=zeros(size(sig1,1),3);
    end
    sig=imf(:,1)+imf(:,2)+imf(:,3);
else
    sig=sig1;
end
    
ons=sig;
sigl=length(sig);
ts=zscore(sig);  % Denoised time series
nosl=find(ts<1); % Non onset locations
ons(nosl)=0;

[nv,FCN]=hmf_nv(sig,ons,fftv,crng);
fftl=FCN(1);
Qc=FCN(3);
nmcc=FCN(4);
[res,del1]=cceps(sig,fftl);
cep_hrf(1:Qc)=res(1:Qc);
cep_hrf(Qc:fftl)=0;
cep_nv(1:Qc)=0;
cep_nv(Qc:fftl)=res(Qc:fftl);
es_hrf=single(icceps(cep_hrf,del1));
es_hrf=es_hrf(1:sigl);
es_nv=single(icceps(cep_nv,del1));
es_nv=es_nv(1:sigl);
qc=res(Qc);


figure(1)
sig2=sig1-sig;
sig2=(sig2-min(sig2))/max(sig2);
sig1=(sig1-min(sig1))/max(sig1);
sig=(sig-min(sig))/max(sig);
ons=sig;
ons(nosl)=0;

   
subplot(3,2,1)
plot(sig1,'Color',[0.8500, 0.3250, 0.0980],'Linewidth',1.5)
text(length(sig1)+1,mean(sig1),'(a)','fontweight','bold','fontsize',11);
ylabel('$\bf{m[n]}$','Interpreter','latex','fontsize',12)
xlabel("TR (Sec)")
xlim([1 length(sig1)])
title(strcat('Preprocessed time course'))

% subplot(5,1,2)
% % plot(sig,'Color',[0, 0.5, 0],'Linewidth',1.5) % Original
% plot(sig2,'Color',[0, 0.5, 0],'Linewidth',1.5) % Original
% text(length(sig1)+1,mean(sig2),'(b)','fontweight','bold','fontsize',12);
% ylabel('$\bf{m_{r}[n]}$','Interpreter','latex','fontsize',14)
% xlabel("TR (Sec)")
% xlim([1 length(sig1)])
% % title('Band passed filtered time course') % Original
% title('Residual (Low frequency drift)')

subplot(3,2,3)
plot(res,'Color',[0 0 0],'Linewidth',1.5);
hold on
X=Qc:1:fftl;
Y=single(res(Qc:fftl,1));
plot(X,Y,'Color',[1 0 1],'Linewidth',1.5);
hold on
line([Qc Qc],[min(res)-0.3 qc],'Color',[0 0 1],'Linewidth',1.5)
title('Optimal cutoff quefrency')
% title(strcat('Qc of the time course', ' FFT length=',num2str(fftl),', Qc point=',num2str(Qc)))
text(length(Y)+4,double(mean(Y)),'(b)','fontweight','bold','fontsize',11);
ylabel('$\bf{\tilde{m}_{d}[n]}$','Interpreter','latex','fontsize',12)
xlabel('Quefrency','fontsize',12)
xlim([0 length(res)])
ylim([min(res)-0.3 max(res)+0.3])
hold off

% % subplot(5,1,3)
% % % es_nv=nv;
% % es_nv=(es_nv-min(es_nv))/max(es_nv);
% % plot(es_nv,'Color',[1 0 1],'Linewidth',1.5)
% % % hold on
% % % plot(ons,'Color',[0 0 1],'Linewidth',1.5)
% % text(length(es_nv)+1,double(mean(es_nv)),'(c)','fontweight','bold','fontsize',12);
% % ylabel('$\bf{\hat{s}[n]}$','Interpreter','latex','fontsize',14)
% % legend('Estimated NV','Probable NV');
% % title(strcat('NAS - HMF',', MNCC = ',num2str(nmcc)))
% % xlabel('TR')
% % xlim([1 length(sig1)])


subplot(3,2,5)
es_hrf=(es_hrf-min(es_hrf))/max(es_hrf);
plot(es_hrf,'Color',[0 0 0],'Linewidth',1.5)
text(length(es_hrf)+1,double(mean(es_hrf)+0.4),'(c)','fontweight','bold','fontsize',11);
ylabel('$\bf{\hat{h}[n]}$','Interpreter','latex','fontsize',12)
xlim([1 length(sig1)])
% xlim([1 20])
title('HRF-HMF')
xlabel('TR (Sec)')

subplot(3,2,4)
% es_hrf=(es_hrf-min(es_hrf))/max(es_hrf);
plot(es_hrf,'Color',[0 0 0],'Linewidth',1.5)
hold on
hrf2=(hrf1-min(hrf1))/max(hrf1);
plot(hrf2,'Color','b','Linewidth',1.5)
text(length(es_hrf(1:30))+0.2,double(mean(es_hrf)+0.1),'(d)','fontweight','bold','fontsize',11);
ylabel('$\bf{\hat{h}[n]}$','Interpreter','latex','fontsize',12)
% xlim([1 length(sig1)])
xlim([1 30])
ylim([0 2])
title('Estimated HRF(Truncated)')
legend('Non-parametric HRF','Canonical HRF')
xlabel('TR (Sec)')


