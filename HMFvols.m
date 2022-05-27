% This code is for taking the preprocessed data and find the NAS using
% Homomorphic filtering and write in volumes

clc
clear all

mskp='/home/sukesh/fMRI/Data/Mask/BB_HealthyControlRestMask.img';
datap='/home/sukesh/fMRI/Data/BB_healthyNDisPre/HealthyControl4D';
sdir='/home/sukesh/fMRI/Data/BB_healthyNDisPre/HC_NAS_HMF_CR_3D/';


v=spm_vol(mskp);
brain=spm_read_vols(v);  % reading atlas
data_tmp = zeros(size(brain));
% m1=min(min(min(brain)));
m1=0;
voxel_ind = find(brain>=m1); %% if I am using greater than equal(zero) to means whole volume(row vector)
num_voxel = length(voxel_ind);

notp=152;%119; 
st=4;% How many scan will b discarded 
tp=notp-st;
% nor=116; % no of ROI
% tnotc=nor; % Total no of time courses
Band=[0 0.1]; % Frequency band(0.01- 0.08)
TR=2;%3;
% hrf1=spm_hrf(TR);
% fftl=115;
% Ql=4;
% Qh=115;
fftv=[148 152 158 162 168 172 182 192 256];%[400 432 452 482 492 512];%;
crng=[3 4 5 6 7 8 9 11 13 15];

subs = dir(datap); 
subs(1:2)=[]; %this removes the "." and ".." 
% subs(2:end)=[];
subs=subs(1:9);
nos=length(subs); %1;%length(sub); % No of subjects
DDD=zeros(nos*num_voxel,4);

for isub=1:1:nos
    disp('Reading data ...');
    sub=subs(isub).name
    file_path=fullfile(datap,sub,strcat(sub,'.nii'));
    %%% Data read %%%%%
    data=spm_read_vols(spm_vol(file_path));
    
    sd=size(data);
    rsig = zeros(tp,num_voxel); % % Every column of it is a BOLD time series
    for iscan=1:1:sd(4)-st   % Upto no of scan
        %%% Reading every scan %%%%
        data_scan=data(:,:,:,iscan+st);
        rsig(iscan,:)=data_scan(voxel_ind);
        clear data_scan;
    end
    clear data;
    hmf_sig=zeros(tp,num_voxel);
    for k=1:1:num_voxel
        tc1=rsig(:,k);
        tc=rest_IdealFilter(tc1, TR, Band);
        imf=emd(tc,'Display',0);
        if(size(imf,2)>=3)
%             imf=zeros(size(tc1,1),3);
           tc=imf(:,1)+imf(:,2)+imf(:,3);
        end  
        mtc=mean(tc);
        ons=tc;
        nosl=find(zscore(tc)<1); % Non onset locations
        ons(nosl)=0;
%         onsl1=find(ons~=0);
%         ons(onsl1)=1;
        if(sum(tc)==0)
            nv=tc;
            FCN=[length(tc) 1 3 0]; % If zero time series
        else
            [nv,FCN]=hmf_nv(tc,ons,fftv,crng); % ons for NAS NMCC % nv=estimated neuronal variable, FCN= FFT point, Cepstrum index, cepstrum point, NMCC
        end
% %         [res,del]=cceps(tc,fftl);
% % %         cep_hrf(1:Ql)=res(1:Ql);
% % %         cep_hrf(Ql:fftl)=0;
% %         cep_nv(1:Ql)=0;
% %         cep_nv(Ql:Qh)=res(Ql:Qh);
% %         cep_nv(Qh:fftl)=0;
% % %         es_hrf=single(icceps(cep_hrf,del));
% % %         es_hrf=es_hrf(1:sigl);
% %         es_nv=single(icceps(cep_nv,del));
% %         es_nv=es_nv(1:tp);
% %         ql=res(Ql);
% %         qh=res(Qh);
        hmf_sig(:,k)=zscore(nv)+mtc;%es_nv;
        DDD((isub-1)*num_voxel+1,:)=FCN;
        clear FCN;clear nv;
    end
%     v.dt=[16,0]; 
    sdir1=fullfile(sdir,sub);  % subject dir
    mkdir(sdir1);
    count=1;
    for j1=1:1:tp
       fp=fullfile(sdir1,strcat(sub,'_',num2str(count),'.nii'));
       vol1=hmf_sig(j1,:);
       v.fname=fp;
       sVol=data_tmp;
       sVol(voxel_ind)=vol1';
       spm_write_vol(v,sVol);
       count=count+1;
    end
    clear hmf_sig;
 end