clc
clear all


fltr=1;
Band=[0.01 0.1]; % Frequency band(0.01- 0.08)
TR=3;
% hrf1=spm_hrf(TR);

% den_sig=csvread('/home/sukesh/fMRI/HomomorphicFiltering/extracted_data/sub01361_roi12vox_nor_tc.csv');  % For Original normalised signal. Though later it is written as denoised
% den_sig=csvread('/home/sukesh/fMRI/HomomorphicFiltering/extracted_data/sub00156_nor_roi_64sig.csv');
den_sig=csvread('/home/sukesh/fMRI/HomomorphicFiltering/extracted_data/sub01361_nor_dmn_sig.csv');

% %%% HMF %%%%%%
roi=125;%125;%290;%30;%3;%15;%%23,5,118;
sig1=den_sig(:,roi);
if(fltr==1)
%     sig=rest_IdealFilter(sig1, TR, Band);
    imf=emd(sig1,'Display',0);
    if(size(imf,2)<3)
       imf=zeros(size(sig1,1),3);
    end
    sig=imf(:,1)+imf(:,2)+imf(:,3);
end

ons=sig;
sigl=length(sig);
ts=zscore(sig);  % Denoised time series(Spatio temporal regularization)
nosl=find(ts<1); % Non onset locations
ons(nosl)=0;
% onsl1=find(ons~=0);
% ons(onsl1)=1;
% ons=hrf1;

[nv,FCN]=hmf_nv(sig,ons);
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
sig1=(sig1-min(sig1))/max(sig1);
sig=(sig-min(sig))/max(sig);
   
subplot(5,1,1)
plot(sig1,'Color',[0.8500, 0.3250, 0.0980],'Linewidth',1.5)
text(length(sig1)+1,mean(sig1),'(a)','fontweight','bold','fontsize',12);
ylabel('$\bf{m[n]}$','Interpreter','latex','fontsize',14)
xlabel("TR (Sec)")
xlim([1 length(sig1)])
title(strcat('Preprocessed time course'))

subplot(5,1,2)
plot(sig,'Color',[0, 0.5, 0],'Linewidth',1.5)
text(length(sig)+1,mean(sig),'(b)','fontweight','bold','fontsize',12);
ylabel('$\bf{m_{den}[n]}$','Interpreter','latex','fontsize',14)
xlabel("TR (Sec)")
xlim([1 length(sig1)])
title('Band passed filtered time course')

subplot(5,1,3)
plot(res,'Color',[0 0 0],'Linewidth',1.5);
hold on
X=Qc:1:fftl;
Y=single(res(Qc:fftl,1));
plot(X,Y,'Color',[1 0 1],'Linewidth',1.5);
hold on
line([Qc Qc],[min(res)-0.3 qc],'Color',[0 0 1],'Linewidth',1.5)
title(strcat('Qc of time course', 'FFT length=',num2str(fftl),', Qc point=',num2str(Qc)))
text(length(Y)+4,double(mean(Y)),'(c)','fontweight','bold','fontsize',12);
ylabel('$\bf{\tilde{m}_{den}[n]}$','Interpreter','latex','fontsize',14)
xlabel('Quefrency(TR)','fontsize',12)
xlim([0 length(res)])
ylim([min(res)-0.3 max(res)+0.3])
hold off

subplot(5,1,4)
es_hrf=(es_hrf-min(es_hrf))/max(es_hrf);
plot(es_hrf,'Color',[0 0 0],'Linewidth',1.5)
text(length(es_hrf)+1,double(mean(es_hrf)+0.4),'(d)','fontweight','bold','fontsize',12);
ylabel('$\bf{\hat{h}[n]}$','Interpreter','latex','fontsize',14)
xlim([1 length(sig1)])
title('HRF-HMF')

subplot(5,1,5)
es_nv=(es_nv-min(es_nv))/max(es_nv);
plot(es_nv,'Color',[1 0 1],'Linewidth',1.5)
hold on
plot(ons,'Color',[0 0 1],'Linewidth',1.5)
text(length(es_nv)+1,double(mean(es_nv)),'(e)','fontweight','bold','fontsize',12);
ylabel('$\bf{\hat{s}[n]}$','Interpreter','latex','fontsize',14)
legend('Estimated NV','Probable NV');
title(strcat('NAS - HMF',', NMCC = ',num2str(nmcc)))
xlabel('TR (Sec)')
xlim([1 length(sig1)])

