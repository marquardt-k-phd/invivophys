
% set directories
clear all; clc

cd('E:\Documents\MATLAB')
addpath(genpath('E:\Documents\MATLAB\Cavanagh Class (Lec+Scripts)'));
addpath(genpath('E:\Documents\MATLAB\eeglab_addons'))


%% WAVELET FILTER

srate=1000; %start here and then downstample to 500Hz
numcycles=4; %lower = better temporal resolution; higher = better frequency resolution
%start_cycles=4;
%end_cycles=10; 
num_wavelets=80; %since going all the way to 100Hz this makes since -> this is what book 12.4 uses
start_freq=1; %epoch length should cover at least 3 cycle of the smallest frequency (ie: 1sec epoch=4Hz lowest frequency)
%an epoch from -1.5 to 7sec (correct) =8.5sec so can do 1Hz activity
end_freq=80;
%frex=linspace(start_freq,end_freq, num_wavelets);%linear increment from 2 to 100 for number of wavelets
frex=logspace(log10(start_freq),log10(end_freq),num_wavelets); %this is a scaling s which starts out time specific and beocmes more freq specific at higher frex
s=numcycles./(2*pi.*frex); %sigma for gaussian wave
%s=logspace(log10(start_cycles),log10(end_cycles),num_wavelets)./(2*pi.*frex);
t=-2.5:1/srate:2.499; %% time, from -2 to 2 second in steps of 1/sampling-rate
clear i w
for fq=1:num_wavelets
    w(fq,:) = exp(2*1i*pi*frex(fq).*t) .* exp(-t.^2./(2*s(fq)^2)); %sine by gaussian
    
      % time res: 2sigma_t *1000 for ms
    timeres(fq)=2*s(fq) * 1000;
    % freq res: 2sigma_f
    frexres(fq)=2* (frex(fq) / numcycles);
end

    tx=(-1000:1:3000)'; %total length of data
    t1=find(tx==-1000); %beginning of trial --> use whole trial to get power to graph
    t2=find(tx==3000); %end of trial
    tf_tx=-1000:1:3000;  %total trial time
    
    
%% Phase Synchrony per mouse  
Numses=5;
Numtrt=2;
NumTriTypes=2; 

ncrit=12;

%LFP.SYNCtcticMouse=cell(1,2);
warning('off')
for ses=1:5
    for trt=1:2
        Num=size(LFP.epochstctic,1);
        for mou=1:Num
             if size(var1{ses,mou,trt},2)==0 %done based on TC in DS, but OFC will follow same mouse pattern
                    continue %don't bother to do functions below if there is no data/ mouse
             end
             
   for ch=1:7
            
        for condi=1:4
            clear DATA        
            if condi==1,         data=squeeze(LFP.refepochstctic{mou,ses,1,2,trt}(:,ch,:));% OFC t1   (mou x ses x type x reg x trt)
            elseif condi==2 ,    data=squeeze(LFP.refepochstctic{mou,ses,1,1,trt}(:,ch,:));% DS t1
                
            elseif condi==3,     data=squeeze(LFP.refepochstctic{mou,ses,2,2,trt}(:,ch,:));% OFC t2
            elseif condi==4 ,    data=squeeze(LFP.refepochstctic{mou,ses,2,1,trt}(:,ch,:));% DS t2
            end
            

            if size(data,2) < ncrit %skip mouse if has lesss than crit num of trials
                if condi==1 || condi== 2
                PHASE_SYNC1=double.empty(80,4001,7,0); %correct
                %ANGavg{condi}=zeros(80,4001);  %correct
                elseif condi== 3|| condi ==4
                PHASE_SYNC2=double.empty(80,4001,7,0); %incorrects
                %ANGavg{condi}=zeros(80,4001); %incorrects
                end
                
                continue
            end
            
            for perms=1:250
                DATA=datasample(data,ncrit,2,'replace',false); %permutations - control for trial number
                dims=size(DATA);
                N=size(DATA,2);

                for fq=1:num_wavelets
                     w_con=fconv_JFC(reshape(DATA,1,dims(1)*dims(2)),w(fq,:));  %reshape data into one dimension -> all trials are lined up end to end (this makes it faster)
                     w_con=w_con((size(w,2)-1)/2+1:end-(size(w,2)-1)/2); %remove 1/2 of the wavelet length from each of the data b/c of edge effects
                     w_con=reshape(w_con,dims(1),dims(2)); %reform the data back into two dimensions

                   if condi==1 
                     SEEDANGLE1(fq,:,:)=angle(w_con(t1:t2,:)); %seedangle is from OFC
                   elseif condi==2
                     TARGANGLE1(fq,:,:)=angle(w_con(t1:t2,:)); %targangle is from DS
                     PHASE_SYNC1(fq,:,ch,perms)    =  abs( mean(  exp(1i*( squeeze(TARGANGLE1(fq,:,:)-SEEDANGLE1(fq,:,:)) ))  ,2)); % phase synch between seed & target
                     %PHASE_LAG1(fq,:,perms)= (1000 .* ( mean(squeeze(TARGANGLE1(fq,:,:)))-abs(squeeze(SEEDANGLE1(fq,:,:)),2)   ) ) ./ 2*pi*frex(fq);                                     
                   end
                   
                   if condi==3
                     SEEDANGLE2(fq,:,:)=angle(w_con(t1:t2,:));
                   elseif condi==4
                     TARGANGLE2(fq,:,:)=angle(w_con(t1:t2,:));
                     PHASE_SYNC2(fq,:,ch,perms)    =  abs(mean(exp(1i*( squeeze(TARGANGLE2(fq,:,:)-SEEDANGLE2(fq,:,:)) )),2)); % phase synch between seed & target
                     %PHASE_LAG2(fq,:,perms)= (1000 .* ( mean(squeeze(TARGANGLE2(fq,:,:,perms)))-abs(squeeze(SEEDANGLE2(fq,:,:,perms)),2)   ) ) ./ 2*pi*frex(fq);
                   end
                clear w_con
                end
%                 if condi==1
%                 TEMPSEEDAVG1(:,:,perms)=mean(SEEDANGLE1,3); %MEAN OVER TRIALS W/IN THE PERMUTATION
%                     if perms==250
%                      ANGavg{1}=mean(TEMPSEEDAVG1,3); clear TEMPSEEDAVG1%after last permutation --> mean x-perms
%                     end
%                 elseif condi==2
%                 TEMPTARGAVG1(:,:,perms)=mean(TARGANGLE1,3); %MEAN OVER TRIALS W/IN THE PERMUTATION
%                     if perms==250
%                      ANGavg{2}=mean(TEMPTARGAVG1,3); clear TEMPTARGAVG1 %after last permutation --> mean x-perms
%                     end
%                 elseif condi==3
%                 TEMPSEEDAVG2(:,:,perms)=mean(SEEDANGLE2,3); %MEAN OVER TRIALS W/IN THE PERMUTATION
%                     if perms==250
%                      ANGavg{3}=mean(TEMPSEEDAVG2,3); clear TEMPTARGAVG1 %after last permutation --> mean x-perms
%                     end
%                 elseif condi==4
%                 TEMPTARGAVG2(:,:,perms)=mean(TARGANGLE2,3); %MEAN OVER TRIALS W/IN THE PERMUTATION
%                     if perms==250
%                      ANGavg{4}=mean(TEMPTARGAVG2,3); clear TEMPTARGAVG1 %after last permutation --> mean x-perms
%                     end
%                 end
            end
        end
        
   end
        
      PHASE_SYNC(:,:,:,1)=mean(PHASE_SYNC1,4); %correct
      PHASE_SYNC(:,:,:,2)=mean(PHASE_SYNC2,4); %incorrects 
      clear PHASE_SYNC1 PHASE_SYNC2
      
    % PHASE_LAG(:,:,1)=(1000 .* (  bsxfun(@minus,ANGavg{2}, ANGavg{1})  )) ./ repmat((2*pi*frex)',1,size(ANGavg{2},2))  ;  %MEAN EACH SET OF ANGLES IN TIME AND PERM DIMENSIONS
     %PHASE_LAG(:,:,2)=(1000 .* (  bsxfun(@minus,ANGavg{4}, ANGavg{3})  )) ./ repmat((2*pi*frex)',1,size(ANGavg{3},2))  ; 
 
     SYNC{mou,ses,trt}=PHASE_SYNC; clear PHASE_SYNC* SEEDANGLE* TARGANGLE*
    % LAG{mou,ses,trt}=PHASE_LAG; clear PHASE_LAG* ANGavg 

    
        end
       
    end
end
LFP.refISPC=SYNC;
    save('DRBR_allchanLFP','LFP','-v7.3');  
warning('on')



LFPmou.lag2=LAG;
%% individual session graphs
CONT=zeros(80,4001);
  frexlow=1; %1Hz
  frexhigh=63; %31Hz
            
for ses=1:5
    for trt=1:1
        for beh=1:1
      CONT=zeros(80,4001);
       Num=size(LFP.epochs,1);
 for j=1:8
     if size(LFP.refISPC{j,ses,trt},2)==0 %if there is no data
        continue 
     end
     temp=LFP.refISPC{j,ses,trt}(:,:,:,beh); 
     temp(isnan(temp))=0;
       if sum(sum(temp,2),1)==0
           continue
       elseif sum(sum(temp,2),1)>0
           CONT=cat(3,CONT,squeeze(mean(temp,3)));
       end

 end
 
 CONT(:,:,1)=[];%clears first set of zeros
 gCONT=mean(CONT,3);
 

figure
hold on
contourf(tf_tx,frex(frexlow:frexhigh),gCONT(frexlow:frexhigh,:),40,'linecolor','none'); 
set(gca,'clim',[.24 .4],  'yscale', 'log','YTick', logspace(log10(frex(frexlow)),log10(frex(frexhigh)),6))
if ses==4 && trt==1
contour(tf_tx,frex(frexlow:frexhigh),LFP.syncROI05{beh}(frexlow:frexhigh,:),1,'linecolor','k','linewidth',2)
end
%title('L vs. R Hemisphere TC Phase Synchrony '); xlabel('Time (ms)');ylabel('Frequency (Hz)');
colorbar
colormap(jet)
title(['Session' num2str(ses) ' Behavior' num2str(beh) ' Treatment' num2str(trt)])
        end
    end
end
%  CONT(:,:,1)=[];%clears first set of zeros
%  gCONT=mean(CONT,3);
%  
% 
% figure
% hold on
% contourf(tf_tx,frex,gCONT,40,'linecolor','none'); 
% set(gca,'clim',[.24 .4],  'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6))
% %contour(tf_tx,frex,SIGBOUND{ses,1},1,'linecolor','k','linewidth',3)
% %title('L vs. R Hemisphere TC Phase Synchrony '); xlabel('Time (ms)');ylabel('Frequency (Hz)');
% colorbar

%% find ROI region of sig - must do for combined ROI across trts and regions (but not behaviors)

%define things
    Numses=5;
    Numtrts=2;
    Numbeh=2;
    nummice=size(LFP.epochs,1);

for B=1:Numbeh % have average within beh, but not between 
    
    for T=1:Numtrts
        for S=2:Numses
            for M=1:nummice
                if size(LFP.refISPC{M,S,T},2)==0
                 continue
                end
             temp= squeeze(mean(LFP.refISPC{M,S,T}(:,:,:,B),3));
             temp(isnan(temp))=0;
              
               if sum(sum(temp,2),1)==0
                   continue
               elseif sum(sum(temp,2),1)>0
                   roi(S,M,:,:)=temp; %*** must change dimensions for another measure
                  
               end
            end
        end 
        avg1(:,:,T)=squeeze(mean(mean(roi,1),2));
    end
    
    avg3(:,:,B)=mean(avg1,3); % time,freq x beh
end


 cluster_pval = 0.05;
 pval=0.05;
for beh=1:2

      DATA=avg3(:,:,beh);
     
         N=1;
          while N<5001
          ranvals=sign(randn(size(DATA)));
          randomdata=DATA(logical(ranvals==1));
            ranavg=mean(randomdata,1);
            ransd=2*std(randomdata,[],1);
            sig(N,1)=ranavg+ransd;
            N=N+1;
          end
          sig=abs(mean(sig,1)); clear N
        
          sigpix=logical(abs(DATA)>= sig);
      
      
%       
%       N=;          %sumBehTri(beh); %total trials per beh
% 
%       nullZ=atanh(DATA);
%       critval=sqrt( -log(pval) ./ N);
%       sigpix=logical(abs(DATA)>=critval); %find regions that are greater than or less than the critical value 

             %cluster correction for sigpix...
             clustinfo = bwconncomp(sigpix);
             clust_info = cellfun(@numel,clustinfo.PixelIdxList);
             
             cluster_pval = 0.05;
 
              mean_clust_info=mean(clust_info);
              med_clust_info=median(clust_info);
              clust_info(:,find(clust_info < med_clust_info))=[];
              
             clust_threshold=prctile(clust_info,100-cluster_pval*100);
                clust_info = cellfun(@numel,clustinfo.PixelIdxList);
             whichclusters2remove = find(clust_info<clust_threshold); %identify which clusters to remove
            
            % remove clusters
            for i=1:length(whichclusters2remove)
                sigpix(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
            end

        SIGPIX{beh}=sigpix; clear sigpix clust* CONT nullZ N EXP
    

end

LFP.syncROI05=SIGPIX;
%was finding TIC ROI, but was b/c all background and SD very low

%Test graphs

SIGPIX=LFPmou.syncROI05;
 figure
 hold on
   contourf(tf_tx,frex, avg3(:,:,1),40,'linecolor','none');
   contour(tf_tx,frex,SIGPIX{1},1,'linecolor','k','linewidth',3)
   %plot(tf_tx,2,'k','linewidth',3)
   %plot(tf_tx,4,'k','linewidth',3)
   %plot(tf_tx,7,'k','linewidth',3)
    %set(gca, 'clim', [.146 .18], 'yscale', 'log' ,'YTick', [],'XTick',[])
   set(gca, 'clim',[.2 .5],'yscale', 'log' ,'YTick', logspace(log10(start_freq),log10(end_freq),6),'XTick',-1000:1000:3000)%logspace(log10(start_freq),log10(end_freq),6))
   colormap(jet)   
   colorbar
      
 % Write ROI values to new excel file --> to transfer to statview for
% analysis

% measures={'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';...
%           'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'ph';...
%           
%           
%           'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';...
%           'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'ph'};
%        
% Treatment={'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';...
%            'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'ph';...
%            
%            'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';...
%            'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'ph'};
%        
% MouNum={'2.0';'3.1';'11.0';'16.1';'17.0';'17.1';'17.2';'20.1';'22.2';'25.0';'26.0';'27.0';...
%         '4.1';'14.1';'18.0';'19.0';'19.1';'21.9';'21.1';'22.1';'22.2';'24.0';'28.0';'ph';...
%         
%         '2.0';'3.1';'11.0';'16.1';'17.0';'17.1';'17.2';'20.1';'22.2';'25.0';'26.0';'27.0';...
%         '4.1';'14.1';'18.0';'19.0';'19.1';'21.9';'21.1';'22.1';'22.2';'24.0';'28.0';'ph'};
%     
   
    
measures={'ROI';'ROI';'ROI';'ROI';'ROI';'ROI';'ROI';'ph';...
          'ROI';'ROI';'ROI';'ROI';'ROI';'ROI';'ROI';'ROI'};
       
Treatment={'FLOX';'FLOX';'FLOX';'FLOX';'FLOX';'FLOX';'FLOX';'ph';...
           'NULL';'NULL';'NULL';'NULL';'NULL';'NULL';'NULL';'NULL'};
       
MouNum={'2.2';'3.2';'4.0';'4.1';'5.0';'7.1';'7.2';'ph';...
        '2.0';'3.0';'3.1';'5.1';'6.0';'6.1';'7.0';'8.0'};
       

% no TIC synchrony --- measuresTIC={'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';...
                                    %'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'ph'};

     clear control                               
for beh=1:1 %finding ROI x-y limits for both behaviors
        if beh==1, [row, col]=find(LFP.syncROI05{1}==1); %row freq, col are time
%            lowbound=find(row<=11);
%             xtcL=row(lowbound);
%             ytcL=col(lowbound); 
%            highbound=find(row>11);
%             xtcH=row(highbound);
%             ytcH=col(highbound);
%             clear row col 
%          elseif beh==2,[row, col]=find(SIGPIX{2}'==1); %row freq, col are time
           ytic=col; %freq
           xtic=row; %time
%            clear row col
        end

  %picking out ROI limits from phasedata    
   for trt=1:2
         for ses=1:5 
             mouserowstart=1;
             for mou=1:8
                 
        if size(LFP.refISPC{mou,ses,trt},2)==0 %if there is no data
            continue 
        end
  
        DATA =  squeeze(LFP.refISPC{mou,ses,trt}(:,:,:,beh)); DATA(isnan(DATA))=0;
       % DATAlag=LFPmou.lagExc{mou,ses,trt}(:,:,beh); DATAlag(isnan(DATAlag))=0;

        if sum(sum(DATA,2),1)==0 %IF ALL DATA IS NANS - also means no data
            control(mouserowstart:mouserowstart+6,ses,trt,beh)=zeros(7,1);
            mouserowstart=mouserowstart+7;
        elseif sum(sum(DATA,2),1)~=0
             control(mouserowstart:mouserowstart+6,ses,trt,beh)=mean(mean(DATA(ytc,xtc,:),2),1); %control, tc
            mouserowstart=mouserowstart+7;
        end
        
            clear DATA
             end
         end  

   end

end
clear D85 RS1 RS3 RS4 R50

        D85=[squeeze(control(:,1,1,1)); squeeze(control(:,1,2,1))]; 
        RS1=[squeeze(control(:,2,1,1)); squeeze(control(:,2,2,1))]; 
        RS3=[squeeze(control(:,3,1,1)); squeeze(control(:,3,2,1))]; 
        RS4=[squeeze(control(:,4,1,1)); squeeze(control(:,4,2,1))]; 
        R50=[squeeze(control(:,5,1,1)); squeeze(control(:,5,2,1))]; 
               %    sac  tc                    cre+     TC            
% 
%         D85=[squeeze(control(1,:,1))'; squeeze(experimental(1,:,1))'];% squeeze(control(2,:,1))'; squeeze(experimental(2,:,1))']; 
%         RS1=[squeeze(control(1,:,2))'; squeeze(experimental(1,:,2))']; %squeeze(control(2,:,2))'; squeeze(experimental(2,:,2))']; 
%         RS3=[squeeze(control(1,:,3))'; squeeze(experimental(1,:,3))']; %squeeze(control(2,:,3))'; squeeze(experimental(2,:,3))']; 
%         RS4=[squeeze(control(1,:,4))'; squeeze(experimental(1,:,4))']; %squeeze(control(2,:,4))'; squeeze(experimental(2,:,4))']; 
%         R50=[squeeze(control(1,:,5))'; squeeze(experimental(1,:,5))']; %squeeze(control(2,:,5))'; squeeze(experimental(2,:,5))']; 
%         %       sac TC LOW                        pae  TC LOW               SAC TC HIGH                  PAE TC HIGH
%          
%         
%         D85TIC=[squeeze(control(3,:,1))';  squeeze(experimental(3,:,1))'];
%         RS1TIC=[squeeze(control(3,:,2))';  squeeze(experimental(3,:,2))'];
%         RS3TIC=[squeeze(control(3,:,3))';  squeeze(experimental(3,:,3))'];
%         RS4TIC=[squeeze(control(3,:,4))';  squeeze(experimental(3,:,4))'];
%         R50TIC=[squeeze(control(3,:,5,1))';  squeeze(experimental(3,:,5))'];


T=table(measures,MouNum,Channel,Treatment,D85,RS1,RS3,RS4,R50); LFP.syncTC_chan=T;
%S=table(measuresTIC,MouNum,Treatment,D85TIC,RS1TIC,RS3TIC,RS4TIC,R50TIC); LFPmou.syncTIC=S;
    writetable(T,'PhaseSync_DRBR_ref_180117.xls','Sheet',1);    
    %writetable(S,'PhaseSync_ROIval.xls','Sheet',2);     


    
    writetable(LFP.cohROI_TC,'DRBR_PhaseMeasures.xls','Sheet',1); 
    writetable(LFP.cohROI_TIC, 'DRBR_PhaseMeasures.xls','Sheet',2);
    writetable(LFP.syncTC,  'DRBR_PhaseMeasures.xls','Sheet',3);

%% PHASE LAG HISTOGRAM X- FREQ

        
phlag=LFPmou.lag2;
SIGPIX=LFPmou.syncROI05;

CONT=zeros(80,4001);
for ses=5:5
    for trt=2:2
        for beh=1:1
 for j=1:12
     if size(phlag{j,ses,trt},1)==0 %if there is no data
        continue 
     end
     temp=phlag{j,ses,trt}(:,:,beh); 
     temp(isnan(temp))=0;
       if sum(sum(temp,2),1)==0
           continue
       elseif sum(sum(temp,2),1)~=0
           CONT=cat(3,CONT,temp);
       end
    
     %CONT=cat(3,CONT,temp); 
     
 end
        end
    end
end
 CONT(:,:,1)=[];%clears first set of zeros
 gCONT=mean(CONT,3);


 figure
 hold on
   contourf(tf_tx,frex, gCONT,40,'linecolor','none');
   contour(tf_tx,frex,SIGPIX{1},1,'linecolor','k','linewidth',3)
   plot(tf_tx,1.75,'k','linewidth',3)
    %set(gca, 'clim', [0 .4],  'yscale', 'log' ,'YTick', [],'XTick',[])
    set(gca, 'clim', [-80 80], 'yscale', 'log' ,'YTick', logspace(log10(start_freq),log10(end_freq),6),'XTick',-1000:1000:3000)%logspace(log10(start_freq),log10(end_freq),6))
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')
    colorbar
      

%histograph of lag at fq x-ROI
  %frequency bars for phase lag
       %1-1.5: 1-9
       %1.5-2: 10-13
       %2: 13-20
       %3: 21-25
       %4: 26-30
       %5: 31-33
       %6: 34-36
       
       freq=[1 1.5 2 3 4 5];
   
%1st - define fq and time row/ colums for ROI for each frequency
[row, col]=find(LFPmou.syncROI05{1}==1); %row freq, col are time
for range=1:5
    if range==1, ulimit=13; llimit=1;  %define range limits based on frex
   % elseif range==2,ulimit=13; llimit=10;
    elseif range==3,ulimit=20; llimit=14;
    elseif range==4,ulimit=25; llimit=21;
    elseif range==5,ulimit=30; llimit=26;
    elseif range==6,ulimit=33; llimit=31;
    %elseif range==7,ulimit=34; llimit=36; %freq = 6
    end
        
    avglag_fq{:,range}=row(find(row < ulimit & row > llimit));
    avglag_time{:,range}=col(find(row < ulimit & row > llimit));
end
for beh=1:1
for mou=1:12
    for ses=1:5
        for trt=1:2
    
        temp=phlag{mou,ses,trt}(:,:,beh);
        if size(temp,2)==0
            continue
        end
        
        for range=1:5
        avglag(mou,ses,trt,range)=mean(mean(temp(avglag_fq{:,range},avglag_time{:,range}),1),2);
        end
    
        end
    end
end
end

for ses=1:5
    freq=[1 2 3 4 5];
    for range=1:5 %average at each range across mice to graph
    tdata=avglag(:,ses,1,range);%for now only controls z=1;
        tdata(isnan(tdata))=[]; %clear nans
        tdatascatter(:,range)=tdata;
        tdata(tdata==0)=[]; %clear zeros
    data(1,range)=squeeze(mean(tdata,1));
    sem(1,range)=std(tdata,0,1) / sqrt(size(tdata,1)); clear tdata
    
    tdata_exp=avglag(:,ses,2,range);%for now only controls z=1;
        tdata_exp(isnan(tdata_exp))=[]; %clear nans
        tdatascatter_exp(:,range)=tdata_exp;
        tdata_exp(tdata_exp==0)=[]; %clear zeros
    data_exp(1,range)=squeeze(mean(tdata_exp,1)); 
    sem_exp(1,range)=std(tdata_exp,0,1) / sqrt(size(tdata_exp,1));clear tdata_exp
    end
    
figure
subplot(1,2,1)
hold on
bar(freq,data,'k')
%errorbar(freq,data,sem,'k')
ylim([-2 6])

subplot(1,2,2)
hold on
bar(freq,data_exp,'r')
%errorbar(freq,data_exp,sem_exp,'k')
ylim([-2 6])
%bar(freq,squeeze(mean(avglag(:,ses,2,1:6),1)),'r')
end

%question: are they sig different by session? --> if not collapse across
%session
clear data* tdata* tdatascatter* sem*

%phase lag for upper and lower fq @ 1sec (tone off) - 1000 to 1050ms
%(average over 50ms) - do this so not as variale
[row, col]=find(LFPmou.syncROI05{1}==1); %row freq, col are time
maxfreq=max(row); %find max frequency row #
lowfqrange=1:11; %1-1.75 Hz
highfqrange=12:maxfreq; %use this as upper range 1.75 to 5.58
time=2000:2050; %tone and 50ms following
for ses=1:5
    for trt=1:2
        for beh=1:1
        Num=size(phlag,1);
        for mou=1:Num

        temp=phlag{mou,ses,trt}(:,:,beh);
        if size(temp,2)==0
            continue
        end
        
        laginstL(mou,ses,trt)=mean(mean(temp(lowfqrange,time),2),1);
        laginstH(mou,ses,trt)=mean(mean(temp(highfqrange,time),2),1);

        end
        end
    end
end

%prelim graph
% for ses=4:5
%     
%     low=laginstL(:,ses,1);
%     low(low==0)=[];
%     lowsd=std(low,0,1);
%     lowsem=lowsd ./ sqrt(size(low,1));
%     lowmean=mean(low,1);
%         crit=lowmean + 2*lowsd;
%         outliers=find(abs(low)>crit);
%         low(outliers)=[];
%     
%     %lowexp=mean(laginstL(:,ses,1),1);
%     %lowexp(low==0)=[];
%     %lowexpsem=std(low,0,1) ./ sqrt(size(low,1));
%     
%     figure
%     boxplot(low)
%     
%     
% end

%make table for excel export --> do stats in R LMM (lmr4)
frequency={'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';...
           'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'PH';...
          
          
          'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';...
          'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'PH'};
       
Treatment={'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';...
           'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'ph';...
           
           'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';...
           'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'ph'};
       
MouNum={'2.0';'3.1';'11.0';'16.1';'17.0';'17.1';'17.2';'20.1';'22.2';'25.0';'26.0';'27.0';...
        '4.1';'14.1';'18.0';'19.0';'19.1';'21.9';'21.1';'22.1';'22.2';'24.0';'28.0';'ph';...
        
        '2.0';'3.1';'11.0';'16.1';'17.0';'17.1';'17.2';'20.1';'22.2';'25.0';'26.0';'27.0';...
        '4.1';'14.1';'18.0';'19.0';'19.1';'21.9';'21.1';'22.1';'22.2';'24.0';'28.0';'ph'};

        D85=[laginstL(:,1,1); laginstL(:,1,2);laginstH(:,1,1);laginstH(:,1,2)]; 
        RS1=[laginstL(:,2,1); laginstL(:,2,2);laginstH(:,2,1);laginstH(:,2,2)];
        RS3=[laginstL(:,3,1); laginstL(:,3,2);laginstH(:,3,1);laginstH(:,3,2)];
        RS4=[laginstL(:,4,1); laginstL(:,4,2);laginstH(:,4,1);laginstH(:,4,2)];
        R50=[laginstL(:,5,1); laginstL(:,5,2);laginstH(:,5,1);laginstH(:,5,2)];


PL=table(frequency,MouNum,Treatment,D85,RS1,RS3,RS4,R50); %LFPmou.syncTC=T;
%S=table(measuresTIC,MouNum,Treatment,D85TIC,RS1TIC,RS3TIC,RS4TIC,R50TIC); LFPmou.syncTIC=S;
    writetable(PL,'PhaseLag_Inst_2ROI.xls','Sheet',1);    


%% graphing w/ each subject displayed as point --SYNC

for fqr=1:1

       for ses=1:5
            
            tempdata=control(fqr,:,ses);
            tempdata(tempdata==0)=[];
              means_cont(1,ses)=mean(tempdata);
              errors_cont(1,ses)=std(tempdata,0,2)/ sqrt(size(tempdata,2)); clear tempdata
              
            tempdata=experimental(fqr,:,ses);
            tempdata(tempdata==0)=[];
              means_exp(1,ses)=mean(tempdata);
              errors_exp(1,ses)=std(tempdata,0,2)/ sqrt(size(tempdata,2)); clear tempdata

       end       

    figure %plot the individual points
%     hold on
%     %control
%     scatter(ones(length(control(fqr,:,1)),1),control(fqr,:,1),90,'k','LineWidth',1); %^ = triangle
%     scatter(repmat(2,length(control(fqr,:,2)),1),control(fqr,:,2),90,'k','LineWidth',1);
%     scatter(repmat(3,length(control(fqr,:,3)),1),control(fqr,:,3),90,'k','LineWidth',1);
%     scatter(repmat(4,length(control(fqr,:,4)),1),control(fqr,:,4),90,'k','LineWidth',1);
%     scatter(repmat(5,length(control(fqr,:,5)),1),control(fqr,:,5),90,'k','LineWidth',1);
%     
%     %experimental
%     scatter(ones(length(experimental(fqr,:,1)),1),experimental(fqr,:,1),90,'r','LineWidth',1);
%     scatter(repmat(2,length(experimental(fqr,:,2)),1),experimental(fqr,:,2),90,'r','LineWidth',1);
%     scatter(repmat(3,length(experimental(fqr,:,3)),1),experimental(fqr,:,3),90,'r','LineWidth',1);
%     scatter(repmat(4,length(experimental(fqr,:,4)),1),experimental(fqr,:,4),90,'r','LineWidth',1);
%     scatter(repmat(5,length(experimental(fqr,:,5)),1),experimental(fqr,:,5),90,'r','LineWidth',1);
%     xlim([0 6])
%     
%   if fqr==1 || fqr==2
%         ylim([.15 .45])
%      elseif fqr==2
%          ylim([0 8])
%      elseif fqr==3
         ylim([.1 .55])
%   end

    hold on %plot the average + standard error
    errorbar(means_cont,errors_cont,'-sk','MarkerFaceColor','k','markers',12)
    errorbar(means_exp,errors_exp,'-sr','MarkerFaceColor','r','markers',12)
    
    title(['Measure' num2str(fqr) ]);


end

























