%% Compute Phase Consistency + Synchrony  w/in mouse
%exclude any mouse with < 20 trials
%if less than three mice with >20 trials for any session -> cannot compute
%session
load DRBR_allchanLFP

% set directories
clear all; clc

cd('C:\Users\Kristen\Documents\MATLAB')
addpath(genpath('C:\Users\Public\Desktop\Plexon Inc Documents'));
addpath(genpath('E:\Documents\MATLAB\Cavanagh Class (Lec+Scripts)'));
addpath(genpath('E:\MATLAB\eeglab_addons'));
addpath(genpath('E:\Documents\MATLAB\eeglab12_0_2_1b'));
addpath(genpath('E:\Documents\MATLAB\kakearney-boundedline-pkg-2112a2b'));

addpath('E:\Documents\MATLAB\DATA_PAE');

addpath(genpath('/Volumes/Kris_SeaDri/MATLAB'))

%% Per region reference
%16 for OFC
%1 for DS
for ses=1:5
    for trt=1:2
  Nummice=size(var1,2);
for mou=1:Nummice
    for reg=1:2
        if reg==1
            ref=LFP.allFP{mou,1,ses,trt}; %ds reference
            chann=2:8;
        elseif reg==2
            ref=LFP.allFP{mou,16,ses,trt};
            chann=9:15;
        end
       
       referencedLFP{mou,ses,reg,trt}=cell2mat(LFP.allFP(mou,chann,ses,trt))- ref;

    end
end
    end
end

LFP.refLFP=referencedLFP; %mou x ses x reg x trt

%% SORTING INTO 4 TRIAL TYPE EPOCHS - ALL CHANNELS DUAL 
Numses=5;
Numtrt=2;
% (1-8=DS) (9-16=OFC) 
for ses=1:Numses
    for trt=1:2
        for reg=1:2
%             if reg==1
%                Chann=1:8;
%             elseif reg==2
%                 Chann=9:16;
%             end
    Numfiles=size(var1,2);
for j=1:Numfiles
        if size(var1{ses,j,trt},1)==0
           continue
        end
        

    %set stuff for correct trials
        tone=cell2mat(LFP.event(j,1,ses,trt));
            tc=bsxfun(@minus,tone,1);
            if tone==-1
                tc=[];
            end


    %set stuff for incorrect trials    
        punish=cell2mat(LFP.event(j,2,ses,trt));
            tic=bsxfun(@minus,punish,10);
            if punish==-1 %this is for if there are no incorrects -> events will be -1
               tic=[]; 
            end
          

         clear  punish tone
            
     %set up count of trials -> start of defining trial types
     trials=union(tc,tic);
        trialtypes=zeros(size(trials,1),3);
        trialtypes(:,1)=[1:numel(trials)]';%numbers the trials in the order they occured
        for t=1:size(trials,1)
            if size(find(trials(t,1)==tc)) >0
                trialtypes(t,2)=1; %corrects
            elseif size(find(trials(t,1)==tic)) >0
                trialtypes(t,2)=2; %incorrects
            end
        end
        
        %define trial types
         for t=2:size(trials,1)
                k=t-1;
                if trialtypes(k,2) == 1 && trialtypes(t,2) == 1
                    trialtypes(t,3)=1;      % correct -> correct
                elseif trialtypes(k,2) == 1 && trialtypes(t,2) == 2
                    trialtypes(t,3)=2;      %correct -> incorrect
                elseif trialtypes(k,2) == 2 && trialtypes(t,2) == 2
                    trialtypes(t,3)=3;      %incorrect -> incorrect
                elseif trialtypes(k,2) == 2 && trialtypes(t,2) == 1
                    trialtypes(t,3)=4;       %incorrect -> correct
                end
         end

        %find lfp times & create trial epochs 
          lfp=cell2mat(LFP.refLFP(j,ses,reg,trt)); %note chanw w/ referencing
          %lfp=cell2mat(LFP.allFP(j,Chann,ses,trt)); 
          int=cell2mat(LFP.interval(j,ses,trt));

           for b=1:size(trialtypes,1) %if everything is perfect
               temp=find(int <= trials(b,1)+.001 & int >= trials(b,1)-.001); %find where the LFP is closest to the spike in time
                    if trials(b,1) > max(int) %some trials are beyond max lfp timepoint
                        templfp(b,:)=zeros(1,4001);
                    elseif trials(b,1)<= max(int) %for those that are not
                        startrow=temp-1000; %1sec prior
                        endrow=temp+3000; %3 sec post
                        templfp(:,:,b)=lfp(startrow:endrow,:); 
                    end
               clear temp startrow endrow
           end
          
        %clear out zeros -> in both LFP and trial type arrays
            tempzero=find(templfp(:,all(templfp==0,1)));
            trialtypes(tempzero,:)=[];
            templfp(tempzero,:)=[];
        
        %split bins into trial types % store     
        LFPEPOCHS{j,ses,1,reg,trt}=templfp(:,:,logical(trialtypes(:,3)==1)); %lfp bins that are type 1 trials 
        LFPEPOCHS{j,ses,2,reg,trt}=templfp(:,:,logical(trialtypes(:,3)==2));
        LFPEPOCHS{j,ses,3,reg,trt}=templfp(:,:,logical(trialtypes(:,3)==3));
        LFPEPOCHS{j,ses,4,reg,trt}=templfp(:,:,logical(trialtypes(:,3)==4));
             

   clear tone tc punish tic trial* type* x y z w lfp int temp*
end
    clear Num
        end
    end
end
LFP.epochs=LFPEPOCHS; 
LFP.refepochs=LFPEPOCHS; 

clear LFPEPOCHS
% SORTING INTO TC and TIC TRIAL TYPE EPOCHS - ALL CHANNELS DUAL 
for ses=1:Numses
    for trt=1:2
        for reg=1:2
%             if reg==1
%                Chann=1:8;
%             elseif reg==2
%                 Chann=9:16;
%             end
    Numfiles=size(var1,2);
for j=1:Numfiles
        if size(var1{ses,j,trt},1)==0
           continue
        end
        
    %set stuff for correct trials
        tone=cell2mat(LFP.event(j,1,ses,trt));
            tc=bsxfun(@minus,tone,1);
            if tone==-1
                tc=[];
            end


    %set stuff for incorrect trials    
        punish=cell2mat(LFP.event(j,2,ses,trt));
            tic=bsxfun(@minus,punish,10);      
            if punish==-1 %if events are -1, means there were none of this trial 
                tic=[];
            end   

         clear  punish tone
            
     %set up count of trials -> start of defining trial types
     trials=union(tc,tic);
        trialtypes=zeros(size(trials,1),3);
        trialtypes(:,1)=[1:numel(trials)]';%numbers the trials in the order they occured
        for t=1:size(trials,1)
            if size(find(trials(t,1)==tc)) >0
                trialtypes(t,2)=1; %corrects
            elseif size(find(trials(t,1)==tic)) >0
                trialtypes(t,2)=2; %incorrects
            end
        end
        

        %find lfp times & create trial epochs 
          lfp=cell2mat(LFP.refLFP(j,ses,reg,trt)); %note chanw w/ referencing
          %lfp=cell2mat(LFP.allFP(j,Chann,ses,trt));
          int=cell2mat(LFP.interval(j,ses,trt));

           for b=1:size(trialtypes,1) %if everything is perfect
               temp=find(int <= trials(b,1)+.001 & int >= trials(b,1)-.001); %find where the LFP is closest to the spike in time
                    if trials(b,1) > max(int) %some trials are beyond max lfp timepoint
                        templfp(b,:)=zeros(1,4001);
                    elseif trials(b,1)<= max(int) %for those that are not
                        startrow=temp-1000; %1sec prior
                        endrow=temp+3000; %3 sec post
                        templfp(:,:,b)=lfp(startrow:endrow,:);
                        tempint(:,:,b)=int(startrow:endrow,:);
                    end
               clear temp startrow endrow
           end
          
        %clear out zeros -> in both LFP and trial type arrays
            tempzero=find(templfp(:,all(templfp==0,1)));
            trialtypes(tempzero,:)=[];
            templfp(tempzero,:)=[];
            tempint(tempzero,:)=[];
        
        %split bins into trial types % store     
        LFPEPOCHS{j,ses,1,reg,trt}=templfp(:,:,logical(trialtypes(:,2)==1)); %lfp bins that are correct trials  
        LFPEPOCHS{j,ses,2,reg,trt}=templfp(:,:,logical(trialtypes(:,2)==2)); %lfp bins that are incorrect trials  
        LFPInt{j,ses,1,trt}=tempint(:,:,logical(trialtypes(:,2)==1));
        LFPInt{j,ses,2,trt}=tempint(:,:,logical(trialtypes(:,2)==2));
             

   clear tone tc punish tic trial* type* x y z w lfp int temp*
end
    clear Num
        end
    end
end
 LFP.epochstctic=LFPEPOCHS; 
  LFP.refepochstctic=LFPEPOCHS; 
 
 save('DRBR_allchanLFP','LFP','-v7.3');  



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
    b1=find(tx==-1000); %beginning of baseline
    b2=find(tx==0); %end of baseline
    t1=find(tx==-1000); %beginning of trial--> at very start so get whole timespan in power graphs
    t2=find(tx==3000); %end of trial
    tf_tx=-1000:1:3000;  %total trial time

  
    
%% Baseline - excluding mice w/ <ncrit trials
Numses=5;
Numtrt=2;
Numreg=2;
Numbeh=2;

ncrit=12;

% Combine epochs across mice - only with >ncrit trials
for ses=1:Numses
    for trt=1:Numtrt
        for reg=1:Numreg

for type=1:Numbeh
       Numfiles=size(var1,2);
       initmou=zeros(1,4001);
        
       for j=1:Numfiles
            if size(var1{ses,j,trt},1)==0 %skip if no data for mouse
                continue
            end
           temp=LFPmou.epochs{j,ses,type,reg,trt};
           if size(temp,1)<ncrit %skip mouse if less than n-crit trials
               continue
           end
           
           initmou=[initmou; temp]; % add all bins from each mouse together
       end
           initmou(all(initmou==0,2),:)=[]; %delete off that first row of all zeros
           
            
          LFPEPOCHScom{ses,type,reg,trt}= initmou;  
          
          
    clear initmou temp*

end
        end
    end
end
 LFPmou.comepochstctic=LFPEPOCHScom; 
 clear LFPEPOCHScom LFPEPOCHS
 
LFP.basetrt2=cell(1,1);
for trt=1:2 %only separate out trts if baselines are too uneven for
   % for reg=1:2
%comparison--> like in DRBR
    for Permi=1:50
        
         N=12; %*change to smallest number of trials*
  
       ALL=zeros(4001,7,1);
       for type=1:1 %correct choices only
           for ses=1:5
               Num=size(var1,2);
               for mou=1:Num
                     if size(var1{ses,mou,trt},1)==0 %skip if no data for mouse
                        continue
                     end
                     if mou==4 && trt==2 %mouse 4 seems to be an outlier
                         continue
                     end
                   
               for reg=1:2 %combine btwn regions 
                %for trt=1:2  %combine between trts
                
                temp=LFP.refepochstctic{mou,ses,type,reg,trt};
                
                if size(temp,3)<ncrit %skip mouse if less than n-crit trials
                    continue
                end

                if ses==1 && type==2
                    continue
                end

               temp=datasample(temp,N,3,'replace', false); %chose N random trials from each ses & type w/out repeating
               ALL=cat(3,ALL, temp); % add all bins from each mouse together
               clear temp
           
               %end
               end
               end
           end
       end    
       ALL(:,:,1)=[]; %delete off that first set of zeros in the third dimension 
       
  
        DATA=squeeze(mean(ALL,2)); %clear ALL
        dims=size(DATA);
        
        warning('off')
         for fq=1:num_wavelets
              w_con=fconv_JFC(reshape(DATA,1,dims(1)*dims(2)),w(fq,:));  %reshape data into one dimension -> all trials are lined up end to end (this makes it faster)
              w_con=w_con((size(w,2)-1)/2+1:end-(size(w,2)-1)/2); %remove 1/2 of the wavelet length from each of the data b/c of edge effects
              w_con=reshape(w_con,dims(1),dims(2)); %reform the data back into two dimensions

              BASE(:,fq,Permi)= mean(  abs(w_con(b1:b2,:)).^2,     2); %baseline power

          clear w_con     
         end
         warning('on')

 clear DATA 
 

    end
    %end
    
   %if by seperate trts
   BASEtrt=LFP.basetrt2;
   BASEtrt{1,trt}=mean(BASE,3);
   LFP.basetrt2=BASEtrt; clear BASEtrt
end
save('DRBR_allchanLFP','LFP','-v7.3');  

figure
hold on
bar(1:80, mean(LFP.basetrt2{1,2}(:,:),1),'g')
bar(1:80,mean(LFP.basetrt2{1,1}(:,:),1))


%if by seperate trts
   BASEtrt{trt}=mean(BASE,3);
   LFP.basetrttctic=BASEtrt;
    
%average together base treatments
    temp=cat(3,LFP.basetrt2{1}, LFP.basetrt2{2});
    BASEPER=mean(temp,3);
    LFP.basecombined=BASEPER;
    
clear BASE



    
    
%% Power and Phase Coherence Analysis Per mouse
 
Numses=5;
Numtrt=2;
numTriTypes=2; 

%BASE=LFPmou.basetrt;

ncrit=12;

warning('off')
%analysis - permouse
 for ses=1:5
    for type=1:1
        for trt=1:Numtrt
            %BASE=LFPmou.basetrt{1,trt}; %define what baseline is

            for reg=1:2
                
           Num=size(LFP.refepochstctic,1);
           
           
  
            for mou=1:Num  
                 if size(LFP.refepochstctic{mou,ses,type,reg,trt},2)==0 %done based on TC in DS, but OFC will follow same mouse pattern
                    continue %don't bother to do functions below if there is no data/ mouse
                 end
                 if mou==4 && trt==2 %b/c an outlier on power
                     continue
                 end
                 
          data=squeeze(LFP.refepochstctic{mou,ses,type,reg,trt}(:,1,:));
           if size(data,2) < ncrit
             continue %skip mouse if has less than critical number of trials
           end
           
     for ch=1:7 
 
            data=squeeze(LFP.refepochstctic{mou,ses,type,reg,trt}(:,ch,:));
            
            for perms=1:150 %still must control for trial number
            DATA=datasample(data,ncrit,2,'replace',false);
            dims=size(DATA);
            N=size(DATA,2);

            
                for fq=1:num_wavelets
                     w_con=fconv_JFC(reshape(DATA,1,dims(1)*dims(2)),w(fq,:));  %reshape data into one dimension -> all trials are lined up end to end (this makes it faster)
                     w_con=w_con((size(w,2)-1)/2+1:end-(size(w,2)-1)/2); %remove 1/2 of the wavelet length from each of the data b/c of edge effects
                     w_con=reshape(w_con,dims(1),dims(2)); %reform the data back into two dimensions
                         
                        POWER(:,fq,ch,perms) = mean(abs(w_con(t1:t2,:)).^2,2); %power over time t1 to t2 for each trial
                        POWER_dB(:,fq,ch,perms) = 10*(  log10(POWER(:,fq,ch,perms)) - log10(repmat( mean(LFP.basecombined(:,fq),1)   ,size(tf_tx,2),1))  );  % dB Conversion -->plot this

                        %PHASE_CONST(:,fq,ch,perms) = abs(mean(exp(1i*(  angle(w_con(t1:t2,:))  )),2));  % Derives to Phase Locking Value (Lachaux et al.)        
                        %PHASEz(:,fq,perms)= N*(PHASE_CONST(:,fq,perms).^2);
                end
           clear w_con DATA
            end

     end
            
        dBpow{mou,ses,type,reg,trt}=  mean(POWER_dB,4); 
            clear POWER_dB POWER
        %phase{mou,ses,type,reg,trt}=mean(PHASE_CONST,4);
       % phaseZ{mou,ses,type,reg,trt}=mean(PHASEz,3); 
            %clear PHASE_CONST PHASEz 
            end 

            end
    clear data BASE
        end
    end
 end
 warning('on')
LFP.refpowComBase=dBpow;    
 save('DRBR_allchanLFP','LFP','-v7.3');  
 

 
  figure
   contourf(tf_tx,frex,squeeze(mean(dBpow{1,2,1,1,1},3))',40,'linecolor','none'); 
        %set(gca, 'clim',[.5 4], 'yscale', 'log' ,'YTick', [],'XTick',[])
        set(gca, 'clim', [-10 10], 'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6))
        colorbar;
        
figure
   contourf(tf_tx,frex,POWER',40,'linecolor','none'); 
        %set(gca, 'clim',[.5 4], 'yscale', 'log' ,'YTick', [],'XTick',[])
       % set(gca, 'clim', [-10 10], 'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6))
        colorbar;
        colormap(jet)


 clear phase
LFPmou.phaseZExc=phaseZ; %phase{mou,ses,type,reg,trt}
LFP.phase=phase; %phase{mou,ses,type,reg,trt}

LFPmou.dBpowtrt=dBpow; %phase{mou,ses,type,reg,trt}

clear dBpow

%% 

  figure
   contourf(tf_tx,frex,avg3,40,'linecolor','none'); 
        %set(gca, 'clim',[.5 4], 'yscale', 'log' ,'YTick', [],'XTick',[])
        set(gca, 'clim', [0 2.5], 'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6))
        colormap(jet)
        colorbar;
        
        % average all stages together - TC only, or TC and TIC? 
        %mou,ses,type,reg,trt
        for reg=1:2
            for trt=1:2
        combined=zeros(4001,80,7);
        for ses=1:5
            for type=1:1
                for mou=1:8
                  if size(LFP.refphase{mou,ses,type,reg,trt},2)==0 %done based on TC in DS, but OFC will follow same mouse pattern
                    continue %don't bother to do functions below if there is no data/ mouse
                 end
        temp=LFP.refphase{mou,ses,type,reg,trt};
        combined=cat(4,combined, temp); clear temp
                end
            end
        end
        %clear out first set of all zeros
        size(combined)
        combined(:,:,:,1)=[];
        size(combined)
        TotAvg{reg,trt}=squeeze(mean(combined,4)); clear combined
            end
        end
        LFP.totrefphaseavgtc=TotAvg;
        
   
        
        
        
 %for by electrode graph in topographic order  
 time=-500:500;
 freq=frex(1,26:38);
figure
subplot(4,2,1)
%    contourf(time,freq,LFP.totrefphaseavgtc{1,1}(1500:2500,26:38,7)',40,'linecolor','none'); 
%       set(gca, 'clim', [0 .45],'YTick', [freq(1,1),freq(1,13)])
%       colormap(jet)
subplot(4,2,2)
   contourf(time,freq,LFP.totrefphaseavgtc{2,2}(1500:2500,26:38,4)',40,'linecolor','none'); 
       set(gca, 'clim', [0 .45],'YTick', [freq(1,1),freq(1,13)])
      colormap(jet)
subplot(4,2,3)
   contourf(time,freq,LFP.totrefphaseavgtc{2,2}(1500:2500,26:38,7)',40,'linecolor','none'); 
       set(gca, 'clim', [0 .45],'YTick', [freq(1,1),freq(1,13)])
      colormap(jet)
subplot(4,2,4)
   contourf(time,freq,LFP.totrefphaseavgtc{2,2}(1500:2500,26:38,3)',40,'linecolor','none'); 
       set(gca, 'clim', [0 .45],'YTick', [freq(1,1),freq(1,13)])
      colormap(jet)
subplot(4,2,5)
   contourf(time,freq,LFP.totrefphaseavgtc{2,2}(1500:2500,26:38,6)',40,'linecolor','none'); 
       set(gca, 'clim', [0 .45],'YTick', [freq(1,1),freq(1,13)])
      colormap(jet)
subplot(4,2,6)
   contourf(time,freq,LFP.totrefphaseavgtc{2,2}(1500:2500,26:38,2)',40,'linecolor','none'); 
       set(gca, 'clim', [0 .45],'YTick', [freq(1,1),freq(1,13)])
      colormap(jet)
subplot(4,2,7)
   contourf(time,freq,LFP.totrefphaseavgtc{2,2}(1500:2500,26:38,5)',40,'linecolor','none'); 
     set(gca, 'clim', [0 .45],'YTick', [freq(1,1),freq(1,13)])
      colormap(jet)
subplot(4,2,8)
    contourf(time,freq,LFP.totrefphaseavgtc{2,2}(1500:2500,26:38,1)',40,'linecolor','none'); 
     set(gca, 'clim', [0 .45],'YTick', [freq(1,1),freq(1,13)])
      colormap(jet)

      %colorbar;
        
 

%% find ROI region of significance
%phase{mou,ses,type,reg,trt}

    %define other things
    Numses=5;
    Numtrts=2;
    Numbeh=1;
    Numreg=2;
    nummice=size(LFPmou.phaseZExc,1);


for B=1:Numbeh % have average within beh, but not between 
    
    for R=1:Numreg  
    for T=1:Numtrts
        
        for S=1:Numses
            for M=1:nummice
             if size(LFPmou.phaseZExc{M,S,B,R,T},2)==0
                 continue
             end
              %*** must change dimensions for another measure
                if size(LFPmou.phaseZExc{M,S,B,R,T},1) > 80
                   roi(S,M,:,:)=mean(LFPmou.phaseZExc{M,S,B,R,T},3)'; 
                elseif size(LFPmou.phaseZExc{M,S,B,R,T},1) == 80
                   roi(S,M,:,:)=mean(LFPmou.phaseZExc{M,S,B,R,T},3);
                end

            
            end
        end 
        
        
        avg1(:,:,T)=squeeze(mean(mean(roi,1),2));
    end
       avg2(:,:,R)=mean(avg1,3);
    end

    avg3(:,:,B)=mean(avg1,3); clear roi avg1 %avg2 % time,freq x beh
end


%% graph per session
LFP.mouseavg=avg3;
LFP.sig05=LFPmou.cohZROI05;
 
CONT=zeros(4001,80);
for ses=2:4
    for beh=1:1
    for reg=2:2
        for trt=2:2
     CONT=zeros(4001,80);
     Num=size(LFPmou.epochs,1);
 for j=1:Num
    
     temp=LFPmou.phaseZExc{j,ses,beh,reg,trt}; 
     if size(temp,2)==0
        continue
     end
     CONT=cat(3,CONT,mean(temp,3)); 
     
 end
 
 CONT(:,:,1)=[];%clears first set of zeros
 CONT=mean(CONT,3); %average across mice
 
 SIGPIX=LFP.sig05{beh}(:,12:36);
 
 figure
 hold on
   contourf(tf_tx,frex,CONT',40,'linecolor','none'); 
   contour(tf_tx,frex(12:36),SIGPIX',1,'linecolor','k','linewidth',3); 
   %line([750,1500],[6,6]);line([750,1500],[2,2]);line([750,750],[2,6]);line([1500,1500],[2,6]);
   line([1000,1550],[1.75,1.75]);line([1000,1550],[0,0]);
       % set(gca, 'clim',[.0 1], 'yscale', 'log' ,'YTick', [],'XTick',[])
        set(gca, 'clim', [0 10],'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6))
       colormap(jet)
        %colorbar;
        %title(['Session' num2str(ses) ' Behavior' num2str(beh) ' Treatment' num2str(trt)])
        
 
        end
    end
    end
end
        %CONT(:,:,1)=[];%clears first set of zeros
        %CONT=mean(CONT,3); %average across mice
 
%  figure
%     contourf(tf_tx,frex,CONT',40,'linecolor','none'); 
%         set(gca, 'clim', [.5 5],'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6))
%         colorbar; 
        
        
        %DRBR LIMITS:
       %for ds ylim [.2 .62],
       %for ofc ylim [.15 .5]
       
        %DRBR LIMITS Z-phase:
       %for ds ylim [.5 5],
       %for ofc ylim [.5 4]

%% find ROI
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

LFP.powSig05=SIGPIX;

%Test graphs
SIGPIX=LFPmou.cohZROI05;

 figure
 hold on
   contourf(tf_tx,frex,avg3(:,:,2),40,'linecolor','none');
   contour(tf_tx,frex,SIGPIX{2},1,'linecolor','k','linewidth',3)
   %plot(tf_tx,2,'k','linewidth',3)
   %plot(tf_tx,4,'k','linewidth',3)
   %plot(tf_tx,7,'k','linewidth',3)
   %plot(tf_tx,1.75,'-k','linewidth',3)
    %set(gca, 'clim',[0 2.5], 'yscale', 'log' ,'YTick', [],'XTick',[])
    set(gca, 'yscale', 'log' ,'YTick', logspace(log10(start_freq),log10(end_freq),6),'XTick',-1000:1000:3000)%logspace(log10(start_freq),log10(end_freq),6))
     colorbar
 
 
% Write ROI values to new excel file --> to transfer to statview for
% analysis

% measures={'TC';'TC';'TC';'TC';'TC';'TC';'TC';'TC';'TC';'TC';'TC';'TC';'ph';'ph';...
%           'TC';'TC';'TC';'TC';'TC';'TC';'TC';'TC';'TC';'TC';'TC';'TC';'ph';'ph';... 
%     
%           'TC';'TC';'TC';'TC';'TC';'TC';'TC';'TC';'TC';'TC';'TC';'TC';'TC';'TC';... 
%           'TC';'TC';'TC';'TC';'TC';'TC';'TC';'TC';'TC';'TC';'TC';'TC';'TC';'TC'};
% 
% Region=   {'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'ph';'ph';...      
%            'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'ph';'ph';...
%            
%            'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';...
%            'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC'}; 
%        
% Treatment={'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'ph';'ph';...
%             'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'ph';'ph';...
%             
%            'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';...
%            'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE'};
%        
% MouNum={'2.0';'3.1';'11.0';'16.1';'17.0';'17.1';'17.2';'20.1';'22.2';'25.0';'26.0';'27.0';'ph';'ph';...
%         '2.0';'3.1';'11.0';'16.1';'17.0';'17.1';'17.2';'20.1';'22.2';'25.0';'26.0';'27.0';'ph';'ph';...
%         
%         '4.1';'14.1';'18.0';'19.0';'19.1';'21.9';'21.1';'22.1';'22.2';'24.0';'28.0';'29.0';'30.0';'31.0';...
%         '4.1';'14.1';'18.0';'19.0';'19.1';'21.9';'21.1';'22.1';'22.2';'24.0';'28.0';'29.0';'30.0';'31.0'};
%        
% 
% measuresTIC={'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'ph';'ph';...
%           'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'ph';'ph';... 
%     
%           'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';... 
%           'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC'};


      SIGPIX=LFP.sig05;
    clear control
   
%control=zeros(56,5,2,2,2);
for beh=1:1 %finding ROI x-y limits for both behaviors
        if beh==1,% [row, col]=find(SIGPIX{1}==1); %col are time, row are fq
           ylow=1:1.75;
           yhigh=2:6;
           xlow=500:1750;
           xhigh=750:1500;
           %ytc=row; %freq
           %xtc=col; %time
         elseif beh==2,%[row, col]=find(SIGPIX{2}==1);  %col are time, row are fq
           %ytic=row;
           %xtic=col;
        end

  %picking out ROI limits from phasedata    
   for trt=1:2
      for reg=1:2
         for ses=1:5 
             mouserowstart=1;
             for mou=1:size(LFP.refphase,1)
  
        DATA = LFP.refphase{mou,ses,beh,reg,trt};
        if size(DATA,2)==0
            control(mouserowstart:mouserowstart+6,ses,reg,trt,beh,1)=zeros(7,1);
            control(mouserowstart:mouserowstart+6,ses,reg,trt,beh,2)=zeros(7,1);
            mouserowstart=mouserowstart+7;
        elseif size(DATA,2) ~= 0
             control(mouserowstart:mouserowstart+6,ses,reg,trt,beh,1)=mean(mean(DATA(xlow,ylow,:),2),1); %control, tc
             control(mouserowstart:mouserowstart+6,ses,reg,trt,beh,2)=mean(mean(DATA(xhigh,yhigh,:),2),1); 
            mouserowstart=mouserowstart+7;
        end

%             control(mouserowstart:mouserowstart+6,ses,reg,trt,beh)=mean(mean(DATA(xtc,ytc,:),2),1); %control, tc
%             mouserowstart=mouserowstart+7;
%        

             end
         end  
      end
   end

end
clear control

        D85L=[squeeze(control(:,1,1,1,1,1)); squeeze(control(:,1,2,1,1,1)); squeeze(control(:,1,1,2,1,1)); squeeze(control(:,1,2,2,1,1))]; 
        RS1L=[squeeze(control(:,2,1,1,1,1)); squeeze(control(:,2,2,1,1,1)); squeeze(control(:,2,1,2,1,1)); squeeze(control(:,2,2,2,1,1))]; 
        RS3L=[squeeze(control(:,3,1,1,1,1)); squeeze(control(:,3,2,1,1,1)); squeeze(control(:,3,1,2,1,1)); squeeze(control(:,3,2,2,1,1))]; 
        RS4L=[squeeze(control(:,4,1,1,1,1)); squeeze(control(:,4,2,1,1,1)); squeeze(control(:,4,1,2,1,1)); squeeze(control(:,4,2,2,1,1))]; 
        R50L=[squeeze(control(:,5,1,1,1,1)); squeeze(control(:,5,2,1,1,1)); squeeze(control(:,5,1,2,1,1)); squeeze(control(:,5,2,2,1,1))]; 
               %    sac ds tc                     sac ofc tc                    pae ds tc                           pae ofc tc
               
        D85H=[squeeze(control(:,1,1,1,1,2)); squeeze(control(:,1,2,1,1,2)); squeeze(control(:,1,1,2,1,2)); squeeze(control(:,1,2,2,1,2))]; 
        RS1H=[squeeze(control(:,2,1,1,1,2)); squeeze(control(:,2,2,1,1,2)); squeeze(control(:,2,1,2,1,2)); squeeze(control(:,2,2,2,1,2))]; 
        RS3H=[squeeze(control(:,3,1,1,1,2)); squeeze(control(:,3,2,1,1,2)); squeeze(control(:,3,1,2,1,2)); squeeze(control(:,3,2,2,1,2))]; 
        RS4H=[squeeze(control(:,4,1,1,1,2)); squeeze(control(:,4,2,1,1,2)); squeeze(control(:,4,1,2,1,2)); squeeze(control(:,4,2,2,1,2))]; 
        R50H=[squeeze(control(:,5,1,1,1,2)); squeeze(control(:,5,2,1,1,2)); squeeze(control(:,5,1,2,1,2)); squeeze(control(:,5,2,2,1,2))]; 
        
        
%         D85TIC=[squeeze(control(2,:,1,1))'; squeeze(control(2,:,1,2))'; squeeze(experimental(2,:,1,1))'; squeeze(experimental(2,:,1,2))'];
%         RS1TIC=[squeeze(control(2,:,2,1))'; squeeze(control(2,:,2,2))'; squeeze(experimental(2,:,2,1))'; squeeze(experimental(2,:,2,2))'];
%         RS3TIC=[squeeze(control(2,:,3,1))'; squeeze(control(2,:,3,2))'; squeeze(experimental(2,:,3,1))'; squeeze(experimental(2,:,3,2))'];
%         RS4TIC=[squeeze(control(2,:,4,1))'; squeeze(control(2,:,4,2))'; squeeze(experimental(2,:,4,1))'; squeeze(experimental(2,:,4,2))'];
%         R50TIC=[squeeze(control(2,:,5,1))'; squeeze(control(2,:,5,2))'; squeeze(experimental(2,:,5,1))'; squeeze(experimental(2,:,5,2))'];


T=table(measures,MouNum,Channel,Region,Treatment,D85L,RS1L,RS3L,RS4L,R50L); LFP.cohROI_TC_chann=T;
S=table(measures,MouNum,Channel,Region,Treatment,D85H,RS1H,RS3H,RS4H,R50H); 
%S=table(measuresTIC,MouNum,Region,Treatment,D85TIC,RS1TIC,RS3TIC,RS4TIC,R50TIC); LFP.cohROI_TIC=S;
    writetable(T,'PhaseCoh_DRBR_ref_180119.xls','Sheet',1);    
    writetable(S,'PhaseCoh_DRBR_ref_180119.xls','Sheet',2);

 %T=LFPmou.cohTC;
 %S=LFPmou.cohTIC;
    
 t=table2array(T(:,5:9)); 
 s=table2array(S(:,5:9));
for mes=3:5

    if mes==1,
        y=[1 11 6 16];
    elseif mes==2,
        y=[2 12 7 17];
%     elseif mes==3,
%         y=[3 13 8 18];
    end

    data=t(y,:);
    sd=s(y,:) ./ sqrt(250); %number of perms
    
    figure
    hold on
        bar(data); legend('d85','rs1','rs3','rs4','r50')
        errorbar(data,sd);
        ylim([0 .5]); %use .5 for tc and .3 for tic
        xlabel('ds sac, ds pae, ofc sac, ofc pae')
 
end



%% Write ROI values to new excel file --> to transfer to statview for analysis
%THIS SECTION RUNS FOR 2 DIFFERENT CLUSTERS OF LOW FREQUENCY TC DATA

measures={'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'ph';'ph';...
          'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'ph';'ph';... 
    
          'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';... 
          'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';'TC-L';...

          'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'ph';'ph';...
          'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'ph';'ph';... 
    
          'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';... 
          'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H';'TC-H'};

Region=   {'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS';'ph';'ph'; ...      
           'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'ph';'ph';...
           
           'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';...
           'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';...
       
           'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS';'ph';'ph'; ...      
           'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'ph';'ph';...
           
           'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';...
           'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC'};
       
Treatment={'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'ph';'ph';...
            'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'ph';'ph';...
            
           'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';...
           'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';...
       
            'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'ph';'ph';...
            'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'ph';'ph';...
            
           'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';...
           'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE'};
       
MouNum={'2.0';'3.1';'11.0';'16.1';'17.0';'17.1';'17.2';'20.1';'22.2';'25.0';'26.0';'27.0';'ph';'ph';...
        '2.0';'3.1';'11.0';'16.1';'17.0';'17.1';'17.2';'20.1';'22.2';'25.0';'26.0';'27.0';'ph';'ph';...
        
        '4.1';'14.1';'18.0';'19.0';'19.1';'21.9';'21.1';'22.1';'22.2';'24.0';'28.0';'29.0';'30.0';'31.0';...
        '4.1';'14.1';'18.0';'19.0';'19.1';'21.9';'21.1';'22.1';'22.2';'24.0';'28.0';'29.0';'30.0';'31.0';...
        
        '2.0';'3.1';'11.0';'16.1';'17.0';'17.1';'17.2';'20.1';'22.2';'25.0';'26.0';'27.0';'ph';'ph';...
        '2.0';'3.1';'11.0';'16.1';'17.0';'17.1';'17.2';'20.1';'22.2';'25.0';'26.0';'27.0';'ph';'ph';...
        
        '4.1';'14.1';'18.0';'19.0';'19.1';'21.9';'21.1';'22.1';'22.2';'24.0';'28.0';'29.0';'30.0';'31.0';...
        '4.1';'14.1';'18.0';'19.0';'19.1';'21.9';'21.1';'22.1';'22.2';'24.0';'28.0';'29.0';'30.0';'31.0'};
       

%measuresTIC={'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'ph';'ph';...
          %'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'ph';'ph';... 
    
          %'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';... 
          %'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC';'TIC'};
      
%RegionTIC={'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS';'ph';'ph'; ...      
           %'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'ph';'ph';...
           
          % 'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';...
           %'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';};
       
%TreatmentTIC={'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'ph';'ph';...
            %'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'ph';'ph';...
            
           %'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';...
           %'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE'};
       
%MouNumTIC={'2.0';'3.1';'11.0';'16.1';'17.0';'17.1';'17.2';'20.1';'22.2';'25.0';'26.0';'27.0';'ph';'ph';...
        %'2.0';'3.1';'11.0';'16.1';'17.0';'17.1';'17.2';'20.1';'22.2';'25.0';'26.0';'27.0';'ph';'ph';...
        
        %'4.1';'14.1';'18.0';'19.0';'19.1';'21.9';'21.1';'22.1';'22.2';'24.0';'28.0';'29.0';'30.0';'31.0';...
        %'4.1';'14.1';'18.0';'19.0';'19.1';'21.9';'21.1';'22.1';'22.2';'24.0';'28.0';'29.0';'30.0';'31.0'};
        
        
        
measures={'ROI-TCL';'ROI-TCL';'ROI-TCL';'ROI-TCL';'ROI-TCL';'ROI-TCL';'ROI-TCL';'ph';...
          'ROI-TCL';'ROI-TCL';'ROI-TCL';'ROI-TCL';'ROI-TCL';'ROI-TCL';'ROI-TCL';'ph';...
          'ROI-TCL';'ROI-TCL';'ROI-TCL';'ROI-TCL';'ROI-TCL';'ROI-TCL';'ROI-TCL';'ROI-TCL';...
          'ROI-TCL';'ROI-TCL';'ROI-TCL';'ROI-TCL';'ROI-TCL';'ROI-TCL';'ROI-TCL';'ROI-TCL';...
          
          'ROI-TCH';'ROI-TCH';'ROI-TCH';'ROI-TCH';'ROI-TCH';'ROI-TCH';'ROI-TCH';'ph';...
          'ROI-TCH';'ROI-TCH';'ROI-TCH';'ROI-TCH';'ROI-TCH';'ROI-TCH';'ROI-TCH';'ph';...
          'ROI-TCH';'ROI-TCH';'ROI-TCH';'ROI-TCH';'ROI-TCH';'ROI-TCH';'ROI-TCH';'ROI-TCH';...
          'ROI-TCH';'ROI-TCH';'ROI-TCH';'ROI-TCH';'ROI-TCH';'ROI-TCH';'ROI-TCH';'ROI-TCH';...
          
          'ROI-TIC';'ROI-TIC';'ROI-TIC';'ROI-TIC';'ROI-TIC';'ROI-TIC';'ROI-TIC';'ph';...
          'ROI-TIC';'ROI-TIC';'ROI-TIC';'ROI-TIC';'ROI-TIC';'ROI-TIC';'ROI-TIC';'ph';...
          'ROI-TIC';'ROI-TIC';'ROI-TIC';'ROI-TIC';'ROI-TIC';'ROI-TIC';'ROI-TIC';'ROI-TIC';...
          'ROI-TIC';'ROI-TIC';'ROI-TIC';'ROI-TIC';'ROI-TIC';'ROI-TIC';'ROI-TIC';'ROI-TIC'};
      
Region=   {'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'ph'; ...      
           'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'ph';...
           'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';...
           'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';...
           
           'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'ph'; ...      
           'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'ph';...
           'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';...
           'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';...
           
           'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'ph'; ...      
           'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'ph';...
           'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';...
           'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC'};
       
Treatment={'FLOX';'FLOX';'FLOX';'FLOX';'FLOX';'FLOX';'FLOX';'ph';...
           'FLOX';'FLOX';'FLOX';'FLOX';'FLOX';'FLOX';'FLOX';'ph';...
           'NULL';'NULL';'NULL';'NULL';'NULL';'NULL';'NULL';'NULL';...
           'NULL';'NULL';'NULL';'NULL';'NULL';'NULL';'NULL';'NULL';...
           
           'FLOX';'FLOX';'FLOX';'FLOX';'FLOX';'FLOX';'FLOX';'ph';...
           'FLOX';'FLOX';'FLOX';'FLOX';'FLOX';'FLOX';'FLOX';'ph';...
           'NULL';'NULL';'NULL';'NULL';'NULL';'NULL';'NULL';'NULL';...
           'NULL';'NULL';'NULL';'NULL';'NULL';'NULL';'NULL';'NULL';...
           
           'FLOX';'FLOX';'FLOX';'FLOX';'FLOX';'FLOX';'FLOX';'ph';...
           'FLOX';'FLOX';'FLOX';'FLOX';'FLOX';'FLOX';'FLOX';'ph';...
           'NULL';'NULL';'NULL';'NULL';'NULL';'NULL';'NULL';'NULL';...
           'NULL';'NULL';'NULL';'NULL';'NULL';'NULL';'NULL';'NULL'};
       
MouNum={'2.2';'3.2';'4.0';'4.1';'5.0';'7.1';'7.2';'ph';...
        '2.2';'3.2';'4.0';'4.1';'5.0';'7.1';'7.2';'ph';...
        '2.0';'3.0';'3.1';'5.1';'6.0';'6.1';'7.0';'8.0';...
        '2.0';'3.0';'3.1';'5.1';'6.0';'6.1';'7.0';'8.0';...
        
        '2.2';'3.2';'4.0';'4.1';'5.0';'7.1';'7.2';'ph';...
        '2.2';'3.2';'4.0';'4.1';'5.0';'7.1';'7.2';'ph';...
        '2.0';'3.0';'3.1';'5.1';'6.0';'6.1';'7.0';'8.0';...
        '2.0';'3.0';'3.1';'5.1';'6.0';'6.1';'7.0';'8.0';...
        
        '2.2';'3.2';'4.0';'4.1';'5.0';'7.1';'7.2';'ph';...
        '2.2';'3.2';'4.0';'4.1';'5.0';'7.1';'7.2';'ph';...
        '2.0';'3.0';'3.1';'5.1';'6.0';'6.1';'7.0';'8.0';...
        '2.0';'3.0';'3.1';'5.1';'6.0';'6.1';'7.0';'8.0'};

    %SIGPIX=LFPmou.cohZROI05;
 %make sure dimensions agree !!!! 
control=zeros(2,8,5,1);
experimental=zeros(2,8,5,1);
for beh=1:1 %finding ROI x-y limits for both behaviors
        if beh==1, [row, col]=find(LFPmou.cohZROI05{1}'==1); %row are FQ, col are TIME
           lowbound=find(row<=11); %COL 11=1.75 CUTOFF; COL 9=1.5HZ CUTOFF
            ytcL=row(lowbound);
            xtcL=col(lowbound); 
           highbound=find(row>11); %CHANGE COL HERE TOO
            ytcH=row(highbound);
            xtcH=col(highbound);
            clear row col
         elseif beh==2,[row, col]=find(LFPmou.cohZROI05{2}'==1); %row are FQ, col are TIME
           ytic=row;
           xtic=col;
           clear row col
        end

  %picking out ROI limits from phasedata    
   for trt=1:2
      for reg=1:2
         for ses=1:5 
             Num=size(LFPmou.epochs,1);
             for mou=1:Num
  
        DATA = LFPmou.phaseZExc{mou,ses,beh,reg,trt}; %time is row fq is col
        if size(DATA,1)==0
            continue
        end
        %DATA=LFP.phaseSDtctic{ses,beh,reg,trt}; %standard deviation between permutations
       
        if trt==1 && beh==1,
            control(1,mou,ses,reg)=mean(mean(DATA(xtcL,ytcL),2),1); %low fq control tc
            control(2,mou,ses,reg)=mean(mean(DATA(xtcH,ytcH),2),1); %hi fq control tc
        elseif trt==1 && beh==2,
            control(3,mou,ses,reg)=mean(mean(DATA(xtic,ytic),2),1); %one fq control tic - same as above
            
        elseif trt==2 && beh==1,
            experimental(1,mou,ses,reg)=mean(mean(DATA(xtcL,ytcL),2),1);%low fq exp tc
            experimental(2,mou,ses,reg)=mean(mean(DATA(xtcH,ytcH),2),1);%hi fq exp tc
        elseif trt==2 && beh==2,
            experimental(3,mou,ses,reg)=mean(mean(DATA(xtic,ytic),2),1);
        end

             end
         end  
      end
   end

end

        D85=[squeeze(control(1,:,1,1))'; squeeze(control(1,:,1,2))'; squeeze(experimental(1,:,1,1))'; squeeze(experimental(1,:,1,2))';...
            squeeze(control(2,:,1,1))'; squeeze(control(2,:,1,2))'; squeeze(experimental(2,:,1,1))'; squeeze(experimental(2,:,1,2))';...
            squeeze(control(3,:,1,1))'; squeeze(control(3,:,1,2))'; squeeze(experimental(3,:,1,1))'; squeeze(experimental(3,:,1,2))']; 
        
        RS1=[squeeze(control(1,:,2,1))'; squeeze(control(1,:,2,2))'; squeeze(experimental(1,:,2,1))'; squeeze(experimental(1,:,2,2))';...
            squeeze(control(2,:,2,1))'; squeeze(control(2,:,2,2))'; squeeze(experimental(2,:,2,1))'; squeeze(experimental(2,:,2,2))';...
            squeeze(control(3,:,2,1))'; squeeze(control(3,:,2,2))'; squeeze(experimental(3,:,2,1))'; squeeze(experimental(3,:,2,2))'];
        
        RS3=[squeeze(control(1,:,3,1))'; squeeze(control(1,:,3,2))'; squeeze(experimental(1,:,3,1))'; squeeze(experimental(1,:,3,2))';...
            squeeze(control(2,:,3,1))'; squeeze(control(2,:,3,2))'; squeeze(experimental(2,:,3,1))'; squeeze(experimental(2,:,3,2))';...
            squeeze(control(3,:,3,1))'; squeeze(control(3,:,3,2))'; squeeze(experimental(3,:,3,1))'; squeeze(experimental(3,:,3,2))'];
        
        RS4=[squeeze(control(1,:,4,1))'; squeeze(control(1,:,4,2))'; squeeze(experimental(1,:,4,1))'; squeeze(experimental(1,:,4,2))';...
            squeeze(control(2,:,4,1))'; squeeze(control(2,:,4,2))'; squeeze(experimental(2,:,4,1))'; squeeze(experimental(2,:,4,2))';...
            squeeze(control(3,:,4,1))'; squeeze(control(3,:,4,2))'; squeeze(experimental(3,:,4,1))'; squeeze(experimental(3,:,4,2))']; 
        
        R50=[squeeze(control(1,:,5,1))'; squeeze(control(1,:,5,2))'; squeeze(experimental(1,:,5,1))'; squeeze(experimental(1,:,5,2))';...
            squeeze(control(2,:,5,1))'; squeeze(control(2,:,5,2))'; squeeze(experimental(2,:,5,1))'; squeeze(experimental(2,:,5,2))';...
            squeeze(control(3,:,5,1))'; squeeze(control(3,:,5,2))'; squeeze(experimental(3,:,5,1))'; squeeze(experimental(3,:,5,2))']; 
        %            sac ds tcL                  sac ofc tcL                 pae ds tcL                         pae ofc tcL
        %           sac ds tcH                  sac ofc tcH                 pae ds tcH                         pae ofc tcH
        
        
%         D85TIC=[squeeze(control(3,:,1,1))'; squeeze(control(3,:,1,2))'; squeeze(experimental(3,:,1,1))'; squeeze(experimental(3,:,1,2))'];
%         RS1TIC=[squeeze(control(3,:,2,1))'; squeeze(control(3,:,2,2))'; squeeze(experimental(3,:,2,1))'; squeeze(experimental(3,:,2,2))'];
%         RS3TIC=[squeeze(control(3,:,3,1))'; squeeze(control(3,:,3,2))'; squeeze(experimental(3,:,3,1))'; squeeze(experimental(3,:,3,2))'];
%         RS4TIC=[squeeze(control(3,:,4,1))'; squeeze(control(3,:,4,2))'; squeeze(experimental(3,:,4,1))'; squeeze(experimental(3,:,4,2))'];
%         R50TIC=[squeeze(control(3,:,5,1))'; squeeze(control(3,:,5,2))'; squeeze(experimental(3,:,5,1))'; squeeze(experimental(3,:,5,2))'];


T=table(measures,MouNum,Region,Treatment,D85,RS1,RS3,RS4,R50); LFPmou.table_cohTCTICroi=T;
%S=table(measuresTIC,MouNumTIC,RegionTIC,TreatmentTIC,D85TIC,RS1TIC,RS3TIC,RS4TIC,R50TIC); LFPmou.cohZTIC2=S;
    writetable(T,'PhaseCoh_DRBR_ROI.xls','Sheet',1);    
    %writetable(S,'PhaseCoh_ROIval_2low.xls','Sheet',2);

 
    
    
%% Line Graphs w/ individual dots per subject   
%plot mean w/ overlayed subject dots
 
%Ltc, Htc, tic 
%control(1,mou,ses,reg)
%experimental(1,mou,ses,reg)


for reg=1:2
    for fqr=1:2
       for ses=1:5
            tempdata=control(fqr,:,ses,reg);
            tempdata(tempdata==0)=[];
              means_cont(1,ses)=mean(tempdata);
              errors_cont(1,ses)=std(tempdata,0,2)/ sqrt(size(tempdata,2)); clear tempdata
              
            tempdata=experimental(fqr,:,ses,reg);
            tempdata(tempdata==0)=[];
              means_exp(1,ses)=mean(tempdata);
              errors_exp(1,ses)=std(tempdata,0,2)/ sqrt(size(tempdata,2)); clear tempdata
        end       

    figure %plot the individual points
    hold on
    %control
   % scatter(ones(length(control(fqr,:,1,reg)),1),control(fqr,:,1,reg),90,'sk','LineWidth',1); %^ = triangle
  %  scatter(repmat(2,length(control(fqr,:,2,reg)),1),control(fqr,:,2,reg),90,'sk','LineWidth',1);
  %  scatter(repmat(3,length(control(fqr,:,3,reg)),1),control(fqr,:,3,reg),90,'sk','LineWidth',1);
   % scatter(repmat(4,length(control(fqr,:,4,reg)),1),control(fqr,:,4,reg),90,'sk','LineWidth',1);
  %  scatter(repmat(5,length(control(fqr,:,5,reg)),1),control(fqr,:,5,reg),90,'sk','LineWidth',1);
    
    %experimental
   % scatter(ones(length(experimental(fqr,:,1,reg)),1),experimental(fqr,:,1,reg),90,'sr','LineWidth',1);
   % scatter(repmat(2,length(experimental(fqr,:,2,reg)),1),experimental(fqr,:,2,reg),90,'sr','LineWidth',1);
   % scatter(repmat(3,length(experimental(fqr,:,3,reg)),1),experimental(fqr,:,3,reg),90,'sr','LineWidth',1);
   % scatter(repmat(4,length(experimental(fqr,:,4,reg)),1),experimental(fqr,:,4,reg),90,'sr','LineWidth',1);
   % scatter(repmat(5,length(experimental(fqr,:,5,reg)),1),experimental(fqr,:,5,reg),90,'sr','LineWidth',1);
   % xlim([0 6])
    
    if fqr==1 || fqr==2
        ylim([0 7]) %10,10, 5
    elseif fqr==2
        ylim([0 7])
    elseif fqr==3
        ylim([0 5])
    end

   
    hold on %plot the average + standard error
    errorbar(means_cont,errors_cont,'-sk','MarkerFaceColor','k','markers',14,'LineWidth',2) %12
    errorbar(means_exp,errors_exp,'-sr','MarkerFaceColor','r','markers',14,'LineWidth',2)
    
  
    
    title(['Measure' num2str(fqr) 'Region' num2str(reg)]);


    end
end







%%



for reg=1:1
    for fqr=1:1
       for ses=1:5
            tempdata=control(fqr,:,ses,reg);
            tempdata(tempdata==0)=[];
              means_cont(1,ses)=mean(tempdata);
              errors_cont(1,ses)=std(tempdata,0,2)/ sqrt(size(tempdata,2)); clear tempdata
              
            tempdata=experimental(fqr,:,ses,reg);
            tempdata(tempdata==0)=[];
              means_exp(1,ses)=mean(tempdata);
              errors_exp(1,ses)=std(tempdata,0,2)/ sqrt(size(tempdata,2)); clear tempdata
        end       

    figure %plot the individual points
    hold on
    %control
    scatter(ones(length(control(fqr,:,1,reg)),1),control(fqr,:,1,reg),90,'^k','LineWidth',1); %^ = triangle
    scatter(repmat(2,length(control(fqr,:,2,reg)),1),control(fqr,:,2,reg),90,'^k','LineWidth',1);
    scatter(repmat(3,length(control(fqr,:,3,reg)),1),control(fqr,:,3,reg),90,'^k','LineWidth',1);
    scatter(repmat(4,length(control(fqr,:,4,reg)),1),control(fqr,:,4,reg),90,'^k','LineWidth',1);
    scatter(repmat(5,length(control(fqr,:,5,reg)),1),control(fqr,:,5,reg),90,'^k','LineWidth',1);
    
    %experimental
    scatter(ones(length(experimental(fqr,:,1,reg)),1),experimental(fqr,:,1,reg),90,'^r','LineWidth',1);
    scatter(repmat(2,length(experimental(fqr,:,2,reg)),1),experimental(fqr,:,2,reg),90,'^r','LineWidth',1);
    scatter(repmat(3,length(experimental(fqr,:,3,reg)),1),experimental(fqr,:,3,reg),90,'^r','LineWidth',1);
    scatter(repmat(4,length(experimental(fqr,:,4,reg)),1),experimental(fqr,:,4,reg),90,'^r','LineWidth',1);
    scatter(repmat(5,length(experimental(fqr,:,5,reg)),1),experimental(fqr,:,5,reg),90,'^r','LineWidth',1);
    xlim([0 6])
    
    if fqr==1 || fqr==2
        ylim([0 6]) %10,10, 5
    elseif fqr==2
        ylim([0 6])
    elseif fqr==3
        ylim([0 2.5])
    end

    hold on %plot the average + standard error
    errorbar(means_cont,errors_cont,'-sk','MarkerFaceColor','k','markers',3,'LineWidth',2) %12
    errorbar(means_exp,errors_exp,'-sr','MarkerFaceColor','r','markers',3,'LineWidth',2)
    
  
    
    title(['Measure' num2str(fqr) 'Region' num2str(reg)]);


    end
end
 