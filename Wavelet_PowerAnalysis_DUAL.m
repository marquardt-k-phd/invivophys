%% Set the current directory and add Plexon code folders and data 

clear all; clc

cd('E:\Documents\MATLAB')
addpath(genpath('C:\Users\Public\Desktop\Plexon Inc Documents'));
addpath(genpath('E:\Documents\MATLAB\Cavanagh Class (Lec+Scripts)'));
addpath(genpath('E:\Documents\MATLAB\eeglab_addons'));
addpath(genpath('E:\Documents\MATLAB\eeglab12_0_2_1b'));
addpath(genpath('E:\Documents\MATLAB\kakearney-boundedline-pkg-2112a2b'));

addpath('E:\Documents\MATLAB\DATA_PAE');



%% Files to open

%Behavioral timestamp files
%Event/ timestamp datafiles
%epoch'ed data files

load DRBR_FP

load PAE_FP_DUAL_all
 
save('PAE_FP_DUAL_all','LFP','-v7.3')


 
%% SORT INTO EPOCHS BASED ON TRIAL TYPES
Numses=5;
Numtrt=2;

for ses=1:Numses
    for trt=1:2
        for reg=1:2
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
          lfp=cell2mat(LFP.allFP(j,chan,ses,reg,trt));
          int=cell2mat(LFP.interval(j,chan,ses,trt));

           for b=1:size(trialtypes,1) %if everything is perfect
               temp=find(int <= trials(b,1)+.001 & int >= trials(b,1)-.001); %find where the LFP is closest to the spike in time
                    if trials(b,1) > max(int) %some trials are beyond max lfp timepoint
                        templfp(b,:)=zeros(1,4001);
                    elseif trials(b,1)<= max(int) %for those that are not
                        startrow=temp-1000; %1sec prior
                        endrow=temp+3000; %3 sec post
                        templfp(b,:)=lfp(startrow:endrow,:); 
                    end
               clear temp startrow endrow
           end
          
        %clear out zeros -> in both LFP and trial type arrays
            tempzero=find(templfp(all(templfp==0,2),:));
            trialtypes(tempzero,:)=[];
            templfp(tempzero,:)=[];
        
        %split bins into trial types % store     
        LFPEPOCHS{j,ses,1,reg,trt}=templfp(logical(trialtypes(:,3)==1),:); %lfp bins that are type 1 trials 
        LFPEPOCHS{j,ses,2,reg,trt}=templfp(logical(trialtypes(:,3)==2),:);
        LFPEPOCHS{j,ses,3,reg,trt}=templfp(logical(trialtypes(:,3)==3),:);
        LFPEPOCHS{j,ses,4,reg,trt}=templfp(logical(trialtypes(:,3)==4),:);
             

   clear tone tc punish tic trial* type* x y z w lfp int temp*
end
    clear Num
        end
    end
end

LFP.epochs=LFPEPOCHS; 

% Combine epochs across mice
for ses=1:Numses
    for trt=1:Numtrt
        for reg=1:2

for type=1:4
       Numfiles=size(var1,2);
        
       initmou=zeros(1,4001);
        
       for j=1:Numfiles
            if size(var1{ses,j,trt},1)==0
                continue
            end
           temp=LFP.epochs{j,ses,type,reg,trt};
           initmou=[initmou; temp]; % add all bins from each mouse together
       end
           initmou(all(initmou==0,2),:)=[]; %delete off that first row of all zeros
           
            
          LFPEPOCHScom{ses,type,reg,trt}= initmou;  
          
          
    clear initmou temp*

end
        end
    end
end

LFP.comepochs=LFPEPOCHScom; clear LFPEPOCHScom LFPEPOCHS

 %save('DRBR_FP','LFP','-v7.3');
 
%% SORT INTO EPOCHS - TC AND TIC ONLY
Numses=5;
Numtrt=2;

for ses=1:Numses
    for trt=1:2
        for reg=1:2
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
        
%         %define trial types
%          for t=2:size(trials,1)
%                 k=t-1;
%                 if trialtypes(k,2) == 1 && trialtypes(t,2) == 1
%                     trialtypes(t,3)=1;      % correct -> correct
%                 elseif trialtypes(k,2) == 1 && trialtypes(t,2) == 2
%                     trialtypes(t,3)=2;      %correct -> incorrect
%                 elseif trialtypes(k,2) == 2 && trialtypes(t,2) == 2
%                     trialtypes(t,3)=3;      %incorrect -> incorrect
%                 elseif trialtypes(k,2) == 2 && trialtypes(t,2) == 1
%                     trialtypes(t,3)=4;       %incorrect -> correct
%                 end
%          end

        %find lfp times & create trial epochs 
          lfp=cell2mat(LFP.avgFP(j,ses,reg,trt));
          int=cell2mat(LFP.interval(j,ses,trt));

           for b=1:size(trialtypes,1) %if everything is perfect
               temp=find(int <= trials(b,1)+.001 & int >= trials(b,1)-.001); %find where the LFP is closest to the spike in time
                    if trials(b,1) > max(int) %some trials are beyond max lfp timepoint
                        templfp(b,:)=zeros(1,4001);
                    elseif trials(b,1)<= max(int) %for those that are not
                        startrow=temp-1000; %1sec prior
                        endrow=temp+3000; %3 sec post
                        templfp(b,:)=lfp(startrow:endrow,:);
                        tempint(b,:)=int(startrow:endrow,:);
                    end
               clear temp startrow endrow
           end
          
        %clear out zeros -> in both LFP and trial type arrays
            tempzero=find(templfp(all(templfp==0,2),:));
            trialtypes(tempzero,:)=[];
            templfp(tempzero,:)=[];
            tempint(tempzero,:)=[];
        
        %split bins into trial types % store     
        LFPEPOCHS{j,ses,1,reg,trt}=templfp(logical(trialtypes(:,2)==1),:); %lfp bins that are correct trials  
        LFPEPOCHS{j,ses,2,reg,trt}=templfp(logical(trialtypes(:,2)==2),:); %lfp bins that are incorrect trials  
        LFPInt{j,ses,1,trt}=tempint(logical(trialtypes(:,2)==1),:);
        LFPInt{j,ses,2,trt}=tempint(logical(trialtypes(:,2)==2),:);
             

   clear tone tc punish tic trial* type* x y z w lfp int temp*
end
    clear Num
        end
    end
end
 LFP.epochstctic=LFPEPOCHS; 
 SFC.LFPepochint=LFPInt;
 SFC.LFPepochs=LFPEPOCHS;

% Combine epochs across mice
for ses=1:Numses
    for trt=1:Numtrt
        for reg=1:2

for type=1:2
       Numfiles=size(var1,2);
        
       initmou=zeros(1,4001);
        
       for j=1:Numfiles
            if size(var1{ses,j,trt},1)==0
                continue
            end
           temp=LFP.epochstctic{j,ses,type,reg,trt};
           initmou=[initmou; temp]; % add all bins from each mouse together
       end
           initmou(all(initmou==0,2),:)=[]; %delete off that first row of all zeros
           
            
          LFPEPOCHScom{ses,type,reg,trt}= initmou;  
          
          
    clear initmou temp*

end
        end
    end
end
 LFP.comepochstctic=LFPEPOCHScom; 
 
 clear LFPEPOCHScom LFPEPOCHS

 %save('PAE_FP_DUAL','LFP','-v7.3');

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

    % Graph... if curious
    %     figure;
    %     subplot(3,1,1)
    %     plot(t,real(w(4,:))); 
    %     subplot(3,1,2)
    %     plot(t,real(w(10,:))); 
    %     subplot(3,1,3)
    %     plot(t,real(w(35,:))); 
    
 
 %% Convlolve wavelet family with epoched data --> COMBINED BASELINE
  %clear  POWER POWER_dB  PHASE PHASEz DATA PHASE_CONST
  %clear BASE BASEP
NumTriTypes=2;  
Numses=5;
  
   %set some time variables
    tx=(-1000:1:3000)'; %total length of data
    b1=find(tx==-1000); %beginning of baseline
    b2=find(tx==0); %end of baseline
    t1=find(tx==-1000); %beginning of trial--> at very start so get whole timespan in power graphs
    t2=find(tx==3000); %end of trial
    tf_tx=-1000:1:3000;  %total trial time
    
  
%Combine trials from all mice w/in a session
%comepochs = session (5) x type x ds/ofc x trt


%for trt=1:2 %only separate out trts if baselines are too uneven for
   % for reg=1:2
%comparison--> like in DRBR
    for Permi=1:500
        
         N=61; %*change to smallest number of trials*
  
       ALL=zeros(1,4001);
       for type=1:NumTriTypes 
           for ses=1:Numses
               for reg=1:2  %combine btwn regions 
                for trt=1:2  %combine between trts
                  
               temp=LFP.comepochstctic{ses,type,reg,trt};
               temp=datasample(temp,N,1,'replace', false); %chose N random trials from each ses & type w/out repeating
               ALL=[ALL; temp]; % add all bins from each mouse together
               clear temp
           
                end
               end
           end
       end    
       ALL(all(ALL==0,2),:)=[]; %delete off that first row of all zeros
       
       
        DATA=ALL'; clear All
        dims=size(DATA);
    
         for fq=1:num_wavelets
              w_con=fconv_JFC(reshape(DATA,1,dims(1)*dims(2)),w(fq,:));  %reshape data into one dimension -> all trials are lined up end to end (this makes it faster)
              w_con=w_con((size(w,2)-1)/2+1:end-(size(w,2)-1)/2); %remove 1/2 of the wavelet length from each of the data b/c of edge effects
              w_con=reshape(w_con,dims(1),dims(2)); %reform the data back into two dimensions

              BASE(:,fq,Permi)= mean(  abs(w_con(b1:b2,:)).^2,     2); %baseline power

          clear w_con     
         end

 clear DATA 

    end
    %end
%end

%if by seperate trts
   % BASEtrt{trt,reg}=mean(BASE,3);
    %LFP.basetrttctic=BASEtrt;
    
%if all combined
    BASEPER=mean(BASE,3);
    LFP.basetctic=BASEPER;
    
clear BASE

 %% Convlolve wavelet family with epoched data --> Permutation 
   clear POWER PHASE_CONST POWER_dB
    clear dBpow dBpowsd phase phasesd

 %set some time variables
    tx=(-1000:1:3000)'; %total length of data
    t1=find(tx==-1000); %beginning of trial --> use whole trial to get power to graph
    t2=find(tx==3000); %end of trial
    tf_tx=-1000:1:3000;  %total trial time

Numses=5;
Numtrt=2;
numTriTypes=2; 

BASE=LFP.basetctic;

 %Power Analysis
 %run the data
    n=100; %set to lowest trial number of all sessions
for ses=1:Numses
    for type=1:numTriTypes
        for trt=1:Numtrt
            for reg=1:2
                data=LFP.comepochstctic{ses,type,reg,trt}';
               
        for Permi=1:250 %tested with 50
            DATA=datasample(data,n,2); %,'replace', false);
            dims=size(DATA);
            %N=size(DATA,2);
            
                for fq=1:num_wavelets
                     w_con=fconv_JFC(reshape(DATA,1,dims(1)*dims(2)),w(fq,:));  %reshape data into one dimension -> all trials are lined up end to end (this makes it faster)
                     w_con=w_con((size(w,2)-1)/2+1:end-(size(w,2)-1)/2); %remove 1/2 of the wavelet length from each of the data b/c of edge effects
                     w_con=reshape(w_con,dims(1),dims(2)); %reform the data back into two dimensions
                         
                        POWER(:,fq,Permi) = mean(abs(w_con(t1:t2,:)).^2,2); %power over time t1 to t2 for each trial
                        POWER_dB(:,fq,Permi) = 10*(  log10(POWER(:,fq)) - log10(repmat(mean(BASE(:,fq),1),size(tf_tx,2),1))  );  % dB Conversion -->plot this

                        PHASE_CONST(:,fq,Permi) = abs(mean(exp(1i*(  angle(w_con(t1:t2,:))  )),2));  % Derives to Phase Locking Value (Lachaux et al.)        
                        %PHASEz(:,fq,Permi)= N*(PHASE_CONST(:,fq,Permi).^2);
                end

           clear w_con DATA
        end
        
        %pow{ses,type,reg,trt}= mean(POWER,3); 
            %powSD{ses,type,reg,trt}=std(POWER,0,3);
            clear POWER
        dBpow{ses,type,reg,trt}= mean(POWER_dB,3); 
            dBpowSD{ses,type,reg,trt}=std(POWER_dB,0,3);
            clear POWER_dB
        phase{ses,type,reg,trt}=mean(PHASE_CONST,3); 
            phaseSD{ses,type,reg,trt}=std(PHASE_CONST,0,3);
            clear PHASE_CONST
        %zphase{sesi,condi}=mean(PHASEz,3); clear PHASEz
           end
    clear data
        end
    end
end

LFP.powtctic=pow;
LFP.powSDtctic=powSD;

LFP.dBpowtctic=dBpow;
LFP.dBpowSDtctic=dBpowSD;

LFP.phasetctic=phase;
LFP.phaseSDtctic=phaseSD;

save('PAE_FP_DUAL_ver','LFP','-v7.3');

%Phase analysis - permouse
BASE=LFP.basetctic;
 for ses=1:Numses
    for type=1:numTriTypes
        for trt=1:Numtrt
            for reg=1:2
           Num=size(LFP.epochstctic,1);
            for mou=1:Num
                if size(LFP.epochstctic{mou,ses,type,reg,trt},1)==0 %done based on TC in DS, but OFC will follow same mouse pattern
                    continue
                end

            DATA=LFP.epochstctic{mou,ses,type,reg,trt}';
            dims=size(DATA);
            N=size(DATA,2);
            
                for fq=1:num_wavelets
                     w_con=fconv_JFC(reshape(DATA,1,dims(1)*dims(2)),w(fq,:));  %reshape data into one dimension -> all trials are lined up end to end (this makes it faster)
                     w_con=w_con((size(w,2)-1)/2+1:end-(size(w,2)-1)/2); %remove 1/2 of the wavelet length from each of the data b/c of edge effects
                     w_con=reshape(w_con,dims(1),dims(2)); %reform the data back into two dimensions
                         
                        POWER(:,fq) = mean(abs(w_con(t1:t2,:)).^2,2); %power over time t1 to t2 for each trial
                        POWER_dB(:,fq) = 10*(  log10(POWER(:,fq)) - log10(repmat(mean(BASE(:,fq),1),size(tf_tx,2),1))  );  % dB Conversion -->plot this

                        PHASE_CONST(:,fq) = abs(mean(exp(1i*(  angle(w_con(t1:t2,:))  )),2));  % Derives to Phase Locking Value (Lachaux et al.)        
                        PHASEz(:,fq)= N*(PHASE_CONST(:,fq).^2);
                end

           clear w_con DATA
   
        dBpow{mou,ses,type,reg,trt}=  POWER_dB; 
            clear POWER_dB
        phase{mou,ses,type,reg,trt}=PHASE_CONST;
        phaseZ{mou,ses,type,reg,trt}=PHASEz; 
            clear PHASE_CONST PHASEz
            
            end 
            end
    clear data
        end
    end
 end
 
 LFP.phaseZtcticMouse=phaseZ;
 
CONT=zeros(4001,80);
for ses=4:4
    for beh=1:1
    for reg=1:1
        for trt=2:2
     CONT=zeros(4001,80);
 for j=1:12
    
     temp=phaseZ{j,ses,beh,reg,trt}; 
     if size(temp,2)==0;
        continue
     end
     CONT=cat(3,CONT,temp); 
     
 end
        end
    end
    end
end
        CONT(:,:,1)=[];%clears first set of zeros
        CONT=mean(CONT,3); %average across mice
 
%  figure
%     contourf(tf_tx,frex,mean(CONT,3)',40,'linecolor','none'); 
%         set(gca, 'clim', [0 6],  'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6))
%         colorbar; 
        
    figure
    hold on
      contourf(tf_tx,frex,CONT',40,'linecolor','none');
      contour(tf_tx,frex,SIGPIX{ses,beh,reg}',1,'linecolor','k','linewidth',3)
      set(gca, 'clim', [0 6],  'yscale', 'log' ,'YTick', [],'XTick',[])%logspace(log10(start_freq),log10(end_freq),6))
      colorbar
        

%% GRAPHING - SPECTRA 

 %set some time variables
    tx=(-1000:1:3000)'; %total length of data
    t1=find(tx==-1000); %beginning of trial --> use whole trial to get power to graph
    t2=find(tx==3000); %end of trial
    tf_tx=-1000:1:3000;  %total trial time
    
for ses=1:1
    for type=1:1
        for reg=1:2
            %for trt=1:2
                
        POWcont=LFP.dBpowtrt{ses,type,reg,1}';
        POWexp=LFP.dBpowtrt{ses,type,reg,2}';
        
        %PHASEcont=LFP.phasetctic{ses,type,reg,1};
        %PHASEexp=LFP.phasetctic{ses,type,reg,2};
        
    figure
    subplot(1,2,1)
    contourf(tf_tx,frex,POWcont,40,'linecolor','none'); 
        set(gca, 'clim', [-6 6],  'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6))
        xlabel('Time (ms)');  ylabel('Frequency (Hz)');
        colorbar; 
    subplot(1,2,2)
    contourf(tf_tx,frex,POWexp,40,'linecolor','none'); 
        set(gca, 'clim', [-6 6],  'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6))
        %title(['POWER:    Session:' num2str(ses) 'Region:' num2str(reg) 'Type:' num2str(type)  ]); 
        xlabel('Time (ms)');  ylabel('Frequency (Hz)');
        colorbar; 

%     figure
%     subplot(1,2,1)
%     contourf(tf_tx,frex,PHASEcont',40,'linecolor','none'); 
%        set(gca,'clim', [0 .6], 'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6))
%        xlabel('Time (ms)');ylabel('Frequency (Hz)');
%        colorbar; 
%     subplot(1,2,2)
%     contourf(tf_tx,frex,PHASEexp',40,'linecolor','none'); 
%        set(gca,'clim', [0 .6], 'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6))
%        title(['PHASE:    Session:' num2str(ses) 'Region:' num2str(reg) 'Type:' num2str(type) ] ); 
%        xlabel('Time (ms)');ylabel('Frequency (Hz)');
%        colorbar; 

  clear POWcont POWexp PHASEcont PHASEexp
            %end
        end
    end
end



%% GRAPHING - LINES


%pick out phase frequency bands and then graph
for condi=1:2
    %for geno=1:2
       % for ses=5:5
            for reg=1:2
                
                phase=DRBR(1).PHASE;
            
   % POWERdeltaWT=mean(dBpow{1,reg,condi}(:,1:24),2);
    %POWERdeltaCKD=mean(dBpow{2,reg,condi}(:,1:24),2);
        PHASEdeltaWT=mean(phase{1,reg,condi}(:,1:24),2); %for logspace
        PHASEdeltaCKD=mean(phase{2,reg,condi}(:,1:24),2); %for logspace
    
   %POWERthetaWT=mean(dBpow{1,reg,condi}(:,25:40),2);
   %POWERthetaCKD=mean(dBpow{2,reg,condi}(:,25:40),2);
        %PHASEthetaWT=mean(phase{1,reg,condi}(:,25:40),2);
        %PHASEthetaCKD=mean(phase{2,reg,condi}(:,25:40),2);

   % POWERgammaWT=mean(dBpow{1,reg,condi}(:,59:76),2);
    %POWERgammaCKD=mean(dBpow{2,reg,condi}(:,59:76),2);
        %PHASEgammaWT=mean(phase{1,reg,condi}(:,59:76),2);
        %PHASEgammaCKD=mean(phase{2,reg,condi}(:,59:76),2);
%     
%     figure
%     hold on
%         plot(tf_tx,POWERdeltaWT,'k');
%         plot(tf_tx,POWERdeltaCKD,'b');
%     
%     figure
%     hold on
%         plot(tf_tx,POWERthetaWT,'k');
%         plot(tf_tx,POWERthetaCKD,'b');
%     
%     figure 
%     hold on
%         plot(tf_tx,POWERgammaWT,'k');
%         plot(tf_tx,POWERgammaCKD,'b');
    %title('Power'); ylim([-2 3.5]);
    
    figure
    hold on
        plot(tf_tx,PHASEdeltaWT,'k','linewidth',2);
        plot(tf_tx,PHASEdeltaCKD,'b','linewidth',2);
        ylim([.05 .3]);
%     figure 
%     hold on
%         plot(tf_tx,PHASEthetaWT,'k','linewidth',2);
%         plot(tf_tx,PHASEthetaCKD,'b','linewidth',2);
%         ylim([0 .3]);
   % plot(tf_tx,PHASEgamma(:,geno,reg),'g');
    %title('Phase');%legend('delta', 'theta');
    
    disp([ ' Condition' num2str(condi) '; Region' num2str(reg)]);
   
            end
      % end
    %end
end

for sesi=1:5
    for condi=1:2
        POWER_dB=cell2mat(dBpow(sesi,condi))';
        PHASE_CONST=cell2mat(phase(sesi,condi))';
        %PHASEz=cell2mat(zphase(sesi,condi))';
        
    figure
    subplot(1,3,1)
    contourf(tf_tx,frex,POWER_dB(:,:),40,'linecolor','none'); 
        set(gca, 'clim', [-5 5],  'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6))
        title('POWER'); xlabel('Time (ms)');ylabel('Frequency (Hz)');
        cbar; 

    subplot(1,3,2)
    contourf(tf_tx,frex,PHASE_CONST(:,:),40,'linecolor','none'); 
        set(gca,'clim', [0 .3], 'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6))
        title('PHASE'); xlabel('Time (ms)');ylabel('Frequency (Hz)');
        cbar; 
     
   % subplot(1,3,3)
    %contourf(tf_tx,frex,PHASEz(:,:),40,'linecolor','none'); 
        %set(gca, 'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6)); %'clim', [0 30],
        %title('PHASE Z '); xlabel('Time (ms)');ylabel('Frequency (Hz)');
        %cbar; 

    end
end

%% Phase synchrony - Average
 clear BASE POWER POWER_dB BASEavg PHASE_CONST BASEang PHASE_SYNC POWER_CORR RANGLE RIGHTPWRT2T LANGLE LANG LEFPWRT2T LEFTANGLE

 %set some time variables
    tx=(-1000:1:3000)'; %total length of data
    t1=find(tx==-1000); %beginning of trial --> use whole trial to get power to graph
    t2=find(tx==3000); %end of trial
    tf_tx=-1000:1:3000;  %total trial time
    
     %gvtest = inline('n.*(icpcmag*exp((-(val).^2)./(4.*pi./n)).*(sqrt(2./n)))');
    %gvtest = @(icpcmag,n,val) n.*(icpcmag*exp((-(val).^2)./(4.*pi./n)).*(sqrt(2./n)));

Numses=5;
Numtrt=2;
NumTriTypes=2; 


LFP.SYNCtctic100=cell(1,2);
%LFP.SYNCPHASEtctic=cell(1,2);
%LFP.SYNCphaselag=cell(1,2);
n=100;
for ses=1:5
    for trt=1:2
        %for type=1:4
        for permi=1:250
        for condi=1:4
            clear DATA
                    
            if condi==1,         data=LFP.comepochstctic{ses,1,1,trt}';% DS t1
                                    DATA=datasample(data,n,2);%,'replace', false);
            elseif condi==2 ,    data=LFP.comepochstctic{ses,1,2,trt}';% OFC t1
                                    DATA=datasample(data,n,2);%,'replace', false);
                
            elseif condi==3,     data=LFP.comepochstctic{ses,2,1,trt}';% DS t2
                                    DATA=datasample(data,n,2);%,'replace', false);
            elseif condi==4 ,    data=LFP.comepochstctic{ses,2,2,trt}';% OFC t2
                                    DATA=datasample(data,n,2);%,'replace', false);
                
%             elseif condi==5,     DATA=LFP.comepochs{ses,3,1,trt}';% DS t3
%             elseif condi==6 ,    DATA=LFP.comepochs{ses,3,2,trt}';% OFC t3
%                 
%             elseif condi==7,     DATA=LFP.comepochs{ses,4,1,trt}';% DS t4
%             elseif condi==8,     DATA=LFP.comepochs{ses,4,2,trt}';% OFC t4

            end
  
            dims=size(DATA);
            N=size(DATA,1);
            

                for fq=1:num_wavelets
                     w_con=fconv_JFC(reshape(DATA,1,dims(1)*dims(2)),w(fq,:));  %reshape data into one dimension -> all trials are lined up end to end (this makes it faster)
                     w_con=w_con((size(w,2)-1)/2+1:end-(size(w,2)-1)/2); %remove 1/2 of the wavelet length from each of the data b/c of edge effects
                     w_con=reshape(w_con,dims(1),dims(2)); %reform the data back into two dimensions


                   %corrects 
                   if condi==1 
                     SEEDANGLE1(fq,:,:)=angle(w_con(t1:t2,:));
                   elseif condi==2
                     TARGANGLE1(fq,:,:)=angle(w_con(t1:t2,:));
                     PHASE_SYNC1(fq,:)    =  abs(mean(exp(1i*( squeeze(TARGANGLE1(fq,:,:)-SEEDANGLE1(fq,:,:)) )),2)); % phase synch between seed & target
                     %PHASE_SYNC1ang(fq,:) =  angle(mean(exp(1i*( squeeze(TARGANGLE1(fq,:,:)-SEEDANGLE1(fq,:,:)) )),2)); % phase synch between seed & target
                     %PHASE_LAG1(fq,:)= mean(((1000 .* (squeeze(TARGANGLE1(fq,:,:)-SEEDANGLE1(fq,:,:))) ) ./ 2*pi*fq), 2);
                   end
                   
                   %incorrects
                   if condi==3
                     SEEDANGLE2(fq,:,:)=angle(w_con(t1:t2,:));
                   elseif condi==4
                     TARGANGLE2(fq,:,:)=angle(w_con(t1:t2,:));
                     PHASE_SYNC2(fq,:)    =  abs(mean(exp(1i*( squeeze(TARGANGLE2(fq,:,:)-SEEDANGLE2(fq,:,:)) )),2)); % phase synch between seed & target
                     %PHASE_SYNC2ang(fq,:) =  angle(mean(exp(1i*( squeeze(TARGANGLE2(fq,:,:)-SEEDANGLE2(fq,:,:)) )),2));
                     %PHASE_LAG2(fq,:)= mean(((1000 .* (squeeze(TARGANGLE2(fq,:,:)-SEEDANGLE2(fq,:,:))) ) ./ 2*pi*fq), 2);
                   end
                   
%                    if condi==5
%                      SEEDANGLE3(fq,:,:)=angle(w_con(t1:t2,:));
%                    elseif condi==6
%                      TARGANGLE3(fq,:,:)=angle(w_con(t1:t2,:));
%                      PHASE_SYNC3(fq,:) =  abs(mean(exp(1i*( squeeze(TARGANGLE3(fq,:,:)-SEEDANGLE3(fq,:,:)) )),2)); % phase synch between seed & target
%                    end
%                    
%                         
%                    if condi==7
%                      SEEDANGLE4(fq,:,:)=angle(w_con(t1:t2,:));
%                    elseif condi==8
%                      TARGANGLE4(fq,:,:)=angle(w_con(t1:t2,:));
%                      PHASE_SYNC4(fq,:) =  abs(mean(exp(1i*( squeeze(TARGANGLE4(fq,:,:)-SEEDANGLE4(fq,:,:)) )),2)); % phase synch between seed & target
%                    end
%                    
                clear w_con
                end
                

        end
      
       
        PHASE_SYNC(:,:,1,permi)=PHASE_SYNC1;
        PHASE_SYNC(:,:,2,permi)=PHASE_SYNC2;
%         PHASE_SYNC(:,:,3)=PHASE_SYNC3;
%         PHASE_SYNC(:,:,4)=PHASE_SYNC4;
%        
        %PHASE_SYNC_ANG(:,:,1)=PHASE_SYNC1ang;
       % PHASE_SYNC_ANG(:,:,2)=PHASE_SYNC2ang;
        
        %CONN_DIR(:,:,1,permi)=PHASE_LAG1;
        %CONN_DIR(:,:,2,permi)=PHASE_LAG2;
        end
        
     SYNC=LFP.SYNCtctic100;   
     %ANG=LFP.SYNCPHASEtctic;
     %CONN=LFP.SYNCphaselag;
     
      SYNC{ses,trt}=mean(PHASE_SYNC,4); 
      SYNC_SD{ses,trt}=std(PHASE_SYNC,0,4);
     %ANG{ses,trt}=PHASE_SYNC_ANG;
     %CONN{ses,trt}=CONN_DIR; 
     clear CONN_DIR SEEDANGLE* TARGANGLE* PHASE_LAG* PHASE_SYNC*
     
     
     LFP.SYNCtctic100=SYNC; clear SYNC
     %LFP.SYNCPHASEtctic=ANG; clear ANG
     %LFP.SYNCphaselag=CONN; clear CONN
    
      % end
    end
end


control=LFP.SYNCtctic100{5,2}(:,:,1);
exp=LFP.SYNCtctic{5,2}(:,:,1);

figure
contourf(tf_tx,frex,control,40,'linecolor','none'); 
set(gca, 'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6))
%title('L vs. R Hemisphere TC Phase Synchrony '); xlabel('Time (ms)');ylabel('Frequency (Hz)');
colorbar


figure
contourf(tf_tx,frex,exp,40,'linecolor','none'); 
set(gca,'clim', [0 1],  'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6))
title('L vs. R Hemisphere TC Phase Synchrony '); xlabel('Time (ms)');ylabel('Frequency (Hz)');
cbar

% pick out phase-sync frequency bands and then graph
mouNcont=13;
mouNexp=15;
for type=1:1
      for ses=1:5
          %for trt=1:2

        PHASEdeltaWT=mean(LFP.SYNCtctic{ses,1}(1:25,:,type),1); %for logspace
           % PHASEdeltaSDWT=(mean(LFP.SYNCtctic{ses,trt}(1:25,:,type),1)) ./ sqrt(mouNcont);
        PHASEdeltaCKD=mean(LFP.SYNCtctic{ses,2}(1:25,:,type),1); %for logspace
            %PHASEdeltaSDCKD=mean(LFP.SYNCtctic{ses,2}(1:25,:,type),1)  ./ sqrt(mouNexp);
    
        PHASEthetaWT=mean(LFP.SYNCtctic{ses,1}(44:55,:,type),1);
          % PHASEthetaSDWT=mean(LFP.SYNCtctic{ses,1}(44:55,:,type),1) ./ sqrt(mouNcont);
        PHASEthetaCKD=mean(LFP.SYNCtctic{ses,2}(44:55,:,type),1);
          % PHASEthetaSDCKD=mean(LFP.SYNCtctic{ses,2}(44:55,:,type),1)  ./ sqrt(mouNexp);

        %PHASEgammaWT=mean(phase{1,reg,condi}(:,59:76),2);
        %PHASEgammaCKD=mean(phase{2,reg,condi}(:,59:76),2);

%     figure
%     subplot(1,2,1)
%     hold on
%         boundedline(tf_tx,PHASEdeltaWT',PHASEdeltaSDWT','k','alpha');
%         boundedline(tf_tx,PHASEdeltaCKD',PHASEdeltaSDCKD','b','alpha');
%         ylim([0 1.2]);
%         title('Delta')
%     subplot(1,2,2)
%     hold on
%         boundedline(tf_tx,PHASEthetaWT',PHASEthetaSDWT','k','alpha');
%         boundedline(tf_tx,PHASEthetaCKD',PHASEthetaSDCKD','b','alpha');
%         ylim([0 1.2]);
%         title('Alpha')

 figure
    subplot(1,2,1)
    hold on
        plot(tf_tx,PHASEdeltaWT','k','linewidth',2);
        plot(tf_tx,PHASEdeltaCKD','b','linewidth',2);
        ylim([.5 .9]);
        title('Delta')
    subplot(1,2,2)
    hold on
        plot(tf_tx,PHASEthetaWT','k','linewidth',2);
        plot(tf_tx,PHASEthetaCKD','b','linewidth',2);
        ylim([.3 .8]);
        title('beta (low) 10-20hz')

   
   %       end
      end
end



%% Phase Synchrony - Individual mice
 clear PHASE*

 %set some time variables
    tx=(-1000:1:3000)'; %total length of data
    t1=find(tx==-1000); %beginning of trial --> use whole trial to get power to graph
    t2=find(tx==3000); %end of trial
    tf_tx=-1000:1:3000;  %total trial time

Numses=5;
Numtrt=2;
NumTriTypes=2; 

%LFP.SYNCtcticMouse=cell(1,2);
for ses=1:5
    for trt=1:2
        Num=size(LFP.epochstctic,1);
        for mou=1:Num
            if size(LFP.epochstctic{mou,ses,1,1,trt},1)==0 %done based on TC in DS, but OFC will follow same mouse pattern
                continue
            end
        for condi=1:4
            clear DATA        
            if condi==1,         DATA=LFP.epochstctic{mou,ses,1,1,trt}';% DS t1   (mou x ses x type x reg x trt)
            elseif condi==2 ,    DATA=LFP.epochstctic{mou,ses,1,2,trt}';% OFC t1
                
            elseif condi==3,     DATA=LFP.epochstctic{mou,ses,2,1,trt}';% DS t2
            elseif condi==4 ,    DATA=LFP.epochstctic{mou,ses,2,2,trt}';% OFC t2
                
%             elseif condi==5,     DATA=LFP.comepochs{ses,3,1,trt}';% DS t3
%             elseif condi==6 ,    DATA=LFP.comepochs{ses,3,2,trt}';% OFC t3
%                 
%             elseif condi==7,     DATA=LFP.comepochs{ses,4,1,trt}';% DS t4
%             elseif condi==8,     DATA=LFP.comepochs{ses,4,2,trt}';% OFC t4
            end
            
  
            dims=size(DATA);
            N=size(DATA,2);

                for fq=1:num_wavelets
                     w_con=fconv_JFC(reshape(DATA,1,dims(1)*dims(2)),w(fq,:));  %reshape data into one dimension -> all trials are lined up end to end (this makes it faster)
                     w_con=w_con((size(w,2)-1)/2+1:end-(size(w,2)-1)/2); %remove 1/2 of the wavelet length from each of the data b/c of edge effects
                     w_con=reshape(w_con,dims(1),dims(2)); %reform the data back into two dimensions

                   if condi==1 
                     SEEDANGLE1(fq,:,:)=angle(w_con(t1:t2,:));
                   elseif condi==2
                     TARGANGLE1(fq,:,:)=angle(w_con(t1:t2,:));
                     
                     if size(TARGANGLE1,3)==1
                          PHASE_SYNC1(fq,:)    =  abs( mean(exp(1i*( squeeze(TARGANGLE1(fq,:,:)-SEEDANGLE1(fq,:,:))' )),2)); 
                          %PHASE_SYNC1ang(fq,:) =  angle( mean(exp(1i*( squeeze(TARGANGLE1(fq,:,:)-SEEDANGLE1(fq,:,:))' )),2)); %for gv test
                     elseif size(TARGANGLE1,3)>1
                          PHASE_SYNC1(fq,:)    =  abs( mean(exp(1i*( squeeze(TARGANGLE1(fq,:,:)-SEEDANGLE1(fq,:,:)) )),2)); % phase synch between seed & target
                          %PHASE_SYNC1ang(fq,:) =  angle( mean(exp(1i*( squeeze(TARGANGLE1(fq,:,:)-SEEDANGLE1(fq,:,:)) )),2)); %for gv test for similar angles
                     end
                   end
                   
                   
                   if condi==3
                     SEEDANGLE2(fq,:,:)=angle(w_con(t1:t2,:));
                   elseif condi==4
                     TARGANGLE2(fq,:,:)=angle(w_con(t1:t2,:));
                     
                     if size(TARGANGLE2,3)==1 %if only one trial... does strange things
                         PHASE_SYNC2(fq,:)    =  abs(mean(exp(1i*( squeeze(TARGANGLE2(fq,:,:)-SEEDANGLE2(fq,:,:))' )),2)); 
                         %PHASE_SYNC2ang(fq,:) =  angle(mean(exp(1i*( squeeze(TARGANGLE2(fq,:,:)-SEEDANGLE2(fq,:,:))' )),2)); 
                     elseif size(TARGANGLE2,3)>1
                         PHASE_SYNC2(fq,:)    =  abs(mean(exp(1i*( squeeze(TARGANGLE2(fq,:,:)-SEEDANGLE2(fq,:,:)) )),2)); % phase synch between seed & target
                         %PHASE_SYNC2ang(fq,:) =  ang(mean(exp(1i*( squeeze(TARGANGLE2(fq,:,:)-SEEDANGLE2(fq,:,:)) )),2));
                     end
                   end
                   
%                    if condi==5
%                      SEEDANGLE3(fq,:,:)=angle(w_con(t1:t2,:));
%                    elseif condi==6
%                      TARGANGLE3(fq,:,:)=angle(w_con(t1:t2,:));
% 
%                      PHASE_SYNC3(fq,:) =  abs(mean(exp(1i*( squeeze(TARGANGLE3(fq,:,:)-SEEDANGLE3(fq,:,:)) )),2)); % phase synch between seed & target
%                    end
%                    
%                         
%                    if condi==7
%                      SEEDANGLE4(fq,:,:)=angle(w_con(t1:t2,:));
%                    elseif condi==8
%                      TARGANGLE4(fq,:,:)=angle(w_con(t1:t2,:));
% 
%                      PHASE_SYNC4(fq,:) =  abs(mean(exp(1i*( squeeze(TARGANGLE4(fq,:,:)-SEEDANGLE4(fq,:,:)) )),2)); % phase synch between seed & target
%                    end
%                    
                clear w_con
                end
                

        end
        
        
       
        PHASE_SYNC(:,:,1)=PHASE_SYNC1; %correct
        PHASE_SYNC(:,:,2)=PHASE_SYNC2; %incorrects

%         PHASE_SYNC(:,:,3)=PHASE_SYNC3;
%         PHASE_SYNC(:,:,4)=PHASE_SYNC4;
%         
     SYNC=LFP.SYNCtcticMouse;   
     SYNC{mou,ses,trt}=PHASE_SYNC; clear PHASE_SYNC* SEEDANGLE* TARGANGLE*
     LFP.SYNCtcticMouse=SYNC; clear SYNC
    
       end
    end
end



CONT=zeros(80,4001);
for ses=5:5
    for trt=1:1
        for beh=1:1
 for j=1:12
     temp=SYNC{j,ses,trt};  
       if size(temp,2)==0
           continue
       elseif size(temp,2)>0
           temp=temp(:,:,beh);
       end
    
     CONT=cat(3,CONT,temp); 
     
 end
        end
    end
end
 CONT(:,:,1)=[];%clears first set of zeros
 gCONT=mean(CONT,3);
 
% control=LFP.SYNCtctic{5,1}(:,:,1);
% exp=LFP.SYNCtctic{5,2}(:,:,1);

figure
hold on
contourf(tf_tx,frex,gCONT,40,'linecolor','none'); 
set(gca,'clim', [0 1],  'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6))
%contour(tf_tx,frex,SIGBOUND{ses,1},1,'linecolor','k','linewidth',3)
%title('L vs. R Hemisphere TC Phase Synchrony '); xlabel('Time (ms)');ylabel('Frequency (Hz)');
colorbar


figure
contourf(1:3002,frex,SYNCbased{2,1,1,1},40,'linecolor','none'); 
set(gca,'clim', [-.6 .6],  'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6))
%title('L vs. R Hemisphere TC Phase Synchrony '); xlabel('Time (ms)');ylabel('Frequency (Hz)');
colorbar

%% pick out phase frequency bands and then graph
for type=1:1
      for ses=1:5
          for reg=1:2

        PHASEdeltaWT=mean(LFP.phasetrt{ses,type,reg,1}(:,1:25),2); %for logspace
            PHASEdeltaSDWT=(mean(LFP.phaseSDtrt{ses,type,reg,1}(:,1:25),2)) ./ sqrt(size(LFP.comepochstctic{ses,type,reg,1},1));
        PHASEdeltaCKD=mean(LFP.phasetrt{ses,type,reg,2}(:,1:25),2); %for logspace
            PHASEdeltaSDCKD=mean(LFP.phaseSDtrt{ses,type,reg,2}(:,1:25),2)  ./ sqrt(size(LFP.comepochstctic{ses,type,reg,2},1));
    
        PHASEthetaWT=mean(LFP.phasetrt{ses,type,reg,1}(:,26:42),2);
           PHASEthetaSDWT=mean(LFP.phaseSDtrt{ses,type,reg,1}(:,26:42),2)  ./ sqrt(size(LFP.comepochstctic{ses,type,reg,1},1));
        PHASEthetaCKD=mean(LFP.phasetrt{ses,type,reg,2}(:,26:38),2);
           PHASEthetaSDCKD=mean(LFP.phaseSDtrt{ses,type,reg,2}(:,26:42),2)  ./ sqrt(size(LFP.comepochstctic{ses,type,reg,2},1));

        %PHASEgammaWT=mean(phase{1,reg,condi}(:,59:76),2);
        %PHASEgammaCKD=mean(phase{2,reg,condi}(:,59:76),2);

    figure
    subplot(1,2,1)
    hold on
        boundedline(tf_tx,PHASEdeltaWT',PHASEdeltaSDWT','k','alpha');
        boundedline(tf_tx,PHASEdeltaCKD',PHASEdeltaSDCKD','b','alpha');
        ylim([0 .6]);
        title('Delta')
    subplot(1,2,2)
    hold on
        boundedline(tf_tx,PHASEthetaWT',PHASEthetaSDWT','k','alpha');
        boundedline(tf_tx,PHASEthetaCKD',PHASEthetaSDCKD','b','alpha');
        ylim([0 .6]);
        title('Theta')

    
    disp([ 'Session ' num2str(ses) '   Region' num2str(reg)]);
   
          end
      end
end


%bar graph average 
choice=850:1150;
tone=1850:2150;
for ses=1:5
    for reg=1:2
        for type=1:1
            
            phasedeltaWT=mean(LFP.phasetrt{ses,type,reg,1}(:,1:25),2);
                phasedeltaSDWT=mean(LFP.phaseSDtrt{ses,type,reg,1}(:,1:25),2) ./ sqrt(size(LFP.comepochstctic{ses,type,reg,1},1));
            phasethetaWT=mean(LFP.phasetrt{ses,type,reg,1}(:,26:42),2);
                phasethetaSDWT=mean(LFP.phaseSDtrt{ses,type,reg,1}(:,26:42),2) ./ sqrt(size(LFP.comepochstctic{ses,type,reg,1},1));
            phasebetaWT=mean(LFP.phasetrt{ses,type,reg,1}(:,43:62),2);
                phasebetaSDWT=mean(LFP.phaseSDtrt{ses,type,reg,1}(:,43:62),2) ./ sqrt(size(LFP.comepochstctic{ses,type,reg,1},1));
            
            phasedeltaEXP=mean(LFP.phasetrt{ses,type,reg,2}(:,1:25),2);
                phasedeltaSDEXP=mean(LFP.phaseSDtrt{ses,type,reg,2}(:,1:25),2) ./ sqrt(size(LFP.comepochstctic{ses,type,reg,2},1));
            phasethetaEXP=mean(LFP.phasetrt{ses,type,reg,2}(:,26:42),2);
                phasethetaSDEXP=mean(LFP.phaseSDtrt{ses,type,reg,2}(:,26:42),2) ./ sqrt(size(LFP.comepochstctic{ses,type,reg,2},1));
            phasebetaEXP=mean(LFP.phasetrt{ses,type,reg,2}(:,43:62),2);
                phasebetaSDEXP=mean(LFP.phaseSDtrt{ses,type,reg,2}(:,43:62),2) ./ sqrt(size(LFP.comepochstctic{ses,type,reg,2},1));
            
            figure
            subplot(1,2,1) %for time of choice -500:500ms
            %control
            hold on
            bar(1, mean(phasedeltaWT(choice,:),1),'k')
               errorbar(1,mean(phasedeltaWT(choice,:),1),mean(phasedeltaSDWT(choice,:),1),'k')
            bar(2,mean(phasethetaWT(choice,:),1),'k')
                errorbar(2,mean(phasethetaWT(choice,:),1),mean(phasethetaSDWT(choice,:),1),'k')
            bar(3,mean(phasebetaWT(choice,:),1),'k')
                errorbar(3,mean(phasebetaWT(choice,:),1),mean(phasebetaSDWT(choice,:),1),'k')
             %experimental
            bar(5, mean(phasedeltaEXP(choice,:),1),'b')
               errorbar(5,mean(phasedeltaEXP(choice,:),1),mean(phasedeltaSDEXP(choice,:),1),'k')
            bar(6,mean(phasethetaEXP(choice,:),1),'b')
                errorbar(6,mean(phasethetaEXP(choice,:),1),mean(phasethetaSDEXP(choice,:),1),'k')
            bar(7,mean(phasebetaEXP(choice,:),1),'b')
                errorbar(7,mean(phasebetaEXP(choice,:),1),mean(phasebetaSDEXP(choice,:),1),'k')
            set(gca,'ylim', [0 .4]);

            subplot(1,2,2) %for time of tone 500: 1500ms
            hold on
            bar(1,mean(phasedeltaWT(tone,:),1),'k')
                errorbar(1,mean(phasedeltaWT(tone,:),1),mean(phasedeltaSDWT(tone,:),1),'k')
            bar(2,mean(phasethetaWT(tone,:),1),'k')
                errorbar(2,mean(phasethetaWT(tone,:),1),mean(phasethetaSDWT(tone,:),1),'k')
            bar(3,mean(phasebetaWT(tone,:),1),'k')
                errorbar(3,mean(phasebetaWT(tone,:),1),mean(phasebetaSDWT(tone,:),1),'k') 
           %experimental
           bar(5,mean(phasedeltaEXP(tone,:),1),'b')
                errorbar(5,mean(phasedeltaEXP(tone,:),1),mean(phasedeltaSDEXP(tone,:),1),'k')
            bar(6,mean(phasethetaEXP(tone,:),1),'b')
                errorbar(6,mean(phasethetaEXP(tone,:),1),mean(phasethetaSDEXP(tone,:),1),'k')
            bar(7,mean(phasebetaEXP(tone,:),1),'b')
                errorbar(7,mean(phasebetaEXP(tone,:),1),mean(phasebetaSDEXP(tone,:),1),'k') 
            set(gca,'ylim', [0 .4]);

        end
    end
end



%% Finding Sig. Phase corrolation differences btwn controls and
% experimentals in average *phase synchrony* --> this is to then chose
% comparison to make with behavior


%Create the Null-Hypothesis Distribution: shuffle data between genotypes of
%mice

    voxel_pval   = 0.01;
    cluster_pval = 0.05;
    pixel_pval=.05;
    pval=0.05;

    
for ses=1:5
   for beh=1:2
     
%COMPARE BETWEEN controls and mutant mice
      CONT=LFP.SYNCtctic{ses,1}(:,:,beh);
      EXP =LFP.SYNCtctic{ses,2}(:,:,beh);
      n=size(LFP.comepochstctic{ses,beh,1,1},1);
      n2=size(LFP.comepochstctic{ses,beh,1,2},1);
      N=mean([n n2],2);
      
      DIFF=bsxfun(@minus,CONT,EXP); %N=size(DIFF,2);
      nullZ=atanh(DIFF);
      
%       pvals=exp(sqrt(1+4*N+4*(N.^2-(bsxfun(@times,N,DIFF)).^2)-(1+2*N)));
%       pvals=exp(-N * DIFF.^2);
 
       critval=sqrt( -log(pval) ./ N);
       sigpix=logical(abs(nullZ)>=critval); %find regions that are greater than or less than the critical value 
       
       %pixel based correction
%        maxsig=max(nullZ,[],1);
%         pix_max_threshold = prctile(maxsig,100-pixel_pval*100);
%        minsig=min(nullZ,[],1);
%         pix_min_threshold = prctile(minsig,100-pixel_pval*100);
%         
%         sigpix=logical(nullZ>= pix_max_threshold | nullZ<= pix_min_threshold);  
       
             %cluster correction for sigpix...
             clustinfo = bwconncomp(sigpix);
             clust_info = cellfun(@numel,clustinfo.PixelIdxList);
            
                cluster_pval=.05; 
                
              mean_clust_info=mean(clust_info);
              med_clust_info=median(clust_info);
              clust_info(:,find(clust_info < mean_clust_info))=[];
              
             clust_threshold=prctile(clust_info,100-cluster_pval*100);
                clust_info = cellfun(@numel,clustinfo.PixelIdxList);
             whichclusters2remove = find(clust_info<clust_threshold); %identify which clusters to remove
            % remove clusters
            for i=1:length(whichclusters2remove)
                sigpix(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
            end
            
  
%             otherstuff=bwboundaries(sigpix); %find boundaries of sigpix for graphing
%                 boundaries=cell2mat(otherstuff(1,1));
%             boundary(:,1)=tf_tx(boundaries(:,2)); %find those boundaries on the x-y axis of time-freq
%             boundary(:,2)=frex(boundaries(:,1));   
             
        %SIGBOUND{ses,beh}=boundary; clear boundary DIFF nullZ N
        SIGPIX{ses,beh}=sigpix; clear sigpix clust* CONT nullZ N EXP
    
   end
end

    LFP.sigcontvsexp=SIGPIX;

%GRAPH
for ses=1:5
    for beh=1:1
        
        CONT=LFP.SYNCtctic{ses,1}(:,:,beh);
        EXP =LFP.SYNCtctic{ses,2}(:,:,beh);
        sig=LFP.sigcontvsexp{ses,beh};
        
        figure
        subplot(1,2,1)
        contourf(tf_tx,frex,CONT,40,'linecolor','none');
        box on
            %set(gca, 'clim', [0 1],  'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6),'XTickLabel','')
            set(gca, 'clim', [0 1],  'yscale', 'log','YTickLabel', '','XTickLabel','')

            % title('POWER'); xlabel('Time (ms)');ylabel('Frequency (Hz)');

        %figure     
        subplot(1,2,2)
        hold on
        contourf(tf_tx,frex,EXP,40,'linecolor','none'); 
        contour(tf_tx,frex,sig,1,'linecolor','k','linewidth',2)
        box on
        %colorbar
        %set(gca, 'clim', [0 1],  'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6),'XTickLabel','')
        set(gca, 'clim', [0 1],  'yscale', 'log','YTickLabel', '','XTickLabel','')

    disp([ 'Session ' num2str(ses) '; Behavior ' num2str(beh)]);
    
    end
end


   
    %gv test for consistency of phase angles 
    
    
    gvtest = inline('n.*(icpcmag*exp((-(val).^2)./(4.*pi./n)).*(sqrt(2./n)))');
    %gvtest = @(icpcmag,n,val) n.*(icpcmag*exp((-(val).^2)./(4.*pi./n)).*(sqrt(2./n)));
    
% number of simulated data points
N = 10000;

u = zeros(2,numUsims);

for i=1:numUsims
    angleDIFF=;
    % make some noise
   phaseangconst1(fq,:) = N.*(  PHASE_SYNC1(fq,:).*exp((-(angleDIFF(fq,:)).^2)./(4.*pi./N)) .* (sqrt(2./N))  ); %gvtest in cohen 2014 pg353
end

    




%% Phase Coherence (intra-region) Significance

    voxel_pval   = 0.01;
    cluster_pval = 0.05;
    pixel_pval=.05;
    pval=0.05;

    
for ses=1:5
   for beh=1:2
     for reg=1:2
%COMPARE BETWEEN controls and mutant mice
      CONT=LFP.phasetrt{ses,beh,reg,1};
      EXP =LFP.phasetrt{ses,beh,reg,2};
      n=size(LFP.comepochstctic{ses,beh,1,1},1);
      n2=size(LFP.comepochstctic{ses,beh,1,2},1);
      N=mean([n n2],2);
      
      DIFF=bsxfun(@minus,CONT,EXP); %N=size(DIFF,2);
      nullZ=atanh(DIFF);
 
       critval=sqrt( -log(pval) ./ N);
       sigpix=logical(abs(nullZ)>=critval); %find regions that are greater than or less than the critical value 

             %cluster correction for sigpix...
             clustinfo = bwconncomp(sigpix);
             clust_info = cellfun(@numel,clustinfo.PixelIdxList);
            
                cluster_pval=.05; 
                
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

            
        SIGPIX{ses,beh,reg}=sigpix; clear sigpix clust* CONT nullZ N EXP
    
   
     end
   end
end

LFP.phasetcticSIG=SIGPIX;   

%GRAPH
for ses=5:5
    for beh=1:1
        for reg=2:2
        
        CONT=LFP.phasetctic{ses,beh,reg,1};
        EXP =LFP.phasetctic{ses,beh,reg,2};
        sig=LFP.phasetcticSIG{ses,beh,reg};
        
        figure
        hold on
        %subplot(1,2,1)
        contourf(tf_tx,frex,CONT',40,'linecolor','none');
        box on
        %set(gca,'clim', [0 .6], 'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6))
        set(gca, 'clim', [0 .6],  'yscale', 'log','YTickLabel', '','XTickLabel','')
        
       figure
        %subplot(1,2,2)
        hold on
        contourf(tf_tx,frex,EXP',40,'linecolor','none'); 
        contour(tf_tx,frex,sig',1,'linecolor','k','linewidth',2)
        box on
         % set(gca,'clim', [0 .6], 'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6))
          set(gca, 'clim', [0 .6],  'yscale', 'log','YTickLabel', '','XTickLabel','')
          %cbar
       
    disp([ 'Session ' num2str(ses) '; Behavior ' num2str(beh)]);
    
        end
    end
end



%% Correlation btwn phase and behavior


%find num first presentation corrects and incorrects...can get from
%combined file, but can do here too b/c data already imported. 
for ses=1:5
    for trt=1:2
    Num=size(SINGLEUNIT.ML(:,ses,trt),1);  %number of mice
    for s=1:Num
        if size(SINGLEUNIT.ML{s,ses,trt},1)==0
            continue
        end
        
%trial numbers for TC and TIC
   code=SINGLEUNIT.ML{s,ses,trt};
    
%Types 1= correct -> correct = first presentation correct
%      2= correct -> incorrect = first presentation incorrect
%      3= incorrect -> incorrect = perseverative correction trial
%      4= incorrect-> correct = correct correction trial 
        for j=1:size(code,1) %use size of combined file ML b/c sometime x-tra timestamps in rec file
            if j==1 && code(j,1)>0
                TYPE(j,1)=1;
            elseif j==1 && code(j,1)==0
                TYPE(j,1)=2;
            elseif j>1
                k=j-1;
                if code(k,1)>0 && code(j,1)> 0
                    TYPE(j,1)=1; % correct -> correct

                elseif code(k,1)>0 && code(j,1)== 0
                     TYPE(j,1)=2; %correct -> incorrect

                elseif code(k,1)==0 && code(j,1)==0
                    TYPE(j,1)=3; %incorrect -> incorrect

                elseif code(k,1)==0 && code(j,1)>0
                    TYPE(j,1)=4; % incorrect -> correct
                end
            end
        end
        
            TotErrors(s,ses,trt)=size(find(code==0),1);
            numTC(s,ses,trt)=size(find(TYPE==1),1);
            numTIC(s,ses,trt)=size(find(TYPE==2),1);
   clear code TYPE

    
    end
    end
end 


%Graph behavior by average theta, delta or beta
    %TC by theta during TC  

choice=850:1150;
tone=1850:2150;

for ses=1:5
    for reg=1:2
        %for trial=1:2
            
%      if trial==1, DATA=numTC; 
%      elseif trial==2, DATA=numTIC; 
%      end
%     
    X=numTC(:,ses,1); %control TC
    X2=numTC(:,ses,2); %experimental TC
    
    X3=TotErrors(:,ses,1); %control TIC
    X4=TotErrors(:,ses,2); %experimental TIC

   Num=size(X,1);
   for mou=1:Num  
       if size(LFP.mousephasetctic{mou,ses,1,reg,1},1)==0 
           Y(mou,1)=0;
       elseif size(LFP.mousephasetctic{mou,ses,1,reg,1},1)~=0 
            temp_Y=mean(LFP.mousephasetctic{mou,ses,1,reg,1}(:,26:42),2); 
            Y(mou,1)=mean(temp_Y(choice,:),1);% control
       end
       
       if size(LFP.mousephasetctic{mou,ses,1,reg,2},1)==0 
           Y2(mou,1)=0;
       elseif size(LFP.mousephasetctic{mou,ses,1,reg,2},1)~=0 
            temp_Y2=mean(LFP.mousephasetctic{mou,ses,1,reg,2}(:,26:42),2); % experimental
            Y2(mou,1)=mean(temp_Y2(choice,:),1);
       end
       
       if size(LFP.mousephasetctic{mou,ses,2,reg,1},1)==0 
           Y3(mou,1)=0;
       elseif size(LFP.mousephasetctic{mou,ses,2,reg,1},1)~=0 
           temp_Y3=mean(LFP.mousephasetctic{mou,ses,2,reg,1}(:,26:42),2); % control
           Y3(mou,1)=mean(temp_Y3(choice,:),1);
       end
       
       if size(LFP.mousephasetctic{mou,ses,2,reg,2},1)==0
           Y4(mou,1)=0;
       elseif size(LFP.mousephasetctic{mou,ses,2,reg,2},1)~=0
           temp_Y4=mean(LFP.mousephasetctic{mou,ses,2,reg,2}(:,26:42),2); % experimental
           Y4(mou,1)=mean(temp_Y4(choice,:),1);
       end
   end

    %Delete any (zero, zero) points... 
        doublezero=find(X==0 & Y==0);
            X(doublezero,:)=[];
            Y(doublezero,:)=[];
        doublezero=find(X2==0 & Y2==0);
            X2(doublezero,:)=[];
            Y2(doublezero,:)=[];
            
        doublezero=find(X3==0 & Y3==0);
            X3(doublezero,:)=[];
            Y3(doublezero,:)=[];
        doublezero=find(X4==0 & Y4==0);
            X4(doublezero,:)=[];
            Y4(doublezero,:)=[];

 
    %best fit lines...
    coeffs = polyfit(X, Y, 1);
    coeffs2= polyfit(X2,Y2,1);
    
    coeffs3 = polyfit(X3, Y3, 1);
    coeffs4= polyfit(X4,Y4,1);
    %Get fitted values
    fittedX = linspace(min(X), max(X), 200);
    fittedX2= linspace(min(X2), max(X2), 200);
    fittedX3 = linspace(min(X3), max(X3), 200);
    fittedX4= linspace(min(X4), max(X4), 200);
    
    fittedY = polyval(coeffs, fittedX);
        slope=coeffs(1); slopeY=num2str(slope); SLOPE(1,1,ses,reg)=slope; %cont x tc x ses x reg
        [R,P]=corrcoef(X,Y); R1=num2str(R(1,2)); P1=num2str(P(1,2)); CORR(1,1,ses,reg)=R(1,2); PVAL(1,1,ses,reg)=P(1,2);
    fittedY2 = polyval(coeffs2, fittedX2);
        slope=coeffs2(1); slopeY2=num2str(slope); SLOPE(2,1,ses,reg)=slope; %exp x tc x ses x reg
        [R,P]=corrcoef(X2,Y2); R2=num2str(R(1,2)); P2=num2str(P(1,2)); CORR(2,1,ses,reg)=R(1,2); PVAL(2,1,ses,reg)=P(1,2);
    fittedY3 = polyval(coeffs3, fittedX3);
        slope=coeffs3(1); slopeY3=num2str(slope); SLOPE(1,2,ses,reg)=slope;%cont x tic x ses x reg
        [R,P]=corrcoef(X3,Y3); R3=num2str(R(1,2)); P3=num2str(P(1,2)); CORR(1,2,ses,reg)=R(1,2); PVAL(1,2,ses,reg)=P(1,2);
    fittedY4 = polyval(coeffs4, fittedX4);
        slope=coeffs4(1); slopeY4=num2str(slope); SLOPE(2,2,ses,reg)=slope; %exp x tic x ses x reg
        [R,P]=corrcoef(X4,Y4); R4=num2str(R(1,2)); P4=num2str(P(1,2)); CORR(2,2,ses,reg)=R(1,2); PVAL(2,2,ses,reg)=P(1,2);
  

    figure
    hold on
    %TC
    scatter(X, Y,70,'filled','k', 'd'); %TC with phasic response, control
        plot(fittedX, fittedY, 'k-', 'LineWidth', 2);
       % text(fittedX(1,200),fittedY(1,200),slopeY);
        %text(fittedX(1,200),fittedY(1,200)-3,R1); %pearson's correlation coefficient
       % text(fittedX(1,200),fittedY(1,200)-7,P1); %P-VALUE OF correlation coefficient
    scatter(X2, Y2,70,'filled','g','d'); %TC with phasic responce, experimental
        plot(fittedX2, fittedY2, 'g-', 'LineWidth', 2);
        %text(fittedX2(1,200),fittedY2(1,200),slopeY2);
       % text(fittedX2(1,200),fittedY2(1,200)-3,R2); %pearson's correlation coefficient
        %text(fittedX2(1,200),fittedY2(1,200)-7,P2); %P-VALUE OF correlation coefficient
        
    %TIC   
    %scatter(X3, Y3,70,'k','LineWidth',2); % TIC with phasic response, control
        %plot(fittedX3, fittedY3, 'k--', 'LineWidth', 2);
        %text(fittedX3(1,200),fittedY3(1,200),slopeY3); %slope of the line
        %text(fittedX3(1,200),fittedY3(1,200)-3,R3); %pearson's correlation coefficient
        %text(fittedX3(1,200),fittedY3(1,200)-7,P3); %P-VALUE OF correlation coefficient
    %scatter(X4, Y4,70,'g','LineWidth',2); % TIC with phasic responce, experimental
        %plot(fittedX4, fittedY4, 'g--', 'LineWidth', 2);
        %text(fittedX4(1,200),fittedY4(1,200),slopeY4);
        %text(fittedX4(1,200),fittedY4(1,200)-3,R4); %pearson's correlation coefficient
        %text(fittedX4(1,200),fittedY4(1,200)-7,P4); %P-VALUE OF correlation coefficient
%         if ses==1
%             set(gca, 'ylim', [0 1],'xlim',[0 30]);
%         elseif ses== 5
%             set(gca,'ylim', [0 1],'xlim',[0 80]);
%         elseif ses== 2|| 3|| 4
%             set(gca,'ylim', [0 1],'xlim',[0 140]); %set(gca,'xlim',[0 30],'ylim',[0 80]);
%         end

        set(gca,'ylim',[0 1], 'xlim', [0 30]);

%     
%     scatter(performancecre, phasicDScre,'filled','c');%TC or TIC with phasic DS response, experimental
%     scatter(performancecre, phasicOFCcre,'filled','m');%TC or TIC with phasic DS response, experimental
    %lsline
    
    title(['Session:' num2str(ses) ' Region' num2str(reg)]);
    
    
       % end
    end
end


%correlation with phase synchrony avg (1 sec pre-choice and 1 sec- post
%choice)
theta=26:42;
beta=43:62;
pre_choice=1:1000;
post_choice=1001:2000;

for ses=1:5
    for reg=1:2
        %for trial=1:2
            
%      if trial==1, DATA=numTC; 
%      elseif trial==2, DATA=numTIC; 
%      end
time=post_choice;
freq=theta;
%     
    X=numTC(:,ses,1); %control TC
    X2=numTC(:,ses,2); %experimental TC
    
    X3=TotErrors(:,ses,1); %control TIC
    X4=TotErrors(:,ses,2); %experimental TIC

   Num=size(X,1);
   for mou=1:Num  
       if size(LFP.mouSYNCtctic{mou,ses,1}(:,:,1),1)==0  %cont, tc
           Y(mou,1)=0;
       elseif size(LFP.mouSYNCtctic{mou,ses,1}(:,:,1),1)~=0 
            temp_Y=mean(LFP.mouSYNCtctic{mou,ses,1}(freq,:,1),1); 
            Y(mou,1)=mean(temp_Y(:,time),2);% control
       end
       
       if size(LFP.mouSYNCtctic{mou,ses,2}(:,:,1),1)==0  %exp, tc
           Y2(mou,1)=0;
       elseif size(LFP.mouSYNCtctic{mou,ses,2}(:,:,1),1)~=0 
            temp_Y2=mean(LFP.mouSYNCtctic{mou,ses,2}(freq,:,1),1); % experimental
            Y2(mou,1)=mean(temp_Y2(:,time),2);
       end
       
       if size(LFP.mouSYNCtctic{mou,ses,1}(:,:,1),1)==0 %cont, tic
           Y3(mou,1)=0;
       elseif size(LFP.mouSYNCtctic{mou,ses,1}(:,:,1),1)~=0 
           temp_Y3=mean(LFP.mouSYNCtctic{mou,ses,1}(freq,:,2),1); % control
           Y3(mou,1)=mean(temp_Y3(:,time),2);
       end
       
       if size(LFP.mouSYNCtctic{mou,ses,2}(:,:,1),1)==0 %exp, tic
           Y4(mou,1)=0;
       elseif size(LFP.mouSYNCtctic{mou,ses,2}(:,:,1),1)~=0
           temp_Y4=mean(LFP.mouSYNCtctic{mou,ses,2}(freq,:,2),1); % experimental
           Y4(mou,1)=mean(temp_Y4(:,time),2);
       end
   end
clear temp*
    %Delete any (zero, zero) points... 
        doublezero=find(X==0 & Y==0);
            X(doublezero,:)=[];
            Y(doublezero,:)=[];
        doublezero=find(X2==0 & Y2==0);
            X2(doublezero,:)=[];
            Y2(doublezero,:)=[];
            
        doublezero=find(X3==0 & Y3==0);
            X3(doublezero,:)=[];
            Y3(doublezero,:)=[];
        doublezero=find(X4==0 & Y4==0);
            X4(doublezero,:)=[];
            Y4(doublezero,:)=[];

 
    %best fit lines...
    coeffs = polyfit(X, Y, 1);
    coeffs2= polyfit(X2,Y2,1);
    
    coeffs3 = polyfit(X3, Y3, 1);
    coeffs4= polyfit(X4,Y4,1);
    %Get fitted values
    fittedX = linspace(min(X), max(X), 200);
    fittedX2= linspace(min(X2), max(X2), 200);
    fittedX3 = linspace(min(X3), max(X3), 200);
    fittedX4= linspace(min(X4), max(X4), 200);
    
    fittedY = polyval(coeffs, fittedX);
        slope=coeffs(1); slopeY=num2str(slope); SLOPE(1,1,ses,reg)=slope; %cont x tc x ses x reg
        [R,P]=corrcoef(X,Y); R1=num2str(R(1,2)); P1=num2str(P(1,2)); CORR(1,1,ses,reg)=R(1,2); PVAL(1,1,ses,reg)=P(1,2);
    fittedY2 = polyval(coeffs2, fittedX2);
        slope=coeffs2(1); slopeY2=num2str(slope); SLOPE(2,1,ses,reg)=slope; %exp x tc x ses x reg
        [R,P]=corrcoef(X2,Y2); R2=num2str(R(1,2)); P2=num2str(P(1,2)); CORR(2,1,ses,reg)=R(1,2); PVAL(2,1,ses,reg)=P(1,2);
    fittedY3 = polyval(coeffs3, fittedX3);
        slope=coeffs3(1); slopeY3=num2str(slope); SLOPE(1,2,ses,reg)=slope;%cont x tic x ses x reg
        [R,P]=corrcoef(X3,Y3); R3=num2str(R(1,2)); P3=num2str(P(1,2)); CORR(1,2,ses,reg)=R(1,2); PVAL(1,2,ses,reg)=P(1,2);
    fittedY4 = polyval(coeffs4, fittedX4);
        slope=coeffs4(1); slopeY4=num2str(slope); SLOPE(2,2,ses,reg)=slope; %exp x tic x ses x reg
        [R,P]=corrcoef(X4,Y4); R4=num2str(R(1,2)); P4=num2str(P(1,2)); CORR(2,2,ses,reg)=R(1,2); PVAL(2,2,ses,reg)=P(1,2);
  

    figure
    hold on
    %TC
    scatter(X, Y,70,'filled','k', 'd'); %TC with phasic response, control
        plot(fittedX, fittedY, 'k-', 'LineWidth', 2);
       % text(fittedX(1,200),fittedY(1,200),slopeY);
        %text(fittedX(1,200),fittedY(1,200)-3,R1); %pearson's correlation coefficient
       % text(fittedX(1,200),fittedY(1,200)-7,P1); %P-VALUE OF correlation coefficient
    scatter(X2, Y2,70,'filled','b','d'); %TC with phasic responce, experimental
        plot(fittedX2, fittedY2, 'b-', 'LineWidth', 2);
        %text(fittedX2(1,200),fittedY2(1,200),slopeY2);
       % text(fittedX2(1,200),fittedY2(1,200)-3,R2); %pearson's correlation coefficient
        %text(fittedX2(1,200),fittedY2(1,200)-7,P2); %P-VALUE OF correlation coefficient
        
    %TIC   
    %scatter(X3, Y3,70,'k','LineWidth',2); % TIC with phasic response, control
        %plot(fittedX3, fittedY3, 'k--', 'LineWidth', 2);
        %text(fittedX3(1,200),fittedY3(1,200),slopeY3); %slope of the line
        %text(fittedX3(1,200),fittedY3(1,200)-3,R3); %pearson's correlation coefficient
        %text(fittedX3(1,200),fittedY3(1,200)-7,P3); %P-VALUE OF correlation coefficient
    %scatter(X4, Y4,70,'g','LineWidth',2); % TIC with phasic responce, experimental
        %plot(fittedX4, fittedY4, 'g--', 'LineWidth', 2);
        %text(fittedX4(1,200),fittedY4(1,200),slopeY4);
        %text(fittedX4(1,200),fittedY4(1,200)-3,R4); %pearson's correlation coefficient
        %text(fittedX4(1,200),fittedY4(1,200)-7,P4); %P-VALUE OF correlation coefficient
%         if ses==1
%             set(gca, 'ylim', [0 1],'xlim',[0 30]);
%         elseif ses== 5
%             set(gca,'ylim', [0 1],'xlim',[0 80]);
%         elseif ses== 2|| 3|| 4
%             set(gca,'ylim', [0 1],'xlim',[0 140]); %set(gca,'xlim',[0 30],'ylim',[0 80]);
%         end

        set(gca,'ylim',[0 1], 'xlim', [0 30]);

%     
%     scatter(performancecre, phasicDScre,'filled','c');%TC or TIC with phasic DS response, experimental
%     scatter(performancecre, phasicOFCcre,'filled','m');%TC or TIC with phasic DS response, experimental
    %lsline
    
    title(['Session:' num2str(ses) ' Region' num2str(reg)]);
    
    
       % end
    end
end





%% Other stuff


%pick out phase sync frequency bands and then graph
for ses=1:1
    for condi=1:2
        %for geno=1:2

        SYNC=DRBR(1).SYNC;
        sync_delta_WT=cell2mat(SYNC(ses,1,condi));%ses,geno,condi
            sync_delta_WT=mean(sync_delta_WT(1:24,:),1);
        sync_delta_CKD=cell2mat(SYNC(ses,2,condi));%ses,geno,condi
            sync_delta_CKD=mean(sync_delta_CKD(1:24,:),1);
        %sync_theta(:,condi,geno,ses)=mean(PHASE_SYNC(25:40,:,condi,geno,ses),1);

        figure
        hold on
        plot(tf_tx,sync_delta_WT,'k','linewidth',2);
        plot(tf_tx,sync_delta_CKD,'b','linewidth',2);
        %plot(tf_tx,sync_theta(:,condi,geno,ses),'g');
           ylim([.45 .85]);
           set(gca,'YTick', .45:.1:.85);

        disp([ 'Session ' num2str(ses) '; Condition ' num2str(condi) ]);  %'; Genotype' num2str(geno)

        %end
    end
end

    

DRBR(1).SYNC=SYNC;
    
    
    save('DRBR_Epochs','DRBR','-v7.3');
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    %% Condi= L or R hemisphere
for condi=1:2
   if condi==1, DATA= SACTICR;
   elseif condi==2, DATA= SACTICL;
   end
        dims=size(DATA);
        N=size(DATA,2);
        for fq=1:num_wavelets
             w_con=fconv_JFC(reshape(DATA,1,dims(1)*dims(2)),w(fq,:));  %reshape data into one dimension -> all trials are lined up end to end (this makes it faster)
             w_con=w_con((size(w,2)-1)/2+1:end-(size(w,2)-1)/2); %remove 1/2 of the wavelet length from each of the data b/c of edge effects
             w_con=reshape(w_con,dims(1),dims(2)); %reform the data back into two dimensions
        
                
                BASE= mean(abs(w_con(b1:b2,:)).^2,2); %baseline power
               
                BASEang=abs(mean(exp(1i*(  angle(w_con(b1:b2,:))  )),2));
        

                POWER(fq,:,condi) = mean(abs(w_con(t1:t2,:)).^2,2); %power over time t1 to t2 for each trial
                PHASE_CONST(fq,:,condi) = abs(mean(exp(1i*(  angle(w_con(t1:t2,:))  )),2));  % Derives to Phase Locking Value (Lachaux et al.)

                %BASE(1,3)=mean(BASE(:,1:2)); BASEang(1,3)=mean(BASEang(:,1:2));
                POWER_dB(fq,:,condi) = 10*(  log10(POWER(fq,:,condi)) - log10(repmat(mean(BASE),1,size(tf_tx,2)))  ); % dB Conversion -->plot this
                        %POWER_dB(fq,:,condi,condj)=10*log10( bsxfun(@rdivide,POWER(fq,:,condi,condj),BASE(:,condi)));
                PHASE_dB(fq,:,condi) = 10*(  log10(PHASE_CONST(fq,:,condi)) - log10(repmat(mean(BASEang),1,size(tf_tx,2)))  ); % dB Conversion
                PHASEz(fq,:,condi)= N*(PHASE_CONST(fq,:,condi).^2);
       
         if condi==1
            SEEDANGLE(fq,:,:)=angle(w_con(t1:t2,:));
            SEEDPWRT2T(fq,:,:)=abs(w_con(t1:t2,:)).^2;
        elseif condi==2
            TARGANGLE(fq,:,:)=angle(w_con(t1:t2,:));
            TARGPWRT2T(fq,:,:)=abs(w_con(t1:t2,:)).^2;
            PHASE_SYNC(fq,:) =  abs(mean(exp(1i*( squeeze(TARGANGLE(fq,:,:)-SEEDANGLE(fq,:,:)) )),2)); % phase synch between seed & target
            for corri=1:size(tf_tx,2) % Power correlation.  Adds a nice completion to the relationship between seed & target
              POWER_CORR(fq,corri) =  corr(squeeze(SEEDPWRT2T(fq,corri,:)),squeeze(TARGPWRT2T(fq,corri,:)),'type','Spearman');
            end
        end
         clear w_con BASE BASEang
        end
     clear DATA N
end






