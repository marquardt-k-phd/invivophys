%% Set the current directory and add Plexon code folders and data 

clear all; clc


cd('E:\Documents\MATLAB')

addpath(genpath('F:\Plexon Inc Documents'));%plexon codes
addpath(genpath('E:\Documents\MATLAB\eeglab_addons'));
addpath(genpath('E:\Documents\MATLAB\kakearney-boundedline-pkg-2112a2b'));
addpath(genpath('E:\Documents\MATLAB\Cavanagh Class (Lec+Scripts)'));
addpath('E:\Documents\MATLAB\DATA_PAE');

DataLoc= 'F:\OFC WT'; %put the folder where your data is located here
addpath(genpath(DataLoc));


%% before you run this script you should have imported plexon file data and saved under a structure 'LFP'
    %load this file now
    %it should contain the following variables: interval avgFP event
    

load PAE_FP_OFC

save('PAE_FP_OFC_ver','LFP','-v7.3')


for trt=1:2
    temp=LFP.avgFP(:,:,trt);
    temp2=LFP.interval(:,:,trt);
    temp3=LFP.event(:,:,:,trt);
    if trt==1
    temp(5,:)=[];
    temp2(5,:)=[];
    temp3(5,:,:)=[];
    
    temp(7,:)=cell(1,5);
    temp2(7,:)=cell(1,5);
    temp3(7,:,:)=cell(1,3,5);
    %elseif trt==2
    end
    
    tempfp(:,:,trt)=temp;
    tempint(:,:,trt)=temp2;
    tempevt(:,:,:,trt)=temp3;
    
    
end

LFP.avgFP=tempfp;
LFP.interval=tempint;
LFP.event=tempevt;
clear temp*




%% EPOCHS --> all four trial types: win-stay, lose-shift, perseverative, regressive
for ses=1:5
    for trt=1:2
    Num=size(var1,2);
   
for j=1:Num
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
            if punish==-1
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
          lfp=cell2mat(LFP.avgFP(j,ses,trt));
          int=cell2mat(LFP.interval(j,ses,trt));

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
        LFPEPOCHS{j,ses,1,trt}=templfp(logical(trialtypes(:,3)==1),:); %lfp bins that are type 1 trials 
        LFPEPOCHS{j,ses,2,trt}=templfp(logical(trialtypes(:,3)==2),:);
        LFPEPOCHS{j,ses,3,trt}=templfp(logical(trialtypes(:,3)==3),:);
        LFPEPOCHS{j,ses,4,trt}=templfp(logical(trialtypes(:,3)==4),:);
             

   clear tone tc punish tic trial* type* x y z w lfp int temp*
end
    end
    clear Num
end

LFP.epochs=LFPEPOCHS; 

% Combine epochs across mice
for ses=1:5
   for trt=1:2
for type=1:4
     Num=size(var1,2);
        
       initmou=zeros(1,4001);
        
        for j=1:Num
            if size(var1{ses,j,trt},1)==0
                continue
            end
           temp=LFP.epochs{j,ses,type,trt};
           initmou=[initmou; temp]; % add all bins from each mouse together
       end
           initmou(all(initmou==0,2),:)=[]; %delete off that first row of all zeros
            
          LFPEPOCHScom{ses,type,trt}= initmou;  
          
          
    clear initmou temp*

end
   end
end

LFP.epochscom=LFPEPOCHScom;

%save('PAE_FP_DS','LFP','-v7.3')%*************************
 %% EPOCHS --> TC and TIC only
Numses=4;
Numtrt=1;

for ses=1:Numses
    for trt=1:Numtrt
       
    Numfiles=size(var1,2);
for j=1:Numfiles
        if size(WTOFCFP.interval{j,ses},1)==0
           continue
        end
        
    %set stuff for correct trials
        tone=cell2mat(WTOFCFP.Evt(j,1,ses,trt));
            tc=bsxfun(@minus,tone,1);
            if tone==-1
                tc=[];
            end


    %set stuff for incorrect trials 
    if j==1 || j==2 || j==3 || j==4 || j==5 || j==6 || j==7 || j==8 || j==9 || j==10 && ses==1 || ses==2 || ses==3
        punish=cell2mat(WTOFCFP.Evt(j,4,ses,trt));
           tic=bsxfun(@minus,punish,10); 
           if punish==-1
              tic=[];
           end
    elseif j==11 || j==12 || j==13 && ses==1 || ses==2 || ses==3
        punish=cell2mat(WTOFCFP.Evt(j,2,ses,trt));
           tic=bsxfun(@minus,punish,10); 
           if punish==-1
              tic=[];
           end
    elseif j==1 || j==2 || j==3 || j==4 || j==5 || j==6 || j==7 || j==8 && ses==4
        punish=cell2mat(WTOFCFP.Evt(j,4,ses,trt));
           tic=bsxfun(@minus,punish,10); 
           if punish==-1
              tic=[];
           end
    elseif j==9 || j==10 || j==11 && ses==4
        punish=cell2mat(WTOFCFP.Evt(j,2,ses,trt));
           tic=bsxfun(@minus,punish,10); 
           if punish==-1
              tic=[];
           end
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
          lfp=cell2mat(WTOFCFP.FP(j,ses,trt));
          int=cell2mat(WTOFCFP.interval(j,ses,trt));

           for b=1:size(trialtypes,1) %if everything is perfect
               temp=find(int <= trials(b,1)+.002 & int >= trials(b,1)-.002); %find where the LFP is closest to the spike in time
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
        LFPEPOCHStctic{j,ses,1,trt}=templfp(logical(trialtypes(:,2)==1),:); %lfp bins that are correct trials  
        LFPEPOCHStctic{j,ses,2,trt}=templfp(logical(trialtypes(:,2)==2),:); %lfp bins that are incorrect trials  
             

   clear tone tc punish tic trial* type* x y z w lfp int temp*
end
    clear Num
    end      
end
 WTOFCLFP.epochstctic=LFPEPOCHStctic; 

% Combine epochs across mice
for ses=1:Numses
    for trt=1:Numtrt

for type=1:2
       Numfiles=size(var1,2);
        
       initmou=zeros(1,4001);
        
       for j=1:Numfiles
            if size(WTOFCFP.interval{j,ses},2)==0
                continue
            end
           temp=WTOFCLFP.epochstctic{j,ses,type,trt};
           initmou=[initmou; temp]; % add all bins from each mouse together
       end
           initmou(all(initmou==0,2),:)=[]; %delete off that first row of all zeros
           
            
          LFPEPOCHScomtctic{ses,type,trt}= initmou;  
          
          
    clear initmou temp*

end
    end
end

WTOFCLFP.comepochstctic=LFPEPOCHScomtctic;

save('WTOFCLFP','WTOFCLFP','-v7.3')%*************************
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
    figure;
    subplot(3,1,1)
    plot(t,real(w(4,:))); 
    subplot(3,1,2)
    plot(t,real(w(10,:))); 
    subplot(3,1,3)
    plot(t,real(w(35,:))); 
    
 


 %% Convlolve wavelet family with epoched data --> COMBINED BASELINE
  %clear  POWER POWER_dB  PHASE PHASEz DATA PHASE_CONST
  %clear BASE BASEP
  
 %set some time variables
    tx=(-1000:1:3000)'; %total length of data
    b1=find(tx==-1000); %beginning of baseline
    b2=find(tx==0); %end of baseline
    t1=find(tx==0); %beginning of trial
    t2=find(tx==3000); %end of trial
    tf_tx=-1000:1:3000;  %total trial time
    
    NumTriTypes=2;
    Numses=5;
    
%first get Combined Baseline
%for trt=1:2
for Permi=1:250
    
    N=26; %*change to smallest number of trials*
  
   ALL=zeros(1,4001);
   for type=1:NumTriTypes
       for ses=1:Numses
           for trt=1:2
               temp=LFP.comepochstctic{ses,type,trt}; %**** make sure epochs from correct place
               temp=datasample(temp,N,1,'replace', false); %chose XX random trials from each ses & type
               ALL=[ALL; temp]; % add all bins from each mouse together
               clear temp
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
%if by seperate trts
%     BASEtrt=mean(BASE,3);
%     LFP.baseTRTtctic=BASEtrt;

%if all combined
    BASEper=mean(BASE,3);
    LFP.basetctic=BASEper;


save('PAE_FP_DS','LFP','-v7.3')%*************************

%%  Convlolve wavelet family with epoched data --> Permutation
    clear POWER PHASE_CONST POWER_dB
    clear dBpow dBpowsd phase phasesd

 %set some time variables
    tx=(-1000:1:3000)'; %total length of data
    t1=find(tx==-1000); %beginning of trial --> use whole trial to get power to graph
    t2=find(tx==3000); %end of trial
    tf_tx=-1000:1:3000;  %total trial time

    NumTriTypes=2;
    Numses=4;
   % BASE=LFP.basetctic;
    

    N=68; %set to lowest trial number of all sessions
for ses=1:Numses %we have 5 sessions
    for trt=1:1
    for type=1:NumTriTypes %four trial types
        data=cell2mat(WTOFCLFP.comepochstctic(ses,type,trt))'; %******************
        
        for Permi=1:500 %tested with 50
            DATA=datasample(data,N,2,'replace', false);
            dims=size(DATA);
            %n=size(DATA,2);
                for fq=1:num_wavelets
                     w_con=fconv_JFC(reshape(DATA,1,dims(1)*dims(2)),w(fq,:));  %reshape data into one dimension -> all trials are lined up end to end (this makes it faster)
                     w_con=w_con((size(w,2)-1)/2+1:end-(size(w,2)-1)/2); %remove 1/2 of the wavelet length from each of the data b/c of edge effects
                     w_con=reshape(w_con,dims(1),dims(2)); %reform the data back into two dimensions

                      %  POWER(:,fq,Permi) = mean(abs(w_con(t1:t2,:)).^2,2); %power over time t1 to t2 for each trial
                      %  POWER_dB(:,fq,Permi) = 10*(  log10(POWER(:,fq,Permi)) - log10(repmat(   mean(BASE(:,fq),1),  size(tf_tx,2),   1) )  );  % dB Conversion -->plot this
                            
                        PHASE_CONST(:,fq,Permi) = abs(mean(exp(1i*(  angle(w_con(t1:t2,:))  )),2));  % Derives to Phase Locking Value (Lachaux et al.)        
                        %PHASEz(:,fq,Permi)= n*(PHASE_CONST(:,fq,Permi).^2);
                end

           clear w_con DATA
        end
        
       %pow{ses,type}= mean(POWER,3); 
        %powsd{ses,type}=std(POWER,[],3);
        clear POWER
        
       %dBpow{ses,type,trt}=mean(POWER_dB,3); 
       % dbpowsd{ses,type,trt}=std(POWER_dB,0,3);
        clear POWER_dB
        
       phase{ses,type,trt}=mean(PHASE_CONST,3); 
        phasesd{ses,type,trt}=std(PHASE_CONST,0,3);
        clear PHASE_CONST

        clear data
    end
    end
end

% dBpow=LFP.dBpower;
% dbpowsd=LFP.dBpowerSD;
% phase=LFP.phase;
% phasesd=LFP.phaseSD;


% save it
%LFP.dBpower=dBpow;
%LFP.dBpowerSD=dbpowsd;
WTOFCLFP.phasetctic=phase;
WTOFCLFP.phaseSDtctic=phasesd;

%save('PAE_FP_DS','LFP','-v7.3')%*************************
 %significance within sessions independently    
    cluster_pval = 0.01;
    pval=0.05;
for ses=1:4
   for beh=1:2
     for reg=1:1
         for trt=1:1
             
      temp=WTOFCLFP.phasetctic{ses,beh}; 
          if size(temp,2)==0;
             continue
           end

          N=1;
          while N<2001
          ranvals=sign(randn(size(temp)));
          randomdata=temp(logical(ranvals==1));
            ranavg=mean(randomdata,1);
            ransd=2*std(randomdata,[],1);
            sig(N,1)=ranavg+ransd;
            N=N+1;
          end
          sig=abs(mean(sig,1)); clear N
        
          sigpix=logical(abs(temp)>= sig);
             
                clear randomdata

%                 
%                 %within session only - determine statistical sig. using crit value
%       DATA=AVG;
%       N=size(LFP.comepochstctic{ses,beh,reg,trt},1); %across trials
% 
%       %nullZ=atanh(DATA);
%       critval=sqrt( -log(pval) ./ N);
%       sigpix=logical(abs(DATA)>=critval); %find regions that are greater than or less than the critical value 

             %cluster correction for sigpix...
             clustinfo = bwconncomp(sigpix);
             clust_info = cellfun(@numel,clustinfo.PixelIdxList);
             
             cluster_pval = 0.01;
 
              mean_clust_info=mean(clust_info);
              med_clust_info=median(clust_info);
              clust_info(:,find(clust_info < med_clust_info))=[];
              
             clust_threshold=prctile(clust_info,100-cluster_pval*100);
                clust_info = cellfun(@numel,clustinfo.PixelIdxList);
             whichclusters2remove = find(clust_info<clust_threshold); %identify which clusters to remove
            
%             % remove clusters
            for i=1:length(whichclusters2remove)
                sigpix(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
            end

        SIGPIX{ses,beh}=sigpix; clear sigpix clust* CONT nullZ N EXP
    
         end
     end
   end
end

WTOFCLFP.phaseSIGsestctic=SIGPIX;


for ses=1:4
    for beh=1:1
   figure
    hold on
      contourf(tf_tx,frex,WTOFCLFP.phasetctic{ses,beh}',40,'linecolor','none');
      contour(tf_tx,frex,WTOFCLFP.phaseSIGsestctic{ses,beh}',1,'linecolor','k','linewidth',3)
      set(gca, 'clim', [0 .4],  'yscale', 'log' ,'YTick', [],'XTick',[])%logspace(log10(start_freq),log10(end_freq),6))
      colorbar
      
    end
end

%% calculate Differences in power
% for sesi=2:4
%     CORR=cell2mat(dBpow(sesi,1));
%     INCORR=cell2mat(dBpow(sesi,2));
%     
%     DIFF{sesi,1}=CORR-INCORR;
%     
% end

% sig change from D85
for type=1:2
    for trt=1:2
    
    d85=cell2mat(LFP.dBpower(1,type,trt));
    
    
    % change from d85
    for ses=2:4
        
        session=cell2mat(LFP.dBpower(ses,type,trt));
        
          diff= d85-session;
          DIFF{ses,type}=diff;
          
          N=1;
          while N<1001
          ranvals=sign(randn(size(session)));
          randomdata=diff(logical(ranvals==1));
            ranavg=mean(randomdata,1);
            ransd=2*std(randomdata,[],1);
            sig(N,1)=ranavg+ransd;
            N=N+1;
          end
          sig=mean(sig,1); clear N
        
          sigpix=logical(diff>= sig);
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
%             otherstuff=bwboundaries(sigpix); %find boundaries of sigpix
%                 boundaries=cell2mat(otherstuff(1,1));
%             boundary(:,1)=tf_tx(boundaries(:,1)); %find those boundaries on the x-y axis of time-freq
%             boundary(:,2)=frex(boundaries(:,2));
            
        SIGBOUND{ses,type,trt}=sigpix;
  clear otherstuff boundary clust* max_clust_info whichclusters2remove sigpix diff boundaries
    end
    end
end

   LFP.sigdBpowtctic=SIGBOUND; %differences from d85
   
%GRAPH
for ses=1:4
    for beh=3:3
        
        CONT=WTOFCLFP.dBpower{ses,beh,1}';
       % EXP =WTOFCLFP.dBpower{ses,beh,2}';
        %sig=LFP.sigdBpowtctic{ses,beh};
        
        figure
        %subplot(1,2,1)
        hold on
        contourf(tf_tx,frex,CONT,40,'linecolor','none');
          colorbar
       % box off
             set(gca, 'clim', [-6 6],  'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6))
            % title('POWER'); xlabel('Time (ms)');ylabel('Frequency (Hz)');

              
%         subplot(1,2,2)
%         hold on
%         contourf(tf_tx,frex,EXP,40,'linecolor','none'); 
%         %contour(tf_tx,frex,sig,1,'linecolor','k','linewidth',2)
%        % box off
%         colorbar
%           set(gca, 'clim', [-6 6],  'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6))

    disp([ 'Session ' num2str(ses) '; Behavior ' num2str(beh)]);
    
    end
end

%% SIGNIFICANCE - DIFFERENCE IN PHASE BETWEEEN CONDITIONS
% sig change from D85

for type=1:4
    for trt=1:1
    for ses=1:4
    
        COMPARISON=cell2mat(WTOFCLFP.phase(ses,type));
        
        for sesi=1:4
        session=cell2mat(WTOFCLFP.phase(sesi,type));
        
          diff= COMPARISON-session;
          DIFF{ses,sesi,type}=diff;
          
          N=1;
          while N<1001
          ranvals=sign(randn(size(session)));
          randomdata=diff(logical(ranvals==1));
            ranavg=mean(randomdata,1);
            ransd=2*std(randomdata,[],1);
            sig(N,1)=ranavg+ransd;
            N=N+1;
          end
          sig=mean(sig,1); clear N
        
          sigpix=logical(diff>= sig);
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
%             otherstuff=bwboundaries(sigpix); %find boundaries of sigpix
%                 boundaries=cell2mat(otherstuff(1,1));
%             boundary(:,1)=tf_tx(boundaries(:,1)); %find those boundaries on the x-y axis of time-freq
%             boundary(:,2)=frex(boundaries(:,2));
            
        SIGBOUND{ses,sesi,type}=sigpix;
  clear otherstuff boundary clust* max_clust_info whichclusters2remove sigpix diff boundaries
        end
    end
    end
end

%GRAPH --> Lines
% Line Graphs 
time=-500:1000;
for ses=1:4
    for type=3:3
    N=250; %number of permutations
    
    PHASEdeltaL=mean(WTOFCLFP.phase{ses,type}(500:2000,1:13),2); %low delta 1-2
        PHdeltaLSEM=mean(WTOFCLFP.phaseSD{ses,type}(500:2000,1:13),2) ./ sqrt(N); %low delta 1-2
        
    PHASEdeltaH=mean(WTOFCLFP.phase{ses,type}(500:2000,14:25),2); %high delta 2-4
        PHdeltaHSEM=mean(WTOFCLFP.phaseSD{ses,type}(500:2000,14:25),2) ./ sqrt(N);
    
   % PHASEtheta=mean(PHASE_CONST(:,25:40),2);
       % PHthetasem=DEVdelta./N;


    figure
    %subplot(1,2,1)
    hold on
    boundedline(time,PHASEdeltaL,PHdeltaLSEM,'k','alpha');
    boundedline(time,PHASEdeltaH,PHdeltaHSEM,'b','alpha');
    set(gca,'ylim',[.22 .42])
%     boundedline(tf_tx,mean(POWtheta,3),POWthetasem,'b','alpha');
%     boundedline(tf_tx,mean(POWgamma,3),POWgammaSD,'g','alpha');
        %title('Power'); ylim([-5 2]);
       % legend('delta', 'theta', 'gamma');
    
%     subplot(1,2,2)
%     hold on
%     boundedline(tf_tx,mean(PHASEdelta,3),PHdeltaSD,'k','alpha');
%     boundedline(tf_tx,mean(PHASEtheta,3),PHthetasem,'b','alpha');
%     boundedline(tf_tx,mean(PHASEgamma,3),PHgammaSD,'g','alpha');
%    % title('Phase'); ylim([.1 .32]);
%     legend('delta', 'theta', 'gamma');
    
  %  disp([  'Condition ' num2str(condi)]);


    end
end


%GRAPH --> Bars
for type=1:4
        
        D85=WTOFCLFP.phase{1,type}';
        RD1=WTOFCLFP.phase{2,type}';
        R50=WTOFCLFP.phase{3,type}';
        R85=WTOFCLFP.phase{4,type}';

        sig_d85vrd1=SIGBOUND{1,2,type}';
        sig_d85vr50=SIGBOUND{1,3,type}';
        sig_d85vr85=SIGBOUND{1,4,type}';
        
        sig_rd1vr50=SIGBOUND{2,3,type}';
        sig_rd1vr85=SIGBOUND{2,4,type}';
        
        sig_r50vr85=SIGBOUND{3,4,type}';
        
        figure
        subplot(1,2,1)
        hold on
        contourf(tf_tx,frex,R50,40,'linecolor','none');
          set(gca, 'clim', [0 .5],  'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6))
          
        subplot(1,2,2)
        hold on
        contourf(tf_tx,frex,R85,40,'linecolor','none');
        contour(tf_tx,frex,sig_r50vr85,1,'linecolor','k','linewidth',2)
          set(gca, 'clim', [0 .5],  'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6))
          colorbar
      
       

   % disp([ 'Session ' num2str(ses) '; Behavior ' num2str(beh)]);
    
end


figure
  subplot(1,4,1)
  contourf(tf_tx,frex,D85,40,'linecolor','none');
  set(gca, 'clim', [0 .5],  'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6))
   subplot(1,4,2)
  contourf(tf_tx,frex,RD1,40,'linecolor','none');
  set(gca, 'clim', [0 .5],  'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6))
   subplot(1,4,3)
  contourf(tf_tx,frex,R50,40,'linecolor','none');
  set(gca, 'clim', [0 .5],  'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6))
   subplot(1,4,4)
  contourf(tf_tx,frex,R85,40,'linecolor','none');
          set(gca, 'clim', [0 .5],  'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6))


%% GRAPHING

for ses=1:4
    for type=4:4
        POWER_dB_cont=cell2mat(WTOFCLFP.dBpower(ses,type,1))';
        PHASE_CONST_cont=cell2mat(WTOFCLFP.phase(ses,type,1))';
        
        %POWER_dB_exp=cell2mat(WTOFCLFP.dBpower(ses,type,2))';
        %PHASE_CONST_exp=cell2mat(WTOFCLFP.phase(ses,type,2))';
        %PHASEz=cell2mat(zphase(sesi,condi))';
        
%         figure
%             subplot(1,2,1)
%             contourf(tf_tx,frex,POWER_dB_cont(:,:),40,'linecolor','none'); 
%                 set(gca, 'clim', [-6 6],  'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6))
%                 title('POWER CONT'); xlabel('Time (ms)');ylabel('Frequency (Hz)');
%                 colorbar; 
% 
%             subplot(1,2,2)
%             contourf(tf_tx,frex,POWER_dB_exp(:,:),40,'linecolor','none'); 
%                 set(gca, 'clim', [-6 6],  'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6))
%                 title('POWER EXP'); xlabel('Time (ms)');ylabel('Frequency (Hz)');
%                 colorbar

    figure
            %subplot(1,2,1)
            contourf(tf_tx,frex,PHASE_CONST_cont(:,:),40,'linecolor','none'); 
               set(gca,'clim', [.1 .5], 'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6))
               title('PHASE'); xlabel('Time (ms)');ylabel('Frequency (Hz)');
               colorbar; 
% 
%             subplot(1,2,2)
%             contourf(tf_tx,frex,PHASE_CONST_exp(:,:),40,'linecolor','none'); 
%                set(gca,'clim', [.1 .5], 'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6))
%                title('PHASE'); xlabel('Time (ms)');ylabel('Frequency (Hz)');
%                colorbar;
     
   % subplot(1,3,3)
   % contourf(tf_tx,frex,PHASEz(:,:),40,'linecolor','none'); 
       % set(gca, 'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6)); %'clim', [0 30],
       % title('PHASE Z '); xlabel('Time (ms)');ylabel('Frequency (Hz)');
       % cbar; 
       
  disp([ 'Session ' num2str(ses) '; Type ' num2str(type)]);
  clear POWER_dB PHASE_CONST
    end
end


for ses=3:4
    for type=1
        
        diff=cell2mat(DIFF(ses,type));
        sig=cell2mat(SIGBOUND(ses,type));
        
        figure
        hold on
        contourf(tf_tx,frex,diff',40,'linecolor','none'); 
        %colormap('parula')
        colorbar
        plot(sig(:,1),sig(:,2),'k','linewidth',2.5);
            set(gca,'clim', [-6 6], 'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6))
            
      
            
    disp([ 'Session ' num2str(ses) '; Type ' num2str(type)]);
    
    end
end



% Line Graphs 
for condi=1:2
    N=12; %number of mice
    POWdelta=mean(POWER_dB(:,1:24,:,condi),2);
        DEVdeltaPOW=std(POWdelta(:,:,:),0,3);
        POWdeltasem=DEVdelta./N;
    PHASEdelta=mean(PHASE_CONST(:,1:24),2);
        DEVdeltaPH=std(POWdelta(:,:,:),0,3);
        PHdeltaSD=DEVdelta./N;
    
    POWtheta=mean(POWER_dB(:,25:40),2);
        DEVthetaPOW=std(POWdelta(:,:,:),0,3);
        POWthetasem=DEVdelta./N;
    PHASEtheta=mean(PHASE_CONST(:,25:40),2);
        DEVthetaPH=std(POWdelta(:,:,:),0,3);
        PHthetasem=DEVdelta./N;
    
    POWgamma=mean(POWER_dB(:,59:76),2);
        DEVgammaPOW=std(POWdelta(:,:,:),0,3);
        POWgammaSD=DEVdelta./N;
    PHASEgamma=mean(PHASE_CONST(:,59:76),2);
        DEVgammaPH=std(POWdelta(:,:,:),0,3);
        PHgammaSD=DEVdelta./N;

    figure
    subplot(1,2,1)
    hold on
    boundedline(tf_tx,mean(POWdelta,3),POWdeltasem,'k','alpha');
    boundedline(tf_tx,mean(POWtheta,3),POWthetasem,'b','alpha');
    boundedline(tf_tx,mean(POWgamma,3),POWgammaSD,'g','alpha');
        %title('Power'); ylim([-5 2]);
        legend('delta', 'theta', 'gamma');
    
    subplot(1,2,2)
    hold on
    boundedline(tf_tx,mean(PHASEdelta,3),PHdeltaSD,'k','alpha');
    boundedline(tf_tx,mean(PHASEtheta,3),PHthetasem,'b','alpha');
    boundedline(tf_tx,mean(PHASEgamma,3),PHgammaSD,'g','alpha');
   % title('Phase'); ylim([.1 .32]);
    legend('delta', 'theta', 'gamma');
    
    disp([  'Condition ' num2str(condi)]);

end

clear POWER_dB PHASE_CONST
%% graph
for sesi=2:4
    for condi=1:2
        POWER_dB=cell2mat(dBpow(sesi,condi))';
        %PHASE_CONST=cell2mat(phase(sesi,condi))';
        %PHASEz=cell2mat(zphase(sesi,condi))';
        
    figure
   % subplot(1,3,1)
    contourf(tf_tx,frex,POWER_dB(:,:),40,'linecolor','none'); 
        set(gca, 'clim', [-5 5],  'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6))
        title('POWER'); xlabel('Time (ms)');ylabel('Frequency (Hz)');
        cbar; 

    %subplot(1,3,2)
   % contourf(tf_tx,frex,PHASE_CONST(:,:),40,'linecolor','none'); 
       % set(gca,'clim', [0 .3], 'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6))
       % title('PHASE'); xlabel('Time (ms)');ylabel('Frequency (Hz)');
       % cbar; 
     
   % subplot(1,3,3)
   % contourf(tf_tx,frex,PHASEz(:,:),40,'linecolor','none'); 
       % set(gca, 'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6)); %'clim', [0 30],
       % title('PHASE Z '); xlabel('Time (ms)');ylabel('Frequency (Hz)');
       % cbar; 
       
  disp([ 'Session ' num2str(sesi) '; Condition ' num2str(condi)]);
    end
end

for sesi=2:4
    DATA=cell2mat(DIFF(sesi,1));
    
   figure
   contourf(tf_tx,frex,DATA',40,'linecolor','none'); 
       set(gca, 'clim', [-5 5],  'yscale', 'log','YTick', logspace(log10(start_freq),log10(end_freq),6))
       title('POWER'); xlabel('Time (ms)');ylabel('Frequency (Hz)');
       cbar; 
        
  % POWER=cell2mat(DIFF(sesi,1));
   
  % POWdelta=mean(POWER(:,1:24),2);
  % POWtheta=mean(POWER(:,25:40),2);
  % POWgamma=mean(POWER(:,59:76),2);
   
   % figure
   % hold on
   % plot(tf_tx,POWdelta,'k');
  %  plot(tf_tx,POWtheta,'b');
  %  plot(tf_tx,POWgamma,'g');
      %  title('Power'); ylim([-2.5 2]) ;%set(gca, 'YTick',-4.5:1.3:2);
       % legend('delta', 'theta', 'gamma');
end

% separate out delta, theta and gamma for analysis
%according to frex: delta(1-4hz)=1-24; theta(4-10Hz)=25-40; gamma(40-80Hz)=59-76;


%plot w/ corr + incorr on same graph
    for sesi=2:4
   % N=Permi;%250 =number of permutations
     SEMPOWtc=cell2mat(SDPOW(sesi,1)).*250/sqrt(52);  SEMPOWtic=cell2mat(SDPOW(sesi,2)).*25/sqrt(52);
     SEPHtc=cell2mat(SDPH(sesi,1)).*250/sqrt(52);  SEPHtic=cell2mat(SDPH(sesi,2)).*250/sqrt(52);
    
     POWtc=cell2mat(dBpow(sesi,1)); POWtic=cell2mat(dBpow(sesi,2));
     PHASEtc=cell2mat(phase(sesi,1)); PHASEtic=cell2mat(phase(sesi,2));
     
    POWdeltatc=mean(POWtc(:,1:24),2); POWdeltatic=mean(POWtic(:,1:24),2);
        POWdeltaSDtc=mean(SEMPOWtc(:,1:24),2);  POWdeltaSDtic=mean(SEMPOWtic(:,1:24),2);  
    PHASEdeltatc=mean(PHASEtc(:,1:24),2);PHASEdeltatic=mean(PHASEtic(:,1:24),2);
        PHdeltaSDtc=mean(SEPHtc(:,1:24),2);PHdeltaSDtic=mean(SEPHtic(:,1:24),2);
    
    POWthetatc=mean(POWtc(:,25:40),2);POWthetatic=mean(POWtic(:,25:40),2);
        POWthetaSDtc=mean(SEMPOWtc(:,25:40),2);  POWthetaSDtic=mean(SEMPOWtic(:,25:40),2); 
    PHASEthetatc=mean(PHASEtc(:,25:40),2);PHASEthetatic=mean(PHASEtic(:,25:40),2);
        PHthetaSDtc=mean(SEPHtc(:,25:40),2);PHthetaSDtic=mean(SEPHtic(:,25:40),2);
    
    POWgammatc=mean(POWtc(:,59:76),2);POWgammatic=mean(POWtic(:,59:76),2);
        POWgammaSDtc=mean(SEMPOWtc(:,59:76),2); POWgammaSDtic=mean(SEMPOWtic(:,59:76),2);  
    PHASEgammatc=mean(PHASEtc(:,59:76),2);PHASEgammatic=mean(PHASEtic(:,59:76),2);
        PHgammaSDtc=mean(SEPHtc(:,59:76),2);PHgammaSDtic=mean(SEPHtic(:,59:76),2);

zeroline=zeros(size(tf_tx,2),1);
    figure
    hold on
    boundedline(tf_tx,POWdeltatc,POWdeltaSDtc,'g','alpha');
    boundedline(tf_tx,POWdeltatic,POWdeltaSDtic,'r','alpha');
    plot(tf_tx,zeroline,'k--'); title('Power delta'); ylim([-4.5 2]) ;set(gca, 'YTick',-4:1:2);
       legend('correct', 'incorrect');
     
   figure
   hold on
   boundedline(tf_tx,POWthetatc,POWthetaSDtc,'g','alpha');
   boundedline(tf_tx,POWthetatic,POWthetaSDtic,'r','alpha');
    plot(tf_tx,zeroline,'k--'); title('Power Theta'); ylim([-4.5 2]) ;set(gca, 'YTick',-4:1:2);
       legend('correct', 'incorrect')
        
    figure
    hold on
   boundedline(tf_tx,POWgammatc,POWgammaSDtc,'g','alpha');
   boundedline(tf_tx,POWgammatic,POWgammaSDtic,'r','alpha');
      plot(tf_tx,zeroline,'k--'); title('Power Gamma'); ylim([-4.5 2]) ;set(gca, 'YTick',-4:1:2);
      legend('correct', 'incorrect')
        
        
   figure
   hold on
   boundedline(tf_tx,PHASEdeltatc,PHdeltaSDtc,'g','alpha');
   boundedline(tf_tx,PHASEdeltatic,PHdeltaSDtic,'r','alpha');
    title('Phase Delta');  ylim([.1 .3]); set(gca, 'YTick',.1:.05:.3);
   % legend('correct', 'incorrect')
    
   figure 
   hold on
   boundedline(tf_tx,PHASEthetatc,PHthetaSDtc,'g','alpha');
   boundedline(tf_tx,PHASEthetatic,PHthetaSDtic,'r','alpha');
    title('Phase Theta'); ylim([.1 .3]); set(gca, 'YTick',.1:.05:.3);
   % legend('correct', 'incorrect')
    
   figure
   hold on
   boundedline(tf_tx,PHASEgammatc,PHgammaSDtc,'g','alpha');
   boundedline(tf_tx,PHASEgammatic,PHgammaSDtic,'R','alpha');
        title('Phase Gamma');  ylim([.1 .3]); set(gca, 'YTick',.1:.05:.3);
      %  legend('correct', 'incorrect')
    
    
        disp([ 'Session ' num2str(sesi)]);
    end

    
for sesi=1:4
    for type=1:1
    
    DELTA_TONE=mean(mean(WTOFCLFP.phase{sesi,type}(1750:2250,13:21),2),1);
        DELTA_TONE_sem=mean(mean(WTOFCLFP.phaseSD{sesi,type}(1750:2250,13:21),2),1) ./ sqrt(12); %sem: num of mice =N
        
   % THETA_TONE=mean(mean(WTOFCLFP.phase{sesi,type}(1750:2250,26:40),2),1);
       % THETA_TONE_sem=mean(mean(WTOFCLFP.phaseSD{sesi,type}(1750:2250,26:40),2),1) ./ sqrt(12);
        
     DELTA_CHOICE=mean(mean(WTOFCLFP.phase{sesi,type}(750:1250,13:21),2),1);
        DELTA_CHOICE_sem=mean(mean(WTOFCLFP.phaseSD{sesi,type}(750:1250,13:21),2),1) ./ sqrt(12);
        
     %THETA_CHOICE=mean(mean(WTOFCLFP.phase{sesi,type}(750:1250,26:40),2),1);
      %  THETA_CHOICE_sem=mean(mean(WTOFCLFP.phaseSD{sesi,type}(750:1250,26:40),2),1) ./ sqrt(12);
    
    


    figure
    hold on
        bar(1,DELTA_TONE)
            errorbar(1,DELTA_TONE,DELTA_TONE_sem)   
       % bar(2,THETA_TONE)
           % errorbar(2,THETA_TONE,THETA_TONE_sem)
            
       bar(3,DELTA_CHOICE)
            errorbar(3,DELTA_CHOICE,DELTA_CHOICE_sem)   
       % bar(4,THETA_CHOICE)
           % errorbar(4,THETA_CHOICE,THETA_CHOICE_sem)
            
    set(gca,'ylim', [0 .45])
    
    disp([ 'Session ' num2str(sesi)]);
    
    end
end



%% Time-Frequency Graph Stastics 
    
%Power Testing - DIFF from baseline

%Create the Null-Hypothesis Distribution: here we shuffle data and baseline
%to test sig. differences from baseline period

    load POWER

    voxel_pval   = 0.01;
    cluster_pval = 0.05;

    baset(1)=1; %
    baset(2)=1001;

    nperms           =250;
    
    power=LFP.dBpower;
    base=LFP.basetctic;
    NumTypes=2;

for sesi=1:5
    for type=1:NumTypes
for condi=1:2
     
 POW=power{sesi,type,condi};
      BASE=base;
      %TRIAL=[BASE; POW(501:8001,:)];
 nTime=size(POW,1);
   
  for permi =1:nperms
   cutpoint = randsample(2:nTime-diff(baset)-2,1);
   permvals(permi,:,:)=10*log10(bsxfun(@rdivide,mean(POW([cutpoint:end 1:cutpoint-1],:),3),mean(BASE,1))); %dB change of random shifted trial
  end

%POW=cell2mat(dBpow(sesi,condi));
%     BASE=BASEP;
%         TRIAL=[BASE; POW(501:8001,:)];
%         nTime=size(TRIAL,1);
%     
%     for permi =1:nperms
%         cutpoint = randsample(2:nTime-diff(baset)-2,1);
%         permvals(permi,:,:)=POW([cutpoint:end 1:cutpoint-1],:);
%         %10*log10(bsxfun(@rdivide,mean(POWER([cutpoint:end 1:cutpoint-1],:,1:250),3),mean(realbase,1))); %dB change of random shifted trial
%     end
    
    perm_vals{sesi,type,condi}=permvals; clear permvals
    
    
end
    end
end

save('PerVals_Cond2','perm_vals','-v7.3');

%cluster based correction
for sesi=4 %this is the loop that may take a while to run
    for condi=2
        permvals=cell2mat(perm_vals(sesi,condi));
    
    for permi=1:nperms
        
      % for cluster correction, apply uncorrected threshold and get maximum cluster sizes
            clstcorr = squeeze((permvals(permi,:,:)-mean(permvals,1)) ./ std(permvals,[],1) );
            clstcorr(abs(clstcorr)<norminv(1-voxel_pval))=0;

            % get number of elements in largest supra-threshold cluster
            clustinfo = bwconncomp(clstcorr);
            max_clust_info(permi) = max([ 0 cellfun(@numel,clustinfo.PixelIdxList) ]); % the zero accounts for empty maps
            % using cellfun here eliminates the need for a slower loop over cells

    end
     maxclust{sesi,condi}=max_clust_info; clear max_clust_info permvals TRIAL
     
    end
end
    save('MaxClust_Cond2','maxclust','-v7.3');
    
for sesi=2:4
    for condi=2
    %define things...    
    realavg=cell2mat(dBpow(sesi,condi));
    permvals=cell2mat(perm_vals(sesi,condi));permvals=permvals(:,502:8502,:);
    max_clust_info=cell2mat(maxclust(sesi,condi));
    
    zmap=(realavg-squeeze(mean(permvals,1))) ./ squeeze(std(permvals,0,1));
    threshmean=realavg;
    threshmean(abs(zmap)<norminv(1-voxel_pval))=0; 

     %apply cluster-level corrected threshold
       zmapthresh = zmap;
        % uncorrected pixel-level threshold
        zmapthresh(abs(zmapthresh)<norminv(1-voxel_pval))=0; 
        % find islands and remove those smaller than cluster size threshold
        clustinfo = bwconncomp(zmapthresh);
        clust_info = cellfun(@numel,clustinfo.PixelIdxList);
            %find max_cluster outliers &remove from percentile clculation
                outliers=find(max_clust_info > (2*std(max_clust_info,0,1)+mean(max_clust_info,1)));
                max_clust_info(outliers)=[];
        clust_threshold = prctile(max_clust_info,100-cluster_pval*100);
        % identify clusters to remove
        whichclusters2remove = find(clust_info<clust_threshold);

        % remove clusters
        for i=1:length(whichclusters2remove)
            zmapthresh(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
        end
        
        sigthres{sesi,condi}=zmapthresh; 
        clear zmapthresh outliers max_clust_info clustinfo clust_info clust_threshold zmap threshmean realavg permvals
        clear whichclusters2remove 
    end
end
save('SigThresh_Cond2','sigthres','-v7.3');

for sesi=2:4
    for condi=2

        dBgraph=cell2mat(dBpow(sesi,condi));
        sigmap=cell2mat(sigthres(sesi,condi)); sigmap=logical(sigmap);
        
        figure
        hold on
            contourf(tf_tx,frex,dBgraph',40,'linecolor','none') %the dB power graph
            contour(tf_tx,frex,sigmap','w','linewidth',1) %the significance outlines
            
            set(gca, 'clim', [-5 5],  'yscale', 'log','YTick', round(logspace(log10(start_freq),log10(end_freq),6)))
            cbar;
        
        disp([  'Condition ' num2str(condi) ';Session ' num2str(sesi)]);

    end
end


