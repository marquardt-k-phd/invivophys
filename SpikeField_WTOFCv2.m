%% Spike-Field cohoerence (FieldTrip ToolBox)

clear all; clc


cd('E:\Documents\MATLAB')
addpath(genpath('D:\Kristin\OFC WT'));
addpath(genpath('F:\Plexon Inc Documents'));%plexon codes
addpath(genpath('E:\Documents\MATLAB'));

DataLoc= 'F:\OFC WT'; %put the folder where your data is located here
addpath(genpath(DataLoc));
%% load files
load R50_PPC1_perm_shorttime
load RD1_PPC1_perm_shorttime
load D85_PPC1_perm_shorttime
load D50_FP_Ep+SpB
load D85_FP_Ep+SpB
load R85_FP_Epoch+SpB %need behaivoral data and TS for LFP epoch

%% Concatenate all trials of one type per session into one matrix
     %For WT OFC there is not conditions to compare - only correct vs. incorrect

Num=11;
for j=1:Num
      T1=cell2mat(Epochs_TIC(j,1));
      T4=cell2mat(Epochs_TC(j,1));
      T2=cell2mat(Epochs_TIC(j,2));
      T5=cell2mat(Epochs_TC(j,2));
      T3=cat(3, T1, T2); 
      T6=cat(3, T4, T5); 
      TIC=mean(T3,3);
      TC=mean(T6,3);
        
      TIC_FP{j,1}=TIC; 
      TC_FP{j,1}=TC;
      clear T1 T2 T3 TIC TC T6 T5 T4
end

TC=[cell2mat(TC_FP(1)) cell2mat(TC_FP(2)) cell2mat(TC_FP(3)) cell2mat(TC_FP(4)) cell2mat(TC_FP(5)) cell2mat(TC_FP(6)) cell2mat(TC_FP(7)) cell2mat(TC_FP(8))...
     cell2mat(TC_FP(9)) cell2mat(TC_FP(10)) cell2mat(TC_FP(11)) ];
clear TC_FP
TIC=[cell2mat(TIC_FP(1)) cell2mat(TIC_FP(2)) cell2mat(TIC_FP(3)) cell2mat(TIC_FP(4)) cell2mat(TIC_FP(5)) cell2mat(TIC_FP(6)) cell2mat(TIC_FP(7)) cell2mat(TIC_FP(8))...
     cell2mat(TIC_FP(9)) cell2mat(TIC_FP(10)) cell2mat(TIC_FP(11))  ];
clear TIC_FP


%% Spike Firing in LFP behavioral Epoch times

%SB= spike firing in binary dimensions: trial x EpochTime x waveforms
     %0 = no spike 1= spike
%STrial = what trial each waveform fired on; dimensions: spikefiring during LFP epoch  x trial#
%*NOTE: In the code below I do not separate out by hemisphere!!!

%Set up the Spike Binary Array
Num=11; %CHANGE TO HOW MANY MICE ARE IN THE SESSION
for j=1:Num %num of files/mice
    for condi=1:2
        if condi==1 %correct trials
           Temp_FP=cell2mat(Epochs_TC_TS(j));   % end + beginning of LFP epoch
           tx=(-1500:1:7000)';
        elseif condi==2 %incorrect trials
           Temp_FP=cell2mat(Epochs_TIC_TS(j));   % end + beginning of LFP epoch
           tx=(-1500:1:10000)';
        end
        
      NumWF=size(Spikes{j},1); %number of waveforms
      NumT=size(Temp_FP,2); %num of trials
      SBi=zeros(1,length(tx));% initilize Spike Binary array
      STriali=zeros(1,length(tx));% initilize Spike Binary Trials array
      
      for w=1:NumWF %num of waveforms per file
          
        Temp_Spike=cell2mat(Spikes{j}(w,1)); %call spikes for correct mouse
        NumS=length(Temp_Spike);

        for k=1:NumS
            for n=1:NumT
                if size(Temp_FP,2)>1
                 Temp_FP_TS=Temp_FP(1,n):.001:Temp_FP(2,n); %make a temporary FP_Trial Timestamp array for each trial
                
                 %find if spike falls w/in LFP epoch time
                a=find(Temp_Spike(k,1)<=Temp_FP(2,n) & Temp_Spike(k,1)>=Temp_FP(1,n));
                    if a==1 %if yes...
                        %use to find at what point in trial spike fires
                        [c, index]=min(abs(Temp_FP_TS-Temp_Spike(k))); %use index as row
                        SBi(n,index,w)=1;  %save w/ waveform number as row and 1 where spike fires in trial
                        %STriali(k,w)=n;
                        STriali(n,index,w)=n;%also save the trial where the spike fired
                    clear a index
                    end   
                elseif size(Temp_FP,2)==1
                    SBi(n,:,w)=0;
                end
            end

        end
        clear Temp_Spike NumS
        
      end
      
       if sum(Temp_FP)==0 & condi==1
            SBi=zeros(1,8501);
            STriali=zeros(1,8501);
       elseif sum(Temp_FP)==0 & condi==2
           SBi=zeros(1,11501);
           sTriali=zeros(1,11501);
       end
       
      %save some stuff
      STrial{j,condi}=STriali; clear STriali
      SB{j,condi}=SBi; clear SBi
 
      clear a c Temp_FP_TS NumT Temp_FP NumWF NumFP 
    end
end
clear j k n w 

%Trial is from 0 to 8501: Measured at rate of 1000Hz so 1000points= 1 sec
    % from 0 to 1500 is pre-choice
    %1500 to 2500 is tone on
    %2500 to 8501 is the 6 sec post tone off 
%New Setup - "sliding" timespots of .5 sec
%-1 to -.5;-.5 to 0;0 to .5; .5 to 1; 1 to 1.5; 1.5 to 2

%Epoching into smaller segments
    %clear Temp SB_Ep SB_EpT
%define Time segments
pre=500:1000;%-1 to -.5
pre2=1000:1500;%-5 to 0
tone=1500:2000;%0 to +.5
tone2=2000:2500;%.5 to 1
post=2500:3000;%1 to 1.5
post2=3000:3500;%1.5 to 2
for j=1:Num
    for condi=1:2
    
    Temp=cell2mat(SB(j,condi));
    TempTri=cell2mat(STrial(j,condi)); 
    NumTrial=size(Temp,1);
    NumWav=size(Temp,3);
    

    if condi==1 & sum(Temp)>0
        Temp=reshape(Temp, [NumTrial.*NumWav 8501]); %**reshape to size of whole trial
        TempTri=reshape(TempTri, [NumTrial.*NumWav 8501]); 
    elseif condi==2 & sum(Temp)>0
        Temp=reshape(Temp, [NumTrial.*NumWav 11501]); %**reshape to size of whole trial
        TempTri=reshape(TempTri, [NumTrial.*NumWav 11501]); 
    end

    %select time segments
     SB_Ep{j,1,condi}=Temp(:,pre); %set for 0 to 1sec tone
     SB_EpT{j,1,condi}=TempTri(:,pre); %for trial counts
     
     SB_Ep{j,2,condi}=Temp(:,pre2); %set for 0 to 1sec tone
     SB_EpT{j,2,condi}=TempTri(:,pre2); %for trial counts
     
     SB_Ep{j,3,condi}=Temp(:,tone); %set for 0 to 1sec tone
     SB_EpT{j,3,condi}=TempTri(:,tone); %for trial counts
     
     SB_Ep{j,4,condi}=Temp(:,tone2); %set for 0 to 1sec tone
     SB_EpT{j,4,condi}=TempTri(:,tone2); %for trial counts
     
     SB_Ep{j,5,condi}=Temp(:,post); %set for 0 to 1sec tone
     SB_EpT{j,5,condi}=TempTri(:,post); %for trial counts
     
     SB_Ep{j,6,condi}=Temp(:,post2); %set for 0 to 1sec tone
     SB_EpT{j,6,condi}=TempTri(:,post2); %for trial counts

    clear NumTrial NumWav Temp TempTri
    
    end
end

%Set LFP trials at time zero for spikes - line them up 
Num=11;
for tim=1:6
    col=1;
    col2=1;
    if tim==1; epoch=pre;
    elseif tim==2; epoch=pre2;
    elseif tim==3; epoch=tone;
    elseif tim==4; epoch=tone2;
    elseif tim==5; epoch=post;
    elseif tim==6; epoch=post2;
    end
    
    for j=1:Num
        for condi=1:2

    Temp_SB=cell2mat(SB_Ep(j,tim,condi));
    NumTrials=size(Temp_SB,1);
    TempTri=cell2mat(SB_EpT(j,tim,condi));
    
        for k=1:NumTrials
            for z=1:size(Temp_SB,2)
                if Temp_SB(k,z)==1
                    if condi==1 %correct
                        LFPSpk{tim,condi}(:,col)= TC(epoch,k); %this will find the LFP epoch based on the spike epochs above
                        LFPSpkT{tim,condi}(:,col)=TempTri(k,z);
                        col=col+1;
                    elseif condi==2 && k > size(TIC,2)
                        LFPSpk{tim,condi}(:,col2)=zeros(1,501);
                    elseif condi==2 && k<=size(TIC,2)
                        LFPSpk{tim,condi}(:,col2)= TIC(epoch,k);  
                        LFPSpkT{tim,condi}(:,col2)=TempTri(k,z);
                        col2=col2+1;
                    end
                end
            end
        end
        
        clear Temp_SB NumTrials TempTri
        end
    end
    clear col col2
end
%clear all zero rows
for tim=1:6
    for condi=1:2
        Temp=cell2mat(LFPSpk(tim,condi));
        TempT=cell2mat(LFPSpkT(tim,condi));
        [col]=find(Temp(:,:)==0);
        
        Temp(~any(Temp,2),:)=[];
        TempT(col,:)=[];
        
        LFPSpk{tim,condi}=Temp;
        LFPSpkT{tim,condi}=TempT;
        clear Temp TempT
    end
end

%clear stuff that isn't important
clear Beh* co* E epoch FP* var*

%clear LFPSpk LFPSpkT
%% Wavelet analysis around each spike w/ -500 to +500 sec around every spike
 clear   ANGLE DATA DATASAMPLE
  

%wavelet filter
srate=1000; %start here and then downstample to 500Hz
numcycles=8; %lower = better temporal resolution; higher = better frequency resolution
num_wavelets=40;
start_freq=1; %epoch length should cover at least 3 cycle of the smallest frequency (ie: 1sec epoch=3Hz lowest frequency)
end_freq=80;
%frex=logspace(log10(start_freq),log10(end_freq),num_wavelets); %this is a scaling which starts out time specific and beocmes more freq specific at higher frex
frex=linspace(start_freq,end_freq, num_wavelets);%linear increment from 1 to 100 for number of wavelets
s=numcycles./(2*pi.*frex); %sigma for gaussian wave
t=-2.5:1/srate:2.499; %% time, from -2 to 2 second in steps of 1/sampling-rate
clear i w
for fq=1:num_wavelets
    w(fq,:) = exp(2*1i*pi*frex(fq).*t) .* exp(-t.^2./(2*s(fq)^2)); %sine by gaussian
    
      % time res: 2sigma_t *1000 for ms
    timeres(fq)=2*s(fq) * 1000;
    % freq res: 2sigma_f
    frexres(fq)=2* (frex(fq) / numcycles);
end

 %set some time variables
    tx=(0:1:500)';
    t1=find(tx==0);
    t2=find(tx==500);
    tf_tx=0:1:3006;

%colvolve
for tim=1:6
    for condi=1:2
        DATA=cell2mat(LFPSpk(tim,condi)); 
        dims=size(DATA);
        N=size(DATA,2);
        for fq=1:num_wavelets
                w_con=fconv_JFC(reshape(DATA,1,dims(1)*dims(2)),w(fq,:));  %reshape data into one dimension -> all trials are lined up end to end (this makes it faster)
                w_con=w_con((size(w,2)-1)/2+1:end-(size(w,2)-1)/2); %remove 1/2 of the wavelet length from each of the data b/c of edge effects
                w_con=reshape(w_con,dims(1),dims(2)); %reform the data back into two dimensions

                for t=1:N
                   ANG(fq,:,t)=angle(w_con(t1:t2,t));
                end
        end
        ANGLE{tim,condi}=ANG; clear ANG
     
  clear w_con DATA dims N
    end
end

clear ANGLE ANG

%% PPC

% Here's your number of pairwise comparisons
N=25;
Iterations=.5*N*(N-1);  disp(['You will have ',num2str(Iterations),' iterations']);    
    %ITPC
    ITPC=abs(mean(exp(1i*( ANGLE(:,:,1:N) )),3)) ;


%with permutations
nperms=200;
N=25;
for condi=1:2
    for tim=1:6
for permi=1:nperms
    
    [TESTANG idx]=datasample(cell2mat(ANGLE(tim,condi)),N,3, 'replace',false);
    Trial=cell2mat(LFPSpkT(tim,condi));
        Trial_ID=Trial(1,idx); clear Trial idx
    
    row=1;
    for ai=1:N
        for bi=1:N
            Trial(:,:,ai,bi) = exp(1i*(  TESTANG(:,:,ai)-TESTANG(:,:,bi)  )) ;
            Trial_Matrix(ai,bi) = Trial_ID(ai)==Trial_ID(bi);   %  matrix of logical question: is the comparison in the same trial?
            row=row+1;
        end
    end
    
    for fi=1:size(TESTANG,1)
        for ti=1:size(TESTANG,2)
          TEMP = squeeze(Trial(fi,ti,:,:)) ;
          Lu=tril(TEMP)~=0;  % Logical mask to remove upper
          Ld=tril(TEMP)~=1;  % Logical mask to remove diagonal
          L0=logical(Lu.*Ld);% Logical mask - for ONLY lower triangle

          % PPC_0 from Vinck et al. 2010, NeuroImage
          %TEMP_PPC0=TEMP(L0);    % Vector based on lower triangle mask
          %PPC_0(fi,ti,permi,tim,condi)=abs(mean(TEMP_PPC0)); % Average over all pairs, take modulus
            
          % PPC_1 from Vinck et al. 2011, J Comput Neurosci
          Lt=tril(Trial_Matrix)~=1; % Logical mask to remove same trials
          L1=logical(Lu.*Ld.*Lt); 
          TEMP_PPC1=TEMP(L1);    % Vector based on this mask
          PPC_1(fi,ti,permi,tim,condi)=abs(mean(TEMP_PPC1)); % Average over all pairs, take modulus
            
          % PPC_2 from Vinck et al. 2011, J Comput Neurosci
          %for triali=1:max(Trial_ID)-1
              %L_ti=repmat(Trial_ID==triali,N,1);  % Mask for trial ti (row-based expanded to matrix)
              %L2=logical(L1.*L_ti);
              %TEMP_PPC2(triali)=abs(mean(TEMP(L2))); % Average within that trial
              %TEMP_PPC2(isnan(TEMP_PPC2))=0; %change nans to zero (occur b/c not all trials have spike)
              %clear L_ti L2;
          %end
            %PPC_2(fi,ti,permi,tim,condi)=abs(mean(TEMP_PPC2)); % Average over all trial-averaged pairs, take modulus
              
          for triali=1:max(Trial_ID)-1
              L_ti=repmat(Trial_ID==triali,N,1);  % Mask for trial ti (row-based expanded to matrix)
              L2=logical(L1.*L_ti);
              TEMP_PPC2(triali)=abs(mean(TEMP(L2))); % Average within that trial
              %clear L_ti L2;
          end
          PPC_2(fi,ti)=abs(mean(TEMP_PPC2)); % Average over all trial-averaged pairs, take modulus

            
            
            
          % Get sizes
          %SIZES(1,permi,tim,condi)=length(TEMP_PPC0);
          SIZES(2,permi,tim,condi)=length(TEMP_PPC1);
          SIZES(3,permi,tim,condi)=0;

          clear TEMP* Lu Ld L0 Lt L1 L_ti L2;
        end
    end
     clear ai bi TESTANG THEGOODS Trial_ID Trial_Matrix row 
end
    end
end

%clear PPC*
%NEXT TIME GRAPH W/ PERMUTATIONS AS SEM
%find standard error of the mean
  n=N;
    
PPC_1avg=mean(PPC_1,3);

PPC_1_SD=std(PPC_1,0,3);
   %PPC_1_SDavg=mean(PPC_1_SD,2);
PPC_1_SEM=PPC_1_SD/sqrt(200);%sem between permutations 
PPC_1graph=mean(PPC_1avg,2);

%smoothing function - smooth across frequency
for tim=1:6
    for condi=1:2
    PPCsmoo=PPC_1graph(:,1,1,tim,condi);
    PPCsmoothed=smooth(PPCsmoo); 
     PPCsm(:,tim,condi)=PPCsmoothed; %for line graphs
     clear PPC1grsm PPCsmoothed PPCsmoo
     
    PPCsmoo=PPC_1avg(:,:,1,tim,condi);
    PPCsmoothed=smooth(PPCsmoo); 
    PPC1grsm=reshape(PPCsmoothed,40,501);
     PPCallsm(:,:,tim,condi)=PPC1grsm;
     clear PPC1grsm PPCsmoothed PPCsmoo
     
    end
end

xlswrite('S-FC_D85_smoothed',PPCsm(:,:,1),1);
% 
% figure
% hold on
% plot(frex,PPCsm(:,1,2),'k'); 
% plot(frex,PPCsm(:,2,2),'b'); 
% plot(frex,PPCsm(:,3,2),'g'); 
% plot(frex,PPCsm(:,4,2),'r'); 
% plot(frex,PPCsm(:,5,2),'y'); 
% plot(frex,PPCsm(:,6,2),'m'); 
% title('PPC1- Incorrect'); xlabel('FREQUENCY (HZ)');ylabel('PPC');
% 
% figure
% hold on
% boundedline(frex,PPCsm(:,1,1),PPC_1_SEM(:,:,:,1,1),'k','alpha');
% boundedline(frex,PPCsm(:,2,1),PPC_1_SEM(:,:,:,2,1),'b','alpha');
% boundedline(frex,PPCsm(:,3,1),PPC_1_SEM(:,:,:,3,1),'g','alpha');
% boundedline(frex,PPCsm(:,4,1),PPC_1_SEM(:,:,:,4,1),'r','alpha');
% boundedline(frex,PPCsm(:,5,1),PPC_1_SEM(:,:,:,5,1),'y','alpha');
% boundedline(frex,PPCsm(:,6,1),PPC_1_SEM(:,:,:,6,1),'m','alpha');
% 
% figure
% hold on
% boundedline(frex,mean(PPCsm(:,1:2,2),2),mean(PPC_1_SEM(:,:,:,1:2,2),4),'k','alpha');
% boundedline(frex,mean(PPCsm(:,3:4,2),2),mean(PPC_1_SEM(:,:,:,3:4,2),4),'b','alpha');
% boundedline(frex,mean(PPCsm(:,5:6,2),2),mean(PPC_1_SEM(:,:,:,5:6,2),4),'g','alpha');
% 
% figure
% hold on
% bar(frex,mean(PPCsm(:,1:2,1),2),'k');
% errorbar(frex,mean(PPCsm(:,1:2,1),2),mean(PPC_1_SEM(:,1,1,1:2,1),4),'k');
% 
% 
% bar(frex,mean(PPCsm(:,3:4,1),2),'b');
% bar(frex,mean(PPCsm(:,5:6,1),2),'g');
% 
% %GRAPH IN A TIMELINE OVER TEH SIX 'EPOCHS'
% PPC1=RANPPC(:,:,:,2);
% %PPC1=PPCallsm(:,:,:,1);
% PPC1=reshape(PPC1,[40 3006]);
% figure
% imagesc(tf_tx,frex,PPC1);
%     set(gca, 'clim', [0 .1])
%     title('SPIKE-LFP: CORRECT'); xlabel('Time (ms)');ylabel('Frequency (Hz)');
%     cbar; 

    
%find difference from "random" permutation -> better hilight change
                            %first need to identify "random" within each condition - shuffle
                            %across time (independent variable)
                            for condi=1:2
                                N=1;
                                while N < 501

                                    msize=size(PPC_1avg(:,:,:,:,condi),2);
                                    w=1;
                                    while w < 41
                                        idx(w,:)=randperm(msize)';
                                        w=w+1;
                                    end

                                    for z=1:40
                                        for b=1:501
                                            col=idx(z,b);
                                            ran(z,b,:,:)=PPC_1avg(z,col,1,:,condi);
                                        end
                                    end

                                   RAND(:,:,:,condi,N)=ran; clear ran

                                    clear msize idx w col
                                N=N+1;
                                end
                            end
                            %average between permutations?- use this to calculate differences 
                            RANPPC=mean(RAND,5); 
                            
PPC_1graph=squeeze(PPC_1graph);%remove singleton dimensions
PPC_1avg=squeeze(PPC_1avg); size(PPC_1avg)

%for differences w/in each condition and time
for condi=1:2
    for tim=1:6
    N=1;
    while N < 501
        
        msize=size(PPC_1graph,1); %Num of frex
        idx=randperm(msize);
     
        for z=1:40
             row=idx(1,z);
             ran(z,:)=PPC_1avg(row,:,tim,condi);%take out random frquencies at each condi and time
        end
        
       RAND(:,:,tim,condi,N)=ran ;
        
        clear msize idx  crow ran
        
    N=N+1;
    end
    %RANDSD=std(mean(RAND(:,:,tim,condi,:),4),0,2); %difference between all fq at each ms timepoint of the 500ms window
    %RANDSEM(:,:,condi,:)=RANDSD./sqrt(6); %take the SEM using 6 timepoints as N
    end
end
RANPPC=mean(RAND,5); 


%for random data (combined conditions) -> for use in calculating
%significance
    N=1;
    while N < 501
        
        msize=size(PPC_1graph,1); %Num of frex
        idx=randperm(msize);
     
        for z=1:40
             row=idx(1,z);
             ran(z,:,:,:)=PPC_1avg(row,:,:,:);%take out random frquencies at each condi and time
        end
        
       RANDzo(:,:,:,N)=mean(ran,4);%average between conditions
        
        clear msize idx  crow ran
        
    N=N+1;
    end
    RANDSDzo=std(mean(RANDzo(:,:,:,:),4),0,2); %difference between all fq at each ms timepoint of the 500ms window
    RANDSEMzo(:,:)=RANDSDzo./sqrt(6); %take the SEM using 6 timepoints as N
    RAND2SD=bsxfun(@times, squeeze(RANDSDzo), 2);
%RANzo=mean(RAND,5); 

clear ran RAND RANDSD RANDSEM

%calculate differences - 
for condi=1:2    
    for tim=1:6
   PPCDIFF(:,tim,condi)=bsxfun(@minus, PPCsm(:,tim,condi), mean(RANPPC(:,:,tim,condi),2));
    end
end

% 
% zeroline=zeros(40,1);
% figure
% hold on
% plot(frex,mean(PPCDiffzZ(:,3:4,2),2),'k'); 
% plot(frex,PPCDIFF(:,2,2),'b'); 
% plot(frex,PPCDIFF(:,3,2),'g'); 
% plot(frex,PPCDIFF(:,4,2),'r'); 
% plot(frex,PPCDIFF(:,5,2),'y'); 
% plot(frex,PPCDIFF(:,6,2),'m'); 
% plot(frex,zeroline, 'k-', 'linewidth', 2);
% set(gca,'ylim',[0 .1]);
% 
% for condi=1:2
%     figure 
%     hold on
% 
%     %bar(frex,mean(PPCDIFF(:,1:2,condi),2),'k'); 
%        % errorbar(frex,mean(PPCDIFF(:,1:2,condi),2),mean(PPC_1_SEM(:,1,1,1:2,condi),4),'k');
%     %bar(frex,mean(PPCDIFF(:,3:4,condi),2),'b'); 
%         %errorbar(frex,mean(PPCDIFF(:,3:4,condi),2),mean(PPC_1_SEM(:,1,1,3:4,condi),4),'k');
%     bar(frex,mean(PPCDIFF(:,5:6,condi),2),'g');
%         errorbar(frex,mean(PPCDIFF(:,5:6,condi),2),mean(PPC_1_SEM(:,1,1,5:6,condi),4),'k');
% 
%         
%    set(gca, 'ylim', [-.03 .035]);
%    disp(['Condition' num2str(condi)]);
% end

% for condi=1:2
%     figure
%     hold on
%         boundedline(frex(1,3:40),mean(PPCDIFF(3:40,1:2,condi),2),mean(PPC_1_SEM(3:40,1,1,1:2,condi),4),'k','alpha');
%         boundedline(frex(1,3:40),mean(PPCDIFF(3:40,3:4,condi),2),mean(PPC_1_SEM(3:40,1,1,3:4,condi),4),'b','alpha');
%         boundedline(frex(1,3:40),mean(PPCDIFF(3:40,5:6,condi),2),mean(PPC_1_SEM(3:40,1,1,5:6,condi),4),'g','alpha');
%         boundedline(frex(1,3:40),zeros(1,38),mean(RANDSEMzo(3:40,5:6),2),'y','alpha');
%         %plot(frex,zeros(1,40));
%     set(gca, 'ylim', [-.03 .035]);
%     disp(['Condition' num2str(condi)]);
% end


for condi=1:2
    figure
    hold on
        boundedline(frex(1,3:40),mean(PPCDIFF(3:40,1:2,condi),2),mean(PPC_1_SEM(3:40,1,1,1:2,condi),4),'k','alpha');
        boundedline(frex(1,3:40),mean(PPCDIFF(3:40,3:4,condi),2),mean(PPC_1_SEM(3:40,1,1,3:4,condi),4),'b','alpha');
        boundedline(frex(1,3:40),mean(PPCDIFF(3:40,5:6,condi),2),mean(PPC_1_SEM(3:40,1,1,5:6,condi),4),'g','alpha');
        %boundedline(frex(1,3:40),zeros(1,38),mean(RANDSEMzo(3:40,5:6),2),'y','alpha');
        %plot(frex(1,3:40),zeros(1,38));
        boundedline(frex(1,3:40),zeros(1,38),mean(RAND2SD(3:40,1:2),2),'y','alpha');
        %plot(frex(1,3:40),mean(RAND2SD(3:40,1:2),2),'k--');
        %plot(frex(1,3:40),-mean(RAND2SD(3:40,1:2),2),'k--');
    set(gca, 'ylim', [-.03 .035]);
    disp(['Condition' num2str(condi)]);
end


    figure
    hold on
        boundedline(frex(1,3:40),mean(PPCDIFF(3:40,1:2,1),2),mean(PPC_1_SEM(3:40,1,1,1:2,1),4),'b','alpha');
        boundedline(frex(1,3:40),mean(PPCDIFF(3:40,1:2,2),2),mean(PPC_1_SEM(3:40,1,1,1:2,2),4),'r','alpha');
        boundedline(frex(1,3:40),zeros(1,38),mean(RAND2SD(3:40,1:2),2),'k','alpha');
        set(gca, 'ylim', [-.03 .04]);
  
    figure
    hold on
        boundedline(frex(1,3:40),mean(PPCDIFF(3:40,3:4,1),2),mean(PPC_1_SEM(3:40,1,1,3:4,1),4),'b','alpha');
        boundedline(frex(1,3:40),mean(PPCDIFF(3:40,3:4,2),2),mean(PPC_1_SEM(3:40,1,1,3:4,2),4),'r','alpha');
        boundedline(frex(1,3:40),zeros(1,38),mean(RAND2SD(3:40,3:4),2),'k','alpha');
        set(gca, 'ylim', [-.03 .04]);
        
    figure
    hold on
        boundedline(frex(1,3:40),mean(PPCDIFF(3:40,5:6,1),2),mean(PPC_1_SEM(3:40,1,1,5:6,1),4),'b','alpha');
        boundedline(frex(1,3:40),mean(PPCDIFF(3:40,5:6,2),2),mean(PPC_1_SEM(3:40,1,1,5:6,2),4),'r','alpha');
        boundedline(frex(1,3:40),zeros(1,38),mean(RAND2SD(3:40,5:6),2),'k','alpha');
        set(gca, 'ylim', [-.03 .035]);



% 
% 
% figure 
% hold on
% plot(frex,mean(PPCDIFF(:,1:2,2),2),'k'); 
% plot(frex,mean(PPCDIFF(:,3:4,2),2),'b'); 
% plot(frex,mean(PPCDIFF(:,5:6,2),2),'g'); 
% plot(frex,zeroline, 'k-', 'linewidth', 2);
% 
% 
% PPC1a=mean(PPCDIFFall(:,:,1:2,1),3);
% PPC1b=mean(PPCDIFFall(:,:,3:4,1),3);
% PPC1c=mean(PPCDIFFall(:,:,5:6,1),3);
% PPC1=[PPC1a PPC1b PPC1c];
% figure
% imagesc(1:1503,frex,PPC1);
%     set(gca, 'clim', [-.02 .07])
%     title('SPIKE-LFP: INCORRECT'); xlabel('Time (ms)');ylabel('Frequency (Hz)');
%     cbar; 
% 




%% origional PPC

N=size(ANGLE,3)
Iterations=.5*N*(N-1);
% Damn!

clear Trial THEGOODS 


%data in frex x time x trials
N=100;
Iterations=.5*N*(N-1);
for ai=1:N
    for bi=1:N

            Trial(:,:,ai,bi) = exp(1i*(  ANGLE(:,:,ai)-ANGLE(:,:,bi)  )) ;
            %Trial(:,:,ai,bi) = cos(  ANGLE_TEST(:,:,ai)-ANGLE_TEST(:,:,bi)  ) ;

    end
    clear ai bi
end

size(Trial)
for fi=1:80
    for ti=1:size(Trial,2)
      TEMP = squeeze(Trial(fi,ti,:,:)) ;
      L1=tril(TEMP)~=0;
      L2=tril(TEMP)~=1;
      L3=logical(L1.*L2);
      TEMP2=TEMP(L3);
      THEGOODS(1,fi,2)=abs(mean(TEMP2));  % Average over all pairs, take modulus %save in 3rd dimension what mouse it is
      clear TEMP TEMP2 L1 L2 L3
    end
end

figure
rose(ANGLE(2,2,:));

figure
hold on
plot(frex,mean(THEGOODS,3)); set(gca,'XScale', 'log');
title('PPC-MAYBE'); xlabel('FREQUENCY (HZ)');ylabel('PPC');


%how can I get something with PPC by time of certain freq....

%% Random other stuff...
%Pre-processing

% filter out 60 Hz noise from LFP 

[60-1*(1/10):(1/10):60+1*(1/10) ];

