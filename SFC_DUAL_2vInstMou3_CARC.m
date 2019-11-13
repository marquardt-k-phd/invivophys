% Set up directories
disp('LOG: Opening files...')

cd('/xena/scratch/kmarquar/MATLAB')

addpath('/xena/scratch/kmarquar/MATLAB/DATA_PAE') %data location
addpath(genpath('/xena/scratch/kmarquar/MATLAB')) %add all folders and subfolders in cd

addpath(genpath('/xena/scratch/kmarquar/MATLAB/Plexon_Inc_Documents'));%plexon codes
addpath(genpath('/xena/scratch/kmarquar/MATLAB/eeglab_addons'));

%--> for use on MAC only
%addpath(genpath('/Users/Kristin/Documents/MATLAB/eeglab_addons'))
%addpath(genpath('/Users/Kristin/Documents/MATLAB/CARC')) 

%profile on

        %load the data files will need
        load SFC_PAEDUAL_inst %will have spikes, interval and FP epcochs saved in per mouse format
        %save('SFC_PAEDUAL_inst','SFC','-v7.3'); 
            disp('LOG: Files opened, paths set up, starting code')
            
            %error(javachk('jvm')) % will give error if java not working...need for parpool
            
            %parpool(15)
            %disp('LOG: Parpool(15) loaded, starting PPC' ); datetime
            

N_mouse=size(SFC.spikes,1); %define how many mice have...
N_bin=4;
N_treatment=2;
N_behavior=2;
N_region=2;

% 1. FIND PHASE ANGLES 

%wavelet filter --> ALTERED FROM ORIGINAL NOW 1-80, STILL LINEAR SPACED
    srate=1000; %start here and then downstample to 500Hz
    numcycles=4; %lower = better temporal resolution; higher = better frequency resolution
    num_wavelets=40;
    start_freq=1; %epoch length should cover at least 3 cycle of the smallest frequency (ie: 1sec epoch=3Hz lowest frequency)
    end_freq=40;
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
    
 disp('LOG: parameters set, starting ANGLE calculation')
 

  for mouse=1:N_mouse
     disp(['Starting SFC for- mouse: ' num2str(mouse) ' Session: ' num2str(ses)]); datetime
             
      
  %convolve data and find wavelets
  %define sessions w/in Matlab call code
  %run each session on a different node in parallel  
      
  warning('off','all') %turn off warnings
  ANGLE=cell(2,2,2); %pre-allocate
    for behavior=1:N_behavior
        for treatment=1:N_treatment
            for region=1:N_region
                %Numfiles=size(SFC.spikes,1); %define how many mice have...
             %for j=1:Numfiles
                if size(SFC.spikes{mouse,ses,region,treatment},1)==0
                   disp(['mouse: ' num2str(mouse) ' in session ' num2str(ses) ' does not exist/ have any spikes in regionion:' num2str(region) ' for treatment: ' num2str(treatment)])
                   continue
                end
                
            DATA=SFC.LFPepochs{mouse,ses,behavior,region,treatment}';
            dims=size(DATA);
            ANG=zeros(size(frex,2),dims(1), dims(2));

            for fq=1:num_wavelets
                w_con=fconv_JFC(reshape(DATA,1,dims(1)*dims(2)),w(fq,:));  %reshape data into one dimension -> all trials are lined up end to end (this makes it faster)
                w_con=w_con((size(w,2)-1)/2+1:end-(size(w,2)-1)/2); %remove 1/2 of the wavelet length from each of the data b/c of edge effects
                w_con=reshape(w_con,dims(1),dims(2)); %reform the data back into two dimensions
                         
               ANG(fq,:,:)=angle(w_con); %fq x time x trials
                
            clear w_con
            end
            clear DATA

            ANGLE{behavior,region,treatment}=ANG; clear ANG %no need to save per session --> will define session in ANA_inst file
           
             %end
            end
        end
    end
    warning('on','all') %turn back on warnings, just in case
    %SFC.angles=ANGLE; %saved angles of whole trial LFP - do this to get angles @ fq lower than 4Hz

 disp('LOG: ANGLES calulated, isolating instantaneous phase angles')
 
% Isolate instantaneous phase angles at time of spike occurance 
%split LFPangles, intervals into bins --> use 1 sec bins this time
BIN_INTERVAL=cell(N_bin, N_behavior, N_treatment);
BIN_ANG=cell(N_bin, N_behavior, N_region, N_treatment);
    for treatment=1:N_treatment
        for behavior=1:N_behavior
            for region=1:N_region
                 %Numfiles=size(ANGLE,1);
               %for j=1:Numfiles
        
             epochint=SFC.LFPepochint{mouse,ses,behavior,treatment}; % trials x time
             epochangle=ANGLE{behavior,region,treatment}; % fq x time x trials
              if size(epochangle, 1)==0
                 continue
              end
             
                  for bin=1:N_bin %split into bins
                      if bin==1, time=1:1001;
                      elseif bin==2, time=1001:2001;
                      elseif bin==3, time=2001:3001;
                      elseif bin==4, time=3001:4001;
                      end

                  BIN_INTERVAL{bin,behavior,treatment}(:,:)=epochint(:,time);
                  BIN_ANG{bin,behavior,region,treatment}(:,:,:)=epochangle(:,time,:); %will save bins(ses,j) in final ANA_inst file
                  end
                  
              % end
                
            end
        end
    end


%find times that spike occurs w/in LFP bins & put angle in matrix
SPIKEang=cell(N_bin,N_behavior,N_region,N_treatment);
TriCOUNT=cell(N_bin,N_behavior,N_region,N_treatment);
    for treatment=1:N_treatment
        for region=1:N_region
           % Numfiles=size(ANGLE,1);
            countup=1; %unique identifier for every trial - even between mice
            %for j=1:Numfiles
             spiketimes=SFC.spikes{mouse,ses,region,treatment};
             numwav=size(spiketimes,1);
             for wave=1:numwav
                 stimes=spiketimes{wave,1};
                 
           for behavior=1:N_behavior
                
               for bin=1:N_bin
                   initangles=zeros(40,1);
                   inittri=zeros(1,1);
                   
                    epochint=BIN_INTERVAL{bin,behavior,treatment}; % trials x time
                    epochangle=BIN_ANG{bin,behavior,region,treatment}; % fq x time x trials
                         if size(epochangle, 1)==0 %if no angle data - continue
                            continue
                         end
                         
                    numtri=size(epochint,1);
                    for t=1:numtri 
                        trialint=epochint(t,:);%the interval of a trial epoch
                        %find which spikes fall w/in the epoch bin...
                        temp=find(stimes>=trialint(1,1) & stimes<=trialint(1,1001)); %gives row number of matching spike time values
                        
                        if size(temp,1)>0 %if the waveform spikes w/in the trial bin
                            for s=1:size(temp,1)
                                %need to know the closest timepoint in lfp epoch that the spikes occur --> so can take angle
                                temptime=find(trialint>=stimes(temp(s,1),1)-.0005 & trialint<=stimes(temp(s,1),1)+.0005 );
                                if size(temptime,2)>1 %if more than one lfp epoch time value is .0005 away from stime
                                    lfptime(s)=temptime(1,1); %use the first value b/c really probably doesn't matter if both so close
                                elseif size(temptime,2)==1
                                    lfptime(s)=temptime;
                                end
                            end
                            initangles=[initangles epochangle(:,lfptime,t)]; %find LFP angle ** this is what we want to save!!
                            inittri=[inittri repmat(countup,1,size(lfptime,2))]; clear lfptime temptime temp trialint %trial count vector
                        end 
                  countup=countup+1;
                    end
                    
                    initangles(:,~any(initangles,1))=[]; %delete off that first COL of all zeros
                    inittri(:,~any(inittri,1))=[];
                    tempSPIKEang{bin,behavior,wave}=initangles; clear initangles
                    tempTriCOUNT{bin,behavior,wave}=inittri; clear inittri
                end
                clear numtri

           end
                clear stimes 
            
             end
             %once done w/ all waveforms in files - combine across the
             %waveforms -keeping w/in bin and behavior
            for bin=1:N_bin
                for behavior=1:N_behavior
                 initanglesW=zeros(40,1);
                 inittriW=zeros(1,1);
             for wave=1:numwav
                 initangles=tempSPIKEang{bin,behavior,wave};
                 inittri=tempTriCOUNT{bin,behavior,wave};
                 
                initanglesW=[initanglesW initangles];
                inittriW=[inittriW inittri];
             end
              initanglesW(:,~any(initanglesW,1))=[]; %delete off that first COL of all zeros
              inittriW(:,~any(inittriW,1))=[];
              
              SPIKEang{bin,behavior,region,treatment}=initanglesW; clear initanglesW
              TriCOUNT{bin,behavior,region,treatment}=inittriW; clear inittriW%instantaneous angles at time of spike w/in each bin -combined across trials and waveforms w/in mousese
               %save these final spike angles and trial id in ANA_inst file
                end
            end
            clear temp* epoch*
                
            
           % end
            clear spiketimes spikeangle initc inittri initang initangles
        end
    end
    

 disp('LOG: Instantaneous phase angles calculated' ); datetime
 
 %license checkout Distrib_Computing_Toolbox


%by mouse -
PPC0=cell(N_bin,N_behavior,N_region,N_treatment);
PPC1=cell(N_bin,N_behavior,N_region,N_treatment);
PPC2=cell(N_bin,N_behavior,N_region,N_treatment);
     for behavior=1:N_behavior
         for region=1:N_region
             for treatment=1:N_treatment
                 if size(unique(TriCOUNT{1,behavior,region,treatment}),2) < 18 && size(unique(TriCOUNT{2,behavior,region,treatment}),2) < 18 && ...
                             size(unique(TriCOUNT{3,behavior,region,treatment}),2) < 18 && size(unique(TriCOUNT{4,behavior,region,treatment}),2) < 18
                         %if there are less than 20 trials (unique means 20
                         %different trials - across all waveforms) in all
                         %bins
                         disp(['mouse: ' num2str(mouse) ' in Session: ' num2str(ses) ' ,region: ' num2str(region) ' ,for behavior: ' num2str(behavior)...
                             ' does not have 18 unique trials across all bins- terminating SFC calc'])
                         continue
                 end
                  
                 for bin=1:N_bin  
                     DATA=SPIKEang{bin,behavior,region,treatment};
                     TRIAL=TriCOUNT{bin,behavior,region,treatment};

       N=200; 
       N_permutations=500; %******
       PPC_0=zeros(size(frex,2),N_permutations);
       PPC_1=zeros(size(frex,2),N_permutations);
       PPC_2=zeros(size(frex,2),N_permutations);
for permi=1:N_permutations %do 500 permutations - each with 250 spikes; 1000 takes too long to do per mousese >96hrs
     %preallocate these  
    Trial=zeros(size(frex,2),N,N);
    Trial_Matrix=zeros(N,N); 
    
     %get data for permutation permi 
     [testang, idx]=datasample(DATA,N,2); % >> w/ replace = bootstrapping- to calculate CI
     Trial_ID=TRIAL(:,idx);
     
     %1. Calculate angle differences
            second_angle=1:1:N;%vectorization of bi (2nd angle)
            for first_angle=1:N
                   Trial(:,first_angle,:) = exp(1i*(  bsxfun(@minus,testang(:,first_angle),testang(:,second_angle))  )); 
                   Trial_Matrix(first_angle,:)= Trial_ID(1,first_angle)==Trial_ID(1,second_angle);
            end
            
            for fi=1:size(frex,2)
                TEMP=squeeze(Trial(fi,:,:));
                 Lu=tril(TEMP)~=0;  % Logical mask to remove upper
                 Ld=tril(TEMP)~=1;  % Logical mask to remove diagonal
                 Lo=logical(Lu.*Ld);% Logical mask - for ONLY lower triangle
                
                %PPC_0 from Vinck et al. 2010, NeuroImage
                  Temp_PPC0=TEMP(Lo);  
                  PPC_0(fi,permi)=abs(mean(Temp_PPC0)); % Average over all pfirst_anglers, take modulus
                  
                %PPC_1 from Vinck et al. 2012, J Comput Neurosci
                  Lt=tril(Trial_Matrix)~=1;
                  L1=logical(Lu.*Ld.*Lt);
                  TEMP_PPC1=TEMP(L1);
                  PPC_1(fi,permi)=abs(mean(TEMP_PPC1));
                
                % PPC_2 from Vinck et al. 2012, J Comput Neurosci
                  for triali=min(Trial_ID):max(Trial_ID)-1
                      L_ti=repmat(Trial_ID==triali,N,1);  % Mask for trial ti (row-based expanded to matrix)
                      L2=logical(L1'.* L_ti);
                      TEMP_PPC2(triali)=abs(mean(TEMP(L2))); % Average within that trial - what happens if no two trials are the same?
                      L_ti=[]; L2=[]; %clearing variables w/out clearing
                  end
                  TEMP_PPC2(isnan(TEMP_PPC2))=[]; %delete out Nan's - trials w/out repetition in sample taken
                  PPC_2(fi,permi)=abs(mean(TEMP_PPC2)); % Average over all trial-averaged pfirst_anglers, take modulus
              
              clear TEMP Lu Ld Lo L1 L2 L_ti
            end
  clear Trial_ID  Trial Trial_Matrix clear Temp*
end

        PPC0{bin,behavior,region,treatment}=mean(PPC_0,2);
        PPC1{bin,behavior,region,treatment}=mean(PPC_1,2);
        PPC2{bin,behavior,region,treatment}=mean(PPC_2,2);
        
        clear TEMP* PPC_1 PPC_2 PPC_0
        
disp(['LOG: X-mouse PPC complete for- behavior: ' num2str(behavior) ' region: ' num2str(region) ' treatment: ' num2str(treatment) ' Bin: ' num2str(bin)]); datetime

                     clear DATA TRIAL 
                 end
             end
         end
     end


disp(['LOG: By mouse PPC0-2 is complete for mouse ' num2str(mouse)])
disp('LOG: Saving data...'); datetime


%Save data by mouse + session - use general lock file so any sessions/
%mouse can call if lock file does not exist
while 1    
if exist('lock.session','file') ==  0 %if file does not exist
   xlswrite('lock.session','lock');
     if exist('SFC_ANA_Instmou.mat','file')==0 % does not exist
         %write variables 
           ANA{mouse,ses}.angles=ANGLE;%no need to load ANA b/c has not been created yet
           ANA{mouse,ses}.BIN_INTERVAL=BIN_INTERVAL;
           ANA{mouse,ses}.BIN_ANG=BIN_ANG;
           ANA{mouse,ses}.SPIKEang=SPIKEang;
           ANA{mouse,ses}.SPIKEtri=TriCOUNT;
          
           ANA{mouse,ses}.PPC0=PPC0; 
           ANA{mouse,ses}.PPC1=PPC1;
           ANA{mouse,ses}.PPC2=PPC2;
           save('SFC_ANA_Instmou','ANA','-v7.3'); %save data
           clear ANA ANGLE BIN_INTERVAL BIN_ANG SPIKEang TriCOUNT PPC0 PPC1 PPC2 DATA
           delete('lock.session')
           break %break the loop once have saved
           
     elseif exist('SFC_ANA_Instmou.mat','file')==2 % does exist - as matlab data file
            load SFC_ANA_Instmou % load the file
            %write variables
           ANA{mouse,ses}.angles=ANGLE;
           ANA{mouse,ses}.BIN_INTERVAL=BIN_INTERVAL;
           ANA{mouse,ses}.BIN_ANG=BIN_ANG;
           ANA{mouse,ses}.SPIKEang=SPIKEang;
           ANA{mouse,ses}.SPIKEtri=TriCOUNT;
          
           ANA{mouse,ses}.PPC0=PPC0; 
           ANA{mouse,ses}.PPC1=PPC1;
           ANA{mouse,ses}.PPC2=PPC2;
           save('SFC_ANA_Instmou','ANA','-v7.3'); %save data
           clear ANA ANGLE BIN_INTERVAL BIN_ANG SPIKEang TriCOUNT PPC0 PPC1 PPC2 DATA
           delete('lock.session')
           break %break the loop once have saved
     end
elseif exist('lock.session','file') ~=  0 %if file does exist
    pause(120) %pause for 2min before re-starting loop and looking for lock file
end
end


disp(['LOG: Session Finished for mousese: ' num2str(mouse) ' in session:' num2str(ses) '!!!!']); datetime

  end
  
  disp('All Mice Completed! Analysis Complete'); datetime
  
  %profsave(profile('info'),'Session1_Profile');
  %profile off
%this loop writes data in session order 1,2,3,4,5
% %holds until previous session finishes 
% while 1
%     if size(dir('done.*'),1)==0 && ses==1 %for first session when no done.# files exist
%        %probably don't need to look for lock file here, b/c first session it
%        %shouldn't exist due to other if statements within the while loop
%        %if exist('lock.session','file') ==  0 %if file does not exist
%           xlswrite('lock.session','lock');
%            
%            ANA{j,ses}.angles=ANGLE;%no need to load ANA b/c has not been created yet
%            ANA{j,ses}.BIN_INTERVAL=BIN_INTERVAL;
%            ANA{j,ses}.BIN_ANG=BIN_ANG;
%            ANA{j,ses}.SPIKEang=SPIKEang;
%            ANA{j,ses}.SPIKEtri=TriCOUNT;
%           
%            ANA{j,ses}.PPC0=PPC0; 
%            ANA{j,ses}.PPC1=PPC1;
%            ANA{j,ses}.PPC2=PPC2;
%            
%            save('SFC_ANA_instmouse','ANA','-v7.3'); %save data
               %xlswrite(['done.' num2str(ses)],'done'); %write done file
%                delete('lock.session')
%                break
%        %elseif exist('lock.session','file') ~=  0
%           % pause(120) %pause for two minutes
%        %end       
%     elseif size(dir('done.*'),1)==0 && ses~=1 %if another session finishes before ses 1
%          pause(120) %pause for 120sec - will keep getting paused until there is a done file in cd
% 
%     elseif size(dir('done.*'),1)~=0 %for all other sessions
%             files=dir('done.*');
%             filelast=files(size(files,1)).name; %find the name of the last file added
%             prevses=ses-1; prev=num2str(prevses); %turn into string b/c can only compare str to str
%             lastses=strsplit(filelast,'.');
%             if lastses{1,2}==prev %only if latest file added is previous session
%                if exist('lock.session','file') ==  0 %look for lock file - if not present
%                    xlswrite('lock.session','lock'); %write a lock file
%                     load SFC_ANA_instmouse
%                     
%                     ANA{j,ses}.angles=ANGLE;%no need to load ANA b/c has not been created yet
%                     ANA{j,ses}.BIN_INTERVAL=BIN_INTERVAL;
%                     ANA{j,ses}.BIN_ANG=BIN_ANG;
%                     ANA{j,ses}.SPIKEang=SPIKEang;
%                     ANA{j,ses}.SPIKEtri=TriCOUNT;
% 
%                     ANA{j,ses}.PPC0=PPC0; 
%                     ANA{j,ses}.PPC1=PPC1;
%                     ANA{j,ses}.PPC2=PPC2;
% 
%                         save('SFC_ANA_instmouse','ANA','-v7.3'); %save data
%                         xlswrite(['done.' num2str(ses)],'done'); %write done file
%                         delete('lock.session')
%                         break
%                elseif exist('lock.session','file') ~=  0  %if file is present
%                    pause(120) %pause for 120sec - then will continue through and re-loop through while 
%                end
%                
%             elseif lastses{1,2}~=prev
%                 pause(120) %pause for 120sec - then will continue through and re-loop while until condition is not met
%             end
%        
%    
%     end
% end


%disp('LOG: Session Finished!!!!'); datetime














