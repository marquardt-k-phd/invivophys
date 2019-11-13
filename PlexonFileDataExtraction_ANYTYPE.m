%%set the current directly and add the Plexon code folders

%NOTE: Befor you run this script you will need to have the plexon codes
%from online downloaded. 

clear all; clc

%define the places where the plexon codes are and where your data is
%located...
cd('E:\Documents\MATLAB');
addpath(genpath('/Users/Kristin/Documents/MATLAB'));%plexon codes

dataloc='insert place here';
addpath(genpath(dataloc));

%create a separate file with data location and all file names -> name it StudyID_Variables
    % Variables 
    %       1=files
    %       2=Channels
    %       3=Events
    %       4=Combined Files
    %       5=Mouse Number
%%
%FILL IN THESE VARIABLES
Numses=5;
Numtrt=2;

for ses=1:Numses
    Numfiles=size(var1,2);
    for j=1:Numfiles
        for trt=1:Numtrt
            if size(var1{ses,j,trt},1)==0
                continue
            end
        [OpenedFileName, Version, Freq, Comment, Trodalness, NPW, PreThresh, SpikePeakV, SpikeADResBits, SlowPeakV, SlowADResBits, Duration, DateTime] = plx_information(var1{ses,j,trt});
        Dur{j,ses,trt}=Duration;
        end
    end
end

        %get out some data if needed and display using below or clear
        %unused variables using below clear fxn
%             disp(['Opened File Name: ' OpenedFileName]);
%             disp(['Version: ' num2str(Version)]);
%             disp(['Frequency : ' num2str(Freq)]);
%             disp(['Comment : ' Comment]);
%             disp(['Date/Time : ' DateTime]);
%             disp(['Duration : ' num2str(Duration)]);
%             disp(['Num Pts Per Wave : ' num2str(NPW)]);
%             disp(['Num Pts Pre-Threshold : ' num2str(PreThresh)]);
%             
             clear OpenedFileName Version Comment Trodalness NPW PreThresh SpikePeakV SpikeADResBits SlowPeakV SlowADResBits DateTime Freq Duration
%%  
%This defines the field potential channels 
var2{1}= 'FP01'; var2{2}='FP02'; var2{3}='FP03'; var2{4}='FP04'; var2{5}='FP05'; var2{6}='FP06'; var2{7}='FP07'; var2{8}='FP08';
var2{9}='FP09'; var2{10}='FP10'; var2{11}='FP11'; var2{12}='FP12'; var2{13}='FP13'; var2{14}='FP14'; var2{15}='FP15'; var2{16}='FP16';

%this loop extracts the a/d firing in 1000Hz and puts into a cell array
%NOTE: This loop for 8 files and 16 channels can take 1hour and 30min
Numchan=16;
for ses=1:Numses
    for trt=1:Numtrt
for j=1:Numfiles
    for k=1:Numchan
         if size(var1{ses,j,trt},1)==0
            continue
         end
    [adfreq, n, ts, fn, ad] = plx_ad_v(var1{ses,j,trt},var2{k});%this line extracts the data for j file and k channel
    FP{j,k,ses,trt}=ad; %puts it into a cell array with each row as a file and each column as an electrode channel
    NumAD{j,ses,trt}=n;
    TS{j,ses,trt}=ts(1,1);
    disp([ 'Finished channel ' num2str(k) ' in file ' num2str(j) ' out of ' num2str(Numfiles) ' for Treatment' num2str(trt) '...' ]);
     clear ad n ts
    end 
end
    end
end



            %other things you can display information for
%             disp(['DigitizationFreq: ' num2str(adfreq)]);
%             disp(['TotalDataPoints: ' num2str(adn)]);
%             disp(['DataPointsinFragment: ' num2str(fn)]);
% 
             %clear adfreq fn ad

LFP.allFP=FP;
             
 %%           
%need to make a time stamp array for FP of each file - the time stamp ends
%after neuronal recording and thus the extra points need to be removed
format long g
        %SINT=LFP.interval; %--> for use with multiple runs
for ses=1:Numses
    for trt=1:Numtrt
for j=1:Numfiles
         if size(var1{ses,j,trt},1)==0
             continue
         end
    INT=TS{j,ses,trt}(1,1):.001:Dur{j,ses,trt}; %a cell is made for each file (by j rows) indicating interval of start TS to duration in 1000Hz increments
    INT=INT'; %transpose this to match ad frequencies stored in FP
    
    End_INT=length(INT); %The end number (per cell) for when deletion of extra numbers needs to end
    Start_INT=length(INT)-((length(INT)-1)-length(FP{j,1,ses,trt})); %the beginning number of where neuronal data ends but time continues
    INT(Start_INT:End_INT)=[]; %remove the extra timestamp values
    
    SINT{j,ses,trt}=INT; %put the proper TS array back into your interval
    clear INT END_INT Start_INT %clear your temporary interval for the next run 
end
    end
end
 %clear TS INT z j k g Start_INT End_INT Numchan NumAD Dur Numtrials 

 LFP.interval=SINT;
 
 
 
 %% AVERAGE FP ACROSS ELECTRODES --> THIS IS SO CAN SAVE VARIABLE
 
 %SINGLE REGION--> within hemispheres
 % (1-8=L) (9-16=R)
%     Num=size(var1,1);
%     for j=1:Num
%         Temp_FP_L=cell2mat(FP(j,1:8));  Temp_Avg_L=zeros(length(Temp_FP_L),1); 
%         Temp_FP_R=cell2mat(FP(j,9:16)); Temp_Avg_R=zeros(length(Temp_FP_R),1); 
%          for g=1:length(cell2mat(FP(j,1)))
%               Temp_Avg_L(g,1)=mean(Temp_FP_L(g,1:8));
%               Temp_Avg_R(g,1)=mean(Temp_FP_R(g,9:16));
%          end
%        FP_Avg{j,1}=Temp_Avg_L;
%        FP_Avg{j,2}=Temp_Avg_R;
%        clear Temp_Avg_L Temp_Avg_R Temp_FP_L Temp_FP_R
%     end
%     clear g j
%     
%     LFP.FP=FP_Avg; 
%     save('DRBR_FP','LFP','-v7.3');
    
 

%SINGLE REGION --> ALL ELECTRODES   
for ses=1:5
    for trt=1:2
        Numfiles=size(var1,2);
        for j=1:Numfiles
            if size(var1{ses,j,trt},1)==0
                continue
            end
            
          Temp_Avg=mean(cell2mat(FP(j,1:16,ses,trt)),2);
            FP_Avg{j,ses,trt}=Temp_Avg;
        end
        
    LFP.avgFP=FP_Avg; clear Temp_Avg
    end
end
 save('PAE_FP_DS','LFP','-v7.3');  %**********************
    
    
%DUAL REGION
% (1-8=DS) (9-16=OFC) 
    %FP_Avg=LFP.avgFP; %--> for use with multiple runs
for ses=4:5
    for trt=2:2
        Numfiles=size(var1,2);
for j=1:Numfiles
        if size(var1{ses,j,trt},1)==0
            continue
        end
        
   Temp_Avg_DS=mean(cell2mat(FP(j,1:8,ses,trt)),2);
   Temp_Avg_OFC=mean(cell2mat(FP(j,9:16,ses,trt)),2);
          
   FP_Avg{j,ses,1,trt}=Temp_Avg_DS;
   FP_Avg{j,ses,2,trt}=Temp_Avg_OFC;
   
   clear Temp_Avg_DS Temp_Avg_OFC
end
    end
    
    LFP.avgFP =FP_Avg; % mouse x session x region(ds=1; ofc=2) x trt
end
clear FP_Avg 
 save('PAE_FP_DUAL_all','LFP','-v7.3'); %**********************

%% Make an array of timestamped events recorded in plexon:

            %This can be used to determine all event names if not sure how they are labeled
           % [n, names]=plx_event_names(var1{1,1,1});

%**make a variable defining events**
var3{1}='EVT03'; var3{2}='EVT06';var3{3}='EVT08'; %var3{4}='EVT07';
    %event 3= tone off
    %event 6= house-light-off
    %event 8=magazine-light-off

%E=LFP.event; % for running multiple times

Numfiles=size(var1,2);
Numevt=3; %*make sure is correct number of events
for ses=1:Numses
    for trt=1:Numtrt
for j=1:Numfiles  
      if size(var1{ses,j,trt},1)==0
         continue
      end
    for g=1:Numevt
    [n, ts, sv]=plx_event_ts(var1{ses,j,trt},var3{g});%this line extracts the data for j file and g event
    E{j,g,ses,trt}=ts; %puts it into a cell array with each row as a file and each column as an event
    disp([ 'Finished event ' num2str(g) ' in file ' num2str(j) ' out of ' num2str(Numfiles) ' ...' ]);
    end
    clear n sv ts%clear these useless variables
end
    end
end
    LFP.event=E;

save('PAE_FP_DUAL','LFP','-v7.3');   %**********************

save('DRBR_allchanLFP','LFP','-v7.3');  

%% SAVE the following variables:

%FP_AVG = field potential (this may take a while to save)
%SINT = short time stamp variable
%E = event timestamps


%% Extraction of spike- timestamp data -> I am using this for extraction of spikes for Spike-LFP analysis


%Spikes=SFC.spikes;

%SINGLE REGION ELECTRODE 
Numses=5;
numtrt=2;
for ses=1:Numses
    Numfiles=size(var1,2);
    for j=1:Numfiles
        for trt=1:Numtrt
            if size(var1{ses,j,trt},1)==0
                continue
            end
            
    %find spikes in the file
     [tscounts, wfcounts, evcounts, contcounts]=plx_info(var1{ses,j,trt},0);
                clear contcounts evcounts wfcounts
            [row, col]=find(tscounts>0);
            chan=col-1; 
            unit=row-1; 
            clear row col

    %import the timestamp of each of those spikes
    for k=1:length(chan)
            c=chan(k);
            u=unit(k);
   
            [n,ts]=plx_ts(var1{ses,j,trt},c,u); 
            spike{k,1}=ts; clear ts
    end

    Spikes{j,ses,trt}=spike;
    clear spike c k u col row chan unit tscounts n
    
disp([ 'Finished Importing Spikes for:  Session ' num2str(ses) ' in file ' num2str(j) ' out of ' num2str(Numfiles) ' for trt' num2str(trt) '...' ]);


        end
    end
end
clear j Num

SFC.spikes=Spikes; clear Spikes


save('PAE_FP','SFC','-v7.3'); %**********************




%DUAL ELECTRODE EXTRACITON 
Numses=5;
numtrt=2;
Numtrt=2;

for ses=1:Numses
    Numfiles=size(var1,2);
    for j=1:Numfiles
        for trt=1:Numtrt
            if size(var1{ses,j,trt},1)==0
                continue
            end
            
         %find the spikes in the file
            [tscounts, wfcounts, evcounts, contcounts]=plx_info(var1{ses,j,trt},0);
                clear contcounts evcounts wfcounts
            [row, col]=find(tscounts>0);
            chan=col-1; 
            unit=row-1; 
            clear row col

            %import the timestamp of each of those spikes
            for k=1:length(chan)
                 c=chan(k);
                 u=unit(k);
                 
                if c<=8
                    [n,ts]=plx_ts(var1{ses,j,trt}, c,u); 
                    spikeDS{k,1}=ts; clear ts %DS waveforms
                    Spikes{j,ses,1,trt}=spikeDS;
                elseif c>8
                    [n,ts]=plx_ts(var1{ses,j,trt}, c,u); 
                    spikeOFC{k,1}=ts; clear ts %OFC waveforms 
                    Spikes{j,ses,2,trt}=spikeOFC;
                end
                
            end

            %Spikes{j,ses,1,trt}=spikeDS; 
            %Spikes{j,ses,2,trt}=spikeOFC; 
            
            clear spikeDS spikeOFC c k u col row chan unit tscounts n
            
            
        disp([ 'Finished Importing Spikes for:  Session ' num2str(ses) ' in file ' num2str(j) ' out of ' num2str(Numfiles) ' for trt' num2str(trt) '...' ]);
        
        end
    end
end

SFC.spikes=Spikes;

save('SFC_PAE_DUAL_inst_all','SFC','-v7.3'); %**********************

% for extra information on waveforms in your file
 %   [n, npw, ts, wave]= plx_waves_v('E:\Desktop\PAE_DataforMatlab\dS D85\PAE 3.0 D85 8-12-13_SortedOnly.plx', 14,2);%number of waveforms, number of points in waveform, array of timestamps arraw of waveforms in mV



