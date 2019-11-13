clear all; clc

%folders from MATLAB folder. Can't use general folder b/c EEGLab/Fieldtrip
%have seperate pca fxn that disrupts matlab programmed one
    addpath(genpath('E:\Documents\MATLAB\Cavanagh Class (Lec+Scripts)'));
    addpath(genpath('E:\Documents\MATLAB\kakearney-boundedline-pkg-2112a2b'));
    addpath(genpath('E:\Documents\MATLAB\ScottsScripts'));
    addpath(genpath('E:\Documents\MATLAB\Shuffle Mat'));
    addpath('E:\Documents\MATLAB\DATA_PAE');

addpath(genpath('D:\Kristin\OFC WT'));
addpath(genpath('F:\Plexon Inc Documents'));%plexon codes

DataLoc= 'F:\OFC WT'; %put the folder where your data is located here
addpath(genpath(DataLoc));

load SINGLEUNIT_PAE_DS_ver

save('SINGLEUNIT_PAE_DS_ver','SINGLEUNIT');

%% Get Data -> import directly from neuroexplorer 

%Create NeuroExploer files as have done previously to make TC and TIC
%timepoints

%under trial bin counts -> MatLab tab and check the box labeled 'Send
%matrix of results to Matlab'

    % you may also wish to send 'All Numberical Results' to Excel as an extra
    % output or check

% NOTE: the mice must be exported in order and organized in following
% manner:

%spikes(mouse, beh, sesion)
%mouse is simply the number you are inputting
%beh indicated TC(1) or TIC(2) 
%Session indicated: 1=D85 2=R25 3=R50 4=R85 etc...
%treatment (trt) is 1=control 2=treatment or genotype --> if all WT leave at 1
    %DO ALL OF ONE TRT FIRST THEN START OVER AGAIN AT MOUSE #1 FOR SECOND
    %TRT
    
    % this order is essential for the following code to work properly
    
 %if re-opening singleunit file then need to re-assign spikes variable so does not become over-written  
   % spikes=SINGLEUNIT.neuroexp;
  
%ctrl+enter will run this selection only
    
        mou=7;
        beh=2;
        ses=5;
        trt=2;

    spikes{mou,beh,ses,trt}=nex; clear nex nexColumnNames
    
    SINGLEUNIT.neuroexp=spikes;

   
    
%% Define Combined Files

%should be defined in your 'Variables' script under var4= location of file
%and var5= mouse number
   

%% read in combined file data... use ML as correct vs. incorrect measure
%Select Stage 9 Entry for each mouse: do not select the last row of zeros if it is there
ML=cell(1,5,2);
for ses=1:5 
    for trt=1:2
        Num=size(SINGLEUNIT.neuroexp(:,:,ses,trt),1);     
        for j=1:Num
            if size(var4{ses,j,trt},1) == 0
                continue
            end
            disp([ 'Select "Stage9:Entry" for mouse ' var5{ses,j,trt} ' , then hit "OK" ...']);
            ML{j,ses,trt}=xlsread (var4{ses,j,trt},-1);
        end

    end
end
SINGLEUNIT.ML=ML;


%count number of trials with ML <3sec
for ses=1:4
    for trt=1:1
        Num=size(SINGLEUNIT.neuroexp(:,:,ses,trt),1); 
        for j=1:Num
             if size(var4{ses,j,trt},1) == 0
                continue
             end
             %first must import all ML data
            disp([ 'Select "Stage9:Entry and Stage9:Exit" for mouse ' var5{ses,j,trt} ' for Session ' num2str(ses) ' , then hit "OK" ...']);
            tempML=xlsread (var4{ses,j,trt},-1);
            tempML(:,3)=tempML(:,2)-tempML(:,1); 
            ML{j,ses,trt}= tempML(:,3); %the ML of all trials
             
            %other things to calculate 
            avgML{j,ses,trt}=mean(tempML(find(tempML(:,3)~=0),3),1);%average ML for each mouse on the stage
            avgbelow(j,ses,trt)=mean(tempML(find(tempML(:,3)<300 & tempML(:,3)~=0),3)); %average for every mouse of trials that are <3sec
            
            num(j,ses,trt)=size(find(tempML(:,3)<300 & tempML(:,3)~=0),1);%number of trials with ML below 3sec (means gets truncated)
            totcorr(j,ses,trt)=size(find(tempML(:,3)~=0),1); %total correct trials
            tunMOU(j,ses,trt)=num(j,ses,trt)/totcorr(j,ses,trt); %proportion of total correct trials are truncated w/in each mouse
            
        end
    end
    tunMOU(isnan(tunMOU)) = 0;%just to be careful... chance all Nans to 0
    tunTOT(ses,trt)=sum(num(:,ses,trt))/sum(totcorr(:,ses,trt)); %proportion TC truncated trials total w/in stage
    
    avgbelow(isnan(avgbelow)) = 0;%must change all Nans to 0 before can calculate
    totavgbelow(ses,trt)=mean(avgbelow(:,ses,trt),1); %total avg w/in stage that is below 3sec 
end

SINGLEUNIT.allML=ML;
SINGLEUNIT.avgML=avgML;
SINGLEUNIT.tunMOU=tunMOU;
SINGLEUNIT.avgbelow3=avgbelow;


%find outliers in ML
for ses=1:4
    for trt=1:1
        Num=size(SINGLEUNIT.neuroexp(:,:,ses,trt),1); 
        for j=1:Num
             if size(var4{ses,j,trt},1) == 0
                continue
             end
             
             tempML=ML{j,ses,trt};
             tempML(tempML==0)=[];
             thresh=2*std(tempML,[],1)+ mean(tempML,1);
             out=find(delta>=thresh);
             
             out=find(tempML>=thresh);
             tempML(out,:)=[];
             
             avgml(j,ses,trt)=mean(tempML,1);
             
        end
    end
end
SINGLEUNIT.avgMLoutrem=avgml;
%% DEFINE TRIAL TYPES based on ML --> could probably do this from just spike data, but I was too lazy to write that program

for ses=1:4
    for trt=1:1
    Num=size(SINGLEUNIT.neuroexp(:,:,ses,trt),1);  %number of mice
    for s=1:Num
        if size(SINGLEUNIT.ML{s,ses,trt},1)==0
            continue
        end
        
  %trial numbers for TC and TIC
   code=SINGLEUNIT.ML{s,ses,trt};
    
%Types 1= correct -> correct
%      2= correct -> incorrect
%      3= incorrect -> incorrect
%      4= incorrect-> correct
        for j=2:size(code,1)

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

%trial numbers for TC and TIC
    x=find(code(:,1)==0);
    y=find(code(:,1)>0);
    

    type{s,1,ses,trt}=TYPE(y,:); clear y %correct
    type{s,2,ses,trt}=TYPE(x,:); clear x %incorrect
    
    clear TYPE code
    
    end
    
    end 
end



%cell with more beh MLs than neuron trials
%     Temp=cell2mat(type(10,2,2));
%     Temp(104,:)=[];
%     type{10,2,2}=Temp;
%     
%     Temp=cell2mat(type(4,2,3));
%     Temp(71:72,:)=[];
%     type{4,2,3}=Temp;

    
  
%split spikes the sort based on trial type
 for beh=1:2
     for ses=1:4
         for trt=1:1
         Num=size(SINGLEUNIT.neuroexp(:,:,ses,trt),1);  
         for k=1:Num
       
     spikes=SINGLEUNIT.neuroexp;
    temp=cell2mat(spikes(k,beh,ses,trt)); %SINGLEUNIT.neuroexp or spikes
    stuff=cell2mat(type(k,beh,ses,trt));
    if size(stuff,1)==0
        continue
    end

    len=size(temp,2);
    num=len/80;
    
    for j=1:num
        for trial=1:size(temp,1)

        last=j*80;
        start=last-79;

        Neurons(trial,:,j)=temp(trial,start:last);
            clear last start
        end
    end
    
    if len==0 
        SPIKES{k,beh,ses,trt}= 0; clear y
    elseif len>0
        SPIKES{k,beh,ses,trt}=Neurons;
    end
        clear  num 
       
                
   %sort into corrects & incorrects only
   if beh==1
       if len==0
           N=0;
       elseif len>0
           TC{k,1,ses,trt}=mean(Neurons,1);
           %TCnum{k,1,ses,trt}=Neurons;
       end
   elseif beh==2
       if len==0
           M=0;
       elseif len>0
           TIC{k,1,ses,trt}=mean(Neurons,1);
           %TICnum{k,1,ses,trt}=Neurons;
       end
   end
   
        
   %sort into trial types
   if beh==1
       if len==0
            Z=0;
            S=0;
       elseif len>0
          if size(stuff,1) > size(Neurons,1)
              over= size(stuff,1) - size(Neurons,1);
              last=size(stuff,1);
              lastover=last-over +1;
              stuff(lastover:last,:)=[];
          end
        Z=find(stuff==1);
            TYPE{k,1,ses,trt}=Neurons(Z,:,:);
        S=find(stuff==4);
            TYPE{k,4,ses,trt}=Neurons(S,:,:);
       end
            
   elseif beh==2
          if size(stuff,1) > size(Neurons,1)
              over= size(stuff,1) - size(Neurons,1);
              last=size(stuff,1);
              lastover=last-over +1;
              stuff(lastover:last,:)=[];
          end
        if len==0
            T=0;
            U=0;
       elseif len>0
        T=find(stuff==3);
            TYPE{k,3,ses,trt}=Neurons(T,:,:);
        U=find(stuff==2);
            TYPE{k,2,ses,trt}=Neurons(U,:,:);
        
        end
   end
    clear Z S T U Neurons len 

         end 
       
         end
     end
 end
       clear stuff type temp over lastover last
 
%  emptysub=double.empty(0,80,4);
%     TYPE{4,1,2}=emptysub;
%     TYPE{4,4,2}=emptysub;

 SINGLEUNIT.trialspikes=TYPE; %mouse x trial type x session
 SINGLEUNIT.TCspikes=TC;
 SINGLEUNIT.TICspikes=TIC;

 save('SINGLEUNIT_PAE_OFC','SINGLEUNIT'); %*****************

%% Average Trials and Concatenate spikes


for ses=1:4
    for ty=1:4
        for trt=1:1
 
     allspikes=zeros(80,1);
     com=zeros(80,1);
     Num=size(SINGLEUNIT.trialspikes(:,:,ses,trt),1);
     
    for j=1:Num
        temp=cell2mat(SINGLEUNIT.trialspikes(j,ty,ses,trt));
        if size(temp,1)==0
            continue
        end
            temp2=squeeze(mean(temp,1));
            if size(temp,3)==1
                temp2=temp2';
            end
        allspikes=[allspikes  temp2]; clear temp
    end
    

        for j=1:Num
            if ty==1, temp=SINGLEUNIT.TCspikes{j,1,ses,trt};
            elseif ty==2,  temp=SINGLEUNIT.TICspikes{j,1,ses,trt};
            elseif ty==3
                continue
            elseif ty==4
                continue
            end

            if size(temp,1)==0
                continue
            end
            if size(temp,3)==1
                temp=temp';
            end
            temp=squeeze(temp);
          
            com=[com temp]; clear temp
        end
  

    
    allspikes(:,1)=[];
    com(:,1)=[];
    ALL{ses,ty,trt}=allspikes; %all types
    COM{ses,ty,trt}=com;  %tc and tic only
    clear allspikes temp*

        end
    end
end

SINGLEUNIT.comspike=ALL; %session x trial type (time x neuron)  x trt
SINGLEUNIT.comspiketctic=COM; % session x type (tc tic) x trt



 save('SINGLEUNIT_OFC','SINGLEUNIT'); %*****************

%total number phasics (within trial type added across) per mouse
for ses=1:4
    for type=1:4
      Num=size(var4(ses,:),2);
    for j=1:Num
        temp=cell2mat(SINGLEUNIT.neuroexp(j,2,ses));
        if size(temp,1)==0
            continue
        end
  
      temp1=SINGLEUNIT.trialspikes{j,type,ses};
      temp2=mean(temp1,1);

        num=size(temp2,3);%finding number of neurons

        
%         last=N*80;
%         start=last-79;
%         Neurons(:,N)=temp1(:,start:last);
%             clear last start

    %define if phasic
    base=[1:20];
    post=[21:80];
    
    for p=1:num
       sig(p)=ttest2(temp2(:,base,p), temp2(:,post,p)); %0=nonphasic 1=phasic
       sig(isnan(sig))=0;
    end

%     if size(find(sig==1),2)==0
%          phasicCount(j,type,ses)=0; clear sig Neurons temp*
%     elseif size(find(sig==1),2)>0
        phasicCount(j,type,ses)=size(find(sig==1),2); clear sig Neurons temp* %count phasics and record
%     end
 
    end
    end
end


%%  THE FOLLOWING ARE DIFFERENT TYPES OF GRAPHS TO MAKE --> need to remove region from all calcs below...

%define baseline and post time bins
base=[1:20];
post=[21:80];

%RAW FIRING RATE PHASICS ONLY w/ separated increasing and decreasing 
        %calculate baseline average
        for j=1:4
           for ses=1:4
                   for trt=1:1

            Temp_DATA=cell2mat(SINGLEUNIT.trialspikes(ses,j,trt));
            %Temp_DATA=cell2mat(SINGLEUNIT.comspike(ses,j,trt));
            %Temp_DATA=SINGLEUNIT.comspiketctic{ses,j,trt};

            sig=ttest2(Temp_DATA(base,:), Temp_DATA(post,:)); %0=nonphasic 1=phasic

                sig(isnan(sig))=0;
            PHASIC{ses,j}=sig;

            z=find(sig==1);
                phasic=Temp_DATA(:,z);
                phasic_base=mean(phasic(base,:),1);
                phasic_post=mean(phasic(post,:),1);

                increasing=find(phasic_post > phasic_base);
                decreasing=find(phasic_post < phasic_base);


            PHASICFIRING{ses,j,trt}=Temp_DATA(:,z);
           % PHASICIN{ses,j,trt}=phasic(:,increasing);
           % PHASICDE{ses,j,trt}=phasic(:,decreasing);

                   end
           end

            clear Temp_DATA phasic increasing decreasing z sig phasic_base phasic_post 
        end

        SINGLEUNIT.phasicrawtctic=PHASICFIRING; %SINGLEUNIT.phasicrawtypes=PHASICFIRING;
        SINGLEUNIT.phasicrawDEC=PHASICDE;
        SINGLEUNIT.phasicrawINC=PHASICIN; 
        
%removal of outliers
        for j=1:2
            for ses=1:5
                    for trt=1:2

                %Temp_DATA=cell2mat(SINGLEUNIT.phasicraw(ses,j,trt));
                Temp_DATA=SINGLEUNIT.comspiketctic{ses,j,trt};
                postavg=mean(Temp_DATA(post,:),1);
                    popavg=mean(Temp_DATA(post,:),2);
                    popstd=std(popavg,0,1);

                posout=mean(popavg,1)  +  2*popstd ;
                negout=mean(popavg,1)  -  2*popstd ;

                outliers=find(postavg>posout | postavg<negout);
                
                Temp_DATA(:,outliers)=[];
                sd=std(Temp_DATA,0,2);
                sem= sd ./ sqrt(size(Temp_DATA,2));

                Rawout{ses,j,trt}=Temp_DATA;
                SEM{ses,j,trt}=sem;

                clear Temp* posout neg out outliers postavg popavg popstd

                    end
            end
        end

            SINGLEUNIT.rawOUTREMOVED=Rawout;
            
% PHASIC NEURONS ONLY Z-SCORED (NO SEPARATION INTO INCREASING AND DECREASING)
%z-score phasic neurons
%Note: commented sections are for sorting out increasing and decreasing
%phasic neurons seperately
        for j=1:2
            for ses=1:5
                    for trt=1:2

                Temp_DATA=cell2mat(SINGLEUNIT.rawOUTREMOVED(ses,j,trt));
                %Temp_DATA=SINGLEUNIT.phasicrawtctic{ses,j,trt};
                %Temp_DATA=SINGLEUNIT.phasicrawrcticOUTREMOVED{ses,j,trt};
                    BASE_AVG=mean(Temp_DATA(base,:),1);
                    BASE_SD=std(Temp_DATA(base,:),0,1);
        %         Temp_IN=cell2mat(SINGLEUNIT.phasicinraw(ses,j));
        %             BASE_IN=mean(Temp_IN(base,:),1);
        %             SD_IN=std(Temp_IN(base,:),0,1);
        %         Temp_DE=cell2mat(SINGLEUNIT.phasicdecraw(ses,j));
        %             BASE_DE=mean(Temp_DE(base,:),1);
        %             SD_DE=std(Temp_DE(base,:),0,1);


               if any(BASE_AVG==0)>0 %if any of the baselines are zero and removes the data
                    [col]=find(BASE_AVG==0);
                    Temp_DATA(:,col)=[];
                    BASE_SD(:,col)=[];
                    BASE_AVG(:,col)=[];
                    clear row col
               end

        %         if any(BASE_IN==0)>0 %if any of the baselines are zero and removes the data
        %             [col]=find(BASE_IN==0);
        %             Temp_IN(:,col)=[];
        %             BASE_IN(:,col)=[];
        %             SD_IN(:,col)=[];
        %             clear row col
        %         end
        %        
        %         if any(BASE_DE==0)>0 %if any of the baselines are zero and removes the data
        %             [col]=find(BASE_DE==0);
        %             Temp_DE(:,col)=[];
        %             BASE_DE(:,col)=[];
        %             SD_DE(:,col)=[];
        %             clear row col
        %         end


                    Temp_Z= bsxfun(@minus, Temp_DATA, BASE_AVG); 
                    Temp_Z= bsxfun(@rdivide, Temp_Z, BASE_SD);
                       % Temp_zSD=std(Temp_Z,0,2);
                       % Temp_ZSE=Temp_zSD ./ sqrt(size(Temp_Z,2));
        %             
        %             Temp_ZIN= bsxfun(@minus, Temp_IN, BASE_IN); 
        %             Temp_ZIN= bsxfun(@rdivide, Temp_ZIN, SD_IN);
        %             
        %             Temp_ZDE= bsxfun(@minus, Temp_DE, BASE_DE); 
        %             Temp_ZDE= bsxfun(@rdivide, Temp_ZDE, SD_DE);


                %smooth the z-scored data
                    sm=smooth(Temp_Z); clear COne
                    sm_reshaped=reshape(sm,80,size(Temp_Z,2));
                        Temp_zSDsm=std(sm_reshaped,0,2);
                        Temp_ZSEsm=Temp_zSDsm ./ sqrt(size(Temp_Z,2));
                            PHASICZsmooth{ses,j,trt}=sm_reshaped; clear sm sm_reshaped
                            SEMZsmooth{ses,j,trt}=Temp_ZSEsm; 


        %             sm=smooth(Temp_ZIN); clear COne
        %             sm_reshaped=reshape(sm,80,size(Temp_ZIN,2));
        %                 Temp_zINSDsm=std(sm_reshaped,0,2);
        %                 Temp_ZINSEsm=Temp_zINSDsm ./ sqrt(size(Temp_ZIN,2));
        %                     PHASICZin{ses,j}=sm_reshaped; clear sm_reshaped
        %                     SEMZin{ses,j}=Temp_ZINSEsm;
        %                     
        %             sm=smooth(Temp_ZDE); clear COne
        %             sm_reshaped=reshape(sm,80,size(Temp_ZDE,2));
        %                 Temp_zDESDsm=std(sm_reshaped,0,2);
        %                 Temp_ZDESEsm=Temp_zDESDsm ./ sqrt(size(Temp_ZDE,2));
        %                     PHASICZde{ses,j}=sm_reshaped; clear sm_reshaped
        %                     SEMZde{ses,j}=Temp_ZDESEsm;



              clear Temp* BASE*
                    end 
            end
        end

        SINGLEUNIT.rawZtcticsmoothOUTREM=PHASICZsmooth;
            SINGLEUNIT.rawZtcticsmoothSEMOUTREM=SEMZsmooth;
        %SINGLEUNIT.phasiczin=PHASICZin;
           % SINGLEUNIT.semzin=SEMZin;
        %SINGLEUNIT.phasiczde=PHASICZde;
            %SINGLEUNIT.semzde=SEMZde;
            
%SESSION X TYPE X TRT
  


        %GRAPH
        for ses=1:5
     %           for trt=1:2

            type1=mean(SINGLEUNIT.rawZtcticsmooth{ses,1,1},2);
                type1sem=SINGLEUNIT.rawZtcticsmoothSEM{ses,1,1};
            type1exp=mean(SINGLEUNIT.rawZtcticsmooth{ses,1,2},2);
                type1expsem=SINGLEUNIT.rawZtcticsmoothSEM{ses,1,2};
                
            type2=mean(SINGLEUNIT.rawZtcticsmooth{ses,2,1},2);
                type2sem=SINGLEUNIT.rawZtcticsmoothSEM{ses,2,1};
            type2exp=mean(SINGLEUNIT.rawZtcticsmooth{ses,2,2},2);
                type2expsem=SINGLEUNIT.rawZtcticsmoothSEM{ses,2,2};
                
%             type3=mean(SINGLEUNIT.phasicZtcticsmooth{ses,3},2);
%                 type3sem=SINGLEUNIT.phasicZtcticsmoothSEM{ses,3};
%             type4=mean(SINGLEUNIT.phasicZtcticsmooth{ses,4},2);
%                 type4sem=SINGLEUNIT.phasicZtcticsmoothSEM{ses,4};


                x=[1:80]';

            figure
            hold on
                boundedline(x,type1,type1sem,'b','alpha'); 
                boundedline(x,type1exp,type1expsem,'g','alpha');
                plot(x,zeros(1,80),'k')
                title('Correct');
                set(gca,'ylim',[-1 5]);
                
%             figure
%             hold on
%                 boundedline(x,type2,type2sem,'b','alpha'); 
%                 boundedline(x,type2exp,type2expsem,'g','alpha');
%                 plot(x,zeros(1,80),'k')
%                 title('Incorrect')
                %set(gca,'ylim',[-.4 1.2]);
%                 boundedline(x,type3,type3sem,'y','alpha'); 
%                 boundedline(x,type4,type4sem,'r','alpha'); 

               % end   
        end
        
  
                                        %SES,BEH,TRT
        temp=SINGLEUNIT.rawZtcticsmooth{  2 , 1 , 1}';
        
          %remove single neuron that is causing huge error --> on D85 TC 
        temp=SINGLEUNIT.rawZtcticsmoothOUTREM{2,1,1}';
          temp(6,:)=[];
          temp=temp';
          SINGLEUNIT.rawZtcticsmoothOUTREM{2,1,1}=temp;
         tempsd=std(temp,0,2);
         tempsem=tempsd ./ sqrt(size(temp,2));
            SINGLEUNIT.rawZtcticsmoothSEMOUTREM{2,1,1}=tempsem;
            clear temp*
            
            
                %remove single neuron that is causing huge error --> on D85 TC 
        temp=SINGLEUNIT.rawZtcticsmooth{2,1,1}';
          temp(21,:)=[];
          temp=temp';
          SINGLEUNIT.rawZtcticsmooth{2,1,1}=temp;
         tempsd=std(temp,0,2);
         tempsem=tempsd ./ sqrt(size(temp,2));
            SINGLEUNIT.rawZtcticsmoothSEM{2,1,1}=tempsem;
            clear temp*
       


% ALL NEURONS Z-SCORE(PHASIC AND NON-PHASIC)
%Note: only neurons with baseline zero are removed from data
        base=[1:20];
        post=[21:80];

        for j=1:2
            for ses=1:5
                    for trt=1:2

                 Temp_DATA=cell2mat(Rawout(ses,j,trt));
                 %Temp_DATA=cell2mat(SINGLEUNIT.comspiketctic(ses,j));
                 %Temp_DATA=SINGLEUNIT.rawOUTREMOVED{ses,j,trt};
                    Temp_DATA(isnan(Temp_DATA))=0;

                    BASE_AVG=mean(Temp_DATA(base,:),1);
                    BASE_SD=std(Temp_DATA(base,:),0,1);


                   if any(BASE_AVG==0)>0 %if any of the baselines are zero and removes the data
                        [col]=find(BASE_AVG==0);
                        Temp_DATA(:,col)=[];
                        BASE_SD(:,col)=[];
                        BASE_AVG(:,col)=[];
                        clear row col
                   end

                    Temp_Z= bsxfun(@minus, Temp_DATA, BASE_AVG); 
                    Temp_Z= bsxfun(@rdivide, Temp_Z, BASE_SD);
        %                 Temp_zSD=std(Temp_Z,0,2);
        %                 Temp_ZSE= Temp_zSD ./ sqrt(size(Temp_Z,2));


                %smooth the z-scored data
                    sm=smooth(Temp_Z); clear COne
                    sm_reshaped=reshape(sm,80,size(Temp_Z,2));
                        Temp_zSDsm=std(sm_reshaped,0,2);
                        Temp_ZSEsm=Temp_zSDsm ./ sqrt(size(Temp_Z,2));

                           %RAWSMOOTH{ses,j}=Temp_Z; clear sm sm_reshaped 
                            %SEMRAWsmooth{ses,j}=Temp_ZSE; clear Temp*

                             RAWSMOOTH{ses,j,trt}=sm_reshaped; clear sm sm_reshaped
                             SEMRAWsmooth{ses,j,trt}=Temp_ZSEsm; clear Temp*
                    end
            end 
        end
        
         SINGLEUNIT.rawZtcticsmoothOUT=RAWSMOOTH;
         SINGLEUNIT.rawZtcticsmoothOUTSEM=SEMRAWsmooth;
        
        %remove single neuron that is causing huge error --> on D85 TC 
%         temp=SINGLEUNIT.rawZtcticsmooth{2,1,1}';
%           temp(22,:)=[];
%           temp(21,:)=[];
%           temp=temp';
%           SINGLEUNIT.rawZtcticsmooth{2,1,1}=temp;
%          tempsd=std(temp,0,2);
%          tempsem=tempsd ./ sqrt(size(temp,2));
%             SINGLEUNIT.rawZtcticsmooth{2,1,1}=tempsem;
%             clear temp*
    
         
 %BASELINE COMPARISON 
         for j=1:2
            for ses=1:5
                    for trt=1:2

                % Temp_DATA=SINGLEUNIT.comspike{ses,j};
                 Temp_DATA=SINGLEUNIT.comspiketctic{ses,j,trt};
                    Temp_DATA(isnan(Temp_DATA))=0;

                    BASE_AVG=mean(Temp_DATA(base,:),1);
                   % BASE_SD=std(BASE_AVG,0,2);


                   if any(BASE_AVG==0)>0 %if any of the baselines are zero and removes the data
                        [col]=find(BASE_AVG==0);
                        Temp_DATA(:,col)=[];
                        %BASE_SD(:,col)=[];
                        BASE_AVG(:,col)=[];
                        clear row col
                   end

                BASE_SD=std(BASE_AVG,0,2);
                    BASE_SEM=BASE_SD./ sqrt(size(BASE_AVG,2));


                BASELINE{ses,j,trt}=BASE_AVG; clear sm sm_reshaped BASE_AVG
                %BASESD(ses,j,trt)=BASE_SD; clear Temp* BASE_SD
                BASESEM(ses,j,trt)=BASE_SEM; clear BASE_SEM
                    end
    
            end
        end
        
      SINGLEUNIT.baseuncorr=BASELINE;
      %SINGLEUNIT.baseuncorrSD=BASESD;
      SINGLEUNIT.baseuncorrSEM=BASESEM;
      
      %ses x type x trt
      
      temp=BASELINE{5,2,2}';
      
%graph of the baselines

%for ses=1:4
    for type=1:4
      
            

%incorrects together
%             BASEcorrcont=mean([SINGLEUNIT.baseuncorr{ses,1,reg,1} SINGLEUNIT.baseuncorr{ses,4,reg,1}],2);
%                 errorcorrcont=mean( [SINGLEUNIT.baseuncorrSEM(ses,1,reg,1) SINGLEUNIT.baseuncorrSEM(ses,2,reg,1)],2); 
%             BASEincont=mean([SINGLEUNIT.baseuncorr{ses,2,reg,1} SINGLEUNIT.baseuncorr{ses,3,reg,1}],2);
%                 errorincont=mean( [SINGLEUNIT.baseuncorrSEM(ses,2,reg,1) SINGLEUNIT.baseuncorrSEM(ses,3,reg,1)],2); 
%                 
%             BASEcorrexp=mean([SINGLEUNIT.baseuncorr{ses,1,reg,2} SINGLEUNIT.baseuncorr{ses,4,reg,2}],2);
%                 errorcorrexp=mean( [SINGLEUNIT.baseuncorrSEM(ses,1,reg,2) SINGLEUNIT.baseuncorrSEM(ses,4,reg,2)],2); 
%             BASEinexp=mean([SINGLEUNIT.baseuncorr{ses,2,reg,2} SINGLEUNIT.baseuncorr{ses,3,reg,2}],2);
%                 errorinexp=mean( [SINGLEUNIT.baseuncorrSEM(ses,2,reg,2) SINGLEUNIT.baseuncorrSEM(ses,3,reg,2)],2); 
%             
%             
%             figure
%             hold on
%             bar(1,BASEcorrcont,'k')
%                 errorbar(1,BASEcorrcont,errorcorrcont,'k')
%             bar(2,BASEcorrexp,'b')
%                 errorbar(2,BASEcorrexp,errorincont,'k')
%             bar(3,BASEincont,'k')
%                 errorbar(3,BASEincont,errorcorrexp,'k')
%             bar(4,BASEinexp,'b')
%                 errorbar(4,BASEinexp,errorinexp, 'k')
%                 
%             set(gca,'ylim',[0 .4]);
%             
%             title(['Session:' num2str(ses)  'Region:' num2str(reg) ])

% by type - comparing across sessions

              BASEcont1=mean(SINGLEUNIT.baseuncorr{1,type},2);
              %BASEexpr1=mean(SINGLEUNIT.baseuncorr{1,type,reg,2},2);
              
              BASEcont2=mean(SINGLEUNIT.baseuncorr{2,type},2);
              %BASEexpr2=mean(SINGLEUNIT.baseuncorr{2,type,reg,2},2);
              
              BASEcont3=mean(SINGLEUNIT.baseuncorr{3,type},2);
              %BASEexpr3=mean(SINGLEUNIT.baseuncorr{3,type,reg,2},2);
              
              BASEcont4=mean(SINGLEUNIT.baseuncorr{4,type},2);
              %BASEexpr4=mean(SINGLEUNIT.baseuncorr{4,type,reg,2},2);
              
              %BASEcont5=mean(SINGLEUNIT.baseuncorr{5,type},2);
              %BASEexpr5=mean(SINGLEUNIT.baseuncorr{5,type,reg,2},2);
              

              figure
              hold on
              %subplot(1,4,1)
              %boxplot(SINGLEUNIT.baseuncorr{1,type})
               %set(gca,'ylim',[0 1.2]);
                bar(1,BASEcont1,'k')
                    errorbar(1,BASEcont1, SINGLEUNIT.baseuncorrSEM(1,type),'k');
               % bar(2,BASEexpr1,'b')
                   % errorbar(2,BASEexpr1 ,SINGLEUNIT.baseuncorrSEM(1,type,reg,2),'k');
                % subplot(1,4,2)
              %  boxplot(SINGLEUNIT.baseuncorr{2,type})
                 %set(gca,'ylim',[0 1.2]);
               bar(3,BASEcont2,'k')
                    errorbar(3,BASEcont2 ,SINGLEUNIT.baseuncorrSEM(2,type),'k');
                %bar(4,BASEexpr2,'b')
                   %  errorbar(4,BASEexpr2 ,SINGLEUNIT.baseuncorrSEM(2,type,reg,2),'k');
                % subplot(1,4,3)
               % boxplot(SINGLEUNIT.baseuncorr{3,type})
                 %set(gca,'ylim',[0 1.2]);
                   bar(5,BASEcont3,'k')
                    errorbar(5,BASEcont3 ,SINGLEUNIT.baseuncorrSEM(3,type),'k');
               % bar(6,BASEexpr3,'b')
                   % errorbar(6, BASEexpr3,SINGLEUNIT.baseuncorrSEM(3,type,reg,2),'k');
                 %subplot(1,4,4)
               % boxplot(SINGLEUNIT.baseuncorr{4,type})
                 %set(gca,'ylim',[0 1.2]);
                   bar(7,BASEcont4,'k')
                    errorbar(7, BASEcont4,SINGLEUNIT.baseuncorrSEM(4,type),'k');
                %bar(8,BASEexpr4,'b')
                  %  errorbar(8,BASEexpr4 ,SINGLEUNIT.baseuncorrSEM(4,type,reg,2),'k');
                
               % bar(9,BASEcont5,'k')
                    %errorbar(9,BASEcont5 ,SINGLEUNIT.baseuncorrSEM(5,type,reg,1),'k');
               % bar(10,BASEexpr5,'b')
                    %errorbar(10,BASEexpr5 ,SINGLEUNIT.baseuncorrSEM(5,type,reg,2),'k');

        set(gca,'ylim',[0 .4]);

                title(['Type:' num2str(type)   ])

%incorrects together
%             BASEcorrcont=mean(SINGLEUNIT.baseuncorrtctic{ses,1,1},2);
%                 errorcorrcont=mean( SINGLEUNIT.baseuncorrSEMtctic(ses,1,1),2); 
%             BASEincont=mean(SINGLEUNIT.baseuncorrtctic{ses,2,1} ,2);
%                 errorincont=mean( SINGLEUNIT.baseuncorrSEMtctic(ses,2,1),2); 
%                 
%             BASEcorrexp=mean(SINGLEUNIT.baseuncorrtctic{ses,1,2},2);
%                 errorcorrexp=mean( SINGLEUNIT.baseuncorrSEMtctic(ses,1,2) ,2); 
%             BASEinexp=mean(SINGLEUNIT.baseuncorrtctic{ses,2,2},2);
%                 errorinexp=mean(SINGLEUNIT.baseuncorrSEMtctic(ses,2,2) ,2); 
%             
%             
%             figure
%             hold on
%             bar(1,BASEcorrcont,'k')
%                 errorbar(1,BASEcorrcont,errorcorrcont,'k')
%             bar(2,BASEcorrexp,'b')
%                 errorbar(2,BASEcorrexp,errorincont,'k')
%             bar(3,BASEincont,'k')
%                 errorbar(3,BASEincont,errorcorrexp,'k')
%             bar(4,BASEinexp,'b')
%                 errorbar(4,BASEinexp,errorinexp, 'k')
%                 
%             set(gca,'ylim',[0 .35]);
%             
%             title(['Session:' num2str(ses)   ])


     
            
   
    end
%end


 SINGLEUNIT.baseuncorrtctic{5,2,2}';     


 % ALL NEURONS NO Z-SCORE --> B/C BASELINES SIG DIFF 
    %SMOOTHED ONLY
%Note: only neurons with baseline zero are removed from data
        base=[1:20];
        post=[21:80];

        for j=1:4
            for ses=1:5
                for reg=1:2
                    for trt=1:2

                 Temp_DATA=SINGLEUNIT.comspike{ses,j,reg,trt};
                    Temp_DATA(isnan(Temp_DATA))=0;

                    BASE_AVG=mean(Temp_DATA(base,:),1);
                    BASE_SD=std(Temp_DATA(base,:),0,1);


                   if any(BASE_AVG==0)>0 %if any of the baselines are zero and removes the data
                        [col]=find(BASE_AVG==0);
                        Temp_DATA(:,col)=[];
                        BASE_SD(:,col)=[];
                        BASE_AVG(:,col)=[];
                        clear row col
                   end

%                     Temp_Z= bsxfun(@minus, Temp_DATA, BASE_AVG); 
%                     Temp_Z= bsxfun(@rdivide, Temp_Z, BASE_SD);         

                %smooth the z-scored data
                    sm=smooth(Temp_DATA); clear COne
                    sm_reshaped=reshape(sm,80,size(Temp_DATA,2));
                        Temp_zSDsm=std(sm_reshaped,0,2);
                        Temp_ZSEsm=Temp_zSDsm ./ sqrt(size(Temp_DATA,2));


                             SMOOTH{ses,j,reg,trt}=sm_reshaped; clear sm sm_reshaped
                             SEMsmooth{ses,j,reg,trt}=Temp_ZSEsm; clear Temp*
                    end
                end
            end
        end
 
      SINGLEUNIT.rawnonz=SMOOTH;
      SINGLEUNIT.rawnonzsem=SEMsmooth;
      
      
   % NO Z-SCORE WITH COMBINED BASELINE W/IN TRIAL TYPE AND SESSION
    %SMOOTHED ONLY
%Note: only neurons with baseline zero are removed from data
        base=[1:20];
        post=[21:80];
        
        for j=1:4
            for ses=1:5
                for reg=1:2
                    for trt=1:2

                    Temp_DATA=SINGLEUNIT.comspike{ses,j,reg,trt};
                    
                    Temp_DATA1=SINGLEUNIT.comspike{ses,j,reg,1};
                    Temp_DATA2=SINGLEUNIT.comspike{ses,j,reg,2};
                    Temp_DATA1(isnan(Temp_DATA1))=0;
                    Temp_DATA2(isnan(Temp_DATA2))=0;

                    BASE_AVG1=mean(Temp_DATA1(base,:),1);
                    BASE_AVG2=mean(Temp_DATA2(base,:),1);
                    
                    BASE_AVG=mean([BASE_AVG1 BASE_AVG2],2);
                    BASE_SD=std([BASE_AVG1 BASE_AVG2],0,2);


                   if any(BASE_AVG==0)>0 %if any of the baselines are zero and removes the data
                        [col]=find(BASE_AVG==0);
                        Temp_DATA(:,col)=[];
                        BASE_SD(:,col)=[];
                        BASE_AVG(:,col)=[];
                        clear row col
                   end

                    Temp_Z= bsxfun(@minus, Temp_DATA, BASE_AVG); 
                    Temp_Z= bsxfun(@rdivide, Temp_Z, BASE_SD);         

                %smooth the z-scored data
                    sm=smooth(Temp_DATA); clear COne
                    sm_reshaped=reshape(sm,80,size(Temp_DATA,2));
                        Temp_zSDsm=std(sm_reshaped,0,2);
                        Temp_ZSEsm=Temp_zSDsm ./ sqrt(size(Temp_DATA,2));


                             SMOOTHcombase{ses,j,reg,trt}=sm_reshaped; clear sm sm_reshaped
                             SEMsmoothcombase{ses,j,reg,trt}=Temp_ZSEsm; clear Temp*
                    end
                end
            end
        end
 
              SINGLEUNIT.rawcombase=SMOOTHcombase;
              SINGLEUNIT.rawcombaseSEM=SEMsmoothcombase;
              


%% GRAPH
%Types 1= correct -> correct
%      2= correct -> incorrect
%      3= incorrect -> incorrect
%      4= incorrect-> correct


for ses=1:4 %in order: D85, RD1, RD3, RD4, R50
   % for trt=1:2 %trt 1 is always control

    type1=mean(cell2mat(SINGLEUNIT.rawZtcticsmooth(ses,1)),2);
        type1sem=cell2mat(SINGLEUNIT.rawZtcticsmoothSEM(ses,1));
    type2=mean(cell2mat(SINGLEUNIT.rawZtcticsmooth(ses,2)),2);
        type2sem=cell2mat(SINGLEUNIT.rawZtcticsmoothSEM(ses,2));
%     type3=mean(cell2mat(SINGLEUNIT.rawcombase(ses,3,1)),2);
%         type3sem=cell2mat(SINGLEUNIT.rawcombaseSEM(ses,3,1));
%     type4=mean(cell2mat(SINGLEUNIT.rawcombase(ses,4,1)),2);
%         type4sem=cell2mat(SINGLEUNIT.rawcombaseSEM(ses,4,1));

%     type1t2=mean(cell2mat(SINGLEUNIT.rawZtcticsmooth(ses,1,2)),2);
%         type1t2sem=cell2mat(SINGLEUNIT.rawZtcticsmoothSEM(ses,1,2));
%     type2t2=mean(cell2mat(SINGLEUNIT.rawZtcticsmooth(ses,2,2)),2);
%         type2t2sem=cell2mat(SINGLEUNIT.rawZtcticsmoothSEM(ses,2,2));
%     type3t2=mean(cell2mat(SINGLEUNIT.rawcombase(ses,3,2)),2);
%         type3t2sem=cell2mat(SINGLEUNIT.rawcombaseSEM(ses,3,2));
%     type4t2=mean(cell2mat(SINGLEUNIT.rawcombase(ses,4,2)),2);
%         type4t2sem=cell2mat(SINGLEUNIT.rawcombaseSEM(ses,4,2));
   
        x=[1:80]';
        
    figure %for corrects
    hold on
 
        boundedline(x,type1,type1sem,'k','alpha'); 
        %boundedline(x,type4,type4sem,'r','alpha');
        
       % boundedline(x,type1t2,type1t2sem,'b','alpha'); 
        %boundedline(x,type4t2,type4t2sem,'m','alpha');
        
        set(gca,'ylim',[-1 2.5]);
        title(['Session:' num2str(ses)  ])
      
    figure %for incorrects
    hold on   
   
        boundedline(x,type2,type2sem,'k','alpha'); 
        %boundedline(x,type3,type3sem,'y','alpha'); 
        
        %boundedline(x,type2t2,type2t2sem,'g','alpha'); 
        %boundedline(x,type3t2,type3t2sem,'r','alpha'); 
        
        set(gca,'ylim',[-.4 1.2]);
        title(['Session:' num2str(ses) ])
        
    %end
 
end


%GRAPH
for j=2:2
    if j==1, DATA=SINGLEUNIT.phasiczin; SEM=SINGLEUNIT.semzin;
    elseif j==2, DATA=SINGLEUNIT.phasiczde; SEM=SINGLEUNIT.semzde;
    end
    for ses=1:4

        type1=mean(cell2mat(DATA(ses,1)),2);
            type1sem=cell2mat(SEM(ses,1));
            
        type2=mean(cell2mat(DATA(ses,2)),2);
            type2sem=cell2mat(SEM(ses,2));
            
        type3=mean(cell2mat(DATA(ses,3)),2);
            type3sem=cell2mat(SEM(ses,3));
            
        type4=mean(cell2mat(DATA(ses,4)),2);
            type4sem=cell2mat(SEM(ses,4));


            x=[1:80]';

        figure
        hold on
            boundedline(x,type1,type1sem,'k','alpha'); 
            boundedline(x,type2,type2sem,'b','alpha'); 
            boundedline(x,type3,type3sem,'r','alpha'); 
            boundedline(x,type4,type4sem,'g','alpha'); 
            %set(gca,'ylim',[-1 5]);

    end
end


%BASELINE Bar Graph Comparison




%%
%z-score phasic neurons
for j=1:4
    for ses=1:4

        Temp_DATA=cell2mat(SINGLEUNIT.rawphasicout(ses,j));
            BASE_AVG=mean(Temp_DATA(base,:),1);
            BASE_SD=std(Temp_DATA(base,:),0,1);
%         Temp_IN=cell2mat(SINGLEUNIT.phasicinraw(ses,j));
%             BASE_IN=mean(Temp_IN(base,:),1);
%             SD_IN=std(Temp_IN(base,:),0,1);
%         Temp_DE=cell2mat(SINGLEUNIT.phasicdecraw(ses,j));
%             BASE_DE=mean(Temp_DE(base,:),1);
%             SD_DE=std(Temp_DE(base,:),0,1);

        
       if any(BASE_AVG==0)>0 %if any of the baselines are zero and removes the data
            [col]=find(BASE_AVG==0);
            Temp_DATA(:,col)=[];
            BASE_SD(:,col)=[];
            BASE_AVG(:,col)=[];
            clear row col
       end
       
%         if any(BASE_IN==0)>0 %if any of the baselines are zero and removes the data
%             [col]=find(BASE_IN==0);
%             Temp_IN(:,col)=[];
%             BASE_IN(:,col)=[];
%             SD_IN(:,col)=[];
%             clear row col
%         end
%        
%         if any(BASE_DE==0)>0 %if any of the baselines are zero and removes the data
%             [col]=find(BASE_DE==0);
%             Temp_DE(:,col)=[];
%             BASE_DE(:,col)=[];
%             SD_DE(:,col)=[];
%             clear row col
%         end

        
            Temp_Z= bsxfun(@minus, Temp_DATA, BASE_AVG); 
            Temp_Z= bsxfun(@rdivide, Temp_Z, BASE_SD);
               % Temp_zSD=std(Temp_Z,0,2);
               % Temp_ZSE=Temp_zSD ./ sqrt(size(Temp_Z,2));
%             
%             Temp_ZIN= bsxfun(@minus, Temp_IN, BASE_IN); 
%             Temp_ZIN= bsxfun(@rdivide, Temp_ZIN, SD_IN);
%             
%             Temp_ZDE= bsxfun(@minus, Temp_DE, BASE_DE); 
%             Temp_ZDE= bsxfun(@rdivide, Temp_ZDE, SD_DE);
            
              
        %smooth the z-scored data
            sm=smooth(Temp_Z); clear COne
            sm_reshaped=reshape(sm,80,size(Temp_Z,2));
                Temp_zSDsm=std(sm_reshaped,0,2);
                Temp_ZSEsm=Temp_zSDsm ./ sqrt(size(Temp_Z,2));
                    PHASICZsmooth{ses,j}=sm_reshaped; clear sm sm_reshaped
                    SEMZsmooth{ses,j}=Temp_ZSEsm; 
            

%             sm=smooth(Temp_ZIN); clear COne
%             sm_reshaped=reshape(sm,80,size(Temp_ZIN,2));
%                 Temp_zINSDsm=std(sm_reshaped,0,2);
%                 Temp_ZINSEsm=Temp_zINSDsm ./ sqrt(size(Temp_ZIN,2));
%                     PHASICZin{ses,j}=sm_reshaped; clear sm_reshaped
%                     SEMZin{ses,j}=Temp_ZINSEsm;
%                     
%             sm=smooth(Temp_ZDE); clear COne
%             sm_reshaped=reshape(sm,80,size(Temp_ZDE,2));
%                 Temp_zDESDsm=std(sm_reshaped,0,2);
%                 Temp_ZDESEsm=Temp_zDESDsm ./ sqrt(size(Temp_ZDE,2));
%                     PHASICZde{ses,j}=sm_reshaped; clear sm_reshaped
%                     SEMZde{ses,j}=Temp_ZDESEsm;
               
            
            
            %PHASICZ{ses,j}=Temp_Z; %z-scored data to each individual baseline
            %STDZ{ses,j}=Temp_zSD;
            %SEMZ{ses,j}=Temp_ZSE;
            
    
      clear Temp* BASE*
    end 

end

SINGLEUNIT.phasic=PHASIC;
SINGLEUNIT.phasicZ=PHASICZ;
SINGLEUNIT.stdz=STDZ;
SINGLEUNIT.semz=SEMZ;
SINGLEUNIT.phasicZsmoothout=PHASICZsmooth;
SINGLEUNIT.semzsmoothout=SEMZsmooth;
SINGLEUNIT.phasiczin=PHASICZin;
    SINGLEUNIT.semzin=SEMZin;
SINGLEUNIT.phasiczde=PHASICZde;
    SINGLEUNIT.semzde=SEMZde;
    
    
  
%GRAPH
for ses=1:1

    type1=mean(cell2mat(SINGLEUNIT.phasicZsmoothout(ses,1)),2);
        type1sem=cell2mat(SINGLEUNIT.semzsmoothout(ses,1));
    type2=mean(cell2mat(SINGLEUNIT.phasicZsmoothout(ses,2)),2);
        type2sem=cell2mat(SINGLEUNIT.semzsmoothout(ses,2));
%     type3=mean(cell2mat(SINGLEUNIT.phasicZsmoothout(ses,3)),2);
%         type3sem=cell2mat(SINGLEUNIT.semzsmoothout(ses,3));
%     type4=mean(cell2mat(SINGLEUNIT.phasicZsmoothout(ses,4)),2);
%         type4sem=cell2mat(SINGLEUNIT.semzsmoothout(ses,4));
%     
    
        x=[1:80]';
        
    figure
    hold on
        boundedline(x,type1,type1sem,'b','alpha'); 
        boundedline(x,type2,type2sem,'g','alpha'); 
        boundedline(x,type3,type3sem,'y','alpha'); 
        boundedline(x,type4,type4sem,'r','alpha'); 
        %set(gca,'ylim',[-1 5]);
 
end

%% for Data containing all phasic and non-phasic neurons (and outliers)


base=[1:20];
post=[21:80];

for j=1:4
    for ses=1:4
        
         Temp_DATA=cell2mat(SINGLEUNIT.comspike(ses,j));
            Temp_DATA(isnan(Temp_DATA))=0;
            
            BASE_AVG=mean(Temp_DATA(base,:),1);
            BASE_SD=std(Temp_DATA(base,:),0,1);
            
            
           if any(BASE_AVG==0)>0 %if any of the baselines are zero and removes the data
                [col]=find(BASE_AVG==0);
                Temp_DATA(:,col)=[];
                BASE_SD(:,col)=[];
                BASE_AVG(:,col)=[];
                clear row col
           end
            
            Temp_Z= bsxfun(@minus, Temp_DATA, BASE_AVG); 
            Temp_Z= bsxfun(@rdivide, Temp_Z, BASE_SD);
%                 Temp_zSD=std(Temp_Z,0,2);
%                 Temp_ZSE= Temp_zSD ./ sqrt(size(Temp_Z,2));
             

        %smooth the z-scored data
            sm=smooth(Temp_Z); clear COne
            sm_reshaped=reshape(sm,80,size(Temp_Z,2));
                Temp_zSDsm=std(sm_reshaped,0,2);
                Temp_ZSEsm=Temp_zSDsm ./ sqrt(size(Temp_Z,2));
    
                   %RAWSMOOTH{ses,j}=Temp_Z; clear sm sm_reshaped 
                    %SEMRAWsmooth{ses,j}=Temp_ZSE; clear Temp*
               
                     RAWSMOOTH{ses,j}=sm_reshaped; clear sm sm_reshaped
                     SEMRAWsmooth{ses,j}=Temp_ZSEsm; clear Temp*
    end
        
end

 %remove single neuron that is causing huge error --> on RD1 type 4 (I-C/ lose-shift)
    
    tempwave=cell2mat(RAWSMOOTH(2,4));
        tempwave(:,31)=[];
        RAWSMOOTH{2,4}=tempwave;
        tempsd=std(tempwave,0,2);
            tempsem=tempsd ./ sqrt(size(tempwave,2));
            SEMRAWsmooth{2,4}=tempsem;
            clear temp*
            

%GRAPH
for ses=1:4

    type1=mean(cell2mat(RAWSMOOTH(ses,1)),2);
        type1sem=cell2mat(SEMRAWsmooth(ses,1));
    type2=mean(cell2mat(RAWSMOOTH(ses,2)),2);
        type2sem=cell2mat(SEMRAWsmooth(ses,2));
    type3=mean(cell2mat(RAWSMOOTH(ses,3)),2);
        type3sem=cell2mat(SEMRAWsmooth(ses,3));
    type4=mean(cell2mat(RAWSMOOTH(ses,4)),2);
        type4sem=cell2mat(SEMRAWsmooth(ses,4));
    
    
        x=[1:80]';
        
    figure
    hold on
%     plot(x,type1,'b')
%     plot(x,type2,'g')
%     plot(x,type3,'y')
%     plot(x,type4,'r')
        boundedline(x,type1,type1sem,'b','alpha'); 
        boundedline(x,type2,type2sem,'g','alpha'); 
        boundedline(x,type3,type3sem,'y','alpha'); 
        boundedline(x,type4,type4sem,'r','alpha'); 
        set(gca,'ylim',[-1 2.5]);
 
end

SINGLEUNIT.phasicZsmoothout=Rawout;
SINGLEUNIT.semzsmoothout=SEM;

%GRAPH
for j=2:2
    if j==1, DATA=SINGLEUNIT.phasiczin; SEM=SINGLEUNIT.semzin;
    elseif j==2, DATA=SINGLEUNIT.phasiczde; SEM=SINGLEUNIT.semzde;
    end
    for ses=1:4

        type1=mean(cell2mat(DATA(ses,1)),2);
            type1sem=cell2mat(SEM(ses,1));
            
        type2=mean(cell2mat(DATA(ses,2)),2);
            type2sem=cell2mat(SEM(ses,2));
            
        type3=mean(cell2mat(DATA(ses,3)),2);
            type3sem=cell2mat(SEM(ses,3));
            
        type4=mean(cell2mat(DATA(ses,4)),2);
            type4sem=cell2mat(SEM(ses,4));


            x=[1:80]';

        figure
        hold on
            boundedline(x,type1,type1sem,'k','alpha'); 
            boundedline(x,type2,type2sem,'b','alpha'); 
            boundedline(x,type3,type3sem,'r','alpha'); 
            boundedline(x,type4,type4sem,'g','alpha'); 
            %set(gca,'ylim',[-1 5]);

    end
end

%% Counting how neurons fires to different types --> glyph graphs

for ses=1:4
    temp=cell2mat(SINGLEUNIT.phasic(ses,1))';  
    allphasic=zeros(size(temp,1),1);
    for j=1:4 %concatenate all 4 types 
        temp=cell2mat(SINGLEUNIT.phasic(ses,j))';  
        allphasic=[allphasic , temp];
    end
     allphasic(:,1)=[];
     comphasic{ses}=allphasic; clear allphasic temp
end

    
    %Types 1= correct -> correct
%      2= correct -> incorrect
%      3= incorrect -> incorrect
%      4= incorrect-> correct

    
% define different types of neurons
Type={'TOTAL';'ALL'; 'CORRECT'; 'C-C'; 'I-C'; 'C-I'; 'I-I'; 'INCORRECT';'SWITCH'; 'NONE'};
for ses=1:4
    allphasic=cell2mat(comphasic(ses));
        
    TOTAL(1,ses)=size(allphasic, 1);
    
    ALLP(1,ses)= sum( all(allphasic==1,2) );   
    NONE(1,ses)=sum( all(allphasic==0,2) );
    
    
    CORR(1,ses)=   size(  find(allphasic(:,1)==1 & allphasic(:,4)==1 & allphasic(:,2)==0 & allphasic(:,3)==0)  ,1);
    INCORR(1,ses)= size(  find(allphasic(:,2)==1 & allphasic(:,3)==1 & allphasic(:,1)==0 & allphasic(:,4)==0)  ,1);
   
    T1(1,ses)= size(  find(allphasic(:,1)==1 & allphasic(:,2)==0 & allphasic(:,3)==0 & allphasic(:,4)==0)  ,1);
    T2(1,ses)= size(  find(allphasic(:,1)==0 & allphasic(:,2)==1 & allphasic(:,3)==0 & allphasic(:,4)==0)  ,1);
    T3(1,ses)= size(  find(allphasic(:,1)==0 & allphasic(:,2)==0 & allphasic(:,3)==1 & allphasic(:,4)==0)  ,1);
    T4(1,ses)= size(  find(allphasic(:,1)==0 & allphasic(:,2)==0 & allphasic(:,3)==0 & allphasic(:,4)==1)  ,1);
    
    SWITCH(1,ses)=size(  find(allphasic(:,1)==0 & allphasic(:,2)==1 & allphasic(:,3)==0 & allphasic(:,4)==1)  ,1);


end

% PUT INTO AN EASY TO READY TABLE
DATA=cat(1,TOTAL, ALLP, CORR,T1, T4, T2, T3, INCORR, SWITCH, NONE);
    clear TOTAL NONE all CORR INCORR T1 T2 T3 T4 SWITCH
tyCOUNT=table(DATA, 'RowNames', Type); clear DATA

SINGLEUNIT.TypeCount=tyCOUNT;

%Make cool graphs 
b=table2array(tyCOUNT);
    prop=b(2:9,:);
    for row=1:8
        for col=1:4
            proportion(row,col)= prop(row,col)/ b(1,col);
        end
    end
    
figure
hold on
 labels={'ses 1'}; %, 'ses2','ses3', 'ses4'};
 vartype={'ALL'; 'CORRECT'; 'C-C'; 'I-C'; 'C-I'; 'I-I'; 'INCORRECT';'SWITCH';};
 figure
glyphplot(proportion(:,1)','standardize','off','varlabels',vartype) % ,'obslabels',labels);




%% PCA -> do on ALL variable --> also need to delete region in all 


%Run PCA
for j=1:2
   for ses=1:5
       for reg=1:2
           for trt=1:2

        Temp_DATA=SINGLEUNIT.rawZtcticsmooth{ses,j, reg,trt}'; %recall data
         [coeff,score,latent,tsquared,explained,mu]=pca(Temp_DATA); %run PCA
        
        COEFF{ses,j,reg,trt}= coeff; %store important data
        EXPLANINED{ses,j,reg,trt}=explained; %store important data
        
        clear coeff score latent tsquared explained mu Temp_DATA
           end
       end
   end
end 


SINGLEUNIT.PCAcoeff=COEFF;
SINGLEUNIT.PCAexp=EXPLANINED;


%Assign units to PCA value - UnitType based on Least Squares fit
for j=1:2
    for ses=1:5
        for reg=1:2
            for trt=1:2
        Temp_DATA=SINGLEUNIT.comspiketctic{ses,j,reg,trt}';
        coeff=SINGLEUNIT.PCAcoeff{ses,j,reg,trt};
        
%         figure 
%         plot(1:80, Temp_DATA);
%         plot(1:80, coeff);
%         
        for t = 1:size(Temp_DATA,1)
            
            % find sum of squares for each PC for this unit
            for g=1:size(coeff,2)
                Diffsq(:,g)=(Temp_DATA(t,:)-coeff(:,g)').^2;
            end
                ssquares=sum(Diffsq,1);

            % find minimum ssquares and set UnitType equal to its number
            minimum=min(ssquares);
            [col]=find(ssquares==minimum);
           
            if size(col,2)>1
             UnitType(t,:)=col(1,1);
            elseif size(col,2)==1
             UnitType(t,:)=col;
            end
            
            Unit{ses,j,reg,trt}=UnitType; %save important thing...
            
            clear minimum col Diffsq ssquares 

        end
        clear Temp_DATA coeff UnitType
            end
        end
    end
end


%Extract units of "type 1 similarity"
for j=1:2
    for ses=1:5
        for reg=1:2
            for trt=1:2
        
        Temp_DATA=SINGLEUNIT.comspiketctic{ses,j,reg,trt}';
        UnitType=Unit{ses,j,reg,trt};
        len=size(UnitType,1);

        row=1;
        for u=1:len
            if UnitType(u,1)==1

                COne(row,:)=Temp_DATA(u,:);
                row=row+1;

            end

        end
        
        %finding standard error of the mean
            N=size(COne,1);
            %deviation=std(COne,0,1);
            %deviation=deviation/sqrt(N);

            %Dev{ses,j}=deviation; clear deviation
            %ComOne{ses,j}=COne; 
        
        %smooth the curve and find smoothed SEM
            smoothCOne=smooth(COne'); clear COne
            COne_sm=reshape(smoothCOne,80,N)';

            deviation=std(COne_sm,0,1);
            deviation=deviation/sqrt(N);

            Devsm{ses,j,reg,trt}=deviation; clear deviation N
            ComOnesm{ses,j,reg,trt}=COne_sm; clear COne_sm smoothCOne

        
      clear UnitType Temp_DATA   
            end
        end
    end
end

SINGLEUNIT.PCAcom1smooth=ComOnesm;
SINGLEUNIT.PCAcom1semsmooth=Devsm;


%Plot

for ses=1:5
    for reg=1:2
        %for trt=1:2
    
    type1=mean(cell2mat(SINGLEUNIT.PCAcom1smooth(ses,1,reg,1)),1);
        type1sem=cell2mat(SINGLEUNIT.PCAcom1semsmooth(ses,1,reg,1));
    type2=mean(cell2mat(SINGLEUNIT.PCAcom1smooth(ses,2,reg,1)),1);
        type2sem=cell2mat(SINGLEUNIT.PCAcom1semsmooth(ses,2,reg,1));
%     type3=mean(cell2mat(SINGLEUNIT.PCAcom1smooth(ses,3)),1);
%         type3sem=cell2mat(SINGLEUNIT.PCAcom1semsmooth(ses,3));
%     type4=mean(cell2mat(SINGLEUNIT.PCAcom1smooth(ses,4)),1);
%         type4sem=cell2mat(SINGLEUNIT.PCAcom1semsmooth(ses,4));


    type1exp=mean(cell2mat(SINGLEUNIT.PCAcom1smooth(ses,1,reg,2)),1);
        type1expsem=cell2mat(SINGLEUNIT.PCAcom1semsmooth(ses,1,reg,2));
    type2exp=mean(cell2mat(SINGLEUNIT.PCAcom1smooth(ses,2,reg,2)),1);
        type2expsem=cell2mat(SINGLEUNIT.PCAcom1semsmooth(ses,2,reg,2));

xaxis=1:80;

figure
hold on
boundedline(xaxis,type1,type1sem,'k','alpha'); 
boundedline(xaxis,type2,type2sem,'b','alpha');
% boundedline(xaxis,type3,type3sem,'r','alpha');
% boundedline(xaxis,type4,type4sem,'g','alpha');
title(['Session:' num2str(ses) 'Region:' num2str(reg) ])
set(gca,'ylim', [0 .7]);

figure
hold on
boundedline(xaxis,type1exp,type1expsem,'k','alpha'); 
boundedline(xaxis,type2exp,type2expsem,'g','alpha');
title(['Session:' num2str(ses) 'Region:' num2str(reg) ])
set(gca,'ylim', [0 .7]);


      %  end
    end
end





%% baseline by time graph --> other stuff
    base=[1:20];
    post=[21:80];
         for j=1:2
            for ses=1:5
                    for trt=1:2

                 %Temp_DATA=SINGLEUNIT.comspike{ses,j,reg,trt};
                 Temp_DATA=SINGLEUNIT.comspiketctic{ses,j,trt};
                    Temp_DATA(isnan(Temp_DATA))=0;

                    BASE_AVG=mean(Temp_DATA(base,:),1);
                   % BASE_SD=std(BASE_AVG,0,2);


                   if any(BASE_AVG==0)>0 %if any of the baselines are zero and removes the data
                        [col]=find(BASE_AVG==0);
                        Temp_DATA(:,col)=[];
                        %BASE_SD(:,col)=[];
                        BASE_AVG(:,col)=[];
                        clear row col
                   end
                   
                   baseline=Temp_DATA(base,:); %baseline only with all zeros and nans removed
                   
                BASE_SD=std(baseline,0,2);
                    BASE_SEM=BASE_SD./ sqrt(size(BASE_AVG,2)); %


                BASELINE{ses,j,trt}=mean(baseline,2); clear sm sm_reshaped BASE_AVG
                %BASESD(ses,j,trt)=BASE_SD; clear Temp* BASE_SD
                BASESEM{ses,j,trt}=BASE_SEM; clear BASE_SEM 
                    end
    
            end
         end
        
         SINGLEUNIT.basetimetctic=BASELINE;
         SINGLEUNIT.basetimeSEMtctic=BASESEM;
         %graph of the baselines

for ses=1:2
    for type=1:2
      
            
            BASEcontr=SINGLEUNIT.basetimetctic{ses,type,1};
                BASEcontrSEM=SINGLEUNIT.basetimeSEMtctic{ses,type,1};
            BASEexpr=SINGLEUNIT.basetimetctic{ses,type,2};
                BASEexprSEM=SINGLEUNIT.basetimeSEMtctic{ses,type,2};
            
        x=[1:20]';

        figure
        hold on
            boundedline(x,BASEcontr,BASEcontrSEM,'k','alpha'); 
            boundedline(x,BASEexpr,BASEexprSEM,'b','alpha'); 
           % boundedline(x,type3,type3sem,'r','alpha'); 
            %boundedline(x,type4,type4sem,'g','alpha'); 
            %set(gca,'ylim',[-1 5]);

            title(['Session:' num2str(ses)  'Type:' num2str(type) ])

    end
end

  
        