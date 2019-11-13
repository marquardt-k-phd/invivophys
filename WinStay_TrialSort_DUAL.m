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

load SINGLEUNIT_PAE_DUAL_all
         save('SINGLEUNIT_PAE_DUAL_all','SINGLEUNIT');
% Input save file name ^^^^  

%% Get Data -> import directly from neuroexplorer 

%Create NeuroExploer files as have done previously to make TC and TIC
%timepoints

%under trial bin counts -> MatLab tab and check the box labeled 'Send
%matrix of results to Matlab'

    % you may also wish to send 'All Numberical Results' to Excel as an extra
    % output or check

% NOTE: the mice must be exported in order and organized in following
% manner

%spikes(mouse, beh, sesion)
%mouse is simply the number you are inputting
%beh indicated TC(1) or TIC(2) 
%Session indicated: 1=D85 2=R25 3=R50 4=R85 etc...
%treatment (trt) is 1=control 2=treatment or genotype --> if all WT leave at 1
    %DO ALL OF ONE TRT FIRST THEN START OVER AGAIN AT MOUSE #1 FOR SECOND
    %TRT
    
  %*****if re-opening singleunit file then need to re-assign spikes variable so does not become over-written  
  %spikes=SINGLEUNIT.neuroexp;
  
%ctrl+enter will run this selection only


    mouse=14;
    beh=2;
    sess=5;
    trt=2;

DS{1,1}='SPK01'; DS{1,2}='SPK02';  DS{1,3}='SPK03';  DS{1,4}='SPK04';  DS{1,5}='SPK05';  DS{1,6}='SPK06';  DS{1,7}='SPK07';  DS{1,8}='SPK08'; 
OFC{1,1}='SPK09'; OFC{1,2}='SPK10'; OFC{1,3}='SPK11'; OFC{1,4}='SPK12'; OFC{1,5}='SPK13'; OFC{1,6}='SPK14'; OFC{1,7}='SPK15'; OFC{1,8}='SPK16';  
    for spk=1:8
        ds=DS{1,spk};
        ofc=OFC{1,spk};
        TF=strfind(nexColumnNames, ds);
        dsneuroncols{spk,:}=find(~cellfun(@isempty,TF)); clear TF
        TF=strfind(nexColumnNames, ofc);
        ofcneuroncols{spk,:}=find(~cellfun(@isempty,TF)); clear TF
    end
    dscols=zeros(1,1);
    ofccols=zeros(1,1);
    for spk=1:8
        ds=cell2mat(dsneuroncols(spk,1));
        dscols=[dscols ds]; clear ds
        
        ofc=cell2mat(ofcneuroncols(spk,1));
        ofccols=[ofccols ofc]; clear ofc
    end
        dscols(:,1)=[];
        ofccols(:,1)=[];

    spikes{mouse,beh,sess,1,trt}=nex(:,dscols);
    spikes{mouse,beh,sess,2,trt}=nex(:,ofccols);
    clear dscols dsneuroncols nex nexColumnNames ofccols ofcneuroncols spk 
    
    SINGLEUNIT.neuroexp=spikes;
    

%% Define Combined Files

%should be defined in your 'Variables' script under var4= location of file
%and var5= mouse number


%% read in combined file data... use ML as correct vs. incorrect measure
%Select Stage 9 Entry for each mouse: do not select the last row of zeros if it is there
ML=cell(1,5,2);
for ses=1:5 
    for trt=1:2
        Num=size(SINGLEUNIT.neuroexp(:,:,ses,:,trt),1);     
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

%% DEFINE TRIAL TYPES based on ML

for ses=1:5
    for trt=1:2
    Num=size(SINGLEUNIT.neuroexp(:,:,ses,:,trt),1);  %number of mice
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
        for j=2:size(code,1) %use size of combined file ML b/c sometime x-tra timestamps in rec file

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


%split spikes based on trial type
 for beh=1:2
     for ses=1:5
         for trt=1:2
         Num=size(SINGLEUNIT.neuroexp(:,:,ses,:,trt),1);  
         for reg=1:2
         for k=1:Num
       
         spikes=SINGLEUNIT.neuroexp;
    temp=cell2mat(spikes(k,beh,ses,reg,trt)); %SINGLEUNIT.neuroexp or spikes
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
    
    if len==0 %to put in when there are no trials
        SPIKES{k,beh,ses,reg,trt}= 0; clear y
    elseif len>0
        SPIKES{k,beh,ses,reg,trt}=Neurons;
    end
        clear  num 
        
   %sort into corrects & incorrects only
   if beh==1
       if len==0
           N=0;
       elseif len>0
           TC{k,1,ses,reg,trt}=mean(Neurons,1);
           %TCnum{k,1,ses,reg,trt}=Neurons;
       end
   elseif beh==2
       if len==0
           M=0;
       elseif len>0
           TIC{k,1,ses,reg,trt}=mean(Neurons,1);
           %TICnum{k,1,ses,reg,trt}=Neurons;
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
            TYPE{k,1,ses,reg,trt}=Neurons(Z,:,:);
        S=find(stuff==4);
            TYPE{k,4,ses,reg,trt}=Neurons(S,:,:);
       end
            
   elseif beh==2
        if len==0
            T=0;
            U=0;
       elseif len>0
           if size(stuff,1) > size(Neurons,1)
              over= size(stuff,1) - size(Neurons,1);
              last=size(stuff,1);
              lastover=last-over +1;
              stuff(lastover:last,:)=[];
          end
        T=find(stuff==3);
            TYPE{k,3,ses,reg,trt}=Neurons(T,:,:);
        U=find(stuff==2);
            TYPE{k,2,ses,reg,trt}=Neurons(U,:,:);
        
        end
   end
    clear Z S T U Neurons len

         end
         end
        end
     end
 end
 
%  emptysub=double.empty(0,80,4);
%     TYPE{4,1,2}=emptysub;
%     TYPE{4,4,2}=emptysub;

 SINGLEUNIT.trialspikes=TYPE; %mouse x trial type x session x region x trt
 SINGLEUNIT.TCspikes=TC;
 SINGLEUNIT.TICspikes=TIC;
 
         save('SINGLEUNIT_PAE_DUAL','SINGLEUNIT'); %**************************

%% Average Trials and Concatenate spikes


for ses=1:5
    for ty=1:4
        for trt=1:2
            for reg=1:2
     allspikes=zeros(80,1);
     com=zeros(80,1);
     Num=size(SINGLEUNIT.trialspikes(:,:,ses,reg,trt),1);
     
    for j=1:Num
        temp=cell2mat(SINGLEUNIT.trialspikes(j,ty,ses,reg,trt));
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
            if ty==1, temp=TC{j,1,ses,reg,trt};
            elseif ty==2,  temp=TIC{j,1,ses,reg,trt};
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
    ALL{ses,ty,reg,trt}=allspikes; %all types
    COM{ses,ty,reg,trt}=com;  %tc and tic only?
    clear allspikes temp*
            end
        end
    end
end

SINGLEUNIT.comspike=ALL; %session x trial type (time x neuron) x region x trt
SINGLEUNIT.comspiketctic=COM; % session x type (tc tic) x region x trt



        save('SINGLEUNIT_PAE_DUAL','SINGLEUNIT'); %************************
%%  THE FOLLOWING ARE DIFFERENT TYPES OF GRAPHS TO MAKE

%define baseline and post time bins
base=[1:20];
post=[21:80];

%stuff=SINGLEUNIT.rawZsmooth{5,4,2,2}';
   %SESSION X TYPE X REGION X TRT

%RAW FIRING RATE PHASICS ONLY w/ separated increasing and decreasing 
        %calculate baseline average
        for j=1:2
           for ses=1:5
               for reg=1:2
                   for trt=1:2

            %Temp_DATA=cell2mat(SINGLEUNIT.comspike(ses,j,reg,trt));
            Temp_DATA=SINGLEUNIT.comspiketctic{ses,j,reg,trt};

            sig=ttest2(Temp_DATA(base,:), Temp_DATA(post,:)); %0=nonphasic 1=phasic

                sig(isnan(sig))=0;
            PHASIC{ses,j}=sig;

            z=find(sig==1);
                phasic=Temp_DATA(:,z);
                phasic_base=mean(phasic(base,:),1);
                phasic_post=mean(phasic(post,:),1);

                increasing=find(phasic_post > phasic_base);
                decreasing=find(phasic_post < phasic_base);


            PHASICFIRING{ses,j,reg,trt}=Temp_DATA(:,z);
            PHASICIN{ses,j,reg,trt}=phasic(:,increasing);
            PHASICDE{ses,j,reg,trt}=phasic(:,decreasing);

                   end
               end
           end

            clear Temp_DATA phasic increasing decreasing z sig phasic_base phasic_post 
        end

        SINGLEUNIT.phasicrawtctic=PHASICFIRING;  SINGLEUNIT.phasicrawtypes=PHASICFIRING;
        SINGLEUNIT.phasicrawDEC=PHASICDE;
        SINGLEUNIT.phasicrawINC=PHASICIN; 
       
        
        
%removal of outliers
        for j=1:2
            for ses=1:5
                for reg=1:2
                    for trt=1:2

                %Temp_DATA=cell2mat(SINGLEUNIT.phasicraw(ses,j,reg,trt));
                Temp_DATA=SINGLEUNIT.comspiketctic{ses,j,reg,trt};
                postavg=mean(Temp_DATA(post,:),1);
                    popavg=mean(Temp_DATA(post,:),2);
                    popstd=std(popavg,0,1);

                posout=mean(popavg,1)  +  2*popstd ;
                negout=mean(popavg,1)  -  2*popstd ;

                outliers=find(postavg>posout | postavg<negout);
                
                Temp_DATA(:,outliers)=[];
                sd=std(Temp_DATA,0,2);
                sem= sd ./ sqrt(size(Temp_DATA,2));

                Rawout{ses,j,reg,trt}=Temp_DATA;
                SEM{ses,j,reg,trt}=sem;

                clear Temp* posout neg out outliers postavg popavg popstd

                    end
                end
            end
        end

            SINGLEUNIT.rawOUTREMOVED=Rawout;
            

    
% PHASIC NEURONS ONLY Z-SCORED (NO SEPARATION INTO INCREASING AND DECREASING)
%z-score phasic neurons
%Note: commented sections are for sorting out increasing and decreasing
%phasic neurons seperately
        for j=1:4
            for ses=1:5
                for reg=1:2
                    for trt=1:2

                %Temp_DATA=cell2mat(SINGLEUNIT.phasicrawOUTREMOVED(ses,j,reg,trt));
                %Temp_DATA=SINGLEUNIT.phasicrawtctic{ses,j,reg,trt};
                Temp_DATA=SINGLEUNIT.phasicrawrcticOUTREMOVED{ses,j,reg,trt};
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
                            PHASICZsmooth{ses,j,reg,trt}=sm_reshaped; clear sm sm_reshaped
                            SEMZsmooth{ses,j,reg,trt}=Temp_ZSEsm; 


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
        end

        SINGLEUNIT.phasicZtcticsmoothOUTREMOVED=PHASICZsmooth;
            SINGLEUNIT.phasicZtcticsmoothSEMOUTREMOVED=SEMZsmooth;
        %SINGLEUNIT.phasiczin=PHASICZin;
           % SINGLEUNIT.semzin=SEMZin;
        %SINGLEUNIT.phasiczde=PHASICZde;
            %SINGLEUNIT.semzde=SEMZde;
            
%SESSION X TYPE X REGION X TRT
  
        %GRAPH
        for ses=1:5
            for reg=2:2
               % for trt=1:2

            type1=mean(SINGLEUNIT.phasicZtcticsmooth{ses,1,reg,1},2);
                type1sem=SINGLEUNIT.phasicZtcticsmoothSEM{ses,1,reg,1};
            type1exp=mean(SINGLEUNIT.phasicZtcticsmooth{ses,1,reg,2},2);
                type1expsem=SINGLEUNIT.phasicZtcticsmoothSEM{ses,1,reg,2};
                
            type2=mean(SINGLEUNIT.phasicZtcticsmooth{ses,2,reg,1},2);
                type2sem=SINGLEUNIT.phasicZtcticsmoothSEM{ses,2,reg,1};
            type2exp=mean(SINGLEUNIT.phasicZtcticsmooth{ses,2,reg,2},2);
                type2expsem=SINGLEUNIT.phasicZtcticsmoothSEM{ses,2,reg,2};
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
                set(gca,'ylim',[-2 7]);
                
            figure
            hold on
                boundedline(x,type2,type2sem,'b','alpha'); 
                boundedline(x,type2exp,type2expsem,'g','alpha');
                plot(x,zeros(1,80),'k')
                title('Incorrect')
                set(gca,'ylim',[-1 2.5]);
%                 boundedline(x,type3,type3sem,'y','alpha'); 
%                 boundedline(x,type4,type4sem,'r','alpha'); 
                %set(gca,'ylim',[-1 5]);

               % end
            end
        end
        

        
% ALL NEURONS Z-SCORE(PHASIC AND NON-PHASIC)
%Note: only neurons with baseline zero are removed from data
        base=[1:20];
        post=[21:80];

        for j=1:2
            for ses=1:5
                for reg=1:2
                    for trt=1:2

                 Temp_DATA=cell2mat(SINGLEUNIT.comspiketctic(ses,j,reg,trt));
                 %Temp_DATA=SINGLEUNIT.rawOUTREMOVED{ses,j,reg,trt};
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

                             RAWSMOOTH{ses,j,reg,trt}=sm_reshaped; clear sm sm_reshaped
                             SEMRAWsmooth{ses,j,reg,trt}=Temp_ZSEsm; clear Temp*
                    end
                end
            end
        end

         SINGLEUNIT.rawZtcticsmooth=RAWSMOOTH;
         SINGLEUNIT.rawZtcticsmoothSEM=SEMRAWsmooth;
         
        %remove single neuron that is causing huge error --> on D85 TC
        temp=SINGLEUNIT.rawZtcticsmooth{5,2,1,2}';
          temp(6,:)=[];
         % temp(21,:)=[];
          temp=temp';
          SINGLEUNIT.rawZtcticsmooth{5,2,1,2}=temp;
         tempsd=std(temp,0,2);
         tempsem=tempsd ./ sqrt(size(temp,2));
            SINGLEUNIT.rawZtcticsmoothSEM{5,2,1,2}=tempsem;
            clear temp*
    
         
        
 %BASELINE COMPARISON 
         for j=1:2
            for ses=1:5
                for reg=1:2
                    for trt=1:2

                 Temp_DATA=SINGLEUNIT.comspiketctic{ses,j,reg,trt};
                 %Temp_DATA=SINGLEUNIT.rawOUTREMOVED{ses,j,reg,trt};
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


                BASELINE{ses,j,reg,trt}=BASE_AVG; clear sm sm_reshaped BASE_AVG
                BASESD(ses,j,reg,trt)=BASE_SD; clear Temp* BASE_SD
                BASESEM(ses,j,reg,trt)=BASE_SEM; clear BASE_SEM
                    end
                end
            end
        end
        
      SINGLEUNIT.baseuncorrtctic=BASELINE;
      %SINGLEUNIT.baseuncorrSD=BASESD;
      SINGLEUNIT.baseuncorrtcticSEM=BASESEM;
      

      %ses x type x reg x trt
      
  temp=SINGLEUNIT.baseuncorrtctic{5,2,2,2}';
      
      
  
%graph of the baselines

for ses=1:5
   % for type=1:4
        for reg=1:2
            
%             BASEcontr1=mean(SINGLEUNIT.baseuncorr{ses,type,1,1},2);
%             BASEexpr1=mean(SINGLEUNIT.baseuncorr{ses,type,1,2},2);
%             
%             BASEcontr2=mean(SINGLEUNIT.baseuncorr{ses,type,2,1},2);
%             BASEexpr2=mean(SINGLEUNIT.baseuncorr{ses,type,2,2},2);
%             
%             
%             figure
%             hold on
%             bar(1,BASEcontr1,'k')
%             bar(2,BASEexpr1,'b')
%             
%             bar(3,BASEcontr2,'k')
%             bar(4,BASEexpr2,'b')

%             title(['Session:' num2str(ses)  'Type:' num2str(type) ])

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

%               BASEcont1=mean(SINGLEUNIT.baseuncorr{1,type,reg,1},2);
%               BASEexpr1=mean(SINGLEUNIT.baseuncorr{1,type,reg,2},2);
%               
%               BASEcont2=mean(SINGLEUNIT.baseuncorr{2,type,reg,1},2);
%               BASEexpr2=mean(SINGLEUNIT.baseuncorr{2,type,reg,2},2);
%               
%               BASEcont3=mean(SINGLEUNIT.baseuncorr{3,type,reg,1},2);
%               BASEexpr3=mean(SINGLEUNIT.baseuncorr{3,type,reg,2},2);
%               
%               BASEcont4=mean(SINGLEUNIT.baseuncorr{4,type,reg,1},2);
%               BASEexpr4=mean(SINGLEUNIT.baseuncorr{4,type,reg,2},2);
%               
%               BASEcont5=mean(SINGLEUNIT.baseuncorr{5,type,reg,1},2);
%               BASEexpr5=mean(SINGLEUNIT.baseuncorr{5,type,reg,2},2);
%               
%               
%               figure
%               hold on
%                 bar(1,BASEcont1,'k')
%                     errorbar(1,BASEcont1, SINGLEUNIT.baseuncorrSEM(1,type,reg,1),'k');
%                 bar(2,BASEexpr1,'b')
%                     errorbar(2,BASEexpr1 ,SINGLEUNIT.baseuncorrSEM(1,type,reg,2),'k');
%                 
%                 bar(3,BASEcont2,'k')
%                     errorbar(3,BASEcont2 ,SINGLEUNIT.baseuncorrSEM(2,type,reg,1),'k');
%                 bar(4,BASEexpr2,'b')
%                      errorbar(4,BASEexpr2 ,SINGLEUNIT.baseuncorrSEM(2,type,reg,2),'k');
%                 
%                 bar(5,BASEcont3,'k')
%                     errorbar(5,BASEcont3 ,SINGLEUNIT.baseuncorrSEM(3,type,reg,1),'k');
%                 bar(6,BASEexpr3,'b')
%                     errorbar(6, BASEexpr3,SINGLEUNIT.baseuncorrSEM(3,type,reg,2),'k');
%                 
%                 bar(7,BASEcont4,'k')
%                     errorbar(7, BASEcont4,SINGLEUNIT.baseuncorrSEM(4,type,reg,1),'k');
%                 bar(8,BASEexpr4,'b')
%                     errorbar(8,BASEexpr4 ,SINGLEUNIT.baseuncorrSEM(4,type,reg,2),'k');
%                 
%                 bar(9,BASEcont5,'k')
%                     errorbar(9,BASEcont5 ,SINGLEUNIT.baseuncorrSEM(5,type,reg,1),'k');
%                 bar(10,BASEexpr5,'b')
%                     errorbar(10,BASEexpr5 ,SINGLEUNIT.baseuncorrSEM(5,type,reg,2),'k');
% 
%         
% 
%                 title(['Type:' num2str(type)  '  Region:' num2str(reg) ])

%incorrects together
            BASEcorrcont=mean(SINGLEUNIT.baseuncorrtcticOUTREMOVED{ses,1,reg,1},2);
                errorcorrcont=mean( SINGLEUNIT.baseuncorrSEMtcticOUTREMOVED(ses,1,reg,1),2); 
            BASEincont=mean(SINGLEUNIT.baseuncorrtcticOUTREMOVED{ses,2,reg,1} ,2);
                errorincont=mean( SINGLEUNIT.baseuncorrSEMtcticOUTREMOVED(ses,2,reg,1),2); 
                
            BASEcorrexp=mean(SINGLEUNIT.baseuncorrtcticOUTREMOVED{ses,1,reg,2},2);
                errorcorrexp=mean( SINGLEUNIT.baseuncorrSEMtcticOUTREMOVED(ses,1,reg,2) ,2); 
            BASEinexp=mean(SINGLEUNIT.baseuncorrtcticOUTREMOVED{ses,2,reg,2},2);
                errorinexp=mean(SINGLEUNIT.baseuncorrSEMtcticOUTREMOVED(ses,2,reg,2) ,2); 
            
            
            figure
            hold on
            bar(1,BASEcorrcont,'k')
                errorbar(1,BASEcorrcont,errorcorrcont,'k')
            bar(2,BASEcorrexp,'b')
                errorbar(2,BASEcorrexp,errorincont,'k')
            bar(3,BASEincont,'k')
                errorbar(3,BASEincont,errorcorrexp,'k')
            bar(4,BASEinexp,'b')
                errorbar(4,BASEinexp,errorinexp, 'k')
                
            set(gca,'ylim',[0 .45]);
            
            title(['Session:' num2str(ses)  'Region:' num2str(reg) ])


     
            
        end
    %end
end

      


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


for ses=1:5 %in order: D85, RD1, RD3, RD4, R50
    for reg=1:2 %DS is region 1 and OFC is region 2
   % for trt=1:2 %trt 1 is always control

    type1=mean(cell2mat(SINGLEUNIT.rawZtcticsmooth(ses,1,reg,1)),2);
        type1sem=cell2mat(SINGLEUNIT.rawZtcticsmoothSEM(ses,1,reg,1));
    type2=mean(cell2mat(SINGLEUNIT.rawZtcticsmooth(ses,2,reg,1)),2);
        type2sem=cell2mat(SINGLEUNIT.rawZtcticsmoothSEM(ses,2,reg,1));
%     type3=mean(cell2mat(SINGLEUNIT.rawcombase(ses,3,reg,1)),2);
%         type3sem=cell2mat(SINGLEUNIT.rawcombaseSEM(ses,3,reg,1));
%     type4=mean(cell2mat(SINGLEUNIT.rawcombase(ses,4,reg,1)),2);
%         type4sem=cell2mat(SINGLEUNIT.rawcombaseSEM(ses,4,reg,1));

    type1t2=mean(cell2mat(SINGLEUNIT.rawZtcticsmooth(ses,1,reg,2)),2);
        type1t2sem=cell2mat(SINGLEUNIT.rawZtcticsmoothSEM(ses,1,reg,2));
    type2t2=mean(cell2mat(SINGLEUNIT.rawZtcticsmooth(ses,2,reg,2)),2);
        type2t2sem=cell2mat(SINGLEUNIT.rawZtcticsmoothSEM(ses,2,reg,2));
%     type3t2=mean(cell2mat(SINGLEUNIT.rawcombase(ses,3,reg,2)),2);
%         type3t2sem=cell2mat(SINGLEUNIT.rawcombaseSEM(ses,3,reg,2));
%     type4t2=mean(cell2mat(SINGLEUNIT.rawcombase(ses,4,reg,2)),2);
%         type4t2sem=cell2mat(SINGLEUNIT.rawcombaseSEM(ses,4,reg,2));
   
        x=[1:80]';
        
    figure %for corrects
    hold on
       if reg==1 
           set(gca,'ylim',[-1 5]);
       elseif reg==2
           set(gca,'ylim',[-1 2]);
       end
        boundedline(x,type1,type1sem,'k','alpha'); 
        %boundedline(x,type4,type4sem,'r','alpha');
        
        boundedline(x,type1t2,type1t2sem,'b','alpha'); 
        %boundedline(x,type4t2,type4t2sem,'m','alpha'); 
        
        title(['Session:' num2str(ses) 'Region:' num2str(reg) ])
      
    figure %for incorrects
    hold on   
       if reg==1 
           set(gca,'ylim',[-1 4]);
       elseif reg==2
           set(gca,'ylim',[-1 4]);
       end
        boundedline(x,type2,type2sem,'k','alpha'); 
        %boundedline(x,type3,type3sem,'y','alpha'); 
        
        boundedline(x,type2t2,type2t2sem,'g','alpha'); 
        %boundedline(x,type3t2,type3t2sem,'r','alpha'); 

        title(['Session:' num2str(ses) 'Region:' num2str(reg) ])
        
    %end
    end
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


%% Counting how neurons fires to different types of behaviors

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




%% PCA -> do on ALL variable


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

%% Other stuff



for ses=1:5
    for type =1:4
        for reg=1:2
            for trt=1:2
                
            %find max post-firing rate
           avgmax(ses,type,reg,trt) =  mean( max(SINGLEUNIT.comspike{ses,type,reg,trt},[],1)  ,2);
           avgmin(ses,type,reg,trt) = mean( min(SINGLEUNIT.comspike{ses,type,reg,trt},[],1)  ,2);
           

            end
        end
    end
end


for ses=1:5
    for type=1:4
        for trt=1:2
            Num=size(SINGLEUNIT.trialspikes,1);
 
 for j=1:Num
     if size(SINGLEUNIT.trialspikes{j,type,ses,reg,trt},1) == 0
        continue
     end
     
     numTC(j,ses,type,trt)=size(SINGLEUNIT.trialspikes{j,type,ses,1,trt},1); %redundant to do both DS and OFC
     
     avgmaxDS(j,ses,type,trt)=max(mean( mean( SINGLEUNIT.trialspikes{j,type,ses,1,trt}, 1), 3));
     avgmaxOFC(j,ses,type,trt)=max(mean( mean( SINGLEUNIT.trialspikes{j,type,ses,2,trt}, 1), 3));
     
     
 end
        end
    end
end

for ses=1:5
        maxDScont=avgmaxDS(:,ses,2,1);
        maxOFCcont=avgmaxOFC(:,ses,2,1);
        
        maxDScre=avgmaxDS(:,ses,2,2);
        maxOFCcre=avgmaxOFC(:,ses,2,2);
        
        numTCcont=numTC(:,ses,2,1);
        numTCcre=numTC(:,ses,2,2);
        
        
        figure
        hold on
        scatter3(maxOFCcont,maxDScont,numTCcont,'filled','k');
        scatter3(maxOFCcre,maxDScre,numTCcre,'filled','b');
        view(40,35)

end




%%  CORROLATE INDIVIDUAL MOUSE BEHAIVOR TO PHASIC RECURITMENT 
base=[1:20];
post=[21:80];
%find phasics per mouse
for type=1:2
    for ses=1:5
       for reg=1:2
           for trt=1:2
               
           if type==1,
              Temp= SINGLEUNIT.TC;
              Num=size(Temp,1);
           elseif type==2, 
               Temp=SINGLEUNIT.TIC;
                Num=size(Temp,1);
           end
           
           
            for mou=1:Num
                DATA=Temp{mou,1,ses,reg,trt};
                if size(DATA,3)==1
                    DATA=DATA';
                elseif size(DATA,3)>1
                    DATA=squeeze(Temp{mou,1,ses,reg,trt});
                end

                if size(DATA,1)==0
                   continue
                end
              
                totnum=size(DATA,2);
                
                sig=ttest2(DATA(base,:), DATA(post,:)); %0=nonphasic 1=phasic
                z=size(find(sig==1),2); %find number that are significant
                
                if z==0
                    propphasic(mou,type,ses,reg,trt)=0;
                elseif z>0
                    propphasic(mou,type,ses,reg,trt)=z/totnum*100;
                end  
            end
             

          end
       end
    end
  clear Temp_DATA phasic increasing decreasing z sig phasic_base phasic_post 
end



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


%Graph behavior by porportion phasic
    %TC by phasic TC  and phasic TIC
    %TIC by phasic TIC  and phasic TC
    %Tot errors by phasic TIC  and phasic TC
    
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
    
    Y=propphasic(:,1,ses,reg,1); % control
    Y2=propphasic(:,1,ses,reg,2); % experimental
    
    Y3=propphasic(:,2,ses,reg,1); % control
    Y4=propphasic(:,2,ses,reg,2); % experimental

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
    scatter(X3, Y3,70,'k','LineWidth',2); % TIC with phasic response, control
        plot(fittedX3, fittedY3, 'k--', 'LineWidth', 2);
        %text(fittedX3(1,200),fittedY3(1,200),slopeY3); %slope of the line
        %text(fittedX3(1,200),fittedY3(1,200)-3,R3); %pearson's correlation coefficient
        %text(fittedX3(1,200),fittedY3(1,200)-7,P3); %P-VALUE OF correlation coefficient
    scatter(X4, Y4,70,'g','LineWidth',2); % TIC with phasic responce, experimental
        plot(fittedX4, fittedY4, 'g--', 'LineWidth', 2);
        %text(fittedX4(1,200),fittedY4(1,200),slopeY4);
        %text(fittedX4(1,200),fittedY4(1,200)-3,R4); %pearson's correlation coefficient
        %text(fittedX4(1,200),fittedY4(1,200)-7,P4); %P-VALUE OF correlation coefficient
        if ses==1
            set(gca, 'ylim', [0 100],'xlim',[0 30]);
        elseif ses== 5
            set(gca,'ylim', [0 100],'xlim',[0 80]);
        elseif ses== 2|| 3|| 4
            set(gca,'ylim', [0 100],'xlim',[0 140]); %set(gca,'xlim',[0 30],'ylim',[0 80]);
        end

%     
%     scatter(performancecre, phasicDScre,'filled','c');%TC or TIC with phasic DS response, experimental
%     scatter(performancecre, phasicOFCcre,'filled','m');%TC or TIC with phasic DS response, experimental
    %lsline
    
    title(['Session:' num2str(ses) ' Region' num2str(reg)]);
    
    
       % end
    end
end












