%% SET DIRECTORIES

cd('E:\Documents\MATLAB')
%addpath(genpath('G:\CARC_test'))
addpath(genpath('E:\Documents\MATLAB\kakearney-boundedline-pkg-2112a2b'))
%addpath(genpath('/Users/Kristin/Documents/MATLAB/CARC'))

load SFC_ANA_Instmou

save('SFC_ANA_DRBR_ANA','ANA','-v7.3')

load AnalyzedSFC_PAE_all_V2


save('AnalyzedSFC_PAE_all_V2','AnalyzedSFC','-v7.3')



%% set some stuff ---these must replicate what did in the actuall analysis

srate=1000; %start here and then downstample to 500Hz
numcycles=4; %lower = better temporal resolution; higher = better frequency resolution
num_wavelets=40;
start_freq=4; %epoch length should cover at least 3 cycle of the smallest frequency (ie: 1sec epoch=3Hz lowest frequency)
end_freq=40;
frex=linspace(start_freq,end_freq, num_wavelets);%linear increment from 1 to 100 for number of wavelets
    
%set some time variables -> for 50ms trials
    tx=(1:1:501)';
    t1=find(tx==1);
    t2=find(tx==501);
    tf_tx=-250:1:250;
    
    
    
%% CALCULATE DIFFERENCES FROM RANDOM CHANCE


% For difference values within condition and trt --> allows for comparison
% (between bins)
% NOT STATISTIC 
for PPCtest=1:3
for ses=1:5 
    for bin=1:4
        for reg=1:2
            for trt=1:2
                for beh=1:2
                       %calculate differences separately for test, session,
                       %behavior, region treatment and time bin
                       NumMice=size(ANA,1);
                     initilize_mouse=zeros(length(frex),1);
                     for mou=1:NumMice
                         if PPCtest==1
                            tempANA=ANA{mou,ses}.PPC0{bin,beh,reg,trt};
                         elseif PPCtest==2
                            tempANA=ANA{mou,ses}.PPC1{bin,beh,reg,trt};
                         elseif PPCtest==3
                            tempANA=ANA{mou,ses}.PPC2{bin,beh,reg,trt};
                         end
                         
                         initilize_mouse=[initilize_mouse tempANA]; %adds col for each mouse
                         

                     end
                     initilize_mouse(:,~any(initilize_mouse,1))=[];
                     
                     tempXmou{bin,reg,trt,beh,PPCtest,ses}=initilize_mouse; clear initilize_mouse
                     
                     
                 end
             end
        end

    end
    
end
    
end

%calcualte difference from random mix within session and type -> mixes by frequency and time                   
for PPCtest=1:3
for ses=1:5
        for reg=1:2
            for trt=1:2
                for beh=1:2

                %for bin=1:4
                   e1=mean(tempXmou{1,reg,trt,beh,PPCtest,ses},2);
                   e2=mean(tempXmou{2,reg,trt,beh,PPCtest,ses},2);
                   e3=mean(tempXmou{3,reg,trt,beh,PPCtest,ses},2);
                   e4=mean(tempXmou{4,reg,trt,beh,PPCtest,ses},2);

                  test=[e1 e2 e3 e4];
                  %test=[e1 e2]; %fq x bin
                 

               if size(test,2)==0 %if come across placeholder
                   continue
               end
               
               N=1;
               while N < 5001 %do the random mixture 5000 times to get a true average
                    msize=size(test); %Num of frex
                    idx_fq=randperm(msize(1,1)); %determine random order   
                   for z=1:msize(1,1) 
                        row=idx_fq(z);
                        col=randperm(msize(1,2),1); %choose a random bin to take from 
                        ran(z,:,N)=test(row,col);%mixes up data in random order random fq and time
                   end           
                 clear msize idx  row 
                 N=N+1;
               end
               
               randomDerZero(:,ses,reg,trt,beh,PPCtest)=mean(mean(mean(ran,3),1),2); %random derived zero at each fq
               RAND_DIFF_AVG{ses,beh,reg,trt,PPCtest}=bsxfun(@minus, test, randomDerZero); clear test  %difference from random of X-mouse (total avg)
               
               %once random derived zero calcualted - need to find difference w/in
               %each mouse
                NumMice=size(ANA,1);
               for mou=1:NumMice
                   
                   if PPCtest==1
                       tempANA_bin1=ANA{mou,ses}.PPC0{1,beh,reg,trt};
                       tempANA_bin2=ANA{mou,ses}.PPC0{2,beh,reg,trt};
                       tempANA_bin3=ANA{mou,ses}.PPC0{3,beh,reg,trt};
                       tempANA_bin4=ANA{mou,ses}.PPC0{4,beh,reg,trt};
                       
                      % tempANA_mou=[tempANA_bin1 tempANA_bin2];
                       tempANA_mou=[tempANA_bin1 tempANA_bin2 tempANA_bin3 tempANA_bin4];
                      
                   elseif PPCtest==2
                       tempANA_bin1=ANA{mou,ses}.PPC1{1,beh,reg,trt};
                       tempANA_bin2=ANA{mou,ses}.PPC1{2,beh,reg,trt};
                       tempANA_bin3=ANA{mou,ses}.PPC1{3,beh,reg,trt};
                       tempANA_bin4=ANA{mou,ses}.PPC1{4,beh,reg,trt};
                       
                       %tempANA_mou=[tempANA_bin1 tempANA_bin2];
                       tempANA_mou=[tempANA_bin1 tempANA_bin2 tempANA_bin3 tempANA_bin4];
                      
                   elseif PPCtest==3
                       tempANA_bin1=ANA{mou,ses}.PPC2{1,beh,reg,trt};
                       tempANA_bin2=ANA{mou,ses}.PPC2{2,beh,reg,trt};
                       tempANA_bin3=ANA{mou,ses}.PPC2{3,beh,reg,trt};
                       tempANA_bin4=ANA{mou,ses}.PPC2{4,beh,reg,trt};
                       
                       %tempANA_mou=[tempANA_bin1 tempANA_bin2];
                       tempANA_mou=[tempANA_bin1 tempANA_bin2 tempANA_bin3 tempANA_bin4];
                      
                   end

                   RAND_DIFF_MOU{mou,ses,beh,reg,trt,PPCtest}=bsxfun(@minus, tempANA_mou , randomDerZero); clear test ran %difference from random
               
                   MOU_NODIFF{mou,ses,beh,reg,trt,PPCtest}=tempANA_mou; clear tempANA* %each mouse w/ bins lined up, but no difference taken
                  
            
               end
               
                 end
             end   
        end   
end
end

AnalyzedSFC.random=randomData;
AnalyzedSFC.RAND_DIFF_AVG=RAND_DIFF_AVG;
AnalyzedSFC.RAND_DIFF_MOU=RAND_DIFF_MOU;
AnalyzedSFC.MOU_NODIFF=MOU_NODIFF;

%put all mice into same matrix
N_bins=size(AnalyzedSFC.RAND_DIFF_MOU{mou,ses,beh,reg,trt,PPCtest},2);
for PPCtest=1:3
for ses=1:5
    for reg=1:2
        for trt=1:2
            for beh=1:2
                    
NumMice=size(ANA,1);
initilize_mouse=zeros(length(frex),N_bins,1);
for mou=1:NumMice

    tempMou=AnalyzedSFC.MOU_NODIFF{mou,ses,beh,reg,trt,PPCtest};
    if size(tempMou,2) < size(initilize_mouse,2)
        continue
    end
    
    initilize_mouse=cat(3,initilize_mouse, tempMou);

end

    initilize_mouse(:,:,1)=[];
    
    
    RAND_NODIFF_AllMOU{ses,beh,reg,trt,PPCtest}=initilize_mouse; clear initilize_mouse

            end
        end 
    end
end
end

AnalyzedSFC.RAND_NODIFF_AllMOU=RAND_NODIFF_AllMOU;

%% GRAPH

ses=1;
beh=1;
reg=1;
trt=2;
PPCtest=2;

freq_theta=2:5; %4.9231 - 7.6923
freq_beta=8:16; %10.4615 - 17.8462

 graph_temp=RAND_NODIFF_AllMOU{ses,beh,reg,trt,PPCtest};
 
 SEM_test= std(graph_temp,[],3) ./ sqrt(size(graph_temp,3));
    graph_SEM=[SEM_test(:,1)' SEM_test(:,2)' SEM_test(:,3)' SEM_test(:,4)'];
 MEAN_test=mean(graph_temp,3);
    graph=[MEAN_test(:,1)' MEAN_test(:,2)' MEAN_test(:,3)' MEAN_test(:,4)'];
 
 figure
 hold on
% plot(zeros(1,40));
% plot(repmat(frex(3:40),1,4),graph(3:40,:))
 errorbar(MEAN_test(3:40,4)', SEM_test(3:40,4)','-sk','MarkerFaceColor','k','markers',12)

 
  
 SEM_test= std(graph_temp,[],3) ./ sqrt(size(graph_temp,3));
    graph_SEM=[mean(SEM_test(freq_theta,1),1) mean(SEM_test(freq_theta,2),1) mean(SEM_test(freq_theta,3),1) mean(SEM_test(freq_theta,4),1)];
 MEAN_test=mean(graph_temp,3);
    graph=[mean( MEAN_test(freq_theta,1),1) mean( MEAN_test(freq_theta,2),1) mean( MEAN_test(freq_theta,3),1) mean( MEAN_test(freq_theta,4),1)];
 
 figure
 hold on
plot(zeros(1,5));
% plot(repmat(frex(3:40),1,4),graph(3:40,:)
 errorbar(graph, graph_SEM,'-sk','MarkerFaceColor','k','markers',12)
 scatter(1:4, mean(graph_temp(freq_theta,:,4),1))
 
 % Line graph - avg fq compared across sessions w/ points per mouse
bin=3;
for PPCtest=2:2
    for reg=2:2
        for beh=1:1
            for frequency=1:1
                if frequency==1,fqr=1:13;%low
                elseif frequency==2,fqr=14:24;% med
                elseif frequency==3,fqr=25:40;%high
                end
 
    for ses=1:5
            
        tempdata=squeeze(mean(AnalyzedSFC.RAND_NODIFF_AllMOU{ses,beh,reg,1,PPCtest}(fqr,:,:),1)); %control mice - mean in fqr dimension
          points_cont{ses}=tempdata(bin,:)'; %bin x mouse x ses
          MEAN_cont(ses,:)=mean(tempdata,2)'; %session x bin
          SEM_cont(ses,:)=(std(tempdata,0,2)/ sqrt(size(tempdata,2)))'; clear tempdata
              
        tempdata=squeeze(mean(AnalyzedSFC.RAND_NODIFF_AllMOU{ses,beh,reg,2,PPCtest}(fqr,:,:),1)); %experimental mice - mean in fqr dimension
          points_exp{ses}=tempdata(bin,:)'; %bin x mouse x ses   
          MEAN_exp(ses,:)=mean(tempdata,2)'; %session x bin
          SEM_exp(ses,:)=(std(tempdata,0,2)/ sqrt(size(tempdata,2)))'; clear tempdata

    end       

    for bin=3:3
    figure %plot the individual points
    hold on
    %control...
    errorbar(MEAN_cont(:,bin), SEM_cont(:,bin),'-ok','MarkerFaceColor','k','markers',12)
    %scatter(ones(length(points_cont{1}),1),points_cont{1},90,'^k','LineWidth',1); %^ = triangle
    %scatter(repmat(2,length(points_cont{2}),1),points_cont{2},90,'^k','LineWidth',1);
    %scatter(repmat(3,length(points_cont{3}),1),points_cont{3},90,'^k','LineWidth',1);
    %scatter(repmat(4,length(points_cont{4}),1),points_cont{4},90,'^k','LineWidth',1);
    %scatter(repmat(5,length(points_cont{5}),1),points_cont{5},90,'^k','LineWidth',1);
    
    %treatment/ experimental...
    errorbar(MEAN_exp(:,bin), SEM_exp(:,bin),'-ok','MarkerFaceColor','r','markers',12)
    %scatter(ones(length(points_exp{1}),1),points_exp{1},90,'^r','LineWidth',1); %^ = triangle
    %scatter(repmat(2,length(points_exp{2}),1),points_exp{2},90,'^r','LineWidth',1);
    %scatter(repmat(3,length(points_exp{3}),1),points_exp{3},90,'^r','LineWidth',1);
    %scatter(repmat(4,length(points_exp{4}),1),points_exp{4},90,'^r','LineWidth',1);
    %scatter(repmat(5,length(points_exp{5}),1),points_exp{5},90,'^r','LineWidth',1);
    ylim([.005 .04])
    
    % if fqr==1 || fqr==2
       % ylim([0 .055])
%     %elseif fqr==2
%        % ylim([0 8])
%     elseif fqr==3
%         ylim([0 5])
 % end
 
    end
            end
    clear points_cont MEAN_cont SEM_cont SEM_exp MEAN_exp points_exp
        end
    end
    
end


%AnalyzedSFC.random(fq,ses,reg,trt,beh,PPCtest)

randomtemp=squeeze(mean(mean(AnalyzedSFC.random,1),4));
%AnalyzedSFC.random(ses,reg,beh,PPCtest)

% BOUNDED LINE GRAPHS - 
for PPCtest=2:2
for ses=1:5
    for beh=1:1
        for reg=1:1
            
        tempdata_cont=AnalyzedSFC.RAND_NODIFF_AllMOU{ses,beh,reg,1,PPCtest}; %bin x mouse x ses
        tempmean_cont= mean(tempdata_cont,3);
        tempsem_cont= (std(tempdata_cont,0,3)/ sqrt(200));
        
        tempdata_exp=AnalyzedSFC.RAND_NODIFF_AllMOU{ses,beh,reg,2,PPCtest};
        tempmean_exp= mean(tempdata_exp,3);
        tempsem_exp= (std(tempdata_exp,0,3)/ sqrt(200)); %number of spikes per permuation iteration
        
        temprandom=randomtemp(ses,reg,beh,PPCtest);
        %figure
       % plot(frex,tempmean(:,1))

    figure
    hold on
    boundedline(frex(1,2:40),tempmean_cont(2:40,3),tempsem_cont(2:40,3),'k','alpha');
    boundedline(frex(1,2:40),tempmean_exp(2:40,3),tempsem_exp(2:40,3),'r','alpha');
    plot(frex(1,2:40),repmat(temprandom,1,size(frex(1,2:40),2)),':k','LineWidth',1.5)
   % boundedline(times2save,GRAPHDATA_experimental,MOUSESEM_experimental,'r','alpha');
    ylim([0.008 .035])
     
    clear temp*
        end
    end
end
end

%% exporting every frequency

%final set up:   ses  beh  regon    M#     fq1      fq2....
alldata=zeros(1,44);
for ses=1:5
    for reg=1:2
        for trt=1:2
    
    NumMice=size(AnalyzedSFC.MOU_NODIFF,1);
    combined_data=zeros(NumMice,size(frex,2)+1);
    for mou=1:NumMice
     %ROI(1,mou,ses,reg,trt,beh)=
     if size(AnalyzedSFC.MOU_NODIFF{mou,ses,beh,reg,trt,PPCtest},2)== 0
        temp=zeros(1,size(frex,2)); 
        combined_data(mou,:)=cat(2,str2double(var5(ses,mou,trt)),temp); clear temp
     elseif size(AnalyzedSFC.MOU_NODIFF{mou,ses,beh,reg,trt,PPCtest},2)> 0
        temp=AnalyzedSFC.MOU_NODIFF{mou,ses,beh,reg,trt,PPCtest}(:,bin,:)'; 
        combined_data(mou,:)=cat(2,str2double(var5(ses,mou,trt)),temp); clear temp
     end
    end
    
    tempses=repmat(ses,NumMice,1);
    tempreg=repmat(reg,NumMice,1);
    temptrt=repmat(trt,NumMice,1);
    
    layout=[tempses,tempreg,temptrt,combined_data]; clear temp*
    alldata=[alldata; layout]; clear layout
    
        end
        
    end
end

alldata(1,:)=[]; %delete first row of zeros

allmice=size(alldata,1);
for j=1:allmice
    
    if sum(alldata(j,5:40))==0
        alldata(j,:)=[];
        allmice=allmice-1;
    end
    
end


 %% save data in table - for ANOVA analysis in 'R'
ROI=zeros(3,14,5,2,2,1);

for PPCtest=2:2
for beh=1:1
    for ses=1:5
        for trt=1:2
            for reg=1:2
                
        NumMice=size(AnalyzedSFC.MOU_NODIFF,1);
        for mou=1:NumMice
            
            bin=3;
            
            frequency1=1:13;%low
            frequency2=14:24;% med
            frequency3=25:40;%high
            
            if size(AnalyzedSFC.MOU_NODIFF{mou,ses,beh,reg,trt,PPCtest},2)>0
                ROI(1,mou,ses,reg,trt,beh)=squeeze(mean(AnalyzedSFC.MOU_NODIFF{mou,ses,beh,reg,trt,PPCtest}(frequency1,bin),1)); %taking average at frequency
                ROI(2,mou,ses,reg,trt,beh)=squeeze(mean(AnalyzedSFC.MOU_NODIFF{mou,ses,beh,reg,trt,PPCtest}(frequency2,bin),1));
                ROI(3,mou,ses,reg,trt,beh)=squeeze(mean(AnalyzedSFC.MOU_NODIFF{mou,ses,beh,reg,trt,PPCtest}(frequency3,bin),1)); 
            elseif size(AnalyzedSFC.MOU_NODIFF{mou,ses,beh,reg,trt,PPCtest},2)==0
                ROI(1,mou,ses,reg,trt,beh)=0; %taking average at frequency
                ROI(2,mou,ses,reg,trt,beh)=0;
                ROI(3,mou,ses,reg,trt,beh)=0; 
            end

        end   
            end
        end
    end     
end
end



        D85=[ROI(1,:,1,1,1)'; ROI(1,:,1,1,2)';ROI(1,:,1,2,1)';ROI(1,:,1,2,2)';...
             ROI(2,:,1,1,1)'; ROI(2,:,1,1,2)';ROI(2,:,1,2,1)';ROI(2,:,1,2,2)';...
             ROI(3,:,1,1,1)'; ROI(3,:,1,1,2)';ROI(3,:,1,2,1)';ROI(3,:,1,2,2)']; 
    
        RS1=[ROI(1,:,2,1,1)'; ROI(1,:,2,1,2)';ROI(1,:,2,2,1)';ROI(1,:,2,2,2)';...
             ROI(2,:,2,1,1)'; ROI(2,:,2,1,2)';ROI(2,:,2,2,1)';ROI(2,:,2,2,2)';...
             ROI(3,:,2,1,1)'; ROI(3,:,2,1,2)';ROI(3,:,2,2,1)';ROI(3,:,2,2,2)']; 
    
        RS3=[ROI(1,:,3,1,1)'; ROI(1,:,3,1,2)';ROI(1,:,3,2,1)';ROI(1,:,3,2,2)';...
             ROI(2,:,3,1,1)'; ROI(2,:,3,1,2)';ROI(2,:,3,2,1)';ROI(2,:,3,2,2)';...
             ROI(3,:,3,1,1)'; ROI(3,:,3,1,2)';ROI(3,:,3,2,1)';ROI(3,:,3,2,2)']; 
    
        RS4=[ROI(1,:,4,1,1)'; ROI(1,:,4,1,2)';ROI(1,:,4,2,1)';ROI(1,:,4,2,2)';...
             ROI(2,:,4,1,1)'; ROI(2,:,4,1,2)';ROI(2,:,4,2,1)';ROI(2,:,4,2,2)';...
             ROI(3,:,4,1,1)'; ROI(3,:,4,1,2)';ROI(3,:,4,2,1)';ROI(3,:,4,2,2)']; 
    
        R50=[ROI(1,:,5,1,1)'; ROI(1,:,5,1,2)';ROI(1,:,5,2,1)';ROI(1,:,5,2,2)';...
             ROI(2,:,5,1,1)'; ROI(2,:,5,1,2)';ROI(2,:,5,2,1)';ROI(2,:,5,2,2)';...
             ROI(3,:,5,1,1)'; ROI(3,:,5,1,2)';ROI(3,:,5,2,1)';ROI(3,:,5,2,2)']; 

              %low ds sac     low ds pae      low ofc sac     low ofc pae
 
 
 
 
Frequency={'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'ph';'ph';...
            'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';...
            'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'ph';'ph';...
            'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';...
     
            'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'ph';'ph';... 
            'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';...
            'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'ph';'ph';... 
            'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';...
          
          'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'ph';'ph';... 
          'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';...
          'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'ph';'ph';... 
          'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH'};

Region=   {'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS';'ph';'ph';...   
           'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';...
           'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'ph';'ph';...
           'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC'
           
           'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS';'ph';'ph';...   
           'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';...
           'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'ph';'ph';...
           'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC'
           
           'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS';'ph';'ph';...   
           'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';...
           'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'ph';'ph';...
           'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC'}; 
       
Treatment={'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'ph';'ph';...
           'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';...
           'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'ph';'ph';...
           'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';
           
           'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'ph';'ph';...
           'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';...
           'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'ph';'ph';...
           'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';

           'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'ph';'ph';...
           'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';...
           'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'ph';'ph';...
           'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE'};
       
MouNum={'2.0';'3.1';'11.0';'16.1';'17.0';'17.1';'17.2';'20.1';'22.2';'25.0';'26.0';'27.0';'ph';'ph';...
        '4.1';'14.1';'18.0';'19.0';'19.1';'21.9';'21.1';'22.1';'22.2';'24.0';'28.0';'29.0';'30.0';'31.0';...
        '2.0';'3.1';'11.0';'16.1';'17.0';'17.1';'17.2';'20.1';'22.2';'25.0';'26.0';'27.0';'ph';'ph';...
        '4.1';'14.1';'18.0';'19.0';'19.1';'21.9';'21.1';'22.1';'22.2';'24.0';'28.0';'29.0';'30.0';'31.0'
        
        '2.0';'3.1';'11.0';'16.1';'17.0';'17.1';'17.2';'20.1';'22.2';'25.0';'26.0';'27.0';'ph';'ph';...
        '4.1';'14.1';'18.0';'19.0';'19.1';'21.9';'21.1';'22.1';'22.2';'24.0';'28.0';'29.0';'30.0';'31.0';...
        '2.0';'3.1';'11.0';'16.1';'17.0';'17.1';'17.2';'20.1';'22.2';'25.0';'26.0';'27.0';'ph';'ph';...
        '4.1';'14.1';'18.0';'19.0';'19.1';'21.9';'21.1';'22.1';'22.2';'24.0';'28.0';'29.0';'30.0';'31.0'
    
        '2.0';'3.1';'11.0';'16.1';'17.0';'17.1';'17.2';'20.1';'22.2';'25.0';'26.0';'27.0';'ph';'ph';...
        '4.1';'14.1';'18.0';'19.0';'19.1';'21.9';'21.1';'22.1';'22.2';'24.0';'28.0';'29.0';'30.0';'31.0';...
        '2.0';'3.1';'11.0';'16.1';'17.0';'17.1';'17.2';'20.1';'22.2';'25.0';'26.0';'27.0';'ph';'ph';...
        '4.1';'14.1';'18.0';'19.0';'19.1';'21.9';'21.1';'22.1';'22.2';'24.0';'28.0';'29.0';'30.0';'31.0'};

    T=table(Frequency,Region, MouNum,Treatment,D85,RS1,RS3,RS4,R50); AnalyzedSFC.table_FQ_PPC1=T;
    writetable(T,'SFC_ALL_PPC1.xls','Sheet',1); 

    %%
     
Frequency={'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'ph';'ph';...
            'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';...
            'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'ph';'ph';...
            'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';'LOW';...
     
            'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'ph';'ph';... 
            'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';...
            'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'ph';'ph';... 
            'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';'MED';...
          
          'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'ph';'ph';... 
          'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';...
          'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'ph';'ph';... 
          'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH';'HIGH'};

Region=   {'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'ph';'ph';...   
           'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';...
           'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'ph';'ph';...
           'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';...
           
           'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'ph';'ph';...   
           'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';...
           'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'ph';'ph';...
           'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';...
           
           'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'DS'; 'ph';'ph';...   
           'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';'DS';...
           'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'ph';'ph';...
           'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC';'OFC'}; 
       
Treatment={'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'ph';'ph';...
           'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';...
           'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'ph';'ph';...
           'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';...
           
           'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'ph';'ph';...
           'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';...
           'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'ph';'ph';...
           'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';...

           'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'ph';'ph';...
           'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';...
           'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'ph';'ph';...
           'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE'};
       
MouNum={'2.0';'3.1';'11.0';'16.1';'17.0';'17.1';'17.2';'20.1';'22.2';'25.0';'26.0';'27.0';'ph';'ph';...
        '4.1';'14.1';'18.0';'19.0';'19.1';'21.9';'21.1';'22.1';'22.2';'24.0';'28.0';'29.0';'30.0';'31.0';...
        '2.0';'3.1';'11.0';'16.1';'17.0';'17.1';'17.2';'20.1';'22.2';'25.0';'26.0';'27.0';'ph';'ph';...
        '4.1';'14.1';'18.0';'19.0';'19.1';'21.9';'21.1';'22.1';'22.2';'24.0';'28.0';'29.0';'30.0';'31.0';...
        
        '2.0';'3.1';'11.0';'16.1';'17.0';'17.1';'17.2';'20.1';'22.2';'25.0';'26.0';'27.0';'ph';'ph';...
        '4.1';'14.1';'18.0';'19.0';'19.1';'21.9';'21.1';'22.1';'22.2';'24.0';'28.0';'29.0';'30.0';'31.0';...
        '2.0';'3.1';'11.0';'16.1';'17.0';'17.1';'17.2';'20.1';'22.2';'25.0';'26.0';'27.0';'ph';'ph';...
        '4.1';'14.1';'18.0';'19.0';'19.1';'21.9';'21.1';'22.1';'22.2';'24.0';'28.0';'29.0';'30.0';'31.0';...
    
        '2.0';'3.1';'11.0';'16.1';'17.0';'17.1';'17.2';'20.1';'22.2';'25.0';'26.0';'27.0';'ph';'ph';...
        '4.1';'14.1';'18.0';'19.0';'19.1';'21.9';'21.1';'22.1';'22.2';'24.0';'28.0';'29.0';'30.0';'31.0';...
        '2.0';'3.1';'11.0';'16.1';'17.0';'17.1';'17.2';'20.1';'22.2';'25.0';'26.0';'27.0';'ph';'ph';...
        '4.1';'14.1';'18.0';'19.0';'19.1';'21.9';'21.1';'22.1';'22.2';'24.0';'28.0';'29.0';'30.0';'31.0'};

 