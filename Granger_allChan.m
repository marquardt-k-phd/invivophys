% GRANGER PREDICTION (CASAULITY) - BIC & Granger sliding time window

addpath(genpath('E:\Documents\MATLAB\BSMART'))
addpath(genpath('E:\Documents\MATLAB\kakearney-boundedline-pkg-2112a2b'))


save('DRBR_allchanGRANG','GRANG','-v7.3');  

%load FP data (need events and avg FP)- only first time run script to make extended epochs
%load LFP_DUAL_PAE_All
load DRBR_allchanGRANG


%% CREATE EXTENDED EPOCHS -2 to +4 sec - to be able to cover granger at ends
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
        
        %find lfp times & create trial epochs 
          lfp=cell2mat(LFP.refLFP(j,ses,reg,trt));
          int=cell2mat(LFP.interval(j,ses,trt));

           for b=1:size(trialtypes,1) %if everything is perfect
               temp=find(int <= trials(b,1)+.001 & int >= trials(b,1)-.001); %find where the LFP is closest to the spike in time
                    if trials(b,1) > max(int) %some trials are beyond max lfp timepoint
                        templfp(b,:)=zeros(1,6001);
                    elseif trials(b,1)<= max(int) %for those that are not
                        startrow=temp-2000; %1sec prior - add buffer second
                        endrow=temp+4000; %3 sec post - add buffer second
                        templfp(:,:,b)=lfp(startrow:endrow,:);
                    end
               clear temp startrow endrow
           end
          
        %clear out zeros -> in both LFP and trial type arrays
            tempzero=find(templfp(:,all(templfp==0,1)));
            trialtypes(tempzero,:)=[];
            templfp(tempzero,:)=[];
        
        %split bins into trial types % store     
        LFPEPOCHS{j,ses,1,reg,trt}=templfp(:,:,logical(trialtypes(:,2)==1)); %lfp bins that are type 1 trials 
        LFPEPOCHS{j,ses,2,reg,trt}=templfp(:,:,logical(trialtypes(:,2)==2)); %lfp bins that are incorrect trials  

   clear tone tc punish tic trial* type* x y z w lfp int temp*
end
    clear Num
        end
    end
end
GRANG.epochs_extended=LFPEPOCHS; 
 
 
 clear LFP

 
%% DOWNSAMPLE
 srate=1000; %actual collection rate (Hz)
 dsrate=200; ratio=srate/dsrate; %the rate to downsample to (Hz)
 down=1:ratio:6001;
 
 
 %Downsample epochs b/c already sorted + z-score along individual trials
for ses=1:5
    for reg=1:2
        for trt=1:2
            for beh=1:2
                Num=size(GRANG.epochs_extended,1);
                for mou=1:Num
                   if size (GRANG.epochs_extended{mou,ses,beh,reg,trt},2)==0
                       continue
                   end
                
                    tempFP=GRANG.epochs_extended{mou,ses,beh,reg,trt};
                    
                    tempdown=tempFP(down,:,:);
                    
                    downsample{mou,ses,beh,reg,trt}=permute(tempdown,[3 1 2]); %downsampled and averaged between electrodes

                end
            end
        end
    end
end
GRANG.epoch_downsample=downsample;

%% BIC CALCULATION - this will take several hours all together

%set window parameters
timewin = 750; % in ms
timewin_points = round(timewin/(1000/dsrate));
LFPtime=-2000:ratio:4000;
times2save = -1000: timewin/2 :3000; %every 125 ms  - timewindow of whole epoch in downsampled time: each point that want to take a timepoint around
    %will take 100(timewin_points around here from data) data is
    %downsampled 200Hz so this is = 500ms
times2saveidx = dsearchn( LFPtime',times2save');% what those times are in downsampled data cols

ncrit=12; %critical number of trials *********
 
%BIC on sliding timescale -- run overnight
for ses=1:5
    for beh=1:2
        for trt=1:2
            Num=size(GRANG.epochs_extended,1);
            for mou=1:Num
               
               if size(GRANG.epoch_downsample{mou,ses,beh,1,trt},1)==0 || size(GRANG.epoch_downsample{mou,ses,beh,2,trt},1)==0
                    %if either region does not have any data cannot do
                    %granger
                    disp(['Mouse: ' num2str(mou) ' Treatment: ' num2str(trt) ' Behavior: ' num2str(beh) ' Session: ' num2str(ses) ' Does not have spikes in region 1 and/ or 2']);
                   continue
                end
                %skips mice that have less than critical number of trials
                if size(GRANG.epoch_downsample{mou,ses,beh,1,trt},1)< ncrit %used region==1 b/c both regions have same #trials
                    disp(['Mouse: ' num2str(mou) ' Treatment: ' num2str(trt) ' Behavior ' num2str(beh) ' Session: ' num2str(ses) ' Does not have > ' num2str(ncrit) ' trials']);
                    continue
                end
                
      for chan =1:7      

       eegdata_OFC=GRANG.epoch_downsample{mou,ses,beh,2,trt}(:,:,chan)'; %{mou,ses,beh,reg,trt};
       eegdata_DS=GRANG.epoch_downsample{mou,ses,beh,1,trt}(:,:,chan)';
       eegdata(1,:,:)=eegdata_OFC; %region x time x trials
       eegdata(2,:,:)=eegdata_DS;
       
       eegdata = bsxfun(@minus,eegdata(:,:,:),mean(eegdata(:,:,:),3)); %subtract ERP (mean EEG)
       
  bic = zeros(length(times2save),25); % Bayes info criteria (go up to MO of 30)
  tempdataZ=zeros(2,timewin_points,size(eegdata,3));
for timei=1:length(times2save) %cant do first or last
    
    % data from all trials in this time window  
        tempdata = squeeze(eegdata(:,   floor(times2saveidx(timei))-floor(timewin_points/2):floor(times2saveidx(timei))+floor(timewin_points/2)-mod(timewin_points+1,2)   ,:));
        
       % -mod(timewin_points+1,2)
    
    % detrend and zscore all data - to improve stationarity
    for triali=1:size(tempdata,3)
        tempdataZ(1,:,triali)=zscore(detrend(squeeze(tempdata(1,:,triali))));
        tempdataZ(2,:,triali)=zscore(detrend(squeeze(tempdata(2,:,triali))));
    end

    % reshape tempdata for armorf
    tempdataRS = reshape(tempdataZ,2,timewin_points*size(tempdata,3));
    
        %for timei=1:segleg+1; %do we need to go through each timepoint...
            % test BIC for optimal model order at each time point
            for bici=1:size(bic,2)
                % run model
                [~,E] = armorf(tempdataRS,size(tempdataZ,3),timewin_points,bici);
                % compute Bayes Information Criteria
                bic(timei,bici) = log(det(E)) + (log(length(tempdataRS))*bici*2^2)/length(tempdataRS);
            end
        %end
        clear temp*


end % end time loop

        [bestbicVal,bestbicIdx]=min(mean(bic,1));
        BIC(mou,chan,ses,beh,trt)=bestbicVal; %actual value
        BICidx(mou,chan,ses,beh,trt)=bestbicIdx; %model order --> to convert to ms bestbicIdx*(1000/dsrate)
        clear bic bici best*
       
      end  
        disp(['Mouse: ' num2str(mou) ' Treatment: ' num2str(trt) ' Behavior ' num2str(beh) ' Session: ' num2str(ses) ' BIC complete']);

    clear eegdata*
            end     
        end
    end
end


GRANG.BICidxChan=BICidx;

BICidx(BICidx==0)=[];
model_order=mean(BICidx);
GRANG.ModelOrderChan=model_order;

% test figure - for use with bic
figure
subplot(121)
plot((1:size(bic,2))*(1000/dsrate),mean(bic,1),'--.')
xlabel('Order (converted to ms)')
ylabel('Mean BIC over all time points')
[bestbicVal,bestbicIdx]=min(mean(bic,1));
hold on
plot(bestbicIdx*(1000/dsrate),bestbicVal,'mo','markersize',15)
title([ 'Optimal order is ' num2str(bestbicIdx) ' (' num2str(bestbicIdx*(1000/dsrate)) ' ms)' ])

subplot(122)
[junk,bic_per_timepoint] = min(bic,[],2);
plot(times2save,bic_per_timepoint*(1000/dsrate),'--.')
xlabel('Time (ms)')
ylabel('Optimal order (converted to ms)')
title('Optimal order (in ms) at each time point')





%% Granger on sliding timescale
model_order=round(GRANG.ModelOrderChan);    %round(GRANG.modelorder);
order_time = model_order*(1000/dsrate); % in ms

%order_points = round(order_time/(1000/dsrate)); %basicaly same as model order - except rounded

min_freq=1; %hz
max_freq=10;%hz

frequencies=logspace(log10(min_freq),log10(max_freq),model_order);
N_sessions=5;
N_behaviors=2;
N_treatments=2;

GC=cell(size(GRANG.epoch_downsample,1),N_sessions,N_behaviors,N_treatments);
for ses=1:N_sessions
    for beh=1:N_behaviors
        for trt=1:N_treatments
            Num=size(GRANG.epochs_extended,1);
            for mou=1:Num
               
               if size(GRANG.epoch_downsample{mou,ses,beh,1,trt},1)==0 || size(GRANG.epoch_downsample{mou,ses,beh,2,trt},1)==0
                    %if either region does not have any data cannot do
                    %granger
                    disp(['Mouse: ' num2str(mou) ' Treatment: ' num2str(trt) ' Behavior: ' num2str(beh) ' Session: ' num2str(ses) ' Does not have spikes in region 1 and/ or 2']);
                   continue
                end
                %skips mice that have less than critical number of trials
                if size(GRANG.epoch_downsample{mou,ses,beh,1,trt},1)< ncrit %used region==1 b/c both regions have same #trials
                    disp(['Mouse: ' num2str(mou) ' Treatment: ' num2str(trt) ' Behavior ' num2str(beh) ' Session: ' num2str(ses) ' Does not have > ' num2str(ncrit) ' trials']);
                    continue
                end
for chan=1:7
       eegdata_OFC=GRANG.epoch_downsample{mou,ses,beh,2,trt}(:,:,chan)'; %{mou,ses,beh,reg,trt};
       eegdata_DS=GRANG.epoch_downsample{mou,ses,beh,1,trt}(:,:,chan)';
       eegdata(1,:,:)=eegdata_OFC; %region x trials x time
       eegdata(2,:,:)=eegdata_DS;
       
       eegdata = bsxfun(@minus,eegdata(:,:,:),mean(eegdata(:,:,:),3)); %subtract ERP (mean EEG)
       
  bic = zeros(length(times2save),25); % Bayes info criteria (go up to MO of 30)
  tempdataZ=zeros(2,timewin_points,size(eegdata,3)); 
 for timei=1:length(times2saveidx) %cant do first or last
    
    % data from all trials in this time window  
        tempdata = squeeze(eegdata(:,   floor(times2saveidx(timei))-floor(timewin_points/2):floor(times2saveidx(timei))+floor(timewin_points/2)-mod(timewin_points+1,2)   ,:));
        
       % -mod(timewin_points+1,2)
    
    % detrend and zscore all data - to improve stationarity
    for triali=1:size(tempdata,3)
        tempdataZ(1,:,triali)=zscore(detrend(squeeze(tempdata(1,:,triali))));
        tempdataZ(2,:,triali)=zscore(detrend(squeeze(tempdata(2,:,triali))));
    end
    
     %test of stationary data - ADF test - unsure what regression
     %coefficient to use- must be <5
%     [tstat,cval]=mvgc_adf(tempdataZ, .05,lags, 5); mean_tstat=mean(tstat,1);
%     Mean_Tstat{mou,ses,beh,trt}(timei,:)=mean_tstat;
%     Cval_t{mou,ses,beh,trt}(timei,:)=cval;

    %reshape tempdata for armorf
    tempdataRS = reshape(tempdataZ,2,timewin_points*size(tempdata,3));
    
     %test of stationary data - KPSS unit root test
     lags=round(sqrt(size(tempdata,3)));
     [ksstat, cval]= mvgc_kpss(tempdataRS,.05,lags); 
     Ksstat(timei,:,mou,ses,beh,trt)=ksstat; clear ksstat 
     Cval_k(timei,mou,ses,beh,trt)=cval; clear cval
    
 
 
    % fit AR models
    [Ax,Ex] = armorf(tempdataRS(1,:),size(tempdataZ,3),timewin_points,model_order); %ofc
    [Ay,Ey] = armorf(tempdataRS(2,:),size(tempdataZ,3),timewin_points,model_order); %ds
    [Axy,E] = armorf(tempdataRS     ,size(tempdataZ,3),timewin_points,model_order); 
    
    % code below is adapted from bsmart toolbox function pwcausal.m
    % corrected covariance
    eyx = E(2,2) - E(1,2)^2/E(1,1); %ofc-->ds
    exy = E(1,1) - E(2,1)^2/E(2,2); %ds-->ofc
    N = size(E,1);
    
    for fi=1:length(frequencies)
        % transfer matrix (note the similarity to Fourier transform)
        H = eye(N);
        for m = 1:model_order
            H = H + Axy(:,(m-1)*N+1:m*N)*exp(-1i*m*2*pi*frequencies(fi)/dsrate);
        end
        Hi = inv(H);
        S  = H\E*Hi'/dsrate;
        
        % granger prediction per frequency
        tf_granger(1,fi,timei,chan) = log( abs(S(2,2))/abs(S(2,2)-(Hi(2,1)*exy*conj(Hi(2,1)))/dsrate) ); %ofc -->ds granger for y-->x(2)
        tf_granger(2,fi,timei,chan) = log( abs(S(1,1))/abs(S(1,1)-(Hi(1,2)*eyx*conj(Hi(1,2)))/dsrate) ); %ds -->ofc granger for y-->x(1)
    end %end freq loop
     clear eyx exy N Axy E Ay Ax tempdataRS tempdataZ tempdata
  end %end of time loop
end
  GC{mou,ses,beh,trt}=tf_granger;  
  clear tf_granger H S Hi eegdata* 
 
            end
        end
    end
end

GRANG.grangerChan=GC; %mou x ses x beh x trt




% - Test plot - %
figure, set(gcf,'Name', 'Granger prediction between electrodes OFC and DS.' );
subplot(211)
contourf(times2save,frequencies,squeeze(tf_granger(1,:,:)),40,'linecolor','none')
set(gca,'clim',[0 .2])
colorbar
title([ 'OFC -> DS' ])

subplot(212)
contourf(times2save,frequencies,squeeze(tf_granger(2,:,:)),40,'linecolor','none')
set(gca,'clim',[0 .2])
colorbar
title( 'DS -> OFC')


%% Combine across mice and average
for ses=1:N_sessions
    for beh=1:N_behaviors
        for trt=1:N_treatments
            Num=size(GRANG.epochs_extended,1);
            
            temp_avgMice=zeros(2,length(frequencies),length(times2saveidx),7 );
            for mou=1:Num

                temp=GRANG.grangerChan{mou,ses,beh,trt}; %(1,fi,timei,chan)
                if size(temp,2)==0
                   continue
                end
                temp_avgMice=cat(5,temp_avgMice,temp);
                
            end
            
            temp_avgMice(:,:,:,:,1)=[];%clears first set of zeros
            avgMice=mean(temp_avgMice,5);
            AVG_MICE{ses,beh,trt}=avgMice;
            MOUSE_COMB{ses,beh,trt}=temp_avgMice; clear temp* avgMice
            
        end
    end
end
GRANG.AVG_MICE_Chan=AVG_MICE; %averaged across all mice
GRANG.COMB_MICE_Chan=MOUSE_COMB; %still has individual mice - but in same dimension now
 




%% FIND SIG ROI

%1)combine across sessions + trts (but not beh)
for beh=1:2
    for direction=1:2
    t=1;
    for trt=1:2
        s=1;
        for ses=1:5
    
        temp=squeeze(GRANG.COMB_MICE{ses,beh,trt}(direction,:,:,:));%(direction,freq,time,mouse);
        if size(temp,4)==0 %means has no mice
            continue
        end
        
        roi(:,:,t,s)=mean(temp,3); %mean across mouse so have same dimensions
            
        s=s+1;
        end
        t=t+1;
    end
    
    avg1(:,:,beh,direction)=mean(mean(roi,4),3); clear roi s t %frequency x time x behavior x direction
    
    end
end
%2) find sig ROI- no cluster correction done here!!
pval=0.05;
for beh=1:2
    for direction=1:2

      DATA=avg1(:,:,beh,direction);
     
         N=1;
          while N<10001
          ranvals=sign(randn(size(DATA)));
          randomdata=DATA(logical(ranvals==1));
            ranavg=mean(randomdata,1);
            ransd=2*std(randomdata,[],1);
            sig(N,1)=ranavg+ransd;
            N=N+1;
          end
          sig=abs(mean(sig,1)); clear N
        
          sigpix=logical(abs(DATA)>= sig);
      
        SIGPIX{direction,beh}=sigpix; clear sigpix clust* CONT nullZ N EXP
    end
end
GRANG.sigROI05=SIGPIX;
%3)graph significant ROI
beh = 1;
direction = 2;
figure
hold on
 contourf(times2save,frequencies,squeeze(avg1(:,:,beh,direction)),40,'linecolor','none')
 %contour(times2save,frequencies,squeeze(SIGPIX{1,beh}(direction,:,:)),1,'linecolor','k','linewidth',3)
    %set(gca,'clim',[0 .14])
    set(gca, 'clim',[0 .14],'YTick', [],'XTick',[])
    %colorbar
   % title( 'OFC -> DS:Significance' )
   
   
SIGPIX=GRANG.sigROI05;


%% PLOTS 

tempzero=logical(zeros(1,7,9));
tempSP=SIGPIX{1,beh}(1,:,1:8);
choiceSP=cat(3,tempSP,tempzero);

tempzero=logical(zeros(1,7,8));
tempSP=SIGPIX{1,beh}(1,:,9:17);
toneSP=cat(3,tempzero,tempSP);

% Actual plot - mouse averaged T-F spectra plot
for ses=1:1
    for beh=1:1
        
        GRAPHDATA_control=mean(GRANG.AVG_MICE_Chan{ses,beh,1},4);
        GRAPHDATA_experimental=mean(GRANG.AVG_MICE{ses,beh,2},4); 
        
    figure, %set(gcf,'Name', 'Granger prediction between electrodes OFC and DS.' );
    %treatment #1
    subplot(221)
    hold on
    contourf(times2save,frequencies,squeeze(GRAPHDATA_control(1,:,:)),40,'linecolor','none')
%     if ses==5
%     line([-250 250], [1 1]); line([-250 250], [8 8]); line([-250 -250], [1 8]);line([250 250], [1 8]);
%     line([750 1250], [1 1]); line([750 1250], [4 4]); line([750 750], [1 4]);line([1250 1250], [1 4]);
%     end
    %contour(times2save,frequencies,squeeze(choiceSP),1,'linecolor','k','linewidth',3)
    %contour(times2save,frequencies,squeeze(toneSP),1,'linecolor','k','linewidth',3,'LineStyle', '-.')
       %contour(times2save,frequencies,squeeze(SIGPIX{1,beh}(1,:,:)),1,'linecolor','w','linewidth',3)
        %contour(times2save,frequencies,squeeze(SIGPIX{1,beh}(1,:,9:17)),1,'linecolor','--w','linewidth',3)
    set(gca,'clim',[0 .035])
     %set(gca, 'clim',[0 .14], 'YTick', [],'XTick',[])
    colorbar
    colormap(jet)
    title( 'OFC -> DS: Control' )

    %figure
    subplot(222)
    hold on
    contourf(times2save,frequencies,squeeze(GRAPHDATA_control(2,:,:)),40,'linecolor','none')
     %contour(times2save,frequencies,squeeze(SIGPIX{2,beh}(:,:)),1,'linecolor','w','linewidth',3)
    set(gca,'clim',[0 .035])
     %set(gca, 'clim',[0 .14], 'YTick', [],'XTick',[])
     colormap(jet)
    %colorbar
    title( 'DS -> OFC: Control')

    %figure
    %treatment #2
    subplot(223)
    hold on
    contourf(times2save,frequencies,squeeze(GRAPHDATA_experimental(1,:,:)),40,'linecolor','none')
%     if ses==2 || ses==5
%     line([-250 250], [1 1]); line([-250 250], [8 8]); line([-250 -250], [1 8]);line([250 250], [1 8]);
%     line([750 1250], [1 1]); line([750 1250], [4 4]); line([750 750], [1 4]);line([1250 1250], [1 4]);
%     end
      %contour(times2save,frequencies,squeeze(toneSP),1,'linecolor','k','linewidth',3,'LineStyle', ':')
    set(gca,'clim',[0 .035])
    %set(gca, 'clim',[0 .14], 'YTick', [],'XTick',[])
    colormap(jet)
    %colorbar
    title( 'OFC -> DS: Experimental' )

    %figure
    subplot(224)
    hold on
    contourf(times2save,frequencies,squeeze(GRAPHDATA_experimental(2,:,:)),40,'linecolor','none')
      %contour(times2save,frequencies,squeeze(SIGPIX{2,beh}(:,:)),1,'linecolor','w','linewidth',3)
    set(gca,'clim',[0 .035])
  % set(gca, 'clim',[0 .14], 'YTick', [],'XTick',[])
     colormap(jet)
    %colorbar
    title( 'DS -> OFC: Experimental')

    end
end

% Actual plot - fq averaged plots x-time (x-session); comparison between trts w/ dots
% per mouse?
freq=1:6; %this corresponds to fq 1-9.5516 Hz 
for ses=1:5
    for beh=1:1
        for direction=2:2
            %dir 1 = OFC --> DS
            %dir 2 = DS --> OFC
            
  GRAPHDATA_control=squeeze(mean(GRANG.AVG_MICE{ses,beh,1}(direction,freq,:),2));
  GRAPHDATA_experimental=squeeze(mean(GRANG.AVG_MICE{ses,beh,2}(direction,freq,:),2));
  
  MOUSEDATA_control=squeeze(mean(GRANG.COMB_MICE{ses,beh,1}(direction,freq,:,:),2)); size(MOUSEDATA_control)
    MOUSESEM_control=squeeze(std( MOUSEDATA_control   ,[],2) ./ sqrt(size(MOUSEDATA_control,2)));
  MOUSEDATA_experimental=squeeze(mean(GRANG.COMB_MICE{ses,beh,2}(direction,freq,:,:),2)); %size(MOUSEDATA_control)
    MOUSESEM_experimental=squeeze(std(  MOUSEDATA_experimental  ,[],2) ./ sqrt(size(MOUSEDATA_experimental,2)));
    
    figure
    hold on
    %plot(times2save,GRAPHDATA_control);
    boundedline(times2save,GRAPHDATA_control,MOUSESEM_control,'k','alpha');
    boundedline(times2save,GRAPHDATA_experimental,MOUSESEM_experimental,'r','alpha');
    ylim([0 .2])
    title([ 'OFC -> DS:' num2str(ses)] )
        
        
        end
    end
end


line([-250 250], [1 1]); line([-250 250], [8 8]); line([-250 -250], [1 8]);line([250 250], [1 8]);
line([750 1250], [1 1]); line([750 1250], [4 4]); line([750 750], [1 4]);line([1250 1250], [1 4]);

%pick out avg connecivity during ROI 0-500 and at 1000ms
%ROIchoice=zeros(14,5,2,1,1);
%ROItone=zeros(14,5,2,1,1);
for beh=1:1
    for direction=1:2
        %dir 1 = OFC --> DS
        %dir 2 = DS --> OFC
        
        if beh==1  %[row, col]=find(squeeze(GRANG.sigROI05{1,1}(direction,:,:))==1); %row freq, col are time
            %choice=find(col<=find(times2save==750));
            xCH=1:4; %fq
            yCH=3:4; %time

           %tone=find(col>=find(times2save==750));
            xT=1:3; %fq
            yT=6:7; %time
              clear row col 
    %    elseif beh==2,[row, col]=find(squeeze(GRANG.sigROI05{1,2}(direction,:,:))==1);%row freq, col are time
    %       ytic=col;
    %       xtic=row;
    %       clear row col

        end
        
    for ses=1:5
        NumMice=size(GRANG.granger,1);
        for mou=1:NumMice
       
            if size(GRANG.granger{mou,ses,beh,1},1)~=0
                MOUSEDATA_control=squeeze(GRANG.granger{mou,ses,beh,1}(direction,:,:)); 
                ROIchoice(mou,ses,1,beh,direction)=mean(mean(MOUSEDATA_control(xCH,yCH,:),2),1); %fq by time
                ROItone(mou,ses,1,beh,direction)=mean(mean(MOUSEDATA_control(xT,yT,:),2),1);
            elseif size(GRANG.granger{mou,ses,beh,1},1)==0
                ROIchoice(mou,ses,1,beh,direction)=0;
                ROItone(mou,ses,1,beh,direction)=0;
            end
            
            if size(GRANG.granger{mou,ses,beh,2},1)~=0
                MOUSEDATA_experimental=squeeze(GRANG.granger{mou,ses,beh,2}(direction,:,:));
                ROIchoice(mou,ses,2,beh,direction)=mean(mean(MOUSEDATA_experimental(xCH,yCH,:),2),1);
                ROItone(mou,ses,2,beh,direction)=mean(mean(MOUSEDATA_experimental(xT,yT,:),2),1);
            elseif size(GRANG.granger{mou,ses,beh,2},1)==0
                ROIchoice(mou,ses,2,beh,direction)=0;
                ROItone(mou,ses,2,beh,direction)=0;
            end
     
           clear  MOUSEDATA_control MOUSEDATA_experimental
        end   
    end     
    end 
end

        D85=[ROIchoice(:,1,1); ROIchoice(:,1,2); ROItone(:,1,1); ROItone(:,1,2);]; 
        RS1=[ROIchoice(:,2,1); ROIchoice(:,2,2); ROItone(:,2,1); ROItone(:,2,2);]; 
        RS3=[ROIchoice(:,3,1); ROIchoice(:,3,2); ROItone(:,3,1); ROItone(:,3,2);]; 
        RS4=[ROIchoice(:,4,1); ROIchoice(:,4,2); ROItone(:,4,1); ROItone(:,4,2);]; 
        R50=[ROIchoice(:,5,1); ROIchoice(:,5,2); ROItone(:,5,1); ROItone(:,5,2);]; 
            %SAC choice         PAE choice        sac tone        pae tone
            
            
measures= {'CHOICE';'CHOICE';'CHOICE';'CHOICE';'CHOICE';'CHOICE';'CHOICE';'CHOICE';'CHOICE';'CHOICE';'CHOICE';'CHOICE';'ph';'ph';...
          'CHOICE';'CHOICE';'CHOICE';'CHOICE';'CHOICE';'CHOICE';'CHOICE';'CHOICE';'CHOICE';'CHOICE';'CHOICE';'CHOICE';'CHOICE';'CHOICE'...
          
          'TONE';'TONE';'TONE';'TONE';'TONE';'TONE';'TONE';'TONE';'TONE';'TONE';'TONE';'TONE';'ph';'ph';...
          'TONE';'TONE';'TONE';'TONE';'TONE';'TONE';'TONE';'TONE';'TONE';'TONE';'TONE';'TONE';'TONE';'TONE'};
       
Treatment={'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'ph';'ph';...
           'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE'...
           
           'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'SAC';'ph';'ph';...
           'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE';'PAE'};
       
MouNum={'2.0';'3.1';'11.0';'16.1';'17.0';'17.1';'17.2';'20.1';'22.2';'25.0';'26.0';'27.0';'ph';'ph';...
        '4.1';'14.1';'18.0';'19.0';'19.1';'21.9';'21.1';'22.1';'22.2';'24.0';'28.0';'29.0';'30.0';'31.0'...
        
        '2.0';'3.1';'11.0';'16.1';'17.0';'17.1';'17.2';'20.1';'22.2';'25.0';'26.0';'27.0';'ph';'ph';...
        '4.1';'14.1';'18.0';'19.0';'19.1';'21.9';'21.1';'22.1';'22.2';'24.0';'28.0';'29.0';'30.0';'31.0'};

    T=table(measures,MouNum,Treatment,D85,RS1,RS3,RS4,R50); GRANG.table_ROItc=T;
    writetable(T,'Granger_ROI.xls','Sheet',1); 



errorbar(2,mean(timeatpeak_exp),peaktime_SEM_exp)
scatter(repmat(2,1,size(MOUSEDATA_experimental,2)),timeatpeak_exp)
view(-90,90)






%% GluN2B - ROI definition

ROI=zeros(8,5,2,1,2); %mou x session x conditions x behavior x direction
for beh=1:1
    for direction=1:2
        %dir 1 = OFC --> DS
        %dir 2 = DS --> OFC
        
        if beh==1, [row, col]=find(squeeze(GRANG.sigROI05{direction,1}(:,:))==1); %row freq, col are time
            xfreq=row;
            ytime=col;
            
              clear row col 
    %    elseif beh==2,[row, col]=find(squeeze(GRANG.sigROI05{1,2}(direction,:,:))==1);%row freq, col are time
    %       ytic=col;
    %       xtic=row;
    %       clear row col

        end
        
    for ses=1:5
        NumMice=size(GRANG.granger,1);
        for mou=1:NumMice
       
            if size(GRANG.granger{mou,ses,beh,1},1)~=0
                MOUSEDATA_control=squeeze(GRANG.granger{mou,ses,beh,1}(direction,:,:)); 
                ROI(mou,ses,1,beh,direction)=mean(mean(MOUSEDATA_control(xfreq, ytime,:),2),1);
            elseif size(GRANG.granger{mou,ses,beh,1},1)==0
                ROI(mou,ses,1,beh,direction)=0;
            end
            
            if size(GRANG.granger{mou,ses,beh,2},1)~=0
                MOUSEDATA_experimental=squeeze(GRANG.granger{mou,ses,beh,2}(direction,:,:));
                ROI(mou,ses,2,beh,direction)=mean(mean(MOUSEDATA_experimental(xfreq, ytime,:),2),1);
            elseif size(GRANG.granger{mou,ses,beh,2},1)==0
                ROI(mou,ses,2,beh,direction)=0;
            end
     
           clear  MOUSEDATA_control MOUSEDATA_experimental
        end   
    end     
    end 
end

        D85=[ROI(:,1,1,1); ROI(:,1,2,1); ROI(:,1,1,2); ROI(:,1,2,2)]; 
        RS1=[ROI(:,2,1,1); ROI(:,2,2,1); ROI(:,2,1,2); ROI(:,2,2,2)]; 
        RS3=[ROI(:,3,1,1); ROI(:,3,2,1); ROI(:,3,1,2); ROI(:,3,2,2)]; 
        RS4=[ROI(:,4,1,1); ROI(:,4,2,1); ROI(:,4,1,2); ROI(:,4,2,2)]; 
        R50=[ROI(:,5,1,1); ROI(:,5,2,1); ROI(:,5,1,2); ROI(:,5,2,2)]; 
            %control o>d        experimental o>d      %control d>o        experimental d>o  
            
for beh=1:1
    for direction=1:2
        %dir 1 = OFC --> DS
        %dir 2 = DS --> OFC
        
        if beh==1  %row freq, col are time
            xCH=1:4; %fq
            yCH=3:4; %time

            xT=1:3; %fq
            yT=6:7; %time
              clear row col 

        end
        
    for ses=1:5
        NumMice=size(GRANG.granger,1);
        for mou=1:NumMice
       
            if size(GRANG.granger{mou,ses,beh,1},1)~=0
                MOUSEDATA_control=squeeze(GRANG.grangerChan{mou,ses,beh,1}(direction,:,:,:)); 
                ROIchoice(mou,ses,:,beh,direction,1)=squeeze(mean(    mean(MOUSEDATA_control(xCH,yCH,:),2)    ,1)); %fq by time
                ROItone(mou,ses,:,beh,direction,1)=squeeze(mean(mean(MOUSEDATA_control(xT,yT,:),2),1));
            elseif size(GRANG.grangerChan{mou,ses,beh,1},2)==0
                ROIchoice(mou,ses,:,beh,direction,1)=0;
                ROItone(mou,ses,:,beh,direction,1)=0;
            end
            
            if size(GRANG.grangerChan{mou,ses,beh,2},1)~=0
                MOUSEDATA_experimental=squeeze(GRANG.grangerChan{mou,ses,beh,2}(direction,:,:,:));
                ROIchoice(mou,ses,:,beh,direction,2)=squeeze(mean(mean(MOUSEDATA_experimental(xCH,yCH,:),2),1));
                ROItone(mou,ses,:,beh,direction,2)=squeeze(mean(mean(MOUSEDATA_experimental(xT,yT,:),2),1));
            elseif size(GRANG.granger{mou,ses,beh,2},1)==0
                ROIchoice(mou,ses,:,beh,direction,2)=0;
                ROItone(mou,ses,:,beh,direction,2)=0;
            end
     
           clear  MOUSEDATA_control MOUSEDATA_experimental
        end   
    end     
    end 
end
            
        D85=[reshape(squeeze(ROIchoice(:,1,:,1,1,1))',7*NumMice,1 ); reshape(squeeze(ROIchoice(:,1,:,1,1,2))',7*NumMice,1 ); reshape(squeeze(ROIchoice(:,1,:,1,2,1))',7*NumMice,1 ); reshape(squeeze(ROIchoice(:,1,:,1,2,2))',7*NumMice,1 );...
             reshape(squeeze(ROItone(:,1,:,1,1,1))',7*NumMice,1 );   reshape(squeeze(ROItone(:,1,:,1,1,2))',7*NumMice,1 );  reshape(squeeze(ROItone(:,1,:,1,2,1))',7*NumMice,1 );   reshape(squeeze(ROItone(:,1,:,1,2,2))',7*NumMice,1 )]; 
        
        RS1=[reshape(squeeze(ROIchoice(:,2,:,1,1,1))',7*NumMice,1 ); reshape(squeeze(ROIchoice(:,2,:,1,1,2))',7*NumMice,1); reshape(squeeze(ROIchoice(:,2,:,1,2,1))',7*NumMice,1); reshape(squeeze(ROIchoice(:,2,:,1,2,2))',7*NumMice,1);...
             reshape(squeeze(ROItone(:,2,:,1,1,1))',7*NumMice,1 );   reshape(squeeze(ROItone(:,2,:,1,1,2))',7*NumMice,1);  reshape(squeeze(ROItone(:,2,:,1,2,1))',7*NumMice,1);   reshape(squeeze(ROItone(:,2,:,1,2,2))',7*NumMice,1)];
         
        RS3=[reshape(squeeze(ROIchoice(:,3,:,1,1,1))',7*NumMice,1); reshape(squeeze(ROIchoice(:,3,:,1,1,2))',7*NumMice,1); reshape(squeeze(ROIchoice(:,3,:,1,2,1))',7*NumMice,1); reshape(squeeze(ROIchoice(:,3,:,1,2,2))',7*NumMice,1);...
             reshape(squeeze(ROItone(:,3,:,1,1,1))',7*NumMice,1);   reshape(squeeze(ROItone(:,3,:,1,1,2))',7*NumMice,1);  reshape(squeeze(ROItone(:,3,:,1,2,1))',7*NumMice,1);   reshape(squeeze(ROItone(:,3,:,1,2,2))',7*NumMice,1)];
         
        RS4=[reshape(squeeze(ROIchoice(:,4,:,1,1,1))',7*NumMice,1); reshape(squeeze(ROIchoice(:,4,:,1,1,2))',7*NumMice,1); reshape(squeeze(ROIchoice(:,4,:,1,2,1))',7*NumMice,1); reshape(squeeze(ROIchoice(:,4,:,1,2,2))',7*NumMice,1);...
             reshape(squeeze(ROItone(:,4,:,1,1,1))',7*NumMice,1);   reshape(squeeze(ROItone(:,4,:,1,1,2))',7*NumMice,1);  reshape(squeeze(ROItone(:,4,:,1,2,1))',7*NumMice,1);   reshape(squeeze(ROItone(:,4,:,1,2,2))',7*NumMice,1)];
         
        R50=[reshape(squeeze(ROIchoice(:,5,:,1,1,1))',7*NumMice,1); reshape(squeeze(ROIchoice(:,5,:,1,1,2))',7*NumMice,1); reshape(squeeze(ROIchoice(:,5,:,1,2,1))',7*NumMice,1); reshape(squeeze(ROIchoice(:,5,:,1,2,2))',7*NumMice,1);...
             reshape(squeeze(ROItone(:,5,:,1,1,1))',7*NumMice,1);   reshape(squeeze(ROItone(:,5,:,1,1,2))',7*NumMice,1);  reshape(squeeze(ROItone(:,5,:,1,2,1))',7*NumMice,1);   reshape(squeeze(ROItone(:,5,:,1,2,2))',7*NumMice,1)];
               %control o>d        experimental o>d      %control d>o        experimental d>o choice then tone 
            
            
Direction=repmat([repmat({'OFC-DS'},7,1); repmat({'OFC-DS'},7,1);repmat({'OFC-DS'},7,1); repmat({'OFC-DS'},7,1);repmat({'OFC-DS'},7,1); repmat({'OFC-DS'},7,1);repmat({'OFC-DS'},7,1); repmat({'ph'},7,1);...  %control
          repmat({'OFC-DS'},7,1); repmat({'OFC-DS'},7,1);repmat({'OFC-DS'},7,1); repmat({'OFC-DS'},7,1);repmat({'OFC-DS'},7,1); repmat({'OFC-DS'},7,1);repmat({'OFC-DS'},7,1); repmat({'OFC-DS'},7,1);... %exp
        
          repmat({'DS-OFC'},7,1); repmat({'DS-OFC'},7,1);repmat({'DS-OFC'},7,1); repmat({'DS-OFC'},7,1);repmat({'DS-OFC'},7,1); repmat({'DS-OFC'},7,1);repmat({'DS-OFC'},7,1); repmat({'ph'},7,1);...
          repmat({'DS-OFC'},7,1); repmat({'DS-OFC'},7,1);repmat({'DS-OFC'},7,1); repmat({'DS-OFC'},7,1);repmat({'DS-OFC'},7,1); repmat({'DS-OFC'},7,1);repmat({'DS-OFC'},7,1); repmat({'DS-OFC'},7,1)],2,1);
      
Time = [repmat({'CHOICE'},7,1); repmat({'CHOICE'},7,1);repmat({'CHOICE'},7,1); repmat({'CHOICE'},7,1);repmat({'CHOICE'},7,1); repmat({'CHOICE'},7,1);repmat({'CHOICE'},7,1); repmat({'ph'},7,1);...  %control
       repmat({'CHOICE'},7,1); repmat({'CHOICE'},7,1);repmat({'CHOICE'},7,1); repmat({'CHOICE'},7,1);repmat({'CHOICE'},7,1); repmat({'CHOICE'},7,1);repmat({'CHOICE'},7,1); repmat({'CHOICE'},7,1);... %exp
        
       repmat({'CHOICE'},7,1); repmat({'CHOICE'},7,1);repmat({'CHOICE'},7,1); repmat({'CHOICE'},7,1);repmat({'CHOICE'},7,1); repmat({'CHOICE'},7,1);repmat({'CHOICE'},7,1); repmat({'ph'},7,1);...  %control
       repmat({'CHOICE'},7,1); repmat({'CHOICE'},7,1);repmat({'CHOICE'},7,1); repmat({'CHOICE'},7,1);repmat({'CHOICE'},7,1); repmat({'CHOICE'},7,1);repmat({'CHOICE'},7,1); repmat({'CHOICE'},7,1);...
          
       repmat({'TONE'},7,1); repmat({'TONE'},7,1);repmat({'TONE'},7,1); repmat({'TONE'},7,1);repmat({'TONE'},7,1); repmat({'TONE'},7,1);repmat({'TONE'},7,1); repmat({'ph'},7,1);...
       repmat({'TONE'},7,1); repmat({'TONE'},7,1);repmat({'TONE'},7,1); repmat({'TONE'},7,1);repmat({'TONE'},7,1); repmat({'TONE'},7,1);repmat({'TONE'},7,1); repmat({'TONE'},7,1);...
          
       repmat({'TONE'},7,1); repmat({'TONE'},7,1);repmat({'TONE'},7,1); repmat({'TONE'},7,1);repmat({'TONE'},7,1); repmat({'TONE'},7,1);repmat({'TONE'},7,1); repmat({'ph'},7,1);...
       repmat({'TONE'},7,1); repmat({'TONE'},7,1);repmat({'TONE'},7,1); repmat({'TONE'},7,1);repmat({'TONE'},7,1); repmat({'TONE'},7,1);repmat({'TONE'},7,1); repmat({'TONE'},7,1)];
       
Treatment=repmat([repmat({'FLOX'},7,1); repmat({'FLOX'},7,1);repmat({'FLOX'},7,1); repmat({'FLOX'},7,1);repmat({'FLOX'},7,1); repmat({'FLOX'},7,1);repmat({'FLOX'},7,1); repmat({'ph'},7,1);...
                 repmat({'CRE+'},7,1); repmat({'CRE+'},7,1);repmat({'CRE+'},7,1); repmat({'CRE+'},7,1);repmat({'CRE+'},7,1); repmat({'CRE+'},7,1);repmat({'CRE+'},7,1); repmat({'CRE+'},7,1);...           

                repmat({'FLOX'},7,1); repmat({'FLOX'},7,1);repmat({'FLOX'},7,1); repmat({'FLOX'},7,1);repmat({'FLOX'},7,1); repmat({'FLOX'},7,1);repmat({'FLOX'},7,1); repmat({'ph'},7,1);...
                repmat({'CRE+'},7,1); repmat({'CRE+'},7,1);repmat({'CRE+'},7,1); repmat({'CRE+'},7,1);repmat({'CRE+'},7,1); repmat({'CRE+'},7,1);repmat({'CRE+'},7,1); repmat({'CRE+'},7,1)],2,1);
       
MouNum=repmat([repmat({'2.2'},7,1); repmat({'3.2'},7,1);repmat({'4.0'},7,1); repmat({'4.1'},7,1);repmat({'5.0'},7,1); repmat({'7.1'},7,1);repmat({'7.2'},7,1); repmat({'ph'},7,1);...  %control
       repmat({'2.0'},7,1); repmat({'3.0'},7,1);repmat({'3.1'},7,1); repmat({'5.1'},7,1);repmat({'6.0'},7,1); repmat({'6.1'},7,1);repmat({'7.0'},7,1); repmat({'8.0'},7,1);... %exp
        
       repmat({'2.2'},7,1); repmat({'3.2'},7,1);repmat({'4.0'},7,1); repmat({'4.1'},7,1);repmat({'5.0'},7,1); repmat({'7.1'},7,1);repmat({'7.2'},7,1); repmat({'ph'},7,1);...
       repmat({'2.0'},7,1); repmat({'3.0'},7,1);repmat({'3.1'},7,1); repmat({'5.1'},7,1);repmat({'6.0'},7,1); repmat({'6.1'},7,1);repmat({'7.0'},7,1); repmat({'8.0'},7,1)],2,1);
   
Channel=repmat(repmat({'1';'2';'3';'4';'5';'6';'7'},32,1),2,1);%       

    T=table(Direction,Time,MouNum,Channel,Treatment,D85,RS1,RS3,RS4,R50); GRANG.table_ROItc_Chan=T;
    writetable(T,'Granger_ROI_DRBR_refChan_allChan.xls','Sheet',1); 
    
 %% Line Graphing of ROI
 
 for type=1:2
     if type==1, data=ROIchoice;
     elseif type==2, data = ROItone;
     end

       for ses=1:5
           
            tempdata=data(:,ses,1); %control
            tempdata(tempdata==0)=[];
              means_cont(1,ses)=mean(tempdata);
              errors_cont(1,ses)=std(tempdata,0,1)/ sqrt(size(tempdata,1)); clear tempdata
              
            tempdata=data(:,ses,2); %experimental
            tempdata(tempdata==0)=[];
              means_exp(1,ses)=mean(tempdata);
              errors_exp(1,ses)=std(tempdata,0,1)/ sqrt(size(tempdata,1)); clear tempdata
        end       

    figure %plot the individual points
    hold on
    %control
    scatter(ones(length(data(:,1,1)),1),data(:,1,1),90,'dk','LineWidth',1); %^ = triangle
    scatter(repmat(2,length(data(:,2,1)),1),data(:,2,1),90,'dk','LineWidth',1);
    scatter(repmat(3,length(data(:,3,1)),1),data(:,3,1),90,'dk','LineWidth',1);
    scatter(repmat(4,length(data(:,4,1)),1),data(:,4,1),90,'dk','LineWidth',1);
    scatter(repmat(5,length(data(:,5,1)),1),data(:,5,1),90,'dk','LineWidth',1);
    
    %experimental
    scatter(ones(length(data(:,1,2)),1),data(:,1,2),90,'dr','LineWidth',1);
    scatter(repmat(2,length(data(:,2,2)),1),data(:,2,2),90,'dr','LineWidth',1);
    scatter(repmat(3,length(data(:,3,2)),1),data(:,3,2),90,'dr','LineWidth',1);
    scatter(repmat(4,length(data(:,4,2)),1),data(:,4,2),90,'dr','LineWidth',1);
    scatter(repmat(5,length(data(:,5,2)),1),data(:,5,2),90,'dr','LineWidth',1);
    xlim([0 6])
    
    ylim([0 .4])
    
%     if fqr==1 || fqr==2
%         ylim([0 6]) %10,10, 5
%     elseif fqr==2
%         ylim([0 6])
%     elseif fqr==3
%         ylim([0 2.5])
%     end

    hold on %plot the average + standard error
    errorbar(means_cont,errors_cont,'-sk','MarkerFaceColor','k','markers',12,'LineWidth',2) %12
    errorbar(means_exp,errors_exp,'-sr','MarkerFaceColor','r','markers',12,'LineWidth',2)


clear data

 end
 
 
 
 
 
 
 %% stationarity test - results

 for mou=1:8
     for ses=1:5
         for beh=1:2
             for trt=1:2
 
 temp_kstat=Ksstat(:,:,mou,ses,beh,trt);
 Kcval=Cval_k(:,1,1,1,1);
 
 
 non_stat_o=0; non_stat_d=0; %total count of all timesegments that are not stationary --> was calculated over all trials
for r=1:2
 for t=1:length(Kcval) %number of time segments
     if r==1,
         if temp_kstat(t,r) > Kcval
         non_stat_d = non_stat_d+1;
         end
     elseif r==2,
         if temp_kstat(t,r) > Kcval 
         non_stat_o = non_stat_o+1;
         end
     end
 end
Non_Stat(1,mou,ses,beh,trt)=non_stat_d; %number of timesegments x- all trials that were non-stationary in ds
Non_Stat(2,mou,ses,beh,trt)=non_stat_o; %number of timesegments x- all trials that were non-stationary in ofc
end

 
             end
         end
     end
 end
 

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 

