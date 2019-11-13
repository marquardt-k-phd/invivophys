
%Run this file for creating/ updating a History File from a Discrimination
%Reversal from the potocol 288100 on the recording side (with lever shadow
%stage)

disp('This MatLab program will help you create or add data to an excel History File.')
disp('****Before you begin you MUST save your combined file as an .xlsx');

disp('Just follow the on screen prompts:')

prompt = 'Are you creating a new History File? Y/N: ';
str = input(prompt,'s');
if str=='N'
    
    disp('Please Chose the K-limbic Combined file to add');
    [filename, pathname]=uigetfile('*.xlsx','Pick a file');
    addpath(pathname);
    [num,txt,raw]=xlsread(filename); %open file to analyze
    
   DAILYDATA=raw;
   DAILYtxt=txt;
   TF=strfind(DAILYtxt,'Subject Id');
     subject=find(~cellfun(@isempty,TF));
        SubjectIDs=DAILYtxt(subject,2); clear subject TF
    
    TF=strfind(DAILYtxt,'Duration');
       time=find(~cellfun(@isempty,TF));
           Time=DAILYDATA(time,2); 
           Time=cell2mat(Time) ./ 6000; %60sec/min and then to convert to sec instead of deci-seconds
                min=fix(Time);
                sec=((Time-min)*60)*.01; %since format is .71 of a min -> must convert to sec 
                DailyTime=min+sec; clear Time min sec TF %finding time
      
     TF=strfind(DAILYtxt,'AC Comment');
        problem=find(~cellfun(@isempty,TF));
           DailyProblem=DAILYtxt(problem,3); clear problem TF %which problem they were on - for record keeping
           
    TF=strfind(DAILYtxt,'Date');
        dates=find(~cellfun(@isempty,TF));
           DateofRun=DAILYtxt(dates,2); clear problem TF %which date the run occured on
        
    TF=strfind(DAILYtxt,'Stage (4)');
     [row col]=find(~cellfun(@isempty,TF));
        stagerow=row(1:size(SubjectIDs,1),:)-3; %finding state row start of each subject
        stage5entrycol=col(1,1); clear TF row col %finding the column of stage 5 - rt
        
     TF=strfind(DAILYtxt,'Stage (9)');
     [row col]=find(~cellfun(@isempty,TF));
        stage10entrycol=col(1,1); clear TF row col%finding column of stage 10 - magazine
        
       [Selection]= listdlg('ListString',SubjectIDs, 'Name','Subject ID List', 'PromptString', 'Select subjects to process'); %selects subjects to process
        SubjectIDs=SubjectIDs(Selection,:);
        stagerow=stagerow(Selection,:);
        DailyProblem=DailyProblem(Selection,:); 
        DailyTime=DailyTime(Selection,:);
        DateofRun=DateofRun(Selection,:);
        
    %find all the data stuff
    for sub=1:size(SubjectIDs, 1)
           startrow=stagerow(sub,1);
       for n=1:200 %for 200 trials --> as far as I know no mouse has done over 150 in one session, including RD1
         if num(startrow,stage5entrycol) <= 0 || isnan(num(startrow,stage5entrycol))==1 %the last row is zero so this will end this loop when it reaches that point or hits a nan
          continue
         end
        tempRT(n,1)=(num(startrow,stage5entrycol+1)-num(startrow,stage5entrycol)) ./ 100;
        tempML(n,1)=(num(startrow,stage10entrycol+1)-num(startrow,stage10entrycol)) ./ 100;
        RT{1,sub}=tempRT;
        ML{1,sub}=tempML;
        startrow=startrow+1;
       end
            %find number of different trial types
             for t=1:size(tempRT,1)
                if t==1 && tempML(t,1)>0
                    ttypes(t,1)=1; %first presentaiton correct - for first trial
                elseif t==1 && tempML(t,1)==0
                    ttypes(t,1)=2; %first presentation incorrect - for first trial
                elseif t>1 && tempML(t,1)>0 && tempML(t-1,1)>0
                    ttypes(t,1)=1; % first presentation correct
                elseif t>1 && tempML(t,1) ==0 && tempML(t-1,1)>0
                    ttypes(t,1)=2; % first presentation incorrect
                elseif t>1 && tempML(t,1)>0 && tempML(t-1,1)==0
                    ttypes(t,1)=3; % correct correction trial error
                elseif t>1 && tempML(t,1)==0 && tempML(t-1,1)==0
                    ttypes(t,1)=4; % incorrect correction trial 

                end 
             end
            
       RTavg(1,sub)=mean(tempRT,1); %find RT
         err=size(find(tempML==0),1);
         tempML(tempML==0)=[]; %remove zeros so don't contaminate the ML
       MLavg(1,sub)=mean(tempML,1); %find ML
       PerCorr(1,sub)=(size(find(ttypes==1),1) ./ 30) .* 100;
       Tc(1,sub)=size(find(ttypes==1),1);
       Tic(1,sub)=size(find(ttypes==2),1);
       Ct(1,sub)=size(find(ttypes==3),1) + size(find(ttypes==4),1);
       
       TotTrials(1,sub)=size(find(ttypes==3),1) + size(find(ttypes==4),1)+size(find(ttypes==2),1)+size(find(ttypes==1),1); %find total trial # 
       TotErr(1,sub)=size(find(ttypes==3),1) + size(find(ttypes==4),1)+size(find(ttypes==2),1); %double checks error count
       
       clear temp*  ttypes
    end

    for sub=1:size(SubjectIDs, 1)

        INFO{sub,1}=DateofRun{sub,1};
        INFO{sub,2}=DailyTime(sub,1);
        INFO{sub,3}=DailyProblem{sub,1};
        INFO{sub,4}=Tc(1,sub); 
        INFO{sub,5}=Tic(1,sub); 
        INFO{sub,6}=Ct(1,sub); 
        INFO{sub,7}=RTavg(1,sub); 
        INFO{sub,8}=MLavg(1,sub); 
        INFO{sub,9}=TotTrials(1,sub); 
        INFO{sub,10}=TotErr(1,sub); 
        INFO{sub,11}=PerCorr(1,sub);
        
    end
    
    
    disp('Open History File'); %open old history file
    [filename, pathname]=uigetfile('*.xls','Choose History File');
        cd(pathname);
        [status, sheets]=xlsfinfo(filename);
        sheets(:,1)=[]; %remove that first sheet, which is always blank
        
        %determine if there are new sheets to add...
        newsubj=zeros(size(SubjectIDs,1),1);
        for sub=1:size(SubjectIDs, 1)
            subinhistoryfile=sheets;
            subinklimbic=SubjectIDs{sub,1};
            ispresent = strmatch(subinklimbic, subinhistoryfile, 'exact');
            if size(ispresent,1)==0
                newsubj(sub,1)=sub;
            end
            clear ispresent
        end
   
        for sub=1:size(SubjectIDs, 1)
            sheet=SubjectIDs{sub,1};
            isnew=newsubj(sub,1);
            if isnew>0
               NewRange{1,sub}='A2:K2'; 
            elseif isnew==0
                OldData =xlsread(filename,sheet);% need to know what subjects where there previously so can add new ones
                day(1,sub)=size(OldData,1); clear OldData %find out which session this is so can create new range
                NewRange{1,sub}=(['A' num2str(day(1,sub)+2) ':K' num2str(day(1,sub)+2)]);
            end
            clear OldData day
        end
        
headerinfo={'DATE', 'Time (min.sec)', 'PROBLEM', 'TC','TIC','CT','RT','ML', 'TOT TRIALS', 'TOT ERRORS', '% CORR'};
    for sub=1:size(SubjectIDs, 1)
        sheet=SubjectIDs{sub,1};
        range=NewRange{1,sub};
        isnew=newsubj(sub,1);
        if isnew>0
            xlswrite(filename, headerinfo, sheet); %for new subjects need to also add the header
            xlswrite(filename, INFO(sub,:),sheet,range);
        elseif isnew==0
            xlswrite(filename, INFO(sub,:),sheet,range);
        end   
    end

    disp( 'You have added data to your History File, woohoo!');
    cd('E:\Documents\MATLAB');
    
    clear all
    
elseif str=='Y'
    disp('Please Chose the K-limbic file to add');
    [filename, pathname]=uigetfile('*.xlsx','Pick a file');
    addpath(pathname);
    [num,txt,raw]=xlsread(filename); %open file to analyze
    
   DAILYDATA=raw;
   DAILYtxt=txt;
   TF=strfind(DAILYtxt,'Subject Id');
     subject=find(~cellfun(@isempty,TF));
        SubjectIDs=DAILYtxt(subject,2); clear subject TF
    
    TF=strfind(DAILYtxt,'Duration');
       time=find(~cellfun(@isempty,TF));
           Time=DAILYDATA(time,2); 
           Time=cell2mat(Time) ./ 6000; %60sec/min and then to convert to sec instead of deci-seconds
                min=fix(Time);
                sec=((Time-min)*60)*.01; %since format is .71 of a min -> must convert to sec 
                DailyTime=min+sec; clear Time min sec TF %finding time
      
     TF=strfind(DAILYtxt,'AC Comment');
        problem=find(~cellfun(@isempty,TF));
           DailyProblem=DAILYtxt(problem,3); clear problem TF %which problem they were on - for record keeping
           
     TF=strfind(DAILYtxt,'Date');
        dates=find(~cellfun(@isempty,TF));
           DateofRun=DAILYtxt(dates,2); clear dates TF %which date the run occured on
        
    TF=strfind(DAILYtxt,'Stage (4)');
     [row col]=find(~cellfun(@isempty,TF));
        stagerow=row(1:size(SubjectIDs,1),:)-3; %finding state row start of each subject
        stage5entrycol=col(1,1); clear TF row col %finding the column of stage 4 - rt
        
     TF=strfind(DAILYtxt,'Stage (9)');
     [row col]=find(~cellfun(@isempty,TF));
        stage10entrycol=col(1,1); clear TF row col%finding column of stage 9 - magazine
        
       [Selection]= listdlg('ListString',SubjectIDs, 'Name','Subject ID List', 'PromptString', 'Select subjects to process'); %selects subjects to process
        SubjectIDs=SubjectIDs(Selection,:);
        stagerow=stagerow(Selection,:);
        DailyProblem=DailyProblem(Selection,:); 
        DailyTime=DailyTime(Selection,:);
        DateofRun=DateofRun(Selection,:);
        
    %find all the data stuff
    for sub=1:size(SubjectIDs, 1)
           startrow=stagerow(sub,1);
      for n=1:200 %for 200 trials --> as far as I know no mouse has done over 150 in one session, including RD1
         if num(startrow,stage5entrycol) <= 0 || isnan(num(startrow,stage5entrycol))==1 %the last row is zero so this will end this loop when it reaches that point or hits a nan
          continue
         end
        tempRT(n,1)=(num(startrow,stage5entrycol+1)-num(startrow,stage5entrycol)) ./ 100;
        tempML(n,1)=(num(startrow,stage10entrycol+1)-num(startrow,stage10entrycol)) ./ 100;
        RT{1,sub}=tempRT;
        ML{1,sub}=tempML;
        startrow=startrow+1;
       end
            %find number of different trial types
             for t=1:size(tempRT,1)
                if t==1 && tempML(t,1)>0
                    ttypes(t,1)=1; %first presentaiton correct - for first trial
                elseif t==1 && tempML(t,1)==0
                    ttypes(t,1)=2; %first presentation incorrect - for first trial
                elseif t>1 && tempML(t,1)>0 && tempML(t-1,1)>0
                    ttypes(t,1)=1; % first presentation correct
                elseif t>1 && tempML(t,1) ==0 && tempML(t-1,1)>0
                    ttypes(t,1)=2; % first presentation incorrect
                elseif t>1 && tempML(t,1)>0 && tempML(t-1,1)==0
                    ttypes(t,1)=3; % correct correction trial error
                elseif t>1 && tempML(t,1)==0 && tempML(t-1,1)==0
                    ttypes(t,1)=4; % incorrect correction trial 

                end 
             end
            
       RTavg(1,sub)=mean(tempRT,1); %find RT
         err=size(find(tempML==0),1);
         tempML(tempML==0)=[]; %remove zeros so don't contaminate the ML
       MLavg(1,sub)=mean(tempML,1); %find ML
       PerCorr(1,sub)=(size(find(ttypes==1),1) ./ 30) .* 100;
       Tc(1,sub)=size(find(ttypes==1),1);
       Tic(1,sub)=size(find(ttypes==2),1);
       Ct(1,sub)=size(find(ttypes==3),1) + size(find(ttypes==4),1);
       
       TotTrials(1,sub)=size(find(ttypes==3),1) + size(find(ttypes==4),1)+size(find(ttypes==2),1)+size(find(ttypes==1),1); %find total trial # 
       TotErr(1,sub)=size(find(ttypes==3),1) + size(find(ttypes==4),1)+size(find(ttypes==2),1); %double checks error count
       
       clear temp*  ttypes
    end

    %set up a cell array for each subject to publish
    headerinfo={'DATE', 'Time (min.sec)', 'PROBLEM', 'TC','TIC','CT','RT','ML', 'TOT TRIALS', 'TOT ERRORS', '% CORR'};

    for sub=1:size(SubjectIDs, 1)

        INFO{sub,1}=DateofRun{sub,1};
        INFO{sub,2}=DailyTime(sub,1);
        INFO{sub,3}=DailyProblem{sub,1};
        INFO{sub,4}=Tc(1,sub); 
        INFO{sub,5}=Tic(1,sub); 
        INFO{sub,6}=Ct(1,sub); 
        INFO{sub,7}=RTavg(1,sub); 
        INFO{sub,8}=MLavg(1,sub); 
        INFO{sub,9}=TotTrials(1,sub); 
        INFO{sub,10}=TotErr(1,sub); 
        INFO{sub,11}=PerCorr(1,sub);
        
    end

    %next step is to assign
    prompt = 'Please name the new history file: ';
    filename=input(prompt,'s');
    
    disp('Where do you want to save this file?');
     folder=uigetdir('C:\Users\Lomas2211\Dropbox\','Where do you want to save this file?');
        cd(folder);

    for sub=1:size(SubjectIDs, 1)
        sheet=SubjectIDs{sub,1};
        xlswrite(filename, headerinfo,sheet);
        xlswrite(filename, INFO(sub,:),sheet,'A2:K2');
    end

    disp( 'You have now made a new History File, congratulations');
    cd('E:\Documents\MATLAB');
    
    clear all
end
 
