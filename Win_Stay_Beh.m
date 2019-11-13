%% win-stay lose-shift trial count

%Types 1= correct -> correct ==> use these as win-stay trials
%      2= correct -> incorrect
%      3= incorrect -> incorrect ==> use these as lose-shift # req to shift beh
%      4= incorrect-> correct 

%before can run this script need:
    %SINGLEUNIT structure with .neuroexp and .ML variables

for ses=1:5
    for trt=1:2
    Num=size(SINGLEUNIT.neuroexp(:,:,ses,:,trt),1);  %number of mice
    for s=1:Num
        if size(SINGLEUNIT.ML{s,ses,trt},1)==0
            continue
        end
        
%trial numbers for TC and TIC
   code=SINGLEUNIT.ML{s,ses,trt};
    
        count=1; string=1;
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

    %count # of Type 1s per mouse (repeated corrects only!!)
     newWS=logical(TYPE(:,:)==1);
     newlose=logical(TYPE(:,:)==2); %when make an error
     AverageWinStay(s,ses,trt)=sum(newWS)./sum(newlose);
     
     TotWinStay(s,ses,trt)=sum(newWS);
     TotWinShift(s,ses,trt)=sum(logical(TYPE(:,:)==2));
     
     %count # of Types 3 (persev only)
     newLS=logical(TYPE(:,:)==4); %when shift back to a win  
     newwin=logical(TYPE(:,:)==2); 
     AverageLoseShift(s,ses,trt)= sum(newLS)./sum(newwin);
     
     TotalLoseShift(s,ses,trt)=sum(logical(TYPE(:,:)==4));
     TotalLoseStay(s,ses,trt)=sum(logical(TYPE(:,:)==3));

    
    clear TYPE code 
    
    end
    end
end 



% export to excel

temp=reshape(AverageWinStay,8,10);
xlswrite('DRBR_WinStay',temp,'Sheet1');

temp=reshape(AverageLoseShift,8,10);
xlswrite('DRBR_WinStay',temp,'Sheet2');



temp=reshape(TotWinStay,8,10);
xlswrite('DRBR_WinStayTot',temp,'Sheet1');

temp=reshape(TotalLoseShift,8,10);
xlswrite('DRBR_WinStayTot',temp,'Sheet2');

temp=reshape(TotalLoseStay,8,10);
xlswrite('DRBR_WinStayTot',temp,'Sheet3');

temp=reshape(TotWinShift,8,10);
xlswrite('DRBR_WinStayTot',temp,'Sheet4');






















