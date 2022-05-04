function [subStr, subNumb, subjectLabel, txtName] = defineSubject(photoname, region,subjectIdx)
% define names of subject and text files
% THIS FUNCTION IS VERY INDIVIDUALIZED

    subStr = [photoname(1:3) '-']; 
    subNumb   = photoname(strfind(photoname, subStr) + length(subStr)); % will get sub
    subjectLabel = photoname(1:5); 
    
    %customize special cases with if statements
    if isempty(str2num(subStr(end))) %isempty(str2num(subNumb))
        idx1 = strfind(photoname, '-') + 1; idx1=idx1(1);
        subNumb   = photoname(idx1); % will get subject number  
        if strcmp(photoname(3),'-')
            subStr = [photoname(1:2) '-'];
        end
        subjectLabel = photoname(1:idx1); 
        %subNumb   = photoname(strfind(photoname, subStr) + length(subStr))-1; % will get sub
    end
    
    %%% Make sure the Timings you use are consistent, or change them every time
    % CUSTOMIZE 
    subName{subjectIdx} = subjectLabel;
    txtName    = [subStr 'Stim_timesM' subNumb '.txt'];    
    
    if ~isfile(txtName)
        txtName = [subStr subjectLabel(end) '_' region '.txt'];
    end
      
end