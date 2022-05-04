function [dat470,dat405,ts,timing] = tdtCONVERT(details, tankdir, tankname, blockname)
%  File to convert all TDT files within
%  static refers to the static system in the Bruchas lab
%  auto refers to the automated stim protocol on the Bruchas cart

fearON = details.fearON;
stimON = details.stimON;
static = details.static;
auto = details.auto;

if static == 0
    storenames = {'470A'}; % name of stores to extract from TDT (usu. 4-letter code) 
else
    storenames = {'465A'}; % name of stores to extract from TDT (usu. 4-letter code) 
end
%LMag is the demodulated data, may also have other timestamps etc

storenames2 = {'405A'};

storenames3 = {'Epo1'}; %stim

storenames4 = {'Epo2'}; % Fear 

if auto
    storenames3 = {'Wi2/'}; %'Per2' or 'Wi2/'
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract

for k = 1:numel(storenames)
  storename = storenames{k};
  S{k} = tdt2mat(tankdir, tankname, blockname, storename);
  
  storename2 = storenames2{k};
  S2{k} = tdt2mat(tankdir, tankname, blockname, storename2);
if stimON
  storename3 = storenames3{k};
  S3{k} = tdt2mat(tankdir, tankname, blockname, storename3);
end
if fearON
  storename4 = storenames4{k};
  S4{k} = tdt2mat(tankdir, tankname, blockname, storename4);
end
end


% Massage data and get time stamps

LMag = S{1}; %add more if you extracted more stores above
% LMag2 = S{2};
% For 2-color rig, LMag data is on channels 1 and 2, channel 1 = 470nm, channel 2 = 405nm
chani1 = LMag.channels==1;
chani2 = LMag.channels==2;

LMag2 = S2{1}; %add more if you extracted more stores above
% LMag2 = S{2};
% For 2-color rig, LMag data is on channels 1 and 2, channel 1 = 470nm, channel 2 = 405nm
chani21 = LMag2.channels==1;
chani22 = LMag2.channels==2;
% chani21 = LMag2.channels==1;
% chani22 = LMag2.channels==2;

if stimON
LMag3 = S3{1}; %add more if you extracted more stores above
% LMag3 = S3{2};
% For 2-color rig, LMag data is on channels 1 and 2, channel 1 = 470nm, channel 2 = 405nm
chani31 = LMag3.channels==1;
chani32 = LMag3.channels==2;
end

if fearON
LMag4 = S4{1}; %add more if you extracted more stores above
% LMag3 = S3{2};
% For 2-color rig, LMag data is on channels 1 and 2, channel 1 = 470nm, channel 2 = 405nm
chani41 = LMag4.channels==1;
chani42 = LMag4.channels==2;
end

% Get LMag data as a vector (repeat for each channel)
dat470 = LMag.data(chani1,:);
dat470 = reshape(dat470', [],1); % unwrap data from m x 256 array
% dat405 = LMag.data(chani21,:);
% dat405 = reshape(dat405', [],1); % unwrap data from m x 256 array

% Get LMag timestamps (use chani1 - timestamps should be the same for all Wpht channels
ts = LMag.timestamps(chani1);
t_rec_start = ts(1);

ts = ts-ts(1); % convert from Unix time to 'seconds from block start'
ts = bsxfun(@plus, ts(:), (0:LMag.npoints-1)*(1./LMag.sampling_rate));
ts = reshape(ts',[],1);

%%%%%%%%%%%%%%%%%%


dat405 = LMag2.data(chani21,:);
dat405 = reshape(dat405', [],1); % unwrap data from m x 256 array
% dat405 = LMag.data(chani21,:);
% dat405 = reshape(dat405', [],1); % unwrap data from m x 256 array

% Get LMag timestamps (use chani1 - timestamps should be the same for all Wpht channels
ts2 = LMag2.timestamps(chani21);
t_rec_start2 = ts2(1);
dt = datetime( t_rec_start2, 'ConvertFrom', 'posixtime' );

ts2 = ts2-ts2(1); % convert from Unix time to 'seconds from block start'
ts2 = bsxfun(@plus, ts2(:), (0:LMag2.npoints-1)*(1./LMag2.sampling_rate));
ts2 = reshape(ts2',[],1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For stim and Fear 

if stimON
dat3 = LMag3.data(chani31,:);
dat3 = reshape(dat3', [],1); % unwrap data from m x 256 array

% Get LMag timestamps (use chani1 - timestamps should be the same for all Wpht channels
ts3 = LMag3.timestamps(chani31);

%multiple stims
if any(diff(ts3)>1)
    dTs3 = find(diff(ts3)>1);
    stimLocs = 1;
    for i = 1:length(dTs3)
        stimLocs = [stimLocs dTs3(i)+1];
    end
end
if exist('stimLocs')
    t_rec_start3 = ts3(1);
    dt3 = datetime( t_rec_start3, 'ConvertFrom', 'posixtime' );
    ts3St = ts3(1); % start of ts3
    ts3 = ts3-ts3(1); % convert from Unix time to 'seconds from block start'
    timing = seconds(dt3 - dt);
    for i = 2:length(stimLocs)
        timing = [timing timing(1)+ts3(stimLocs(i))];
    end
else
    t_rec_start3 = ts3(1);
    dt3 = datetime( t_rec_start3, 'ConvertFrom', 'posixtime' );
    ts3St = ts3(1); % start of ts3
    ts3 = ts3-ts3(1); % convert from Unix time to 'seconds from block start'
    timing = seconds(dt3 - dt);
end
end

if fearON
dat4 = LMag4.data(chani41,:);
dat4 = reshape(dat4', [],1); % unwrap data from m x 256 array

% Get LMag timestamps (use chani1 - timestamps should be the same for all Wpht channels
ts4 = LMag4.timestamps(chani41);
t_rec_start4 = ts4(1:length(ts4)); t_rec_start4=t_rec_start4(logical(dat4));
dt4 = datetime( t_rec_start4, 'ConvertFrom', 'posixtime' ); 

ts4 = ts4-ts4(1); % convert from Unix time to 'seconds from block start'
ts4 = bsxfun(@plus, ts4(:), (0:LMag4.npoints-1)*(1./LMag4.sampling_rate));
ts4 = reshape(ts4',[],1);
timing = seconds(dt4 - dt); 

end

if fearON && stimON
    warning('timing was overwritten, you did fearTTLs and stim!')
    warning('write more code')
end



end