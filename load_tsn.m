function [timetag,data,fR,lR,fs] = load_tsn(Nscan,offset,cfile)
if nargin == 1 % load tsn file
    cfile=Nscan;
    Nscan =100000000000000;
    offset=0;
    [timetag,data,fR,lR,fs]=load_tsn_sub(Nscan,offset,cfile);
elseif ischar(Nscan)&&ischar(offset)% load tsn file by specifical time
    stime=Nscan;
    etime=offset;
    %% load MTU tsn file for given time
    [startime,endtime,fs]=get_tsn_time(cfile);
    [startsp,endsp]=time2sample(startime,endtime,stime,etime,fs);
    sp=(endsp-startsp)+1;
    offsetn=(startsp-1)/fs;
    [timetag,data,fR,lR,fs]= load_tsn_sub(sp,offsetn,cfile);
elseif isnumeric(Nscan)&&isnumeric(offset)% load tsn file by specifical scan and offset.
    [timetag,data,fR,lR,fs] = load_tsn_sub(Nscan,offset,cfile);
else
    return
    errordlg('Please check the input parameters');
end
timetag=timetag';
function  [timetag,data,firstscan_R,lastscan_R,Nr] = load_tsn_sub(Nscan,offset,cfile)
% definition of vars:
% a 'smart' function to scan the TSn files for data
% fid : the file handlclce of the TSn file
% first let's start at the beginning of the file
% Nscan: number of data to load
% offset: time to skip at begining (in seconds)
% cfile: the data file name
% opt:   tell the programme to convert record to physical value (or not)
% Rstart : the index of the starting record
% Rend : the index of the ending record
% Nr : number of records in a single scan (between time tags)/Sampling Rate
% Nch : number of the channels/two electrical and three magnetical
% tsize : teh size (length) of the time tag should be 16 or 32 byte
% ======================================================================= %
% some notes on the ts format of tsn file
% each time record consists of three bytes,  
% let's name them byte 1, byte 2, and byte 3
% the time record should be (+ -) (byte2*256 + byte1)
% please note that the third byte is merely a sigh byte. if read in uint8
% format, it should be 255(-) or 0(+).
% ======================================================================= %
% allocate memory
% data  = zeros(Nc,Nscan); % you don't really nead this
% try to find the length of header
% if nargin<5
%     opt=0;
% end
%  	try opening the ts data file 
fid = fopen(cfile,'r','ieee-le');
if fid>0    
    %   determine tag length (16 or 32)
    fseek(fid,0,-1);
    tag = fread(fid,32,'uint8');
    tsize = tag(14);%tsize:the size of tag,16btyes for old v5,32 bytes for tsn
%     Currently, the TSH and TSL files are 16 bytes long and 32 bytes long in the TSn file. Table C-2 shows the byte allocation in the tag, and the following paragraphs give other details
    Nch = tag(13);
    fseek(fid,0,'eof');
    fsize = ftell(fid);
else
    disp(['YUKI.N> error loading file' cfile]);
end
if (tsize == 0); % we find the 'tag code' here
    tsize = 16;  % set tag length to 16
    disp('YUKI.N> TSH/L file (old v5 type) found') % for debug
% else
%     strind=strfind(cfile,'\');
%     str = cfile(strind(end)+1:end);
%     disp([str ' is loading']); % for debuging
end; 
% get the number of record in a single scan
Nr = tag(11) + tag(12)*256;
% 16 bit number tag(12)=9,tag(11)=96
% first let's start at the beginning of the file
fseek(fid,0,-1);
% allocating memory
while (Nscan+floor(offset*Nr))/Nr*(Nr*3*Nch+tsize) > fsize
     totalscans = fsize/(3*Nch*Nr+tsize)*Nr-floor(offset*Nr);
% warndlg('the load scans overflow!  ','!! Warning !!');
 % uiwait(h);
Nscan = totalscans;
end
Rstart=floor(offset)+1;%bu neng shi zai o second kaishi er shi zai 1 second shi kaishi
Rend=floor((Nscan+floor(offset*Nr)-1)/Nr)+1;
data =zeros(Nch,Nscan);
% determine the location of the first record
first_scan = floor(offset*Nr)+1; % the starting scan
firstscan_R =first_scan-(Rstart-1)*Nr;%
last_scan = first_scan+Nscan-1; % the ending scan
lastscan_R = last_scan-(Rend-1)*Nr;
Ncross = Rend - Rstart-1 ; % see how many scans we have between 
 % the starting and ending records
 timetag = zeros(tsize,Ncross+2);
if Ncross <=0 % we are in the same scan window
    jump_tsn(fid, Rstart, firstscan_R,Nr ,Nch, tsize);
    temp=fread(fid,(last_scan-first_scan+1)*Nch*3);
    Cpos=length(temp)/(Nch*3); % find the current positon
    data(:,1:Cpos)=conv_record(temp,Nch);
    timetag = read_timetag(fid,Rstart,Nr,Nch,tsize);
else % we have to cross over one or more scan windows
    % jump to the first record
    jump_tsn(fid, Rstart, firstscan_R,Nr ,Nch, tsize);
    temp=fread(fid,(Nr-firstscan_R+1)*Nch*3);
    Cpos=length(temp)/(Nch*3); % find the current positon
    data(:,1:Cpos)=conv_record(temp,Nch);
     timetag(:,1)=read_timetag(fid,Rstart,Nr,Nch,tsize);
    if Ncross >=1 % we just have to go over a time tag and scan
        for i=1:Ncross % jump to every scans between
          jump_tsn(fid,Rstart+i,1,Nr ,Nch, tsize);
            temp=fread(fid,Nr*Nch*3);%read the whole scan window
            data(:,Cpos+1:Cpos+Nr)=conv_record(temp,Nch);
            Cpos=Cpos+Nr;
            timetag(:,i+1) = read_timetag(fid,Rstart+i,Nr,Nch,tsize);
%             disp(i)
            % disp(i) for debug
        end
    end
    % jump to the last scan
    jump_tsn(fid, Rend,1,Nr ,Nch, tsize);
    temp=fread(fid,lastscan_R*Nch*3);
    if ~isempty(temp)%isempty yong yu pan duan bian liang shi fou jing guo chu shi hua ru guo shi ze wei jia fan zhi wei zhen ,'~'biao shi 'fei'
        data(:,Cpos+1:Cpos+length(temp)/(Nch*3))=conv_record(temp,Nch);
         timetag(:,Ncross+2) = read_timetag(fid,Rend,Nr,Nch,tsize);
    end
end
fclose(fid);
data=data';




%return
function position = jump_tsn(fid, Rstart,firstscan_R, Nr ,Nch, tsize)
% a 'smart' function to "jump" to a certain record in the TSn files
% fid : the file handle of the TSn file
% idx : the index of the record to "jump to"
% Nr : number of records in a single scan (between time tags)
% Nch : number of the channels
% tsize : teh size (length) of the time tag should be 16 or 32 byte
% position : the location in file 
% ========================================================== %
fseek(fid,(Rstart-1)*(tsize+Nr*Nch*3),-1);%3 means 24bit to store a num 
% The time series data is stored in a 24-bit binary complement format, with each sample being 3 bytes, and at least the first byte is meaningful. A scan is a set of samples, each one at the same time. A complete scan of the sampling time is stored consecutively in the order of the track number. (The track is numbered from 1, not from 0) The scan is stored in the order of the sampling time.
% The first scan in the record always starts exactly one UTC seconds, and the scan rate is always exactly an integer multiple of 1 Hz.
fseek(fid,tsize+(firstscan_R-1)*Nch*3,0);
position = ftell(fid);
return
function data=conv_record(temp,Nc)
% a 'smart' function to convert pheonix TSn file record to a data matrix
Nr=length(temp)/(Nc*3); % number of records we have
data=zeros(Nc,Nr);
for i=1:Nr
    for j=1:Nc
            data(j,i)=(temp((j+(i-1)*Nc)*3)>127)*(256-temp((j+(i-1)*Nc)*3))*(-65536)+... % the last byte (sigh) The time series data is stored in a 24-bit binary complement format(yin wei shi bu ma suo yi zheng fu hao biao shi shang cun zai cha yi)
             (temp((j+(i-1)*Nc)*3)<=127)*(temp((j+(i-1)*Nc)*3))*(65536)+...
            temp((j+(i-1)*Nc)*3-1)*256+...              % the second byte (*256) (temp mei ci zhi neng du yi ge zhi jie er yi ge shu ju bao han san ge zhi jie gu yong  -0 -1 -2 biao shi)
            temp((j+(i-1)*Nc)*3-2);               % the first byte          sheng lue hao biao shi xia hang ji xu   
    end
end
function data = read_timetag(fid,Rstart,Nr,Nch,tsize)
%% a 'smart 'function to read TSn file's tag
% Rstart :the start record
% Rend :the end record
% Nr: sample rate (Hz)
%data :
fseek(fid,(Rstart-1)*(tsize+Nr*Nch*3),-1);
data = fread(fid,tsize,'uint8');




