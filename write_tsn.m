function fname = write_tsn(timetag,data,firstscan_R,lastscan_R,tsfile)
%% function :write binary file ,in the TSn type
% fname : the path and filename of TS file you had writed 
%filename:  filename you want to write
%timetag: time series info for each record,a matrix:rows means recordes
%data : time series data: a matrix:cols means channels
%Nch =timetag(13,1);
% firstscan_R = 1;
%lastscan_R = the sample rate;
%if tsfile is null, write a single ts file,if tsfile has alreadly
%existe,add the data to the end of the tsfile;
narginchk(4,5);
if nargin < 5
    SN = timetag(1,10)*256+timetag(1,9);
    %tsize = timetag(14,1);
    Nr = timetag(1,11)+timetag(1,12)*256;
   if Nr == 2400
        Tsn = 3;
    elseif Nr == 150;
        Tsn =4;
    else Tsn = 5;%if we are dealing with AMT data this sentence should be changed
    end
    surfix = strcat('*.TS',num2str(Tsn));
    defname = num2str(SN);
    [filename,path,ind] = uiputfile(surfix,'save TS binary file',defname);
    if ind == 0
        return
    else
        [tdata,ddata] = choose_data(timetag,data,firstscan_R,lastscan_R,Nr);
        Nred = length(tdata);
        fname = strcat(path,filename);
        fid =  fopen(fname,'w+','ieee-le');
        fseek(fid,0,-1);
        for i = 1:Nred
            fwrite(fid,tdata(i,:),'uint8');
            fwrite(fid,ddata((1:Nr)+(i-1)*Nr,:)','bit24');
        end
        fclose(fid);
    end
else
    Nr = timetag(1,11)+timetag(1,12)*256;
    [tdata,ddata] = choose_data(timetag,data,firstscan_R,lastscan_R,Nr);
    %NNr = size(data,2);
    Nred = length(tdata);
    fid =  fopen(tsfile,'a+','ieee-le');
    %fseek(fid,0,-1);
    for i = 1:Nred
        %fseek(fid,3*Nr*Nch*(i-1),0);
        fwrite(fid,tdata(i,:),'uint8');
        %fseek(fid,tsize*i,0);
        fwrite(fid,ddata((1:Nr)+(i-1)*Nr,:)','bit24');
    end
    fclose(fid);
end


function [tdata,ddata] = choose_data(timetag,data,firstscan_R,lastscan_R,Nr)
%% function : choose data to write  make sure to write entire record data
% timetag:tag data
% data :time series data
% firstscan_R: position of the first scan in the first record
% lastscan_R:position of the last scan in the last record
% tdata: tag data
% ddata: time series data
%% defilition of vars
% Nr:sample rate (Hz)
% Nch = timetag(1,13);
% SN = timetag(1,10)*256+timetag(1,9);
% tsize = timetag(1,14);
% Tsn = timetag(1,13);
% NNr = size(data,2);
% Nred = size(timetag,1);
%% code begining.....
%Nr = timetag(1,12)*256+timetag(1,11);
if  firstscan_R == 1 && lastscan_R ==Nr%%each record has Nr scans;
    tdata =timetag;
    ddata = data;
elseif  firstscan_R == 1 && lastscan_R ~=Nr%%the first record and the last one has not Nr scans
    tdata = timetag(1:end-1,:);
    ddata = data(1:end-lastscan_R,:);
elseif  firstscan_R ~= 1 && lastscan_R ==Nr
    tdata = timetag(2,:);
    ddata = data(Nr-firstscan_R+2:end,:);
else
    tdata = timetag(2:end-1,:);
    ddata = data(Nr-firstscan_R+2:end-lastscan_R,:);
end
