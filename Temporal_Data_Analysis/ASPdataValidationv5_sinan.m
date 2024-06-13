function [ N_SPAD, nframes, lambdam, lambdav,fullData,PARITYout] = ASPdataValidationv5_sinan(B,quadtype, NPixel, shankNumber)
%Validates the data
%   Typecasts the data, removes first few frames
%   Detects any faulty data
%   Sums the data into frames
%   Calculates the poisson variance and mean
% V3 also works for quadprobe2018
% V4 also works for quadprobe2020
% V5 removed numPattern, each Pattern is processed separately

	%dual probe:
		%PROBE 1 | ADDR 6 | PIX 2 | '0' | DATA 6
	%quad probe 2018:
		%PROBE 2 | ADDR 6 | PIX 2 | DATA 6
    %quad probe 2020:
        %PROBE 2 | ADDR 6 | PIX 2 | PARITY | DATA 5
        


%%%%%remove first few frames, find onset.
if quadtype > 0 %all quadprobes
    
    %NPixel = 1024;
    %NPixel = 128;
    
    % if we only measure from 128 pixels for in vivo tests, adjust the expected first pixel address based on which shank is inserted.
    % if we measure all 1024 pixels, start from 0.
    if NPixel == 128
        if shankNumber == 0
            firstPixAdd = 128;
        elseif shankNumber == 1
            firstPixAdd = 384;
        elseif shankNumber == 2
            firstPixAdd = 640;
        elseif shankNumber == 3
            firstPixAdd = 896;
        end
        
    elseif NPixel == 1024
        firstPixAdd = 0;
    end
    
%     cutRows = NPixel*ignoreframes;
    
%    addr = bin2dec(B(:,1:10));
%    data = bin2dec(B(:,end-5:end));
    
    addr = floor(B(1:end)/2^6);   %take 10 highest bits
    data = mod(B(1:end),2^6);     %take 6 lowest bits
    
    if (quadtype == 1)  %2018 quad
        %convert data
% %         parity = ones(size(data)); %take 1 highest bit
% %         
% %         if (mod(addr(1),2)==0)
% %             onset = NPixel-addr(1)/2+1;       %this finds the right onset usually
% %         else
% %             onset = NPixel/2-(addr(1)-1)/2+1; %this finds the right onset usually
% %         end
        
    else %2020 quad
        parity = floor(data/2^5); %take 1 highest bit
        data = mod(data,2^5); %take 5 lowest bits
    end    
else
%   NPixel = 512;
% %     cutRows = NPixel*ignoreframes;
%     addr = floor(B(1:end)/2^7);   %take 9 highest bits
%     data = mod(B(1:end),2^6);     %take 6 lowest bits
% 	parity = ones(size(data)); %take 1 highest bit
% 	onset = NPixel-addr(1);       %this finds the right onset usually
end    

%%%%%extract intact windows, find offset
if quadtype == 1 %the 2018 quadprobe goes addr = [0 2 4 ... 1022 -> 1 3 5 1023]
%     while true
%         if addr(onset)==0 && addr(onset+NPixel/2)==1 && addr(onset+1)==2
%             break
%         end
%         onset = onset + 1;
%     end
% 
%     if (mod(addr(end),2)==0)
%         offset = length(addr)-addr(end)/2-1;
%     else
%         offset = length(addr)-(addr(end)-1)/2-1-NPixel/2;
%     end

else %the dualprobe and 2020 quadprobe goes normally addr = [0 1 2 3 ... 511]

    jj = 0+firstPixAdd+NPixel-addr(1); % 0+NPixel-addr(1);
    if jj==NPixel
        jj = 1;
        onset = jj;
    else
        onset = jj; 
    end
    
    while true
        %addr(jj)==0 && addr(jj+1)==1 && addr(jj+2)==2
        %addr(jj)==127 && addr(jj+1)==128 && addr(jj+2)==129
        %addr(jj)==383 && addr(jj+1)==384 && addr(jj+2)==385
        %addr(jj)==639 && addr(jj+1)==640 && addr(jj+2)==641
        %addr(jj)==895 && addr(jj+1)==896 && addr(jj+2)==897
        
        if addr(jj)==firstPixAdd && addr(jj+1)==firstPixAdd+1 && addr(jj+2)==firstPixAdd+2
            onset = jj;
            break
        end
        jj = jj + 1;
    end
    offset = length(addr)-(addr(end)-firstPixAdd)-1;
end

safeAddr = addr(onset:offset);
safeData = data(onset:offset);
safePARITY = parity(onset:offset);

%% jaebin and adriaan's method to determine the bitflips and extract intact frames.

%%%%%now the array starts at 0 and ends at NPixel-1 (511 or 1023)
    
    
    %check for parity, not functional yet
%     safePARITY = (safePARITY == parityqp(dec2bin(safeData,5)));
%%%%%find locations where an address bit had flipped
    %bitflip sometimes coincides with the -511 or -1023 address shift
% % % % % % % % % % % % % % % deltaAddr = (safeAddr(2:end)-safeAddr(1:end-1));
% % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % if (quadtype == 1) %2018 QP
% % % % % % % % % % % % % % %     bitfliploc = find((deltaAddr~=(-1021)) & (deltaAddr~=(-1023)) & (deltaAddr~=(2)));
% % % % % % % % % % % % % % % elseif (quadtype == 2) %2020QP
% % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % %     bitfliploc = find((deltaAddr~=(-127)) & (deltaAddr~=(1)));
% % % % % % % % % % % % % % %     %bitfliploc = find((deltaAddr~=(-1023)) & (deltaAddr~=(1)));
% % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % else %dualprobe
% % % % % % % % % % % % % % %     bitfliploc = find((deltaAddr~=(-511))  & (deltaAddr~=(1)));
% % % % % % % % % % % % % % % end
% % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % %deltaBitfliploc = (bitfliploc(2:end)-bitfliploc(1:end-1));
% % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % safeAddrOriginal = safeAddr;
% % % % % % % % % % % % % % % safeDataOriginal = safeData;
% % % % % % % % % % % % % % % safePARITYOriginal = safePARITY;
% % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % for ii = 1:length(bitfliploc)
% % % % % % % % % % % % % % %     bitflipaddr = safeAddrOriginal(bitfliploc(ii));
% % % % % % % % % % % % % % %     bitflipmagn = deltaAddr(bitfliploc(ii))-1;
% % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % %     %%%%%remove that entire frame
% % % % % % % % % % % % % % %     if (quadtype == 1) %2018 quad
% % % % % % % % % % % % % % % %         if (bitflipmagn > 0) && (mod(bitflipaddr,2)==0) %skips frames up at an even frame
% % % % % % % % % % % % % % % %             safeData(bitfliploc(ii)-bitflipaddr/2:bitfliploc(ii)-bitflipaddr/2+NPixel-bitflipmagn) = NaN;
% % % % % % % % % % % % % % % %             safeAddr(bitfliploc(ii)-bitflipaddr/2:bitfliploc(ii)-bitflipaddr/2+NPixel-bitflipmagn) = NaN;
% % % % % % % % % % % % % % % %             safePARITY(bitfliploc(ii)-bitflipaddr/2:bitfliploc(ii)-bitflipaddr/2+NPixel-bitflipmagn) = NaN;
% % % % % % % % % % % % % % % %         elseif (bitflipmagn > 0) && (mod(bitflipaddr,2)==1) %skips frames up at an odd frame
% % % % % % % % % % % % % % % %             safeData(bitfliploc(ii)-(bitflipaddr-1)/2-NPixel/2:bitfliploc(ii)-(bitflipaddr-1)/2+NPixel/2-(bitflipmagn-1)/2-1) = NaN;
% % % % % % % % % % % % % % % %             safeAddr(bitfliploc(ii)-(bitflipaddr-1)/2-NPixel/2:bitfliploc(ii)-(bitflipaddr-1)/2+NPixel/2-(bitflipmagn-1)/2-1) = NaN;
% % % % % % % % % % % % % % % %             safePARITY(bitfliploc(ii)-(bitflipaddr-1)/2-NPixel/2:bitfliploc(ii)-(bitflipaddr-1)/2+NPixel/2-(bitflipmagn-1)/2-1) = NaN;
% % % % % % % % % % % % % % % %         else    %skips frames down
% % % % % % % % % % % % % % % %             safeData(bitfliploc(ii)+1:bitfliploc(ii)-bitflipmagn) = NaN;
% % % % % % % % % % % % % % % %             safeAddr(bitfliploc(ii)+1:bitfliploc(ii)-bitflipmagn) = NaN;
% % % % % % % % % % % % % % % %             safePARITY(bitfliploc(ii)+1:bitfliploc(ii)-bitflipmagn) = NaN;
% % % % % % % % % % % % % % %              
% % % % % % % % % % % % % % %             %not implemented
% % % % % % % % % % % % % % % %         end            
% % % % % % % % % % % % % % %     else %dualprobe or 2020QP
% % % % % % % % % % % % % % %         if bitflipmagn > 0 %skips frames up 
% % % % % % % % % % % % % % %             safeData(bitfliploc(ii)-bitflipaddr+128:bitfliploc(ii)-bitflipaddr+128-bitflipmagn+(NPixel-1)) = NaN;
% % % % % % % % % % % % % % %             safeAddr(bitfliploc(ii)-bitflipaddr+128:bitfliploc(ii)-bitflipaddr+128-bitflipmagn+(NPixel-1)) = NaN;
% % % % % % % % % % % % % % %             safePARITY(bitfliploc(ii)-bitflipaddr+128:bitfliploc(ii)-bitflipaddr+128-bitflipmagn+(NPixel-1)) = NaN;
% % % % % % % % % % % % % % %         else    %skips frames down
% % % % % % % % % % % % % % %             safeData(bitfliploc(ii)+1:bitfliploc(ii)-bitflipmagn) = NaN;
% % % % % % % % % % % % % % %             safeAddr(bitfliploc(ii)+1:bitfliploc(ii)-bitflipmagn) = NaN;
% % % % % % % % % % % % % % %             safePARITY(bitfliploc(ii)+1:bitfliploc(ii)-bitflipmagn) = NaN;
% % % % % % % % % % % % % % %         end
% % % % % % % % % % % % % % %     end
% % % % % % % % % % % % % % % end
% % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % safeData(isnan(safeData)) = [];
% % % % % % % % % % % % % % % safeAddr(isnan(safeAddr)) = [];
% % % % % % % % % % % % % % % safePARITY(isnan(safePARITY)) = [];

%% there are bitflips in the safeAddr. Extract the intact frames that count continuously from 0 to 1023. Or 128 to 255.

if NPixel==128
    beginFlags = find(safeAddr == firstPixAdd);
    endFlags = find(safeAddr == firstPixAdd+127);
elseif NPixel==1024
    beginFlags = find(safeAddr == 0);
    endFlags = find(safeAddr == 1023);
end

safeAddrCleaned = [];
safeDataCleaned = [];
safePARITYCleaned = [];

for bb = 1:length(beginFlags)
    remainingEndFlags = endFlags(endFlags > beginFlags(bb));
    
    if remainingEndFlags(1) - beginFlags(bb) == NPixel-1 % 1023 or 127 
        safeAddrCleaned = [safeAddrCleaned, safeAddr(beginFlags(bb):remainingEndFlags(1))];
        safeDataCleaned = [safeDataCleaned, safeData(beginFlags(bb):remainingEndFlags(1))];
        safePARITYCleaned = [safePARITYCleaned, safePARITY(beginFlags(bb):remainingEndFlags(1))];
    else
        continue
    end
end

%prepare outputs
fullData = squeeze(reshape(safeDataCleaned,NPixel,[]));
PARITYout = squeeze(reshape(safePARITYCleaned,NPixel,[]));
N_SPAD = squeeze(sum(fullData,2,'omitnan'));

%calculate the statistics for all pixels
nframes = size(fullData,2);
lambdam = N_SPAD/nframes;
lambdav = squeeze(var(fullData,0,2,'omitnan'));

end

