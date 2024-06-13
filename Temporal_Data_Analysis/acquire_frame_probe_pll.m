function [ N_SPAD, nframes , lambdam, lambdav,fullData] = acquire_frame_probe_pll( outfile, exefile, flashfile, datasize, ignoreframes, clkdiv,spadduty, spadphase, framelength, hotpix,quad,h)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%Acquire data from optical probe
if ~isempty(h)
    h.SC_Enable(0);
end
    system([exefile ...
                 ' -r 1' ...
                 ' -rstc 0007' ... %reset all but fsm
                 ' -f 1' ...
                 ' -ff ' flashfile ...
                 ' -get 1' ...
                 ' -ds ' num2str(datasize) ...
                 ' -o ' outfile ...
                 ' -ep00 ' num2str(hex2dec('10130')) ... %[enable_FSM enable_datapath dont_reset] [MSB LSB] [0150]
                 ' -ep02 ' num2str(framelength) ...
                 ' -clkdiv ' num2str(clkdiv) ... % 
                 ' -rpp 1' ... % reprogram of PLL
                 ' -c0d ' num2str(spadduty*1000) ... 
                 ' -c0p ' num2str(spadphase*1000) ...
                 ' -ns 1']); %... 

if ~isempty(h)
    h.SC_Disable(0);
end
    
%     pause(0.5)
%     global keysight_gpib;
%     I(1) = str2double(query(keysight_gpib,':MEASure:CURRent? (@1)'));
%     I(2) = str2double(query(keysight_gpib,':MEASure:CURRent? (@2)'));
    
    %Open intensity data into matlab 2D-matrix
    fileID = fopen([outfile '_chunk_0']);
    A = fread(fileID, 'uint16=>double');
    fclose(fileID);
    B = dec2bin(A);
	%dual probe:
		%PROBE 1 | ADDR 6 | PIX 2 | '0' | DATA 5
	%quad probe:
		%PROBE 2 | ADDR 6 | PIX 2 | DATA 5
	
    if exist([outfile '_chunk_0'],'file')==2
        delete([outfile '_chunk_0']);            
    end   
    
    if quad
        NPixel = 1024;
    else
        NPixel = 512;
    end
    [ N_SPAD, nframes, lambdam, lambdav, fullData] = ASPdataValidationv2(B,ignoreframes);
    
    %hot pixels
    
    
    
    %Display data
    figure(991)
    set(gcf, 'Position', [50 50 700 700])
    subplot(4,1,1)
    scatter(repmat((1:NPixel)',nframes,1),fullData(:),10)
    axis([1 NPixel 0 63])
    title('Scatterplot of All 5-bit Captures')
    
    subplot(4,1,2)
    scatter(repmat((1:NPixel)',nframes,1),fullData(:),10)
    axis([1 NPixel 0 10])
    title('Scatterplot of All 5-bit Captures')
    
    subplot(4,1,3)
    bar(1:NPixel,N_SPAD,2,'r')
    set(gca,'YScale','log')
    xlim([1 NPixel])
    title('Accumulation, log')
    
%     lambdaplot = lambdam;
%     lambdaplot(find(hotpix)) = 0;
    
%     subplot(5,1,4)
%     bar(1:NPixel,lambdaplot,2,'g')
%     axis([1 NPixel 0 31])
%     title('Poisson lambda parameter')
    
    subplot(4,1,4)
    plot(1:NPixel,N_SPAD)
    xlim([1 NPixel])
    title('Accumulation, linear')  

end

