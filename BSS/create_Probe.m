function [ XProbe ] = create_Probe(NPixel,xpitch,ypitch,probesep,quad)
%   Create probe locations [xloc yloc] x NPixel
%
%   X is along the shank. Y is only taking two values: 0 or 73.22 um.
%
%   [XY] = [0 0] is pixel 3
%   so it actually goes like
%
%   BASE [1 2] [5 6] ============ [253 254]
%   BASE [3 4] [7 8] ============ [255 256]
%   
%   BASE    |
%   BASE    Y
%   BASE    |
%   BASE    ----------------- X -----------------

xseq = [0 1 0 1];
yseq = [1 1 0 0];
Xoff =[];
Yoff = [];


if probesep > 0
    if quad %4 shanks, and separation in the middle 
        cutoff = NPixel/4;
        for nn = 1:cutoff/4
            Xoff(((nn-1)*4+1):(nn*4)) = (xseq+2*(nn-1))*xpitch;
        end

        %repeat 4 times for 4 shanks
        Xoff = [Xoff Xoff Xoff Xoff];

        for nn = 1:cutoff
            Yoff(((nn-1)*4+1):(nn*4)) = yseq*ypitch + floor((4*(nn-1))/cutoff)*probesep;
        end
        
    else
        cutoff = NPixel/4;
        for nn = 1:cutoff/4
            Xoff(((nn-1)*4+1):(nn*4)) = (xseq+2*(nn-1))*xpitch;
        end

        %repeat 2 times for 2 shanks
        Xoff = [Xoff Xoff];

        for nn = 1:cutoff
            Yoff(((nn-1)*4+1):(nn*4)) = yseq*ypitch + floor((2*(nn-1))/cutoff)*probesep;
        end
    end
else   % for single shank probe. no gratings. dualShank2018.
    for nn = 1:NPixel/4
        Xoff(((nn-1)*4+1):(nn*4)) = (xseq+2*(nn-1))*xpitch;
        Yoff(((nn-1)*4+1):(nn*4)) = yseq*ypitch;
    end
end
%space for the large oled pad in the center
%Xoff(Xoff>64*xpitch) = Xoff(Xoff>64*xpitch)+188.04; 
XProbe = [Xoff' Yoff'];
end