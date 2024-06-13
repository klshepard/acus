function [xlocmat,xloc,x_int] = create_random_beads_3D(X,Y,Z,L,xoff,yoff,zoff)
    % Create x-vector, L random fluorophore locations with given intensity.
    % Everything will be stored in micrometer domain
    % only move to rounded indices for display, otherwise reconstruction fails
    xlocmat =[];         %in indices 0 to sizes
    xloc = [];

    Xdist = X(end)-X(1);
    Ydist = Y(end)-Y(1);
    Zdist = Z(end)-Z(1);

    Xsize = length(X);
    Ysize = length(Y);
    Zsize = length(Z);

    Xmid = mean(X);
    Ymid = mean(Y);

    sizes = [Xsize Ysize Zsize];
    x_int   = zeros(prod(sizes),1);
    % fprintf('We place %d fluorophores at uniformly random locations with mean middle of volume\n', mx);
    ii = 1;
    while size(xlocmat,1) < L

        xyz = rand(size(sizes)).*[Xdist Ydist Zdist]+[xoff yoff zoff];
        if ( (xyz(1)>=X(1)) && (xyz(1)<=X(end)) && (xyz(2)>=Y(1)) && (xyz(2)<=Y(end)) && (xyz(3)>=Z(1)) && (xyz(3)<=Z(end)) )
        %if ((xyz(1)>(Xmid-xoff-xyz(3)*tan(deadangle))) && (xyz(1)<(Xmid+xoff+xyz(3)*tan(deadangle))) && (xyz(2)>(Ymid-yoff-xyz(3)*tan(deadangle))) && (xyz(2)<(Ymid+yoff+xyz(3)*tan(deadangle))) && (xyz(3)>=Z(1)) && (xyz(3)<=Z(end)))
            
            xyzround = [X(findnearest(xyz(1),X)) Y(findnearest(xyz(2),Y)) Z(findnearest(xyz(3),Z))]; 
            xlocmat = [xlocmat; xyzround];
            xloc(ii) = sub2ind(sizes,findnearest(xyz(1),X),findnearest(xyz(2),Y),findnearest(xyz(3),Z));
            ii = ii + 1;
        else
            fprintf('Bead is out of defined volume. Skipping this and creating a new one.');
            continue
        end
    end    

    x_int(xloc) = 1;

    mx = sum(x_int);
    fprintf('Total amount of fluorophores created is %i.\n', mx)
end

