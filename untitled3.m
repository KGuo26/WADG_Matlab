[VZ VY VX] = meshgrid(linspace(-1,1,Ky+1),linspace(-1,1,Ky+1),linspace(-1,1,Kx+1));

% plot(VX,VY,'o')
% text(VX(:)+.1,VY(:),num2str((1:length(VX(:)))'))

sk = 1;
for ez = 1:Kz
    for ey = 1:Ky
        for ex = 1:Kx
            id = @(ex,ey,ez) ex + (ey-1)*(Kx+1) + (ez-1)*(Kx+1)*(Ky+1);
            id1 = id(ex,ey,ez);
            id2 = id(ex+1,ey,ez);
            id3 = id(ex+1,ey+1,ez);
            id4 = id(ex,ey+1,ez);
            id5 = id(ex,ey,ez+1);
            id6 = id(ex+1,ey,ez+1);
            id7 = id(ex+1,ey+1,ez+1);
            id8 = id(ex,ey+1,ez+1);
            VXe = VX([id1 id2 id3 id4 id5 id6 id7 id8]);
            VYe = VY([id1 id2 id3 id4 id5 id6 id7 id8]);
            EToV(6*sk-5,:) = [id1 id2 id3 id6];
            EToV(6*sk-4,:) = [id1 id3 id4 id6];
            EToV(6*sk-3,:) = [id1 id2 id3 id];
            EToV(6*sk-2,:) = [id3 id4 id1];
            EToV(6*sk-1,:) = [id1 id2 id3];
            EToV(6*sk,:)   = [id3 id4 id1];
            sk = sk + 1;                               
%         pause
        end
    end
end