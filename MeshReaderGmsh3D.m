function [Nv, VX, VY, VZ, K, EToV] = MeshReaderGmsh3D(FileName)

% function [Nv, VX, VY, VZ, K, EToV] = MeshReaderGmsh3D(FileName)
% Purpose  : Read in basic grid information to build grid

file = fopen(FileName, 'rt');
line = fgetl(file);

while(line ~= -1)
    
    if strfind(line,'$Nodes')
        
        line = fgetl(file);
        Nv = str2num(line);
        VX = zeros(1,Nv);
        VY = zeros(1,Nv);
        VZ = zeros(1,Nv);
        
        for i=1:Nv
            line = fgetl(file);
            out = sscanf(line,'%i%f%f%f');
            VX(i) = out(2);
            VY(i) = out(3);
            VZ(i) = out(4);
        end
    end
    
    if strfind(line,'$Elements')
        
        line = fgetl(file);
        K = str2num(line);
        EToV = zeros(K,4);
        
        kTet = 0;
        for i=1:K
            line = fgetl(file);
            out = sscanf(line,'%i%i');
            if out(2) == 4
                kTet = kTet+1;
                out = sscanf(line,'%i%i%i%i%i%i%i%i%i');
                EToV(kTet,1) = out(6);
                EToV(kTet,2) = out(7);
                EToV(kTet,3) = out(8);
                EToV(kTet,4) = out(9);
            end
        end
        
        K = kTet;
        EToV = EToV(1:K,1:4);
    end
    
    line = fgetl(file);
end
