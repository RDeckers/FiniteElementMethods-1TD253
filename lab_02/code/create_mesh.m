function [ p,e,t ] = create_mesh( dx )
 R1 = [3,4, 0,1,1,0, 0,0,1,1]'; %unit square
 dl = decsg(R1); %make geometry
 [p, e, t] = initmesh(dl,'Hmax', dx); %create a mesh with max segment size dx.
end

