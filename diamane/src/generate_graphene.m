function [filename,len,wid,stacking] = generate_graphene(len,wid,stacking,layers)
%   creates single/multi-layer graphene
%   This code generates aa- or ab- stacked
%   graphene. input = length (nm), width (nm), stacking (ab or aa), layers
%   (1, 2,....)

format long;

%% initialize
cc  = 1.4;                          %% C-C sp2 bond distance
per_arm = 4.26;                     %% C-C distance in armchair direction
per_zig = 2.46;                     %% C-C distance in zigzag direction

buffer = 0.2;
len_w_buffer = len + buffer;             %% add buffer of 0.2
wid_w_buffer = wid + buffer;             %% add buffer of 0.2

nx = ceil((len_w_buffer*10)/per_arm);    %% repetitions in the x direction
ny = ceil((wid_w_buffer*10)/per_zig);    %% repetitions in the y direction
             
cc_height = 3.35;                        %% c-c height in z direction (A)

write = true;
show = false;

%% Size of the unit cell
lx = 3*cc;
ly = sqrt(3)*cc;

%% Coordinates of the 4 basis atoms in the unit cell
base = [ 0.0 , 0.0 , 0.0 ;
    cc/2 , ly/2 , 0.0 ;
    lx/2 , ly/2 , 0.0 ;
    2*cc , 0.0 , 0.0 ];

%% Total number of atoms
N = length(base)*nx*ny;

%% Calculate the coordinates of the atoms in the layer
coords = zeros(N,3);
id = 0;

for ix=1:nx
    for iy=1:ny
        for iatom=1:length(base)
            id = id + 1;
            coords(id,:) = base(iatom,:)+[(ix-1)*lx,(iy-1)*ly,0];
        end
    end
end

total_atoms = N*layers;
filename1 = 'graphene_%d_layers_%s_%dx%d.data';
filename = sprintf(filename1,layers,stacking,floor(len),floor(wid));

%% showing the structure as a figure
if show
    hold on
    plot(coords(:,1),coords(:,2),'o')
    plot(base(:,1),base(:,2),'.r','markersize',20)
    axis equal
end

%% writing output to file
if write
    fid = fopen(filename,'w');
    fprintf(fid,'#graphene %dx%d, a=%g, lx=%g, ly=%g\n',len,wid,cc,nx,ny);
    fprintf(fid,'%g atoms\n\n',total_atoms);
    fprintf(fid,'%g atom types\n\n',layers);
    fprintf(fid,'0 %g xlo xhi\n',lx*nx);
    fprintf(fid,'0 %g ylo yhi\n',ly*ny);
    fprintf(fid,'%g %g zlo zhi\n\n',-2*floor(cc_height*6*layers),2*floor(cc_height*6*layers));
    fprintf(fid,'Masses\n\n');
    
    for u=1:layers
        fprintf(fid,'%g 12.0107\n',u);
    end
    
    fprintf(fid, '\n');
    fprintf(fid,'Atoms\n\n');
    
    atom_layer = total_atoms/layers;
    p=1;
    for i = 1:layers
        for j = 1:atom_layer
            if strcmp(stacking, 'ab')
                fprintf(fid,'%g %g %g %g %g\n',(j+(p-1)*atom_layer),i,(coords(j,1)+(cc*abs(1-rem(i,2)))),coords(j,2),coords(j,3)+(cc_height*i));
            else
                fprintf(fid,'%g %g %g %g %g\n',(j+(p-1)*atom_layer),i,coords(j,1),coords(j,2),coords(j,3)+(cc_height*i));
            end
        end
        p=p+1;
    end
    fclose(fid);
end

end

