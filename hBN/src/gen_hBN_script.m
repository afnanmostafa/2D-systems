%%% Afnan Mostafa
%%% 02/25/2023

%%% ### CODE TO CREATE SINGLE- OR MULTI-LAYERED hBN (AA or AB) ### %%%

%%% This code generates aa- or ab- stacked
%%% hBN. The number of layers is defined by the variable 'layers' and
%%% stacking is defined by the variable 'stacking'
clear

%% parameters
format long;
cc  = 1.45;                 %% interatomic distance
per_arm = 4.35;
per_zig = 2.51;

%% change here
length_sheet = 20.2;        %% in nm with buffer of 0.2
width_sheet = 10.2;         %% in nm with buffer of 0.2

%%
nx = ceil((length_sheet*10)/per_arm);
ny = ceil((width_sheet*10)/per_zig);
cc_height = 3.3;

%% write or show final results
write = true;
show = true;

%% unit cell size
lx = 3*cc;
ly = sqrt(3)*cc;

%% coordinates of the 4 basis atoms in the unit cell
base = [ 0.0 , 0.0 , 0.0 ;
    cc/2 , ly/2 , 0.0 ;
    lx/2 , ly/2 , 0.0 ;
    2*cc , 0.0 , 0.0 ];

%% total number of atoms
N = length(base)*nx*ny;

%% coordinates of the atoms in the layer
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

%% hBN layers, atom types 
atomtype = 2;
stacking = 'ab';
layers = 2;         %% hBN layers
total_atoms = N*layers;

%% show as a figure
if show
    plot(coords(:,1),coords(:,2),'o')
    hold on
    plot(base(:,1),base(:,2),'.r','markersize',20)
    axis equal
end

%% write to file
filename1 = 'hBN_%d_layers_%s_%dnmx%dnm.data';
filename = sprintf(filename1,layers,stacking,floor(length_sheet),floor(width_sheet));

if write
    fid = fopen(filename,'w');
    fprintf(fid,'#hBN %dx%d, a=%g, lx=%g, ly=%g\n',length_sheet,width_sheet,cc,nx,ny);
    fprintf(fid,'%g atoms\n\n',total_atoms);
    fprintf(fid,'%g atom types\n\n',atomtype);
    fprintf(fid,'0 %g xlo xhi\n',lx*nx);
    fprintf(fid,'0 %g ylo yhi\n',ly*ny);
    fprintf(fid,'%g %g zlo zhi\n\n',-2*floor(cc_height*6*layers),2*floor(cc_height*6*layers));
    fprintf(fid,'Masses\n\n');
    
    
    fprintf(fid,'1 14.0067\n');  % N mass
    fprintf(fid,'2 10.811\n');   % B mass
    
    fprintf(fid, '\n');
    fprintf(fid,'Atoms\n\n');
    
    atom_layer = total_atoms/layers;
    p=1;
    t = (1:atomtype);
    type = (repmat(t,1,length(coords)/2))';
    
    for i = 1:layers
        for j = 1:atom_layer
            if strcmp(stacking, 'ab')
                fprintf(fid,'%g %g %g %g %g\n',(j+(p-1)*atom_layer),type(j,1),(coords(j,1)+(cc*abs(1-rem(i,2)))),coords(j,2),coords(j,3)+(cc_height*i));
            else
                fprintf(fid,'%g %g %g %g %g\n',(j+(p-1)*atom_layer),type(j,1),coords(j,1),coords(j,2),coords(j,3)+(cc_height*i));
            end
        end
        p=p+1;
    end
    fclose(fid);
end
