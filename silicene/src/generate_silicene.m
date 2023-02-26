function [filename] = generate_silicene(len,wid,stacking,layers,isBuckled)
%   generates silicene (buckled/planar)
%   takes length, width, stacking, layers, and boolean 0 = no buckled
%   (planar) and 1 = buckled as inputs

%% parameters
format long;
SiSi  = 2.28;                 %% interatomic distance

%% unit cell size
lx = 3*SiSi;
ly = sqrt(3)*SiSi;

%% change here
length_sheet = len + 0.2;        %% in nm with buffer of 0.2
width_sheet = wid + 0.2;         %% in nm with buffer of 0.2

%% buckling or not
if isBuckled == 1
    delta = 0.46;
    Buckle = 'Buckle';
    new_Si_Si = sqrt((SiSi)^2+(delta)^2);
else
    delta = 0.0;
    Buckle = 'NoBuckle';
end
%%
nx = ceil((length_sheet*10)/lx);
ny = ceil((width_sheet*10)/ly);
Si_Si_height = 4.0;

%% write or show final results
write = true;
show = true;

%% coordinates of the 4 basis atoms in the unit cell
base = [ 0.0 , 0.0 , 0.0 ;
    SiSi/2 , ly/2 , delta ;
    lx/2 , ly/2 , 0.0 ;
    2*SiSi , 0.0 , delta ];

%% total number of atoms
N = length(base)*nx*ny;

%% coordinates of the atoms
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

%% show as a figure
if show
    plot3(coords(:,1),coords(:,2),coords(:,3),'o')
    hold on
    plot(base(:,1),base(:,2),'.r','markersize',20)
    axis equal
end

%% write to file
if layers > 1
    filename1 = 'silicene_%d_lay_%s_%dnmx%dnm_%s.data';
    filename = sprintf(filename1,layers,stacking,floor(length_sheet),floor(width_sheet),Buckle);
else
    filename1 = 'silicene_%d_lay_%dnmx%dnm_%s.data';
    filename = sprintf(filename1,layers,floor(length_sheet),floor(width_sheet),Buckle);
end

if write
    fid = fopen(filename,'w');
    if isBuckled == 1
        fprintf(fid,'#Silicene %dnmx%dnm, a=%g(%g), lx=%g, ly=%g\n\n',length_sheet,width_sheet,SiSi,new_Si_Si,nx,ny);
    else
        fprintf(fid,'#Silicene %dnm x %dnm, a=%g, lx=%g, ly=%g\n\n',length_sheet,width_sheet,SiSi,nx,ny);
    end
    fprintf(fid,'%g atoms\n\n',total_atoms);
    fprintf(fid,'%g atom types\n\n',layers);
    fprintf(fid,'0 %g xlo xhi\n',lx*nx);
    fprintf(fid,'0 %g ylo yhi\n',ly*ny);
    fprintf(fid,'%g %g zlo zhi\n\n',-2*floor(Si_Si_height*6*layers),2*floor(Si_Si_height*6*layers));
    fprintf(fid,'Masses\n\n');
    
    for u=1:layers
        fprintf(fid,'%g 28.0855\n',u);
    end
    
    fprintf(fid, '\n');
    fprintf(fid,'Atoms\n\n');
    
    atom_layer = total_atoms/layers;
    p=1;
    for i = 1:layers
        for j = 1:atom_layer
            if strcmp(stacking, 'ab')
                fprintf(fid,'%g %g %g %g %g\n',(j+(p-1)*atom_layer),i,(coords(j,1)+(SiSi*abs(1-rem(i,2)))),coords(j,2),coords(j,3)+(Si_Si_height*i));
            else
                fprintf(fid,'%g %g %g %g %g\n',(j+(p-1)*atom_layer),i,coords(j,1),coords(j,2),coords(j,3)+(Si_Si_height*i));
            end
        end
        p=p+1;
    end
    fclose(fid);
end

end

