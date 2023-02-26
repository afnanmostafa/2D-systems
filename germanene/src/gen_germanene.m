function [filename] = gen_germanene(len,wid,stacking,layers,isBuckled)
%   generates germanene (buckled/planar)
%   takes length, width, stacking, layers, and boolean 0 = no buckled
%   (planar) and 1 = buckled as inputs

format long;
length_sheet = len + 0.2;
width_sheet = wid + 0.2;

if isBuckled == 1
    if strcmp(stacking, 'ab')
        ge_ge  = 2.41;                  %% interatomic distance
        gege_height = 3.34;             %% final len 2.68 = 3.34 - delta
    elseif strcmp(stacking, 'aa')
        ge_ge  = 2.47;                  %% interatomic distance
        gege_height = 2.56;
    end
else
    if strcmp(stacking, 'ab')
        ge_ge  = 2.35;                  %% interatomic distance
        gege_height = 2.68;             %% final len 2.68
    elseif strcmp(stacking, 'aa')
        ge_ge  = 2.56;                  %% interatomic distance
        gege_height = 2.56;
    end    
end

%% unit cell size
lx = 3*ge_ge;
ly = sqrt(3)*ge_ge;

%% buckling or not
if isBuckled == 1
    delta = 0.66;
    Buckle = 'Buckle';
    new_ge_ge = sqrt((ge_ge)^2+(delta)^2);
else
    delta = 0.0;
    Buckle = 'NoBuckle';
end

%% repetitions in x and y directions
nx = ceil((length_sheet*10)/lx);
ny = ceil((width_sheet*10)/ly);

%% write or show final results
write = true;
show = true;

%% coordinates of the 4 basis atoms in the unit cell
base = [ 0.0 , 0.0 , delta ;
    ge_ge/2 , ly/2 , 0 ;
    lx/2 , ly/2 , delta ;
    2*ge_ge , 0.0 , 0 ];

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
    filename1 = 'germanene_%d_lay_%s_%dnmx%dnm_%s.data';
    filename = sprintf(filename1,layers,stacking,floor(length_sheet),floor(width_sheet),Buckle);
else
    filename1 = 'germanene_%d_lay_%dnmx%dnm_%s.data';
    filename = sprintf(filename1,layers,floor(length_sheet),floor(width_sheet),Buckle);
end

if write
    fid = fopen(filename,'w');
    if isBuckled == 1
        fprintf(fid,'#germanene %dx%d, a=%g(%g), lx=%g, ly=%g\n\n',floor(length_sheet),floor(width_sheet),ge_ge,new_ge_ge,nx,ny);
    else
        fprintf(fid,'#germanene %dx%d, a=%g, lx=%g, ly=%g\n\n',floor(length_sheet),floor(width_sheet),ge_ge,nx,ny);
    end
    
    fprintf(fid,'%g atoms\n\n',total_atoms);
    fprintf(fid,'%g atom types\n\n',layers);
    fprintf(fid,'0 %g xlo xhi\n',lx*nx);
    fprintf(fid,'0 %g ylo yhi\n',ly*ny);
    fprintf(fid,'%g %g zlo zhi\n\n',-2*floor(gege_height*6*layers),2*floor(gege_height*6*layers));
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
                fprintf(fid,'%g %g %g %g %g\n',(j+(p-1)*atom_layer),i,(coords(j,1)+(ge_ge*abs(1-rem(i,2)))),coords(j,2),coords(j,3)+(gege_height*i));
            else
                fprintf(fid,'%g %g %g %g %g\n',(j+(p-1)*atom_layer),i,coords(j,1),coords(j,2),coords(j,3)+(gege_height*i));
            end
        end
        p=p+1;
    end
    fclose(fid);
end

end

