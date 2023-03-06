function [finalfile] = generate_diamane(lengths,widths,stackings,layers)
%   code to generate diamane (2D diamond) structures from bilayer graphene
% input: length array, width array, stacking(aa or ab), no of layers (2)

%% reading the data file and storing spatial coordinate values
for h=1:length(lengths)
    for g=1:length(widths)
        for idx=1:length(stackings)
            
            %% edge case 1:
            if sum(strcmpi(stackings(idx), "ab") || strcmpi(stackings(idx), "aa")) ~= 1
                error('wrong stacking, use either aa or ab')
            end
            %% edge case 2:
            if layers ~= 2
                error('use 2 layers')
            end
            %% read atom data
            length_sheet = lengths(h);
            width_sheet = widths(g);
            [file,len,wid,~] = generate_graphene(length_sheet,width_sheet,stackings(idx),layers);
            [s] = readtextfile(file,5,15,'','#');
            datafile = file;
            filename = file;
            index=s{1};
            atom_type=s{2};
            x = s{3};
            y = s{4};
            z = s{5};
            zz_atom_oneStrip = ceil((wid*10)/2.45);
            ac_atom_oneStrip = ceil((len*10)/4.2);
            totAtom_zz = zz_atom_oneStrip*4;
            totAtom_ac = ac_atom_oneStrip*4;
            top_end = length(x);
            top_strt = length(x) - totAtom_zz;
            bot_strt = 1;
            bot_end = bot_strt + totAtom_zz;
            
            %% %%% selecting the aligned atoms %%% %%
            x_bottom = x(1:length(x)/2);
            x_top = x(length(x)/2+1:end);
            
            y_bottom = y(1:length(y)/2);
            y_top = y(length(y)/2+1:end);
            
            z_bottom = z(1:length(z)/2);
            z_top = z(length(z)/2+1:end);
            
            bottom_aligned_atoms = zeros(length(x_top),1);
            top_aligned_atoms = zeros(length(x_top),1);
            counter1 = 0;
            
            for i = 1:length(x_bottom)
                for j = 1:length(x_top)
                    if (y_bottom(i) == y_top(j)) && (x_bottom(i) == x_top(j))
                        counter1 = counter1 + 1;
                        bottom_aligned_atoms(counter1,1) = index(i);
                        top_aligned_atoms(counter1,1) = index(j);
                    end
                end
            end
            
            bottom_aligned_atoms = bottom_aligned_atoms(1:counter1);
            top_aligned_atoms = top_aligned_atoms(1:counter1);
            top_aligned_atoms_added = top_aligned_atoms + length(x_bottom);
            
            id_cc_b = index((bottom_aligned_atoms));
            id_cc_t = index(top_aligned_atoms_added);
            
            %% %%% indices of atoms for hydrogenation %%% %%
            
            if strcmpi(stackings(idx), "ab")
                %% ################ AB ############ %%
                %% bottom layer %%
                A = index(1:length(x_bottom));
                set = [id_cc_b; (bot_strt:4:bot_end)'];
                idxx = ismember(A, set);
                A(idxx) = [];
                
                %% top layer %%
                A1 = index(length(x_top)+1:end);
                set1 = [id_cc_t; (top_strt:4:top_end)'];
                idx1 = ismember(A1, set1);
                A1(idx1) = [];
                A1 = sort(A1);
                A = sort(A);
                id_cc_b = sort(id_cc_b);
                id_cc_t = sort(id_cc_t);
                
            elseif strcmpi(stackings(idx), "aa")
                %% ################ AA ############ %%
                for u = 1:length(atom_type)
                    if (rem((index(u,1)), 2) == 0)
                        A(u,1) = u;
                    else
                        A1(u,1) = index(u,1);
                    end
                end
                
                A(A==0) = [];
                A1(A1==0) = [];
                ch_bot=A1(A1 < (length(x_bottom)+1));
                ch_top=A1(A1 > length(x_top));
                cc_bot=A(A < (length(x_bottom)+1));
                cc_top=A(A > (length(x_bottom)));
            end
            %% %%% writing LAMMPS input script for C-C atoms to be moved towards each other in order to be bonded %%% %%
            finalfilename = '2D_diamond_%s_%dlayers_%dx%d.data';
            finalfile = sprintf(finalfilename,stackings(idx),layers,length_sheet,width_sheet);
            
            script = 'in_%s_%dx%d.lmp';
            fid2 = sprintf(script, stackings(idx),length_sheet,width_sheet);
            fid2 = fopen(fid2, 'w');
            fprintf(fid2, '#2D Diamond\n\n');
            fprintf(fid2,'dimension 3\n');
            fprintf(fid2,'units metal\n');
            fprintf(fid2,'processors * * 1\n');
            fprintf(fid2,'boundary p p p\n');
            fprintf(fid2,'neighbor 0.3 bin\n');
            fprintf(fid2,'neigh_modify every 1 delay 0 check yes\n\n');
            fprintf(fid2,'#define variables\n');
            fprintf(fid2,'variable T equal 300.0\n');
            fprintf(fid2,'variable V equal lx*ly*6.7  # 3.35 (1st layer) + gap (0) + 3.35 (2nd layer)\n');
            fprintf(fid2,'variable dt equal 0.0001\n\n');
            fprintf(fid2,'#lattice\n');
            fprintf(fid2,'atom_style atomic\n');
            fprintf(fid2,'read_data %s\n\n',finalfile);
            fprintf(fid2,'#define interatomic potentials\n');
            fprintf(fid2,'pair_style			airebo 3 1 0\n'); %LJ on (1), torsion off (0)
            fprintf(fid2,'pair_coeff 			* * CH.airebo C C H\n\n');
            fprintf(fid2,'dump equil all xyz 1000 equilibration.xyz\n');
            fprintf(fid2,'dump_modify equil sort id\n\n');
            fprintf(fid2,'#define thermo variable settings\n');
            fprintf(fid2,'timestep ${dt}\n\n');
            fprintf(fid2,'#### BOND FORMATION ####\n\n');
            fprintf(fid2,'group atom1 id \t');
            
            if strcmpi(stackings(idx), "ab")
                fprintf(fid2, '%d ',set);
                fprintf(fid2,'\n\n');
                fprintf(fid2,'group atom2 id \t');
                fprintf(fid2, '%d ',set1);
                
            elseif strcmpi(stackings(idx), "aa")
                fprintf(fid2, '%d ',cc_bot);
                fprintf(fid2,'\n\n');
                fprintf(fid2,'group atom2 id \t');
                fprintf(fid2, '%d ',cc_top);
            end
            
            fprintf(fid2,'\n\n');
            fprintf(fid2,'displace_atoms atom1 move 0 0 1.15 units box\n');
            fprintf(fid2,'displace_atoms atom2 move 0 0 -1.15 units box\n\n');
            fprintf(fid2,'run 2000\n\n');
            fprintf(fid2,'group bonded union atom1 atom2\n');
            fprintf(fid2,'group nonbonded subtract all bonded\n\n');
            fprintf(fid2,'thermo_style  custom step time temp press pxx pyy pzz pe ke etotal\n');
            fprintf(fid2,'thermo    1000\n\n');
            fprintf(fid2,'# ======= minimization =======\n\n');
            fprintf(fid2,'min_style cg\n');
            fprintf(fid2,'minimize  1e-10 1e-10 200000 400000\n\n');
            fprintf(fid2,'write_data post_mini.data\n\n');
            fprintf(fid2,'# ======= equilibration =======\n\n');
            fprintf(fid2,'fix 2 all recenter INIT INIT INIT\n\n');
            fprintf(fid2,'fix 1 nonbonded npt temp $T $T 0.1 x 0.0 0.0 1.0 y 0.0 0.0 1.0 couple xy\n');
            fprintf(fid2,'run 150000\n');
            fprintf(fid2,'unfix 1\n\n');
            fprintf(fid2,'fix 1 all npt temp $T $T 0.1 x 0.0 0.0 1.0 y 0.0 0.0 1.0 couple xy\n');
            fprintf(fid2,'run 200000\n');
            fprintf(fid2,'unfix 1\n');
            fprintf(fid2,'unfix 2\n\n');
            fprintf(fid2,'write_data post_equil.data\n');
            fclose(fid2);
            
            %% %%% data file including the H atoms (No C atoms in this data file) %%% %%
            output = '2d_diamond_%s_%dx%d_onlyH.data';
            outfile = sprintf(output,stackings(idx),length_sheet,width_sheet);
            fid = fopen(outfile,'w');
            
            if strcmpi(stackings(idx), "ab")
                
                l=1;
                for m=1:length(A)
                    fprintf(fid,'%g 3 %g %g %g\n',length(atom_type)+l,x_bottom(A(m)), y_bottom(A(m)),z_bottom(A(m))-0.9);
                    l=l+1;
                end
                
                clear l m
                l=length(A) +1;
                
                for m=1:length(A1)
                    fprintf(fid,'%g 3 %g %g %g\n',length(atom_type)+l,x(A1(m)), y(A1(m)),z(A1(m))+0.9);
                    l=l+1;
                end
                
                fclose(fid);
                
                total_atoms = length(index)+length(A1)+length(A);
                types = 3;
            elseif strcmpi(stackings(idx), "aa")
                l=1;
                
                for m=1:length(ch_bot)
                    fprintf(fid,'%g 3 %g %g %g\n',length(atom_type)+l,x_bottom(ch_bot(m)), y_bottom(ch_bot(m)),z_bottom(ch_bot(m))-0.9);
                    l=l+1;
                end
                
                clear l m
                l=length(ch_top) +1;
                
                for m=1:length(A1)/2
                    fprintf(fid,'%g 3 %g %g %g\n',length(atom_type)+l,x(ch_top(m)), y(ch_top(m)),z(ch_top(m))+0.9);
                    l=l+1;
                end
                
                fclose(fid);
                total_atoms = (length(index)/2)+length(index);
                types = 3;
                
            end
            
            %% final modification to output file
            onetimerun= true;
            if onetimerun == true
                str5 = sprintf("rm %s '",finalfile);
                command5 = (str5);
                system(command5);
                
                command = 'cat %s %s > %s';
                command = sprintf(command,filename,outfile, finalfile);
                system(command);
                
                command = 'sed -i "s/2 atom types/3 atom types/" %s';
                command = sprintf(command, finalfile);
                system(command);
                
                command = 'sed -i "s/%d atoms/%d atoms/1" %s';
                command = sprintf(command, length(index), total_atoms, finalfile);
                system(command);
                
                command = 'sed -i "13 a 3 1.00784" %s';
                command = sprintf(command, length(index), total_atoms, finalfile);
                system(command);
                
                clearvars -except g h idx stackings lengths widths layers
            end
        end
    end
end

end

