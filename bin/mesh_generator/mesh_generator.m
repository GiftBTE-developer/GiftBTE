% 2023/07
% This script is used for generating mesh file for 1D lines, 2D square/rectangle or 3D cubic geometries
% The mesh_control file specifies the dimension, length scale and the number of discrete meshes
% Please put the generated mesh file (inputmesh.txt) in the input folder

%%
close all;clear all;clc;

fid=fopen('mesh_control');
if(fid==-1)
    error('Error: need mesh_control file');
end

while ~feof(fid)
    str=fgetl(fid);
    if (contains(str, 'GeometryDimension'))
        Dimension=str2num(fgetl(fid));
    end
end
clear fid;
%%
% 1D Line geometry
if(Dimension==1)
    fid=fopen('mesh_control');
    while ~feof(fid)
        str=fgetl(fid);
        if (contains(str, 'X_MeshNumber'))
            MeshNumX=str2num(fgetl(fid));
        end
    end
    clear fid;

    %
    fid=fopen('inputmesh.txt','w');
    fprintf(fid, '# COMSOL\n');
    fprintf(fid, '1 # sdim\n');
    fprintf(fid, strcat(num2str(MeshNumX+1),' # number of mesh vertices'));
    fprintf(fid, '\n');
    fprintf(fid,'0 # lowest mesh vertex index\n\n');
    fprintf(fid,'# Mesh vertex coordinates\n');
    for i=1:MeshNumX+1
        temp=(i-1)*(1/MeshNumX);
        fprintf(fid, '%4.16f', temp);
        fprintf(fid, '\n');
    end

    fprintf(fid,'\n2 # number of element types\n\n');
    
    %
    fprintf(fid,'# Type #0\n\n');

    fprintf(fid,'1 # number of vertices per element\n');
    fprintf(fid,'2 # number of elements\n');
    fprintf(fid,'# Elements\n');
    fprintf(fid,'%d\n', 0);
    fprintf(fid,'%d\n', MeshNumX);

    fprintf(fid,'\n2 # number of geometric entity indices\n');
    fprintf(fid,'# Geometric entity indices\n');
    fprintf(fid,'%d\n', 0);
    fprintf(fid,'%d\n', 1);

    %
    fprintf(fid,'\n# Type #1\n\n');

    fprintf(fid,'2 # number of vertices per element\n');
    fprintf(fid,strcat(num2str(MeshNumX),' # number of elements\n'));
    fprintf(fid,'# Elements\n');
    for i=1:MeshNumX
        fprintf(fid, '%d %d\n',i-1, i);
    end
    fprintf(fid, '\n');

    fprintf(fid,strcat(num2str(MeshNumX),' # number of geometric entity indices\n'));
    fprintf(fid, '# Geometric entity indices\n');
    for i=1:MeshNumX
        fprintf(fid, '%d\n',1);
    end

    fclose(fid);

    fprintf('Dimension = ');
    fprintf('%d\n', Dimension);
    fprintf('Total number of discrete meshes = ');
    fprintf('%d\n', MeshNumX);

end

%%
% 2D Square/Rectangle geometry
if(Dimension==2)
    fid=fopen('mesh_control');
    while ~feof(fid)
        str=fgetl(fid);
        if (contains(str, 'X_Length'))
            LX=str2num(fgetl(fid));
        end
        if (contains(str, 'Y_Length'))
            LY=str2num(fgetl(fid));
        end
        if (contains(str, 'X_MeshNumber'))
            MeshNumX=str2num(fgetl(fid));
        end
        if (contains(str, 'Y_MeshNumber'))
            MeshNumY=str2num(fgetl(fid));
        end
    end
    clear fid;
    
    %
    fid=fopen('inputmesh.txt','w');
    fprintf(fid, '# COMSOL\n');
    fprintf(fid, '2 # sdim\n');
    fprintf(fid, strcat(num2str((MeshNumX+1)*(MeshNumY+1)),' # number of mesh vertices'));
    fprintf(fid, '\n');
    fprintf(fid,'0 # lowest mesh vertex index\n\n');
    fprintf(fid,'# Mesh vertex coordinates\n');
    for i=1:MeshNumX+1
        for j=1:MeshNumY+1
            temp1=(i-1)*(1/MeshNumX)*LX;
            temp2=(j-1)*(1/MeshNumY)*LY;
            fprintf(fid, '%4.16f %4.16f', temp1, temp2);
            fprintf(fid, '\n');
        end
    end

    fprintf(fid,'\n3 # number of element types\n\n');

    %
    fprintf(fid,'# Type #0\n\n');
    fprintf(fid,'# Type #1\n\n');

    fprintf(fid,'2 # number of vertices per element\n');
    fprintf(fid,strcat(num2str(MeshNumX*2+MeshNumY*2),' # number of elements\n'));
    fprintf(fid,'# Elements\n');
    for i=1:MeshNumY
        fprintf(fid,'%d %d\n', i-1, i); % nodes in boundary 1
    end
    for i=1:MeshNumX
        fprintf(fid,'%d %d\n', (i-1)*(MeshNumY+1), i*(MeshNumY+1)); % nodes in boundary 2
    end
    for i=1:MeshNumX
        fprintf(fid,'%d %d\n', i*(MeshNumY+1)-1, (i+1)*(MeshNumY+1)-1); % nodes in boundary 3
    end
    for i=1:MeshNumY
        fprintf(fid,'%d %d\n', MeshNumX*(MeshNumY+1)+i-1, MeshNumX*(MeshNumY+1)+i); % nodes in boundary 4
    end

    fprintf(fid, '\n');
    fprintf(fid, strcat(num2str(MeshNumX*2+MeshNumY*2),' # number of geometric entity indices\n'));
    fprintf(fid,'# Geometric entity indices\n');
    for i=1:MeshNumY
        fprintf(fid,'%d\n', 1-1); % boundary 1
    end
    for i=1:MeshNumX
        fprintf(fid,'%d\n', 2-1); % boundary 2
    end
    for i=1:MeshNumX
        fprintf(fid,'%d\n', 3-1); % boundary 3
    end
    for i=1:MeshNumY
        fprintf(fid,'%d\n', 4-1); % boundary 4
    end

    %
    fprintf(fid,'\n# Type #2\n\n');

    fprintf(fid,'4 # number of vertices per element\n');
    fprintf(fid,strcat(num2str(MeshNumX*MeshNumY),' # number of elements\n'));
    fprintf(fid,'# Elements\n');
    for i=1:MeshNumX
        for j=1:MeshNumY
            fprintf(fid, '%d %d %d %d\n', (i-1)*(MeshNumY+1)+j-1, (i-1)*(MeshNumY+1)+j, i*(MeshNumY+1)+j-1, i*(MeshNumY+1)+j);
        end
    end
    fprintf(fid, '\n');

    fprintf(fid,strcat(num2str(MeshNumX*MeshNumY),' # number of geometric entity indices\n'));
    fprintf(fid, '# Geometric entity indices\n');
    for i=1:MeshNumX*MeshNumY
        fprintf(fid, '%d\n', 1);
    end

    fclose(fid);

    fprintf('Dimension = ');
    fprintf('%d\n', Dimension);
    fprintf('Total number of discrete meshes = ');
    fprintf('%d', MeshNumX);
    fprintf(' x ');
    fprintf('%d', MeshNumY);
    fprintf(' = %d\n', MeshNumX*MeshNumY);
end

%%
% 3D Cubic geometry
if(Dimension==3)
    fid=fopen('mesh_control');
    while ~feof(fid)
        str=fgetl(fid);
        if (contains(str, 'X_Length'))
            LX=str2num(fgetl(fid));
        end
        if (contains(str, 'Y_Length'))
            LY=str2num(fgetl(fid));
        end
        if (contains(str, 'Z_Length'))
            LZ=str2num(fgetl(fid));
        end
        if (contains(str, 'X_MeshNumber'))
            MeshNumX=str2num(fgetl(fid));
        end
        if (contains(str, 'Y_MeshNumber'))
            MeshNumY=str2num(fgetl(fid));
        end
        if (contains(str, 'Z_MeshNumber'))
            MeshNumZ=str2num(fgetl(fid));
        end
    end

    %
    fid=fopen('inputmesh.txt','w');
    fprintf(fid, '# COMSOL\n');
    fprintf(fid, '3 # sdim\n');
    fprintf(fid, strcat(num2str((MeshNumX+1)*(MeshNumY+1)*(MeshNumZ+1)),' # number of mesh vertices'));
    fprintf(fid, '\n');
    fprintf(fid,'0 # lowest mesh vertex index\n\n');
    fprintf(fid,'# Mesh vertex coordinates\n');
    for i=1:MeshNumX+1
        for j=1:MeshNumY+1
            for k=1:MeshNumZ+1
                temp1=(i-1)*(1/MeshNumX)*LX;
                temp2=(j-1)*(1/MeshNumY)*LY;
                temp3=(k-1)*(1/MeshNumZ)*LZ;
                fprintf(fid, '%4.16f %4.16f %4.16f', temp1, temp2, temp3);
                fprintf(fid, '\n');
            end
        end
    end

    fprintf(fid,'\n4 # number of element types\n\n');

    %
    fprintf(fid,'# Type #0\n\n');
    fprintf(fid,'# Type #1\n\n');
    fprintf(fid,'# Type #2\n\n');

    fprintf(fid,'4 quad # type name\n');
    fprintf(fid,strcat(num2str(MeshNumX*MeshNumY*2+MeshNumX*MeshNumZ*2+MeshNumZ*MeshNumY*2),' # number of elements'));
    fprintf(fid,'\n# Elements\n');
    for i=1:MeshNumY
        for j=1:MeshNumZ
            fprintf(fid,'%d %d %d %d\n', (i-1)*(MeshNumZ+1)+j-1, (i-1)*(MeshNumZ+1)+j, i*(MeshNumZ+1)+j-1, i*(MeshNumZ+1)+j); % faces in boundary 1
        end
    end
    for i=1:MeshNumX
        for j=1:MeshNumZ
            fprintf(fid,'%d %d %d %d\n', (i-1)*(MeshNumZ+1)*(MeshNumY+1)+j-1, (i-1)*(MeshNumZ+1)*(MeshNumY+1)+j, i*(MeshNumZ+1)*(MeshNumY+1)+j-1, i*(MeshNumZ+1)*(MeshNumY+1)+j); % faces in boundary 2
        end
    end
    for i=1:MeshNumX
        for j=1:MeshNumY
            fprintf(fid,'%d %d %d %d\n', (i-1)*(MeshNumZ+1)*(MeshNumY+1)+(j-1)*(MeshNumZ+1), (i-1)*(MeshNumZ+1)*(MeshNumY+1)+j*(MeshNumZ+1), i*(MeshNumZ+1)*(MeshNumY+1)+(j-1)*(MeshNumZ+1), i*(MeshNumZ+1)*(MeshNumY+1)+j*(MeshNumZ+1)); % faces in boundary 3
        end
    end
    for i=1:MeshNumX
        for j=1:MeshNumY
            fprintf(fid,'%d %d %d %d\n', (i-1)*(MeshNumZ+1)*(MeshNumY+1)+(j-1)*(MeshNumZ+1)+MeshNumZ, (i-1)*(MeshNumZ+1)*(MeshNumY+1)+j*(MeshNumZ+1)+MeshNumZ, i*(MeshNumZ+1)*(MeshNumY+1)+(j-1)*(MeshNumZ+1)+MeshNumZ, i*(MeshNumZ+1)*(MeshNumY+1)+j*(MeshNumZ+1)+MeshNumZ); % faces in boundary 4
        end
    end
    for i=1:MeshNumX
        for j=1:MeshNumZ
            fprintf(fid,'%d %d %d %d\n', (i-1)*(MeshNumZ+1)*(MeshNumY+1)+j-1+(MeshNumZ+1)*MeshNumY, (i-1)*(MeshNumZ+1)*(MeshNumY+1)+j+(MeshNumZ+1)*MeshNumY, i*(MeshNumZ+1)*(MeshNumY+1)+j-1+(MeshNumZ+1)*MeshNumY, i*(MeshNumZ+1)*(MeshNumY+1)+j+(MeshNumZ+1)*MeshNumY); % faces in boundary 5
        end
    end
    for i=1:MeshNumY
        for j=1:MeshNumZ
            fprintf(fid,'%d %d %d %d\n', (i-1)*(MeshNumZ+1)+j-1+(MeshNumZ+1)*(MeshNumY+1)*MeshNumX, (i-1)*(MeshNumZ+1)+j+(MeshNumZ+1)*(MeshNumY+1)*MeshNumX, i*(MeshNumZ+1)+j-1+(MeshNumZ+1)*(MeshNumY+1)*MeshNumX, i*(MeshNumZ+1)+j+(MeshNumZ+1)*(MeshNumY+1)*MeshNumX); % faces in boundary 6
        end
    end

    fprintf(fid, '\n');
    fprintf(fid, strcat(num2str(MeshNumX*MeshNumY*2+MeshNumX*MeshNumZ*2+MeshNumZ*MeshNumY*2),' # number of geometric entity indices\n'));
    fprintf(fid,'# Geometric entity indices\n');
    for i=1:MeshNumY
        for j=1:MeshNumZ
            fprintf(fid,'%d\n', 1-1); % faces in boundary 1
        end
    end
    for i=1:MeshNumX
        for j=1:MeshNumZ
            fprintf(fid,'%d\n', 2-1); % faces in boundary 2
        end
    end
    for i=1:MeshNumX
        for j=1:MeshNumY
            fprintf(fid,'%d\n', 3-1); % faces in boundary 3
        end
    end
    for i=1:MeshNumX
        for j=1:MeshNumY
            fprintf(fid,'%d\n', 4-1); % faces in boundary 4
        end
    end
    for i=1:MeshNumX
        for j=1:MeshNumZ
            fprintf(fid,'%d\n', 5-1); % faces in boundary 5
        end
    end
    for i=1:MeshNumY
        for j=1:MeshNumZ
            fprintf(fid,'%d\n', 6-1); % faces in boundary 6
        end
    end

    %
    fprintf(fid,'\n# Type #3\n\n');

    fprintf(fid,'3 hex # type name\n');
    fprintf(fid,'8 # number of vertices per element\n');
    fprintf(fid,strcat(num2str(MeshNumX*MeshNumY*MeshNumZ),' # number of elements\n'));
    fprintf(fid,'# Elements\n');
    for i=1:MeshNumX
        for j=1:MeshNumY
            for k=1:MeshNumZ
                temp1=(i-1)*(MeshNumZ+1)*(MeshNumY+1)+(j-1)*(MeshNumZ+1)+k-1;
                temp2=i*(MeshNumZ+1)*(MeshNumY+1)+(j-1)*(MeshNumZ+1)+k-1;
                temp3=(i-1)*(MeshNumZ+1)*(MeshNumY+1)+j*(MeshNumZ+1)+k-1;
                temp4=i*(MeshNumZ+1)*(MeshNumY+1)+j*(MeshNumZ+1)+k-1;
                temp5=(i-1)*(MeshNumZ+1)*(MeshNumY+1)+(j-1)*(MeshNumZ+1)+k;
                temp6=i*(MeshNumZ+1)*(MeshNumY+1)+(j-1)*(MeshNumZ+1)+k;
                temp7=(i-1)*(MeshNumZ+1)*(MeshNumY+1)+j*(MeshNumZ+1)+k;
                temp8=i*(MeshNumZ+1)*(MeshNumY+1)+j*(MeshNumZ+1)+k;
                fprintf(fid, '%d %d %d %d %d %d %d %d\n', temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8);
            end
        end
    end
    fprintf(fid, '\n');

    fprintf(fid,strcat(num2str(MeshNumX*MeshNumY*MeshNumZ),' # number of geometric entity indices\n'));
    fprintf(fid, '# Geometric entity indices\n');
    for i=1:MeshNumX*MeshNumY*MeshNumZ
        fprintf(fid, '%d\n', 1);
    end

    clear fid;

    fprintf('Dimension = ');
    fprintf('%d\n', Dimension);
    fprintf('Total number of discrete meshes = ');
    fprintf('%d', MeshNumX);
    fprintf(' x ');
    fprintf('%d', MeshNumY);
    fprintf(' x ');
    fprintf('%d', MeshNumZ);
    fprintf(' = %d\n', MeshNumX*MeshNumY*MeshNumZ);
end

%%
fprintf('\n**********************************************************\n');
fprintf('***     Finish Generating Mesh File: inputmesh.txt     ***\n');
fprintf('***    Please put inputmesh.txt in the input folder    ***\n');
fprintf('**********************************************************\n\n');



