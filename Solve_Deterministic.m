function Solve_Deterministic(process)
    [x,y,theta,index_i,index_j,points,center,R] = Import_Data(process);
    G = zeros(length(index_i),length(index_j));
    Gi = zeros(length(index_i),length(index_j));
    Gj = zeros(length(index_i),length(index_j));
    for j = 1:length(theta)
        G(:,j) = Pseudospectral_G(points,length(theta), theta(j), center, R);
        Gi(:,j) = index_i;
        Gj(:,j) = ones(size(index_i))*(index_j(j));
        disp(j)
    end
    B = Pseudospectral_B(points, center, R);
    file_output_G = append(append('Matlab_buffer/G_',string(process)),'.txt');
    fileID = fopen(file_output_G,'w');
    fprintf(fileID,'%1.10f\n',G(:));
    fclose(fileID);
    file_output_G = append(append('Matlab_buffer/Gi_',string(process)),'.txt');
    fileID = fopen(file_output_G,'w');
    fprintf(fileID,'%d\n',Gi(:));
    fclose(fileID);
    file_output_G = append(append('Matlab_buffer/Gj_',string(process)),'.txt');
    fileID = fopen(file_output_G,'w');
    fprintf(fileID,'%d\n',Gj(:));
    fclose(fileID);
    file_output_B = append(append('Matlab_buffer/B_',string(process)),'.txt');
    fileID = fopen(file_output_B,'w');
    fprintf(fileID,'%1.10f\n',B);
    file_output_B = append(append('Matlab_buffer/Bi_',string(process)),'.txt');
    fileID = fopen(file_output_B,'w');
    fprintf(fileID,'%d\n',index_i);
    fclose(fileID);
return