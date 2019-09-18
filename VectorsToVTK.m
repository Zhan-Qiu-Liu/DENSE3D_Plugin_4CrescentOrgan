function [] = VectorsToVTK(Nodes, Scalars, ScalarName, Vectors, VectorName, OutputFileName)

fileID = fopen([OutputFileName,'.vtk'],'w');

fprintf(fileID, '%s\n', '# vtk DataFile Version 4.0');
fprintf(fileID, '%s\n', 'vtk output');
fprintf(fileID, '%s\n', 'ASCII');
fprintf(fileID, '%s\n', 'DATASET POLYDATA');
fprintf(fileID, '%s ' , 'POINTS');
fprintf(fileID, '%d ' , size(Nodes,1));
fprintf(fileID, '%s\n', 'float');

for i=1:size(Nodes,1)
    fprintf(fileID, '%1.6f %1.6f %1.6f \n', Nodes(i,:));
end

fprintf(fileID, '\n%s ' , 'POINT_DATA');
fprintf(fileID, '%d \n' , size(Nodes,1) );

if (size(Scalars, 1) > 0)
    fprintf(fileID, '%s \n' , ['SCALARS   ',ScalarName,'   double']);
    fprintf(fileID, '%s\n', 'LOOKUP_TABLE default');
    for i=1:size(Scalars,1)
        fprintf(fileID, '%1.6f \n', Scalars(i,:));
    end
end

% if (size(Scalars, 1) == 0)
% end

if (size(Vectors, 1) > 0)
%     fprintf(fileID, '\n%s ' , 'POINT_DATA');
    fprintf(fileID, '%s \n' , ['VECTORS   ',VectorName,'   double']);
    for i=1:size(Vectors,1)
        fprintf(fileID, '%1.6f %1.6f %1.6f \n', Vectors(i,:));
    end
end
fclose(fileID);

end


