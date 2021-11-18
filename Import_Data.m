function [x,y,theta,index_i, index_j,points,center,R] = Import_Data(process)
    delimiterIn = ',';
    headerlinesIn = 0;
    parameters_route = "Input/Buffer/parameters_" + string(process)+".csv";
    parameters_table = importdata(parameters_route,delimiterIn, headerlinesIn);
    int_index(1) = int32(parameters_table(1,1));
    int_index(2) = int32(parameters_table(1,2));
    center(1) = parameters_table(1,3);
    center(2) = parameters_table(1,4);
    R = parameters_table(1,5);
    supp_points_route = "Input/Buffer/knots_"+string(process)+".csv";
    supp_points_table = importdata(supp_points_route,delimiterIn, headerlinesIn);
    index_j = supp_points_table(:,1);
    theta = supp_points_table(:,2);
return 