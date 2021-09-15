function Plot_interpolator(inter_file, sten_file)
    inter_table = readtable(inter_file);
    x_i = [];
    y_i = [];
    sol_PDD_i = [];
    for i = 1:10
        x_i(:,i) = inter_table{:,(i-1)*5 + 2};
        y_i(:,i) = inter_table{:,(i-1)*5 + 3};
        sol_PDD_i(:,i) = inter_table{:,(i-1)*5 + 5};
    end
    clear inter_table
    sten_table = readtable(sten_file);
    x_s = [];
    y_s = [];
    sol_PDD_s = [];
    for i = 1:10
        x_s(:,i) = sten_table{:,(i-1)*5 + 2};
        y_s(:,i) = sten_table{:,(i-1)*5 + 3};
        sol_PDD_s(:,i) = sten_table{:,(i-1)*5 + 5};
        sol_PDD_a(:,i) = sten_table{:,(i-1)*5 + 4};
    end
    
    clear sten_table
    figure1 = figure;
    analytic_numerical = subplot(2,2,1);
    analytic_interpolated = subplot(2,2,2);
    numerical_interpolated = subplot(2,2,3);
    numerical_analytic_interpolated = subplot(2,2,4)
    for i = 1:10
        x_p = linspace(x_i(1,i),x_i(end,i),100);
        x_p =x_p';
        y_p = linspace(y_i(1,i),y_i(end,i),100);
        y_p = y_p';
        z_p = solution(x_p,y_p);
        subplot(analytic_numerical)
        hold on
        plot(y_p,z_p,'blue')
        scatter(y_i(:,i),sol_PDD_i(:,i),'red')
        legend('Analytic solution','Numerical solution')
        subplot(analytic_interpolated)
        hold on
        z_analytic_int = RBFInterpolator(y_s(:,i),sol_PDD_a(:,i),y_p);
        %plot(y_p,z_p,'blue')
        %plot(y_p,z_analytic_int,'green')
        plot(y_p,z_p)
        plot(y_p,z_analytic_int)
        legend('Analytic solution','Interpolated value')
        subplot(numerical_interpolated)
        hold on
        z_numerical_int = RBFInterpolator(y_s(:,i),sol_PDD_s(:,i),y_p);
        orange = [0.8500, 0.3250, 0.0980];
        %scatter(y_i(:,i),sol_PDD_i(:,i),'red')
        %plot(y_p,z_numerical_int,'color',orange)
        scatter(y_i(:,i),sol_PDD_i(:,i))
        plot(y_p,z_numerical_int)
        legend('Numerical solution','Interpolated value')
        subplot(numerical_analytic_interpolated)
        hold on
        %plot(y_p,z_p,'blue')
        %plot(y_p,z_numerical_int,'color',orange)
        plot(y_p,z_p)
        plot(y_p,z_numerical_int)
    end
    figure2 = figure
    for i = 3:10 
    y_p = linspace(y_s(1,i),y_s(end,i),100);
    y_p = y_p';
    hold on
    z_numerical_int = RBFInterpolator(y_s(:,i),sol_PDD_s(:,i),y_p);
    scatter(y_i(:,i),sol_PDD_i(:,i))
    plot(y_p,z_numerical_int)
    legend('Numerical solution','Interpolated value')
    end
return
function z = solution(x,y)
C = 2;
kx = 0.47;
ky = 0.89;
z = C + sin(kx*pi*x).*sin(ky*pi*y);
return
