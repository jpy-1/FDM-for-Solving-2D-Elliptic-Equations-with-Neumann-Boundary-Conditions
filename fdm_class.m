%fdm_class.m 封装有限差分法求解的类
%该类定义了二维空间域的边界，并提供了
%解析解和源项函数
%fun_analytical: 解析解函数
%fun_source: 源项函数
%fun_p: 系数函数p(x,y)
%fun_q: 系数函数q(x,y)

%求解的抛物型方程形式为：
% u_t-\nabla \left( p\left( x \right) \nabla u \right) +q\left( x \right) u=f
%边界条件
%u(x_left,y,t)=0
%u(x,y_left,t)=0
%u_x(x_right,y,t)=0  
%u_y(x,y_right,t)=0

classdef fdm_class
    properties(Constant)
        x_left=0;    % Left boundary of the spatial domain
        x_right=1;   % Right boundary of the spatial domain
        y_left=0;  % Left boundary of the spatial domain
        y_right=1;     % Right boundary of the spatial domain
        t_test=2; %测试时间点
    end
    
    methods
        function y=fun_analytical(~,t,x,y)
            y=t^2.*sin(pi/2*x).*sin(pi/2*y);
        end
        
        function y=fun_source(~,t,x,y)
            y= 2*t*sin(pi/2*x).*sin(pi/2*y)+...
               -pi*t^2.*cos(pi/2*x).*sin(pi/2*y)+...
               +pi^2*x.*t^2.*sin(pi/2*x).*sin(pi/2*y)+...
               +3*t^2.*y^2.*sin(pi/2*x).*sin(pi/2*y);
        end

        function z=fun_p(~,x,y)
            z=2*x;
        end

        function z=fun_q(~,x,y)
            z=3*y^2;
        end
    end
end