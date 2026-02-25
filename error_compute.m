function error=error_compute(Nh,u_fdm,fdm_obj)
    %% 空间离散参数
    x_left=fdm_obj.x_left;
    x_right=fdm_obj.x_right;
    y_left=fdm_obj.y_left;
    y_right=fdm_obj.y_right;
    t_test=fdm_obj.t_test; %测试时间点

    Nx=Nh; %x方向网格数
    Ny=Nh; %y方向网格数
    dx=(x_right-x_left)/Nx; %x方向空间步长
    dy=(y_right-y_left)/Ny; %y方向空间步长
    x=x_left:dx:x_right; %x方向网格节点
    y=y_left:dy:y_right; %y方向网格节点
    N=(Nx+1)*(Ny+1); %总网格数
    u_exact=zeros(N,1); %初始化解析解矩阵
    for i=1:Nx+1
        for j=1:Ny+1
            k=(i-1)*(Ny+1)+j; %当前节点在系数矩阵中的行号
            u_exact(k)=fdm_obj.fun_analytical(t_test,x(i),y(j));
        end
    end
    error_vector=u_fdm-u_exact;
    len=length(error_vector);
    error_Linf=max(abs(error_vector));
    error_L2=sqrt(sum(error_vector.^2)/len);
    error_L1=sum(sum(abs(error_vector)))/len;
    error=[error_Linf,error_L2,error_L1];
end