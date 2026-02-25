function u_fdm=fun_fdm(Nh,Nt,fdm_obj)
%% 空间离散参数
x_left=fdm_obj.x_left;
x_right=fdm_obj.x_right;
y_left=fdm_obj.y_left;
y_right=fdm_obj.y_right;
t_test=fdm_obj.t_test; %测试时间点

%% 网格剖分与步长计算
Nx=Nh; %x方向网格数
Ny=Nh; %y方向网格数

dx=(x_right-x_left)/Nx; %x方向空间步长
dy=(y_right-y_left)/Ny; %y方向空间步长
dt=t_test/Nt; %时间步长
x=x_left:dx:x_right; %x方向网格节点
y=y_left:dy:y_right; %y方向网格节点


%% 系数矩阵的构造
% 初始化数值解矩阵
u=zeros(Nx+1,Ny+1); %t=0时刻的数值解矩阵

%计算u_{0,0}到u_{Nx,Ny}的系数矩阵
N=(Nx+1)*(Ny+1); %总节点数
% L=zeros(N,N); %初始化系数矩阵A
% L=sparse(L); %将矩阵L转换为稀疏矩阵以节省内存
% for i=1:Nx+1
%     for j=1:Ny+1
%         k=(i-1)*(Ny+1)+j; %当前节点在系数矩阵中的行号
%         if i==1 ||j==1 %位于左边界和下边界Dirichlet边界条件的点
%             L(k,k)=0;
%             continue;
%         end
        
%         %在x方向上计算p的值
%         if i==Nx+1  %右边界Neumann边界条件
%             pxp=fdm_obj.fun_p(x(i),y(j));%p_{i+1/2,j}，使用Neumann边界条件
%             pxm=(fdm_obj.fun_p(x(i-1),y(j))+fdm_obj.fun_p(x(i),y(j)))/2;%p_{i-1/2,j}
%         else
%             pxp=(fdm_obj.fun_p(x(i+1),y(j))+fdm_obj.fun_p(x(i),y(j)))/2;%p_{i+1/2,j}
%             pxm=(fdm_obj.fun_p(x(i-1),y(j))+fdm_obj.fun_p(x(i),y(j)))/2;%p_{i-1/2,j}
%         end
%         %在y方向上计算p的值
%         if j==Ny+1 %上边界Neumann边界条件
%             pyp=fdm_obj.fun_p(x(i),y(j));%p_{i,j+1/2}，使用Neumann边界条件
%             pym=(fdm_obj.fun_p(x(i),y(j-1))+fdm_obj.fun_p(x(i),y(j)))/2;%p_{i,j-1/2}
%         else
%             pyp=(fdm_obj.fun_p(x(i),y(j+1))+fdm_obj.fun_p(x(i),y(j)))/2;%p_{i,j+1/2}
%             pym=(fdm_obj.fun_p(x(i),y(j-1))+fdm_obj.fun_p(x(i),y(j)))/2;%p_{i,j-1/2}
%         end

%         q_ij=fdm_obj.fun_q(x(i),y(j));

%         %构造L矩阵
%         L(k,k)=q_ij+(pxp+pxm)/dx^2+(pyp+pym)/dy^2;
%         if i==Nx+1 %右边界Neumann边界条件
%             L(k,k-(Ny+1))=-(pxm+pxp)/dx^2; %u_{i-1,j}
%         else
%             L(k,k-(Ny+1))=-pxm/dx^2; %u_{i-1,j}
%             L(k,k+(Ny+1))=-pxp/dx^2; %u_{i+1,j}
%         end

%         if j==Ny+1 %上边界Neumann边界条件
%             L(k,k-1)=-(pym+pyp)/dy^2; %u_{i,j-1}
%         else
%             L(k,k-1)=-pym/dy^2; %u_{i,j-1}
%             L(k,k+1)=-pyp/dy^2; %u_{i,j+1}
%         end     
%     end
% end

% 预分配非零元素存储空间（关键优化）
max_nonzeros = 5 * (Nx * Ny); % 理论最大非零元素数（每个内部点最多5个非零元素）
row_idx = zeros(max_nonzeros, 1);
col_idx = zeros(max_nonzeros, 1);
vals = zeros(max_nonzeros, 1);
count = 0; % 当前非零元素计数

% 优化2: 仅遍历内部点（避免Dirichlet边界点的处理）
for i = 2:Nx+1
    for j = 2:Ny+1
        k = (i-1)*(Ny+1) + j; % 当前节点索引     
        % 计算扩散系数
         %在x方向上计算p的值
        if i==Nx+1  %右边界Neumann边界条件
            pxp=fdm_obj.fun_p(x(i),y(j));%p_{i+1/2,j}，使用Neumann边界条件
            pxm=(fdm_obj.fun_p(x(i-1),y(j))+fdm_obj.fun_p(x(i),y(j)))/2;%p_{i-1/2,j}
        else
            pxp=(fdm_obj.fun_p(x(i+1),y(j))+fdm_obj.fun_p(x(i),y(j)))/2;%p_{i+1/2,j}
            pxm=(fdm_obj.fun_p(x(i-1),y(j))+fdm_obj.fun_p(x(i),y(j)))/2;%p_{i-1/2,j}
        end
        %在y方向上计算p的值
        if j==Ny+1 %上边界Neumann边界条件
            pyp=fdm_obj.fun_p(x(i),y(j));%p_{i,j+1/2}，使用Neumann边界条件
            pym=(fdm_obj.fun_p(x(i),y(j-1))+fdm_obj.fun_p(x(i),y(j)))/2;%p_{i,j-1/2}
        else
            pyp=(fdm_obj.fun_p(x(i),y(j+1))+fdm_obj.fun_p(x(i),y(j)))/2;%p_{i,j+1/2}
            pym=(fdm_obj.fun_p(x(i),y(j-1))+fdm_obj.fun_p(x(i),y(j)))/2;%p_{i,j-1/2}
        end

        q_ij=fdm_obj.fun_q(x(i),y(j));
    
        % 1. 对角线元素
        count = count + 1;
        row_idx(count) = k;
        col_idx(count) = k;
        vals(count) = q_ij + (pxp + pxm)/dx^2 + (pyp + pym)/dy^2;
        
        % 2.  (i-1,j)


        if i==Nx+1
            %u_{i-1,j}
            count = count + 1;
            row_idx(count)=k;
            col_idx(count)=k-(Ny+1);
            vals(count)=-(pxp+pxm)/dx^2;
        else
            %u_{i-1,j}
            count=count+1;
            row_idx(count)=k;
            col_idx(count)=k-(Ny+1);
            vals(count)=-pxm/dx^2;
            %u_{i+1,j}
            count=count+1;
            row_idx(count)=k;
            col_idx(count)=k+(Ny+1);
            vals(count)=-pxp/dx^2;
        end

        if j==Ny+1
            %u_{i,j-1}
            count=count+1;
            row_idx(count)=k;
            col_idx(count)=k-1;
            vals(count)=-(pym+pyp)/dy^2;
        else
            %u_{i,j-1}
            count=count+1;
            row_idx(count)=k;
            col_idx(count)=k-1;
            vals(count)=-pym/dy^2;
            %u_{i,j+1}
            count=count+1;
            row_idx(count)=k;
            col_idx(count)=k+1;
            vals(count)=-pyp/dy^2;
        end
    end
end

% 优化3: 一次性构建稀疏矩阵（关键性能提升）
L = sparse(row_idx(1:count), col_idx(1:count), vals(1:count), N, N);



%Crank-Nicolson时间推进系数矩阵
I=speye(N); %单位矩阵
A=I+dt/2*L; %左侧系数矩阵
B=I-dt/2*L; %右侧系数矩阵

%% 时间推进求解
F=zeros(N,1); %初始化右侧项
% F_n=zeros(N,1); %初始化右侧项
% F_np=zeros(N,1); %初始化下一个时间点的右侧项
u=reshape(u',1,[])'; %将二维矩阵按行转换为列向量
for n=1:Nt
    u_n=u; %上一个时间步的数值解
    % t_n=(n-1)*dt; %当前时间点
    % t_np=n*dt; %下一个时间点
    for i=1:Nx+1
        for j=1:Ny+1
            k=(i-1)*(Ny+1)+j; %当前节点在系数矩阵中的行号
            if i==1 ||j==1 %位于左边界和下边界Dirichlet边界条件的点
                F(k)=0;
                % F_n(k)=0;
                % F_np(k)=0;
            else
                F(k)=fdm_obj.fun_source(dt*(n-1/2),x(i),y(j));
                % F_n(k)=fdm_obj.fun_source(t_n,x(i),y(j));
                % F_np(k)=fdm_obj.fun_source(t_np,x(i),y(j));
            end
        end
    end 
    b=B*u_n+dt*F; %右侧项计算
    % b=B*u_n+dt/2*(F_n+F_np); %右侧项计算
    %求解线性方程组
    u=A\b;%当前时间步的数值解
end
%将列向量u转换为二维矩阵
% u_fdm=reshape(u,Nx+1,Ny+1)';
u_fdm=u;
end