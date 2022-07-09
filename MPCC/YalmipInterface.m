% Copyright (C) 2018, ETH Zurich, D-ITET, Kenneth Kuchera, Alexander Liniger
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ X,U,dU,info ] = YalmipInterface(stage,MPC_vars,ModelParams)
nx = ModelParams.nx;
nu = ModelParams.nu;
tic
yalmip('clear')
x    = sdpvar(nx+nu, MPC_vars.N+1);          % states + initial state; fifth initial state for discretization
u    = sdpvar(nu, MPC_vars.N);           % (front) steering angle
objective = 0;
constraints = [x(:,1) == blkdiag(MPC_vars.Tx,MPC_vars.Tu)*[stage(1).x0;stage(1).u0]];   % initialize initial state

for i = 1:MPC_vars.N
    constraints = [constraints;
                   x(:, i+1) == stage(i).Ak*x(:,i) + stage(i).Bk*u(:,i)+ stage(i).gk ; %dynamics
                   stage(i).lb <= [x(:,i);u(:,i)] <= stage(i).ub; % bounds 注意这里有归一化了
                   stage(i).lg <= stage(i).Ck*x(:,i) <= stage(i).ug]; %track constraints
%                    pos = [-0.24;0.63];
%                    r = 0.2 ;
%                    b = (x([1:2],i)-pos)'*((x([1:2],i)-pos)) - r^2;
%                    b_next = (x([1:2],i+1)-pos)'*((x([1:2],i+1)-pos)) - r^2;
%                    constraints = [constraints; b_next - b + 0.4 * b >=
%                    0];%CLF CBF
%     objective = objective + 0.5*(x(:,i)'*stage(i).Qk*x(:,i) + u(:,i)'*stage(i).Rk*u(:,i)) ...
%                           + stage(i).fk'*x(:,i)+1000*((x(1,i+1)-x(1,i)-0.02*(x(4,i)*cos(x(3,i)) - x(5,i)*sin(x(3,i))))^2 ...
%                       +(x(2,i+1)-x(2,i)-0.02*(x(4,i)*sin(x(3,i)) + x(5,i)*cos(x(3,i))))^2 ...
%                   +(x(3,i+1)-x(3,i)-0.02*x(6,i))^2);  % cost
        cost_fun=abs(x(:, i+1)-stage(i).Ak*x(:,i) -stage(i).Bk*u(:,i)- stage(i).gk);
        objective = objective + 0.5*(x(:,i)'*stage(i).Qk*x(:,i) + u(:,i)'*stage(i).Rk*u(:,i)) ...
                          + stage(i).fk'*x(:,i)+1000*(cost_fun(1)+cost_fun(2)+cost_fun(3)+cost_fun(4)+cost_fun(5)+cost_fun(6));  % cost  
%       objective = objective + 0.5*(x(:,i)'*stage(i).Qk*x(:,i) + u(:,i)'*stage(i).Rk*u(:,i)) ...
%                           + stage(i).fk'*x(:,i);
end
i = MPC_vars.N+1;
objective = objective + 0.5*(x(:,i)'*stage(i).Qk*x(:,i)) + stage(i).fk'*x(:,i);
constraints = [constraints;
               stage(i).lb(1:nx+nu) <= x(:,i) <= stage(i).ub(1:nx+nu); %bounds
               stage(i).lg <= stage(i).Ck*x(:,i) <= stage(i).ug];  % track constraints
yalmipTime = toc;           
ops = sdpsettings();

% solve QP
tic;
exitflag = solvesdp(constraints,objective,ops);
QPtime = toc;

x_opt = double(x);
u_opt = double(u);

% rescale outputs
X= MPC_vars.invTx*x_opt(1:nx,:);
U = MPC_vars.invTu*x_opt(nx+1:end,2:end);
dU = u_opt;

info.exitflag = exitflag;
info.QPtime = QPtime;

end

