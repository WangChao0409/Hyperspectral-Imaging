% Mar 2023 by Wang Chao
% This file gives the description for the hyperspectral imaging.
% This code is consisted of three parts:
% 1. load data: X=[x1, x2, ..., xk] - the signal data containing the information of objects
%               M=diag matrix - containing the information optical system
%               A - matrix should be designed
%               Y=A*M*X
% 2. process the trained data X by using NMF
%    find matrix W and H such that X=W*H
%
% 3. optimize A (GA with Spline Interpolation)
%                
% 4. inverse problem: y=A*M*x=A*M*D*z
%                     z = (A*M*D)^(-1)*y 
% 5. show results

%%
% initialization
n = 1201; 
m = 9; % m=16; 
k = 234;   % n - dimensional of trained data;
                             % m - number of designed channel;
                             % k - number of the trained samples
% n = 721;  % 400nm-760nm
 nH = 50;

%% 1. Data Load                              
% load data for the optical system
% L(\lambda_j) intensity of light source
load light_spectrums.mat
natural_light = data(3,101:1301);   % intensity for natural light: 400nm-1000nm
halogen_lamp = data(6,101:1301);    % intensity for halogen lamp: 400nm-1000nm
LED = data(8,101:1301);             % intensity for LED: 400nm-1000nm
% 
% natural_light = data(3,101:821);   % intensity for natural light: 400nm-760nm
% halogen_lamp = data(6,101:821);    % intensity for halogen lamp: 400nm-760nm
% LED = data(8,101:821);             % intensity for LED: 400nm-760nm

% T_len(\lambda_j) transmittance for lens
load MTV5-IR(3MP).mat
T_len = data;

% T_sensor quantum efficiency of sensor
load IMX273LLR-C.mat
T_sensor = data;

% load data for the object 
load 24_seka.mat
R_24 = data; % 400nm-1000nm
% R_24 = data(:,1:721); % 400nm-760nm

% load data for the 6 colors
%load ColorData
%x_color = ColorData(:,3:8);

load PANTONE_seka.mat
R_pantone = data; % 400nm-1000nm
% R_pantone = data(:,1:721); % 400nm-760nm

% n_24 = 1;  n_pantone = 1;
% construction of M
M = zeros(n,n);
for i=1:n
   M(i,i) = halogen_lamp(i)*T_len(i)*T_sensor(i);
end

% construction of X
x_24 = R_24';
x_pantone = R_pantone';
X = [ x_24 x_pantone ];  % trained data
% x_train = data_interp(1:300,41:1241);
% X = x_train';


%% % 2. based on NMF method

[p,q] = size(X);

[eNMF,W,H] = NMF(X,nH,'maxiter', 40);

%% 3. Genetic algorithm with interpolation
%%9通道优化算法
nvars=m*10;%带求参数个数
lb=0;
ub=1;
PopulationSize_Data=200;
MaxGenerations_Data=300;
[x,fval,exitflag,output,population,score] = ga_op(nvars,lb,ub,PopulationSize_Data,MaxGenerations_Data);

xx = 400:60:940;
xq = 400:0.5:1000;

for i=1:m
 
    B(i,:) = x((i-1)*10+1:i*10);
    A(i,:) = interp1(xx,B(i,:),xq,'spline');

end


%
%%16通道优化算法
% nvars=m*10;%带求参数个数
% lb=0;
% ub=1;
% PopulationSize_Data=200;
% MaxGenerations_Data=300;
% [x,fval,exitflag,output,population,score] = ga_op(nvars,lb,ub,PopulationSize_Data,MaxGenerations_Data);

% for i=1:m
% 
%     B(i,:) = x((i-1)*10+1:i*10);
%     A(i,:) = interp1(xx,B(i,:),xq,'spline');
% 
% end

%% 4 Inverse problem
for i=1:q
    
    x = X(:,i);  % original signal
    y = A*M*x; 
    
    % first method
    h1 = pinv(A*M*W)*y;
    x1_rec(:,i) = W*h1;
%     y1_rec = A*M*x1_rec;
    % second method
    alpha=0.001;
    h2 = ((A*M*W)'*(A*M*W)+alpha*eye(nH))^(-1)*(A*M*W)'*y;
    x2_rec(:,i) = W*h2;
%     y2_rec = A*M*x2_rec;
    
    Res_Err1_x(i) = norm(x-x1_rec(:,i))/norm(x);
    Err1_x(i) = norm(x-x1_rec(:,i));
    Res_Err2_x(i) = norm(x-x2_rec(:,i))/norm(x);
    Err2_x(i) = norm(x-x2_rec);
%     error1_y(i) = norm(y-y1_rec);
%     error2_y(i) = norm(y-y2_rec);
    
end


%% 5 result 
% figure
% plot(1:1201,X(:,235),'b')
% hold on
% plot(1:1201,X(:,236),'g')
% hold on
% plot(1:1201,X(:,237),'y')
% hold on
% plot(1:1201,X(:,238),'color','[0.85 0.3250 0.0980]')
% hold on
% plot(1:1201,X(:,239),'m')
% hold on
% plot(1:1201,X(:,240),'r')

% figure
% y235=A*M*X(:,235);
% plot(1:m,y235,'b')
% hold on
% y236=A*M*X(:,236);
% plot(1:m,y236,'g')
% hold on
% y237=A*M*X(:,237);
% plot(1:m,y237,'y')
% hold on
% y238=A*M*X(:,238);
% plot(1:m,y238,'color','[0.85 0.3250 0.0980]')
% hold on
% y239=A*M*X(:,239);
% plot(1:m,y239,'m')
% hold on
% y240=A*M*X(:,240);
% plot(1:m,y240,'r')

% % 24 color checker
% result_24seka=zeros(24,24);
% for i=1:24
%    for j=i:24
%       temp_cos = (X(:,i)'*X(:,j))/(norm(X(:,i))*norm(X(:,j)));
%       result_24seka(i,j) = temp_cos;
%       result_24seka(j,i) = result_24seka(i,j);
%    end
% end
% figure
% imagesc(result_24seka);
% colorbar


result_xcos = zeros(24,24);
for i = 1:24
    for j = i:24
        temp_cos = (x1_rec(:,i)'*x1_rec(:,j))/(norm(x1_rec(:,i))*norm(x1_rec(:,j)));
        result_xcos(i,j) = temp_cos;
        result_xcos(j,i) = result_xcos(i,j);
    end
end
figure
imagesc(result_xcos);
colorbar
 
Y=zeros(m,q);
for i = 1:q
   Y(:,i) = A*M*X(:,i); 
end
result_ycos = zeros(24,24);
for i = 1:24
    for j = i:24
        temp_cos = (Y(:,i)'*Y(:,j))/(norm(Y(:,i))*norm(Y(:,j)));
        result_ycos(i,j) = temp_cos;    
        result_ycos(j,i) = result_ycos(i,j);
    end
end
figure
imagesc(result_ycos);
colorbar

% result_xcos = zeros((q-1),(q-2));
% for i = 1:q-1
%     
%     for j = i+1:q
% 
%         temp_cos = (x1_rec(:,i)'*x1_rec(:,j))/(norm(x1_rec(:,i))*norm(x1_rec(:,j)));
%         result_xcos(i,j) = temp_cos;
%         
%     end
% 
% end
% 
% Y=zeros(m,q);
% for i = 1:q
%    Y(:,i) = A*M*X(:,i); 
% end
% result_ycos = zeros((q-1),(q-2));
% for i = 1:q-1
%     
%     for j = i+1:q
% 
%         temp_cos = (Y(:,i)'*Y(:,j))/(norm(Y(:,i))*norm(Y(:,j)));
%         result_ycos(i,j) = temp_cos;
%         
%     end
% 
% end
