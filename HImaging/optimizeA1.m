function result = optimizeA( A1 )
% Chao Wang 
% 28 Jun 2023
% This function gives the objective function of the optimal problem.

A1=optimresults.x;
q=0;
chanel=9;
for i=1:chanel
     for j=1:10
          A2(i,j)=A1(j+q*10);
     end
     q=q+1;
end
p2=1:(10/1334):10;
for i=1:chanel
     A(i,:)=interp1(1:10,A2(i,:),p2,'spline');
end
for i=1:9
   for j=1:1201
       if A(i,j)<0
          A(i,j)=0.0001;
        elseif A(i,j)>1
            A(i,j)=0.999;
        end
    end
end

X1=load( 'X.mat');
W1=load( 'W.mat' );
M1=load( 'M.mat' );
q=0;
%m=12;
chanel=9;
for i=1:chanel
    for j=1:20
        A2(i,j)=A1(j+q*20);
    end
    q=q+1;
end
p2=1:(10/1334):120;
for i=1:chanel
    A(i,:)=interp1(1:20,A2(i,:),p2,'spline');
end
X = X1.X;
W = W1.W;
M = M1.M;

[m,n] = size(X);

result_L2 =0;
alpha = 0.1;
beta = 0.01;

grad_A = gradient(A);
NorGradA =0;
for i = 1:chanel
   temp = grad_A(i,:);
   NorGradA = NorGradA + norm(temp);
end

X_cal =zeros(m,n);
for i = 1:n
   x = X(:,i);
   h = calh(A,M,W,x);
   x_cal = W*h;
   X_cal(:,i) = x_cal;
end

for i = 1:n-1

    temp_L2 = norm(X(:,i)-X_cal(:,i)) ;
    
    result_cos =0;
    for j = i+1:n

        temp_cos = (X_cal(:,i)'*X_cal(:,j))/(norm(X_cal(:,i))*norm(X_cal(:,j)));
        result_cos = result_cos + temp_cos;
        
    end
   
    result_L2 = result_L2 + beta*result_cos + temp_L2;
end

temp_L3 = norm(X(:,n)-X_cal(:,n));
result_L2 = result_L2+temp_L3;

% result_gray = 0;
% gamma = 0.01;
% for i=1:n-1
%    
%     x1 = X(:,i);
%     y1 = A*M*x1*0.5;
%     
%     for j=i+1:n
%        x2 = X(:,j);
%        y2 = A*M*x2*0.5;
%        temp_gray = (y1'*y2)/(norm(y1)*norm(y2));
%        result_gray = result_gray + temp_gray;
%     end
%     
% end

result = ( result_L2  + alpha*NorGradA )/n;
% result = ( result_L2 + beta*result_cos + gamma*result_gray + alpha*NorGradA )/n;

end

function h = calh( A, M, W, x )

        y = A*M*x;
        h = pinv(A*M*W)*y;

end