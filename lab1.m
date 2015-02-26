clear all; close all; clc;
a= 1/sqrt(2);
OA=[0 1 0 0 0 -1 0 0 0 0 0 0 0
    0 0 1 0 0 0 0 0 0 0 0 0 0
    a 0 0 -1 -a 0 0 0 0 0 0 0 0
    a 0 1 0 a 0 0 0 0 0 0 0 0
    0 0 0 1 0 0 0 -1 0 0 0 0 0
    0 0 0 0 0 0 1 0 0 0 0 0 0
    0 0 0 0 a 1 0 0 -a -1 0 0 0
    0 0 0 0 a 0 1 0 a 0 0 0 0
    0 0 0 0 0 0 0 0 0 1 0 0 -1
    0 0 0 0 0 0 0 0 0 0 1 0 0
    0 0 0 0 0 0 0 1 a 0 0 -a 0
    0 0 0 0 0 0 0 0 a 0 1 a 0
    0 0 0 0 0 0 0 0 0 0 0 a 1];
OB=[0 10 0 0 0 0 0 15 0 20 0 0 0];
%% Inverse method
A=OA;B=OB;
X_Inverse = inv(A)*transpose(B);
%% Gaussian elimination
A=OA;B=OB;
for k = 1:length(A)-1
    % Pivot
    p = 0;
    while A(k+p,k) == 0
        p=p+1;
        if k+p > length(A)
            p = -1;
            break;
        end
    end
    % Upper triangular
    if p ~= -1
        A = [A(1:k-1,:);A(k+p,:);A(k:k+p-1,:);A(k+p+1:end,:)]; %jank
        B = [B(1:k-1) B(k+p) B(k:k+p-1) B(k+p+1:end)];
        for i = k+1:length(A)
            lambda=A(i,k)/A(k,k);
            if lambda ~=0
                for j = k:length(A)
                    A(i,j) = A(i,j) - A(k,j)*lambda;
                end
                B(i) = B(i) - B(k)*lambda;
            end
        end
    end
end
% Backwards substitution
for k = length(A):-1:1
    for i = k+1:length(A)
        B(k) = B(k) - A(k,i)*B(i);
    end
    B(k) = B(k)/A(k,k);
    
end

X_Gauss = transpose(B);
%% Gaussian-Seidel
A=OA;B=transpose(OB);
x = [0 0 0 0 0 0 0 0 0 0 0 0 0];
for n=1:10
    for i=1:length(x);
        temp(i) = B(i);
        for j=1:length(x);
            if i ~= j 
                temp(i) = temp(i) - A(i,j)*x(j);
            end
        end
        x(i) = temp(i)/A(i,i);
    end
end
X_Sidel = transpose(x);
