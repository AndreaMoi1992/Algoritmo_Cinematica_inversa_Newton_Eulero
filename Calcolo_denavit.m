function [A] = Calcolo_denavit(DH)
%CALCOLO_DENAVIT Summary of this function goes here
%   Formato della matrice in ingresso
%     [
%         i, alfa_i, a_i, d_i, theta_i;
%         1, alfa_1, a_1, d_1, theta_1
%         ..., ... ,... ,... , ...]
I=eye(3,3);
I_4=eye(4,4);
o=[0, 0, 0];

A=sym(zeros(4,4));

alfa=DH(:,1);
a=DH(:,2);
d=DH(:,3);
theta=DH(:,4);

numBracci=length(DH(:,1));

for i=1:numBracci
    
    qZ=[0; 0; d(i)];
    qX=[a(i); 0; 0];
    
    Rz=[
        cos(theta(i)) -sin(theta(i)) 0;
        sin(theta(i)) cos(theta(i)) 0;
        0 0 1;
        ];

    Rx=[
        1 0 0;
        0 cos(alfa(i)) -sin(alfa(i));
        0 sin(alfa(i)) cos(alfa(i));
        ];    
    
    T_Trasl_z=[[I;o],[qZ;1]];
    T_Trasl_x=[[I;o],[qX;1]];
    T_Rot_z=[[Rz;o],[o';1]];
    T_Rot_x=[[Rx;o],[o';1]];
    
    T=I_4*T_Trasl_z*T_Rot_z*T_Trasl_x*T_Rot_x;
    
    A(:,:,i)=T;
    
end

end


