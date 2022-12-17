% four bar mechanism
clear all
close all
clc
l1=0.4; l2=0.1; l3=0.5; l4=0.4; omg2=4; tho1=0; tho2=90*pi/180; tho3 =36.87*pi/180; tho4=-90*pi/180;
dt=0.05;% Step size
t_end=2; t_start=0; 
n_sol=(t_end-t_start)/dt+1; %Number of Steps
error_tol=1.0E-6;

Res_mat=zeros(n_sol,37);

q_num=[0 0 0 0.5*l2*cos(tho2) 0.5*l2*sin(tho2) tho2 l2*cos(tho2)+0.5*l3*cos(tho3) l2*sin(tho2)+0.5*l3*sin(tho3) tho3 l2*cos(tho2)+l3*cos(tho3)+0.5*l4*cos(tho4) l2*sin(tho2)+l3*sin(tho3)+0.5*l4*sin(tho4) tho4 ];
qd_num=zeros(1,12);
qdd_num=zeros(1,12);
syms q1 q2 q3 q4 q5 q6 q7 q8 q9 q10 q11 q12 t
syms qd1 qd2 qd3 qd4 qd5 qd6 qd7 qd8 qd9 qd10 qd11 qd12
q=[q1 q2 q3 q4 q5 q6 q7 q8 q9 q10 q11 q12];
qd=[qd1 qd2 qd3 qd4 qd5 qd6 qd7 qd8 qd9 qd10 qd11 qd12];

C=[ q1;
    q2;
    q3;
    q4-0.5*l2*cos(q6);
    q5-0.5*l2*sin(q6);
    q4+0.5*l2*cos(q6)-q7+0.5*l3*cos(q9);
    q5+0.5*l2*sin(q6)-q8+0.5*l3*sin(q9);
    q7+0.5*l3*cos(q9)-q10+0.5*l4*cos(q12);
    q8+0.5*l3*sin(q9)-q11+0.5*l4*sin(q12);
    q10+0.5*l4*cos(q12)-l1;
    q11+0.5*l4*sin(q12);
    q6-omg2*t-tho2];

% For position
Cq=jacobian(C,q);
% For velocity you need also
Ct=diff(C,t);
% For acceleration you need in addition to Cq the following
Ctt=diff(Ct,t);
Cqt=diff(Cq,t);

Cq_qd=Cq*qd.';
Cq_qdq=jacobian(Cq_qd,q);

Qd=-Cq_qdq*qd.'-2*Cqt*qd.'-Ctt;

% At t=0
% For verification
C_num1=subs(C,q,q_num);
C_num2=subs(C_num1,t,0);
Cq_num=subs(Cq,q,q_num);
qd_num=-(Cq_num\Ct)';

Qd_num1=subs(Qd,q,q_num);
Qd_num2=subs(Qd_num1,qd,qd_num);
qdd_num=Cq_num\Qd_num2;

Res_mat(1,2:13)=q_num;
Res_mat(1,14:25)=qd_num;
Res_mat(1,26:37)=qdd_num;
Res_mat(1:n_sol,1)=[t_start:dt:t_end];



for i_res=2:n_sol
    t_num=Res_mat(i_res,1);
    q_num_n=Res_mat(i_res-1,2:13)+dt*Res_mat(i_res-1,14:25);
    error1=1.0;
    while abs(error1)>error_tol,
    C_num1=subs(C,q,q_num_n);
    C_num2=subs(C_num1,t,t_num);
    Cq_num=subs(Cq,q,q_num_n);
    
     C_num2 = vpa(C_num2);
    C_num2 = simplify(C_num2);
    Cq_num = vpa(Cq_num);
    Cq_num = simplify(Cq_num);
     

    d_q_num_n=-(Cq_num\C_num2)';
    q_num_np1=q_num_n+d_q_num_n;
    error1=eval(norm(C_num2));
    error2=eval(norm(d_q_num_n));
    q_num_n=q_num_np1;
    
    end
    
    Cq_num=subs(Cq,q,q_num_n);
    qd_num=-(Cq_num\Ct)';

    Qd_num1=subs(Qd,q,q_num_n);
    Qd_num2=subs(Qd_num1,qd,qd_num);
    qdd_num=Cq_num\Qd_num2;
    
    Res_mat(i_res,2:13)=q_num_n;
    Res_mat(i_res,14:25)=qd_num;
    Res_mat(i_res,26:37)=qdd_num;
    
end

%required(1)
figure(1);
plot(Res_mat(:,1), Res_mat(:,13));
title('(1)theta 4');
xlabel('t');
ylabel('theta4');

%required(2)
figure(2);
plot(Res_mat(:,1),Res_mat(:,23));
title('(2)theta 4 dot');
xlabel('t');
ylabel('w');

%required(3)
figure(3);
PX= Res_mat(:,8)+ 0.2*cos(Res_mat(:,10))-0.25*sin(Res_mat (:,10));
PY= Res_mat(:,9)+0.2*sin(Res_mat(:,10))+0.25*cos(Res_mat(:,10));
plot(PX,PY);
title('(3) Trace of  point (0.2, 0.25) defined in body 3 coordinate system.');
xlabel('x');
ylabel('y');



