% --------------------SETUP--------------------
% normal A
syms d1 q1 a1 alfa1
syms q2 theta2 a2 alfa2
syms q3 theta3 a3 alfa3
% center A
syms d1_c q1 a1 alfa1
syms q2 theta2 a2_c alfa2
syms q3 theta3 a3 alfa3
% q's
q = [q1,q2,q3];
% q_d's
syms q1_d q2_d q3_d
q_d = [q1_d,q2_d,q3_d].';
% q_dd's
syms q1_dd q2_dd q3_dd
q_dd = [q1_dd,q2_dd,q3_dd].';
% gravedad
syms g
% --------------------ENDSETUP--------------------

% --------------------GET_Hs--------------------
% normal A
A1 = dh(d1,q1,0,0);
A2 = dh(q2,0,-a2,-sym(pi)/2);
A3 = dh(q3,0,0,0);
% center A
A1_c = dh(d1_c,q1,0,0);
A2_c = dh(q2,0,-a2_c,-sym(pi)/2);
A3_c = dh(q3,0,0,0);

H01 = A1_c;
H02 = A1*A2_c;
H03 = A1*A2*A3_c;
% --------------------ENDGET_Hs--------------------

% --------------------GET_JACOBIAN_1--------------------
% z's
z0 = [0,0,1]';
z1 = H01(1:3,3);
% t's
t0 = [0,0,0]';
tn = H01(1:3,4);

% jacobian
J_1 = [cross(z0,tn-t0), cross(z1,tn-tn), cross(z1,tn-tn);
     z0, z1-z1, z1-z1];
J_1 = simplify(J_1)
% --------------------ENDGET_JACOBIAN_1--------------------
% --------------------GET_JACOBIAN_2--------------------
% Z's
z0 = [0,0,1]';
z1 = A1(1:3,3);
z2 = H02(1:3,3);
% t's
t0 = [0,0,0]';
t1 = A1(1:3,4);
tn = H02(1:3,4);

% jacobian
J_2 = [cross(z0,tn-t0),  z1, cross(z1,tn-tn);
     z0, z1-z1, z1-z1];
J_2 = simplify(J_2)
% --------------------ENDGET_JACOBIAN_2--------------------
% --------------------GET_JACOBIAN_3--------------------
% Z's
z0 = [0,0,1]';
z1 = A1(1:3,3);
z2 = A2(1:3,3);
z3 = H03(1:3,3);
% t's
t0 = [0,0,0]';
t1 = A1(1:3,4);
t2 = A2(1:3,4);
tn = H03(1:3,4);

% jacobian
J_3 = [cross(z0,tn-t0), z1, z2;
      z0, z1-z1, z1-z1];
J_3 = simplify(J_3)
% --------------------ENDGET_JACOBIAN_3--------------------
% --------------------GET_D--------------------
% setup
syms I_x1 I_y1 I_z1
syms I_x2 I_y2 I_z2
syms I_x3 I_y3 I_z3
syms m_1 m_2 m_3
% Inertias
I_1 = [I_x1, 0 , 0;
      0, I_y1, 0;
      0, 0, I_z1];
I_2 = [I_x2, 0 , 0;
      0, I_y2, 0;
      0, 0, I_z2];
I_3 = [I_x3, 0 , 0;
      0, I_y3, 0;
      0, 0, I_z3];
% lineal
J_v1 = J_1(1:3,:);
J_v2 = J_2(1:3,:);
J_v3 = J_3(1:3,:);
linear =  m_1*J_v1.'*J_v1 + ...
          m_2*J_v2.'*J_v2 + ...
          m_3*J_v3.'*J_v3;
% angular
J_w1 = J_1(4:6,:);
J_w2 = J_2(4:6,:);
J_w3 = J_3(4:6,:);
angular = J_w1.'*(H01(1:3,1:3).'*I_1*H01(1:3,1:3))*J_w1 + ...
          J_w2.'*(H02(1:3,1:3).'*I_2*H02(1:3,1:3))*J_w2 + ...
          J_w3.'*(H03(1:3,1:3).'*I_3*H03(1:3,1:3))*J_w3;

D = linear + angular;
D = simplify(D)
% --------------------ENDGET_D--------------------
% --------------------GET_C--------------------
C = get_C(D,q,q_d);
C = simplify(C)
% --------------------ENDGET_C--------------------
% --------------------GET_P--------------------
% para este caso se asume gravedad en direccion -z
P1 = m_1*g*H01(3,4);
P2 = m_2*g*H02(3,4);
P3 = m_3*g*H03(3,4);
P = P1 + P2 + P3;

P_1 = diff(P,q1);
P_2 = diff(P,q2);
P_3 = diff(P,q3);
phi = [P_1, P_2, P_3].';
phi = simplify(phi)
% --------------------ENDGET_P--------------------
% --------------------GET_THEM_ALL--------------------
T = D*q_dd + C*q_d + phi;
T = simplify(T)
% --------------------ENDGET_THEM_ALL--------------------
