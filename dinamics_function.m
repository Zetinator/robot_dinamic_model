function out = pregunta_4_5(input)
  % process inputs
  q1 = input(1); q2 = input(2); q3 = input(3);
  q1_d = input(4); q2_d = input(5); q3_d = input(6);

  % set up q_deseada
  q_deseada_1 = pi/4;
  q_deseada_2 = .2;
  q_deseada_3 = .2;
  q = [q1; q2; q3];
  q_d = [q1_d;q2_d; q3_d];
  q_deseada = [q_deseada_1; q_deseada_2; q_deseada_3];
  error_q = q - q_deseada;

  % gains
  k = 5.0 * eye(3);
  k_d = 1.0 * eye(3);

  % set values to parameters
  m_1 = 1; m_2 = .4; m_3 = .4;
  d1 = 1; d1_c = d1/2;
  a2 = .25; a2_c = a2/2;
  I_x1 = .05; I_x2 = .01; I_x3 = .01;
  I_y1 = .05; I_y2 = .01; I_y3 = .01;
  I_z1 = .1; I_z2 = .02; I_z3 = .05;
  g = 9.81;

  % M(q)q_dd + C(q,q_d)q_d + phi
  M = [[ I_z1 + a2^2*m_3 + a2_c^2*m_2 + m_3*q3^2 + I_y2*cos(q1)^2 + I_y3*cos(q1)^2 + I_x2*sin(q1)^2 + I_x3*sin(q1)^2,         0, -m_3*(a2*cos(q1) + q3*sin(q1))]
      [                                                                                                           0, m_2 + m_3,                              0]
      [                                                                              -m_3*(a2*cos(q1) + q3*sin(q1)),         0,                            m_3]];
  C = [[ q1_d*((I_x2*sin(2*q1))/2 + (I_x3*sin(2*q1))/2 - (I_y2*sin(2*q1))/2 - (I_y3*sin(2*q1))/2) + m_3*q3*q3_d, 0, m_3*(q3*q1_d - q3_d*sin(q1))]
      [                                                                                                      0, 0,                            0]
      [                                                               -m_3*q1_d*(q3 - a2*sin(q1) + q3*cos(q1)), 0,                            0]];
  phi = [[0]
      [g*(m_2 + m_3)]
      [0]];

  % control signal
  U = -k*error_q - k_d*q_d + phi;
  q_dd = inv(M)*(-C*q_d - phi + U);

  out = [q_dd; error_q];
