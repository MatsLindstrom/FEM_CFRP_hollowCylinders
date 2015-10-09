function S = cylCompMat(phi, E, v, G)
  % Transforms the compliance matrix of a ply from local coordinate system in the
  % fiber direction to its equivalent in the global cylindrical coordinate system
  % Input:
    % Ply angle phi in radians: phi
    % Elastic modulus: E = [E_r, E_phi, E_z] in GPa
    % Poisson's ratio: v = [v_phi_z, v_r_z, v_r_phi]
    % Shear modulus: G = [G_phi_z, G_r_z, G_r_phi] in GPa
  % Output: Compliance matrix in global coordinate system

  l = cos(phi);
  l2 = cos(2*phi);
  m = sin(phi);
  m2 = sin(2*phi);

  % Compliance matrix in material direction
  S0 = [1/E(1), -v(3)/E(1), -v(2)/E(1), 0, 0, 0;
        -v(3)/E(1), 1/E(2), -v(1)/E(2), 0, 0, 0;
        -v(2)/E(1), -v(1)/E(2), 1/E(3), 0, 0, 0;
        0, 0, 0, 1/G(1), 0 , 0;
        0, 0, 0, 0, 1/G(2) , 0;
        0, 0, 0, 0, 0, 1/G(3)];

  % Transformation matrix
  T = [1,   0,   0,   0,   0,   0;
       0, l^2, m^2, -m2, 0, 0;
       0, m^2,   l^2, m2, 0, 0;
       0, m2/2, -m2/2, l2, 0, 0;
       0, 0, 0, 0, l, m;
       0, 0, 0, 0, -m, l];

% Material compliance matrix in cylinder principal axis direction
S = inv(T)*S0*T;
