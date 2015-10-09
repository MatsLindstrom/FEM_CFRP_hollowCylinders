function S = effCompMatrix(phi,meshDim,E,nu,G)
  % Computes the effective compliance matrix of a series of laminates
  % Based on algorithm given by Sun, et al. in "Stress analysis of multi-layered hollow anisotropic composite cylindrical structures using the homogenization method", Acta Mechanica, 2014;225(6)
  % Input:
    % Ply angle configuration in radians: phi
    % Mesh dimensions: meshdim = [r_o,r_i,dr,h,n]
      % where meshdim(1:4) is in meters and n is the number of equally spaced nodes
    % Elastic modulus: E = [E_r, E_phi, E_z] in GPa
    % Poisson's ratio: v = [v_phi_z, v_r_z, v_r_phi]
    % Shear modulus: G = [G_phi_z, G_r_z, G_r_phi] in GPa
  % Output: Effective compliance matrix

  r_o = meshDim(1);
  r_i = meshDim(2);
  dr = meshDim(3);
  lamina = struct;

  for i = 1:length(phi)
    lamina(i).S = cylCompMat(phi(i),E,nu,G);
    lamina(i).Q = inv(lamina(i).S);
    lamina(i).r1 = r_i + (i-1)*dr;
    lamina(i).r2 = r_i + (i)*dr;
  end

  t = r_o^2 - r_i^2;
  delta1 = 0;
  delta2 = 0;
  delta3 = 0;
  Qeff = zeros(6,6);
  terms = zeros(4,3,4);

  for i = 1:length(lamina)
    Q = lamina(i).Q;
    r1 = lamina(i).r1;
    r2 = lamina(i).r2;
    V = (r2^2 - r1^2)/t;

    deltai = Q(5,5)*Q(6,6) - Q(5,6)*Q(5,6);
    delta1 = delta1 + V*Q(5,5)/deltai;
    delta2 = delta2 + V*Q(6,6)/deltai;
    delta3 = delta3 + V*Q(5,6)/deltai;

    Qeff(1,1) = Qeff(1,1) + V/Q(1,1);

    terms(:,:,1) = terms(:,:,1) + V.*Q(1:4,2:4);
    terms(:,:,2) = terms(:,:,2) + V.*repmat(Q(1,1:4)',[1,3])./Q(1,1);
    terms(:,:,3) = terms(:,:,3) + V.*repmat(Q(1,2:4),[4,1])./Q(1,1);
    terms(:,:,4) = terms(:,:,4) + V.*repmat(Q(1,1:4)',[1,3]).*repmat(Q(1,2:4),[4,1])./Q(1,1);

    Qeff(5:6,5:6) = Qeff(5:6,5:6) + V.*Q(5:6,5:6)./deltai;

  end

  Qeff(1,1) = 1/Q(1,1);
  delta = delta1*delta2 - delta3^2;
  Qeff(5:6,5:6) = Qeff(5:6,5:6)./delta;
  Qeff(1:4,2:4) = terms(:,:,1) + Qeff(1,1).*terms(:,:,2).*terms(:,:,3) - terms(:,:,4);

  Qeff = triu(Qeff)+triu(Qeff,1)';
  Q = Qeff;
  S = inv(Q);