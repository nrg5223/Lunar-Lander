function [Ax, Az] = LanderAccel(x, z, Vx, Vz, Q, m, T)
% Lunar Lander Acceleration Calculator
% Calculates net acceleration of lunar lander module (LM) for final project
% See Lunar Physics Reference document for definition of coordinate system
% Input:  x    - current horiz. position of LM              {m}
%         z    - current vert. position of LM (altitude)    {m}
%         Vx   - current horiz. velocity of LM              {m/s}
%         Vz   - current vert. velocity of LM               {m/s}
%         Q    - current orientation angle of LM
%                (where 0° = vetical & Q is positive clockwise  {degrees}
%         m    - current mass of LM                         {kg}
%         T    - current thrust output of LM rocket         {N}
% Output: Ax   - current horiz. acceleration of LM          {m/s^2}
%         Az   - current vert. acceleration of LM           {m/s^2}
% Const:  gm   - lunar gravity                              {m/s^2}
%         Rm   - radius of moon                             {m}

    % Initialize constants
    gm = 1.62;      % {m/s^2}
    Rm = 1.738E6;   % {m}

    % Calculate forces on lander
    % Weight
    Wz = m * (-gm);         % {N}
    % Thrust, vert. & horiz. components
    Tz = T * cosd(Q);       % {N}
    Tx = T * sind(Q);       % {N}
    % 'Centrifugal force' to correct for orbital mechanics
    Cz = m * Vx^2 / (Rm+z); % {N}

    % Calculate horiz. & vert. acceleration of lander
    Ax = Tx / m;            % {m/s^2}
    Az = (Wz+Tz+Cz) / m;     % {m/s^2}

end