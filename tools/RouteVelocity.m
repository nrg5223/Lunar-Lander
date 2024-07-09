%Routespeed
clear,clc

%Initialize Data
dist = [560, 440, 490, 530, 370];
time = [10.3, 8.2, 9.1, 10.1, 7.5];
Vhigh = 0;

%Divide distance element by respective time element for velocity
for n = 1:5
    v(n) = dist(n)/time(n);
    
    %find highest vel(n) value by comparison
    if v(n) > Vhigh
        Vhigh = v(n);
    end
end

%Output velocities and highest velocity
disp('Velocities of each route: ')
disp(v)
disp('Highest velocity: ')
disp(Vhigh)

%Output formatted velocities and highest velocity
%complete later