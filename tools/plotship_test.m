clear
clc

%GLOBALIZE
global gm x y

%INITIALIZE
gm = -1.62;      %[m/s^2]
Rm = 1.738E6;    %[m]
Vx = 0;          %[m/s]
Vy = 0;          %[m/s]
x = 1000;        %[m]
y = 2000;        %[m]
dt = 0.1;        %[s]
h_moon = 0;      %[m]
mship = 4900;    %[kg]
mfuel = 10300;   %[kg]
m = mship+mfuel; %[kg]
Q = 0;           %[deg]
T = 0;           %[N]

%% SETUP
%CREATE FIGURE WINDOW
%Determine screen size
scrsize = get(0, 'ScreenSize'); % = [1,1,1920,1080]
%Create Figure Window
%Set bottom-left corner location
swidth = scrsize(3); %Set width
sheight = scrsize(4); %Set height
fig1 = figure('Position', [0,0,swidth/2,sheight]); %Create Figure window on
%left half of screen

%Define axes
xmin = 0;    %[m]
xmax = 2500; %[m]
ymin = xmin; %[m]
ymax = xmax; %[m]
axis([xmin,xmax,ymin,ymax]);
axis equal
axis manual
%Define background, bg
bg = [0, xmax, xmax,    0;...
      0,    0, xmax, xmax];

%PATCHES:
%Create background patch
bg_patch = patch(bg(1,:),bg(2,:),'k');

axis equal
axis on %axis off (I'm gonna want to turn it off later)
hold on
%Create moon patch
xmoon = 100*linspace(0,20,21);
ymoon = 100*[0,1,2,5,7,7,6,7,10,10,9,6,5,5,6,11,13,13,10,9,0];
moon_patch = patch(xmoon,ymoon,'w');
%Create ship patch
%xcords:
xship = 10*[0,-1,-2,-2,-3,-3,-4,-1,-2,-2,   2, 2, 1, 4, 3, 3, 2, 2, 1];
%ycords:
yship = 10*[4, 4, 3, 1, 0,-3,-4,-4,-3,-2,  -2,-3,-4,-4,-3, 0, 1, 3, 4];
hship = 40; %height of ship [m]
ship = [xship;yship]; %concatenate x and y cords
%use patch function to create ship patch
ship_patch = patch(ship(1,:)+x,ship(2,:)+y-hship/2,'b'); %ship "spawns" at
%top-center of the plot and is blue
