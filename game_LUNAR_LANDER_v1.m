%FINAL PROJECT: LUNAR LANDER
%clear workspace
clear
clc

%GLOBALIZE
global x y Q T mfuel m Tmax dT start tryagain %(Variables required for
%keyboard or mouse input)

%INITIALIZE
start = 0;          %variable representing whether the start button has
                    %been pushed.  = 0 means it hasn't, = 1 means it has.
point = 0;          %variable that counts player's points
tryagain = 0;       %variable that acts like the start variable but keeps
                    %track of whether the "Try Again!" button has been
                    %pressed
tabutton = 0;       % = 0 if the "Try Again!" button has not been spawned
                    % = 1 if it has
gm = 1.62;          %gravity of the moon       [m/s^2]
Rm = 1.738E6;       %radius of the moon        [m]
Vx = 0;             %velocity in the x-direction of the ship [m/s]
Vy = 0;             %velocity in the y-direction of the ship [m/s]
x = 1000;           %x-position of the ship    [m]
y = 1960;           %y-position of the ship    [m]
t = 0;              %time                      [s]
dt = 0.1;           %change in time            [s]
h_moon = 0;         %y-coordinate of the moon directly beneath the ship [m]
mship = 4900;       %mass of the ship          [kg]
mfuel = 10300;      %mass of the fuel          [kg]
mtot = mship+mfuel; %total mass of the ship (sum of ship and fuel) [kg]
Q = 0;              %angle between ship's vertical & plot's vertical [deg]
T = 0;              %thurst acting on the ship [N]
Tmax = 45000;       %maximum thrust that can act on the ship[N]
dT = Tmax/10;       %change in thrust          [N]

%% SETUP
%Create figure window using function
[fig1] = create_fig;

%Register KeyPressFcn to allow figure window to respond to keypresses
set(fig1, 'KeyPressFcn', @KeyPress);

%Create Title screen with explanation of game controls and mission
[Readybutton,Mission,CONTROLS,ctr1,ctr2,ctr3,ctr4,GL] = tscreen;

%Create infinite delay loop to allow user to remain on Title screen until
%pressing the "Ready for Takeoff?" button
while start ~= 1
    pause(0.1);
end

%Reset start variable value so function handle @start can be reused for the
%"Start Mission!" button
start = 0;

%Delete Title screen text boxes
deletetscreen(Readybutton,Mission,CONTROLS,ctr1,ctr2,ctr3,ctr4,GL);

%Create plot using function
[bg] = create_axis;

%Create "setting" using function
[bg_patch,xmoon,ymoon,moon_patch] = create_setting(bg);

%Create ship coordinates and ship patch using function
[ship,ship_patch,hship] = create_ship(x,y);

%Create text boxes to display title and important variable values using
%function
text_boxes;

%Create editable text boxes that display important variable values using
%function
[FuelnumBox,AltnumBox,VxnumBox,VynumBox,TnumBox] = ...
    editable_boxes(mfuel,y,Vx,Vy,T);

%% GAMEPLAY

while point ~= 1 %Loop to determine if player has won

%Create infinite delay loop to allow user to start game with start button
while start ~= 1
    pause(0.1);
end

while y-hship > h_moon %game plays while base of ship is above moon surface
    
    %Calculate acceleration variables using function
    [Ax, Ay] = LanderAccel(x, y, Vx, Vy, Q, mtot, T);
    
    %Calculate position and velocity variables values using function
    [x,y,Vx,Vy] = kinematics(x,y,Vx,Vy,Ax,Ay,dt);
    
    %Update position of ship
    set(ship_patch,'XData',ship(1,:)+x,'YData',ship(2,:)+y);
    
    %Create a rotation matrix to update ship's orientation
    RtnMtrx = [cosd(-Q), sind(-Q); -sind(-Q), cosd(-Q)];
    %Multiply the ship coordinates by the rotation matrix to reorient it
    shiprotate = RtnMtrx*ship;
    
    %Update orientation of ship
    set(ship_patch,'XData',shiprotate(1,:)+x,'YData',shiprotate(2,:)+y);
    
    %Calculate h_moon at x-cord corresponding to x, and hship using
    %function
    [h_moon] = calc_h_moon(x,ymoon);
    
    %Calculate new fuel and mass values
    mfuel = mfuel-5*10^(-3)*T*dt; %[kg]
    mtot = mship + mfuel; %[kg]
    
    %Turn off Thrust if there is no fuel left
    if mfuel<=0
        T = 0; %[N]
        %Set fuel equal to 0 (sometimes it becomes negative due to
        %precision issues with increments, and that looks weird)
        mfuel = 0; %[kg]
    end
    
    %Update important values using UIControl function
    %Mass of fuel left
    FuelnumBox.String = num2str(mfuel,'%7.2f');
    %Altitude
    AltnumBox.String = num2str(y,'%6.3f');
    %Velocity in x direction
    VxnumBox.String = num2str(Vx,'%5.3e');
    %Velocity in y direction
    VynumBox.String = num2str(Vy,'%5.3e');
    %Thrust magnitude
    TnumBox.String = num2str(T,'%7.2f');

    pause(dt/10); %pause between iterations; dividing dt by 10 makes it
    %more playable
end

%Evaluate player's conditions and output message in text box using function
[point] = evaluate(Vy,mfuel,x,Q);
    
%Set velocity and thrust values equal to zero now that program is finished
%(program doesn't update them outside of loop)
Vx = 0;
Vy = 0;
T = 0;
%Update velocity and thrust text boxes
%Velocity in x direction
VxnumBox.String = num2str(Vx,'%4.2f');
%Velocity in y direction
VynumBox.String = num2str(Vy,'%4.2f');
%Thrust magnitude
TnumBox.String = num2str(T,'%4.2f');

%Spawn try again button if conditions are met
if point ~= 1 && tabutton ~= 1
    %Function to spawn tryagain button
    tryagain_button
    %Change value of tabutton
    tabutton = 1;
end

while tryagain ~= 1 && point ~= 1
    %Reset position, fuel, and angle
    x = 1000;
    y = 1960;
    mfuel = 10300;
    Q = 0;
   %Pause until user presses "Try Again!" button
    pause(0.1);
 end

%Reorient ship and place at starting point
RtnMtrx = [cosd(-Q), sind(-Q); -sind(-Q), cosd(-Q)];
shiprotate = RtnMtrx*ship;
set(ship_patch,'XData',shiprotate(1,:)+x,'YData',shiprotate(2,:)+y);

%Reset tryagain variable
tryagain = 0;
end

%Pause for 5 seconds to allow for appreciation of spectacular programming
%ability
pause(5) %lol
%Clear everything once game is over
clear all
clc

%% FUNCTIONS

function [fig1] = create_fig
%Creates a figure window on the left half of the screen using screen
%dimensions

%Inputs: none
%Outputs: fig1 - figure window in which the game is played

%Determine screen size
scrsize = get(0, 'ScreenSize'); % = [1,1,1920,1080]
%Create Figure Window
%Set bottom-left corner location
swidth = scrsize(3); %Set width
sheight = scrsize(4); %Set height
fig1 = figure('Position', [0,0,swidth/2,sheight]); %Create Figure window on
%left half of screen
end

function [Readybutton,Mission,CONTROLS,ctr1,ctr2,ctr3,ctr4,GL] = tscreen
%Creates a title screen on the figure, which includes the Title of the
%game, instructions on how to use the controls, a good luck message,
%and a button that brings the player to the play screen

%Inputs: none
%Outputs: Readybutton - button that takes user to play screen
%         Mission     - text box that tells the player about their mission
%         CONTROLS    - text box with the text, "CONTROLS:"
%         ctr1        - text box explaining the up arrow control
%         ctr2        - text box explaining the down arrow control
%         ctr3        - text box explaining the right arrow control
%         ctr4        - text box explaining the left arrow control
%         GL          - text box with a good luck message

%Title
uicontrol('Style','text','String','***LUNAR***LANDER***',...
    'Position',[100,750,600,50],'FontSize',30);
%"Ready for Takeoff?" button
Readybutton = uicontrol('Style','pushbutton','String','Ready for Takeoff?',...
    'Position',[325,710,150,30],'Callback',@Start,'FontSize',10);
%Mission details
Mission = uicontrol('Style','text','String',...
    'YOUR MISSION: land on a flat surface oriented upright, with a low vertical velocity.',...
    'Position',[100,525,600,150],'FontSize',24);
%Control title
CONTROLS = uicontrol('Style','text','String','CONTROLS:','Position',...
    [100,475,600,50],'FontSize',28);
%Up arrow
ctr1 = uicontrol('Style','text','String',...
    'up arrow = increase thrust by an increment','Position',...
    [100,425,600,50],'FontSize',18);
%Down arrow
ctr2 = uicontrol('Style','text','String',...
    'down arrow = decrease thrust by an increment','Position',...
    [100,375,600,50],'FontSize',18);
%Right arrow
ctr3 = uicontrol('Style','text','String',...
    'right arrow = rotate ship clockwise by an increment',...
    'Position',[100,325,600,50],'FontSize',18);
%Left arrow
ctr4 = uicontrol('Style','text','String',...
    'left arrow = rotate ship counterclockwise by an increment','Position',...
    [50,275,700,50],'FontSize',18);
%Good luck message
GL = uicontrol('Style','text','String','GOOD LUCK!',...
    'Position',[240,175,300,50],'FontSize',24);
end

function deletetscreen(Readybutton,Mission,CONTROLS,ctr1,ctr2,ctr3,ctr4,GL)
%Deletes all text boxes and buttons on title screen except for the Title

%Inputs:  Readybutton - button that takes user to play screen
%         Mission     - text box that tells the player about their mission
%         CONTROLS    - text box with the text, "CONTROLS:"
%         ctr1        - text box explaining the up arrow control
%         ctr2        - text box explaining the down arrow control
%         ctr3        - text box explaining the right arrow control
%         ctr4        - text box explaining the left arrow control
%         GL          - text box with a good luck message
%Outputs: none

%Use delete function to delete text boxes from title screen
delete(Readybutton);
delete(Mission);
delete(CONTROLS);
delete(ctr1);
delete(ctr2);
delete(ctr3);
delete(ctr4);
delete(GL);
end

function [bg] = create_axis
%Creates a background and axis for the plot

%Inputs: none
%Outputs: bg - array containing coordinates that make up the background of
%              the plot [m]

%Define axes
xmin = 0;    %[m]
xmax = 2000; %[m]
ymin = xmin; %[m]
ymax = xmax; %[m]
axis([xmin,xmax,ymin,ymax]);
axis manual %let axis use defined values
axis equal %use defined values for axis
%Define background, bg
bg = [0, xmax, xmax,    0;...
      0,    0, xmax, xmax]; %[m]
end

function [bg_patch,xmoon,ymoon,moon_patch] = create_setting(bg)
%Creates a black background for the plot so the background looks like
%space, generates arrays to represent the moon's surface, and creates white
%patches on the plot to appear like the moon

%Inputs:  bg         - array containing coordinates that make up the background of
%                      the plot [m]
%Outputs: bg_patch   - patch for black background on the plot
%         xmoon      - array containing x-coordinates of the moon patch [m]
%         ymoon      - array containing y-coordinates of the moon patch [m]
%         moon_patch - patch for the moon's surface

%Create black background patch
bg_patch = patch(bg(1,:),bg(2,:),'k');

axis off %turn off axis
hold on %keep plot on so more can be added or edited later
%Create moon patch
%Generate x-coordinates of moon
xmoon = linspace(0,2000,2001); %[m]
%Generate y-coordinates of moon
[ymoon] = genymoon; %[m]
moon_patch = patch(xmoon,ymoon,'w');

%Use nested function to create ymoon
function [ymoon] = genymoon
%Use linspace function to create linear segments of moon surface
ymoon = 100*[linspace(0,2,201),linspace(2+.01,5,100),...
         linspace(5+.01,7,100),...
         linspace(7,7,100),linspace(7+.01,6,100),...
         linspace(6+.01,7,100),linspace(7+.01,10,100),...
         linspace(10+.01,10,100),linspace(10+.01,9,100),...
         linspace(9+.01,6,100),linspace(6+.01,5,100),...
         linspace(5+.01,5,100),linspace(5+.01,6,100),...
         linspace(6+.01,11,100),linspace(11+.01,13,100),...
         linspace(13+.01,13,100),linspace(13+.01,10,100),...
         linspace(10+.01,9,100),linspace(9+.01,0,100),]; %[m]
end
end

function [ship,ship_patch,hship] = create_ship(x,y)
%Creates a ship patch by creating x and y coordinates, conctenating them, 
%and creating the patch.  Also creates the hship variable, which allows the
%ship to land when the bottom of the patch touches the moon surface, not
%the middle of the patch

%Inputs:  x          - x-position of ship [m]
%         y          - y-position of ship [m]
%Outputs: ship       - array containing x and y coordinates that make up
%                      ship patch [m]
%         ship_patch - patch for ship
%         hship      - height of the ship (distance between center of ship
%                      and bottom of ship) [m]

%xcords:
xship = 1/10*[0,-1,-2,-2,-3,-3,-4,-1,-2,-2,   2, 2, 1, 4, 3, 3, 2, 2, 1];
%ycords:
yship = 1/10*[4, 4, 3, 1, 0,-3,-4,-4,-3,-2,  -2,-3,-4,-4,-3, 0, 1, 3, 4];

hship = 40; %distance between center and bottom of ship [m]
wship = hship; %distance between center and side of ship [m]

ship = 100*[xship;yship]; %concatenate x and y cords and scale it
%use patch function to create ship patch
ship_patch = patch(ship(1,:)+x,ship(2,:)+y,'b'); %ship "spawns" at
%top-center of the plot and is blue
end

function text_boxes
%Creates text boxes and displays them in the figure window, next to the
%plot using the UI Controls

%Inputs: none
%Outputs: none

%Create Start button using UIControls
uicontrol('Style','pushbutton','String','Start Mission!',...
          'Position',[325,710,150,30],'Callback',@Start,'FontSize',14)

%Title of Game
TitleBox = uicontrol('Style','text','String','***LUNAR***LANDER***',...
           'Position',[100,750,600,50],'FontSize',30);
%Mass of fuel left
FuelBox = uicontrol('Style','text','String','Fuel left:',...
          'Position',[0,650,100,50],'FontSize',14);
%Units of mass of fuel
FuelunitsBox = uicontrol('Style','text','String','kg',...
               'Position',[0,600,100,50],'FontSize',12);
%Altitute
AltBox = uicontrol('Style','text','String','Altitude:',...
         'Position',[0,525,100,50],'FontSize',14);
%Units of altitude
AltunitsBox = uicontrol('Style','text','String','m',...
               'Position',[0,475,100,50],'FontSize',12);
%Velocity in x direction
VxBox = uicontrol('Style','text','String','X-velocity:',...
        'Position',[0,400,100,50],'FontSize',14);
%Units of velocity in x direction
VxunitsBox = uicontrol('Style','text','String','m/s',...
             'Position',[0,350,100,50],'FontSize',12);
%Velocity in y direction
VyBox = uicontrol('Style','text','String','Y-velocity:',...
        'Position',[0,275,100,50],'FontSize',14);
%Units of velocity in y direction
Vyunits = uicontrol('Style','text','String','m/s',...
          'Position',[0,225,100,50],'FontSize',12);
%Thrust magnitude
TBox = uicontrol('Style','text','String','Thrust:',...
        'Position',[0,150,100,50],'FontSize',14);
%Units of thrust magnitude
Tunits = uicontrol('Style','text','String','N',...
          'Position',[0,100,100,50],'FontSize',12);
end

function [FuelnumBox,AltnumBox,VxnumBox,VynumBox,TnumBox] = ...
          editable_boxes(mfuel,y,Vx,Vy,T)
%Creates and displays editable text boxes showing important values needed
%to play the game.  They are stored as variables so they can be edited
%later on in the program

%Inputs:  mfuel      - mass of fuel                    [kg]
%         y          - y position of ship              [m]
%         Vx         - velocity in x-direction of ship [m/s]
%         Vy         - velocity in y-direction of ship [m/s]
%         T          - magnitude of thrust             [N]
%Outputs: FuelnumBox - UI text box displaying value of mfuel
%         AltnumBox  - UI text box displaying value of y
%         VxnumBox   - UI text box displaying value of Vx
%         VynumBox   - UI text box displaying value of Vy
%         TnumBox    - UI text box displaying value of T

%Create and display editable boxes showing important values using UIControl
%Mass of fuel
FuelnumBox = uicontrol('Style','text','String',num2str(mfuel,'%7.2f'),...
             'Position',[0,650,100,25],'FontSize',12);
%Altitude
AltnumBox = uicontrol('Style','text','String',num2str(y,'%6.3f'),...
            'Position',[0,525,100,25],'FontSize',12);
%Velocity in x direction
VxnumBox = uicontrol('Style','text','String',num2str(Vx,'%5.3f'),...
           'Position',[0,400,100,25],'FontSize',12);
%Velocity in y direction
VynumBox = uicontrol('Style','text','String',num2str(Vy,'%5.3f'),...
           'Position',[0,275,100,25],'FontSize',12);
%Thrust magnitude
TnumBox = uicontrol('Style','text','String',num2str(T,'%7.2f'),...
           'Position',[0,150,100,25],'FontSize',12);
end

function Start(pushbutton, ~)
%Allows user to change start variable value by pressing the "Start Mission!"
%button

%Inputs: pushbutton - source of push button event
%        ~          - I'm not sure what this means.  Maybe it is the push
%                     button event?
%Outputs: none

%Give access to global variables
global start

%Change start value
start = 1;
end

function [h_moon] = calc_h_moon(x,ymoon)
%Calculates the height of the moon directly below the ship

%Inputs:  x      - x-position [m]
%         ymoon  - array of y-coordinates for moon patch      [m]
%Outputs: h_moon - y coordinate of moon directly beneath ship [m]

%Calculate counter, n, by rounding x value and using scale variable
n = round(x);
%Set h_moon equal to the y-value of moon corresponding to ship's x-value
h_moon = ymoon(n); %[m]
end

function [Ax, Ay] = LanderAccel(x, y, Vx, Vy, Q, mtot, T)
% Lunar Lander Acceleration Calculator
% Calculates net acceleration of lunar lander module (LM) for final project
% See Lunar Physics Reference document for definition of coordinate system
% Input:  x    - current horiz. position of LM              {m}
%         y    - current vert. position of LM (altitude)    {m}
%         Vx   - current horiz. velocity of LM              {m/s}
%         Vy   - current vert. velocity of LM               {m/s}
%         Q    - current orientation angle of LM
%                (where 0° = vertical & Q is positive clockwise  {degrees}
%         m    - current mass of LM                         {kg}
%         T    - current thrust output of LM rocket         {N}
% Output: Ax   - current horiz. acceleration of LM          {m/s^2}
%         Ay   - current vert. acceleration of LM           {m/s^2}
% Const:  gm   - lunar gravity                              {m/s^2}
%         Rm   - radius of moon                             {m}

    % Initialize constants
    gm = 1.62;      % {m/s^2}
    Rm = 1.738E6;   % {m}

    % Calculate forces on lander
    % Weight
    Wy = mtot * (-gm);         % {N}
    % Thrust, vert. & horiz. components
    Ty = T * cosd(-Q);       % {N}
    Tx = T * sind(-Q);       % {N}
    % 'Centrifugal force' to correct for orbital mechanics
    Cy = mtot * Vx^2 / (Rm+y); % {N}

    % Calculate horiz. & vert. acceleration of lander
    Ax = 100* Tx / mtot;            % {m/s^2}
    Ay = 100* (Wy+Ty+Cy) / mtot;     % {m/s^2}
    %Note: I scaled these values to make them work better with the values
    %that I was already using

end

function [x,y,Vx,Vy] = kinematics(x,y,Vx,Vy,Ax,Ay,dt)
%Calculates the position and velocity vectors of the ship using Euler's
%method

%Inputs:  x  - x-position of ship                  [m]
%         y  - y-position of ship                  [m]
%         Vx - velocity in x-direction of ship     [m/s]
%         Vy - velocity in y-direction of ship     [m/s]
%         Ax - acceleration in x-direction of ship [m/s^2]
%         Ay - acceleration in y-direction of ship [m/s^2]
%         dt - change in time                      [s]
%Outputs: x  - x-position of ship                  [m]
%         y  - y-position of ship                  [m]
%         Vx - velocity in x-direction of ship     [m/s]
%         Vy - velocity in y-direction of ship     [m/s]

%Calc position vectors using kinematics formulas and Euler's method
x = x+Vx*dt+1/2*Ax*dt^2; %[m]
y = y+Vy*dt+1/2*Ay*dt^2; %[m]
    
%Calc velocity vectors using kinematics formulas and Euler's method
Vx = Vx*dt+1/2*Ax*dt^2; %[m/s]
Vy = Vy*dt+1/2*Ay*dt^2; %[m/s]
end

function [point] = evaluate(Vy,mfuel,x,Q)
%Uses conditional operators to evaluate how the player has done.
%Depeneding on the conditions of the ship's landing, a message is displayed
%in the figure window using UI Controls

%Inputs:  Vy    - velocity in y-direction of ship [m/s]
%         mfuel - mass of fuel                    [kg]
%         x     - x-position of ship              [m]
%         Q     - angle from vertical of ship to vertical of plot [deg]
%Outputs: none

%Evaluate player's conditions and output message in text box
if (x>400 && x<500) || (x>800 && x<900)...
        || (x>1200 && x<1300) || (x>1600 && x<1700)
    %Player has landed ship on one of four flat surfaces 
    if abs(Q)<=2 && Vy>-5*10^(-1)
        %Player has landed ship oriented upright, with low velocity
        %Player has met all requirements to win, display message
        uicontrol('Style','text','String','Congratulations! You have successfully landed on the moon!',...
        'Position',[100,50,600,50],'FontSize',12);
        uicontrol('Style','text','String','Good luck getting back home with only',...
        'Position',[95,25,400,50],'FontSize',12);
        uicontrol('Style','text','String',num2str(mfuel,'%7.2f'),...
        'Position',[435,25,60,50],'FontSize',12);
        uicontrol('Style','text','String','kg of fuel left! XD',...
        'Position',[495,25,150,50],'FontSize',12);
        %Add to point counter, point
        point = 1;
    else
        %Player has not met requirements for a smooth landing, display
        %message
        uicontrol('Style','text','String','That landing was too harsh, your ship got destroyed! Try again!',...
        'Position',[100,50,600,50],'FontSize',12);
        %Keep the same point value
        point = 0;
    end
else
    %Player has not landed on one of four flat surfaces, display message
    uicontrol('Style','text','String','Your ship cannot land on a slope that steep! Try again!',...
    'Position',[100,50,600,50],'FontSize',12);
    %Keep the same point value
    point = 0;
end
end

function tryagain_button
%Create Try Again button using UIControls
uicontrol('Style','pushbutton','String','Try Again!',...
          'Position',[325,710,150,30],'Callback',@Tryagain,'FontSize',14)
end

function Tryagain(pushbutton, ~,point)
%Allows user to change start variable value by pressing the "Start Mission!"
%button

%Inputs: pushbutton - source of push button event
%        ~          - I'm not sure what this means.  Maybe it is the push
%                     button event?
%Outputs: none

%Give access to global variables
global tryagain

%Change tryagain value
tryagain = 1;

end

function KeyPress(src,event)
%Takes mouse input and allows user to control ship by changing variable
%values depending on the mouse input

%Inputs: src   - source of the keypress event
%        event - the keypress event
%Outputs: none

%Give access to global variables
global x y Q T mfuel dT Tmax

%Reads the name of the pressed key
keyname = event.Key;

%Check which key was pressed and respond
switch keyname
    case 'rightarrow'
        %Rotate right by decreasing angle
        Q = Q-2; %[deg] 
    case 'leftarrow'
        %Rotate left by increasing angle
        Q = Q+2; %[deg]
    case 'uparrow'
        %Increase Thrust if Thrust is below max value & there is fuel left
        if mfuel > 0 && T < Tmax
            T = T+dT; %[N]
        end
    case 'downarrow'
        %Decrease Thrust if Thrust is not already 0
        if T > 0
            T = T-dT; %[N]
        end
end
end