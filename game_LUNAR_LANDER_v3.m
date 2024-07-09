%FINAL PROJECT: LUNAR LANDER
%clear workspace
clear
clc

%GLOBALIZE
global x y Q T mfuel m Tmax dT start tryagain flength s_sel selected choose
global right left button_on mars choose_exp main main_again
%(Variables required for keyboard or mouse input)

%INITIALIZE
right = 0;          %Variable representing whether the go_right pushbutton 
                    %has been pressed
left = 0;           %Variable representing whether the go_left pushbutton
                    %has been pressed
start = 0;          %variable representing whether the start button has
                    %been pushed.  = 0 means it hasn't, = 1 means it has.
point = 0;          %variable that counts player's points
tryagain = 0;       %variable that acts like the start variable but keeps
                    %track of whether the "Try Again!" button has been
                    %pressed
tabutton = 0;       % = 0 if the "Try Again!" button has not been spawned
                    % = 1 if it has
Vx = 0;             %velocity in the x-direction of the ship [m/s]
Vy = 0;             %velocity in the y-direction of the ship [m/s]
x0 = 1000;          %initial x-position of the ship    [m]
y0 = 2060;          %initial y-position of the ship    [m]
x = x0;             %x-position of the ship    [m]
y = y0;             %y-position of the ship    [m]
t = 0;              %time                      [s]
dt = 0.1;           %change in time            [s]
h_moon = 0;         %y-coordinate of the moon directly beneath the ship [m]
Q = 0;              %angle between ship's vertical & plot's vertical [deg]
T = 0;              %thurst acting on the ship [N]
Tmax = 45000;       %maximum thrust that can act on the ship[N]
dT = Tmax/10;       %change in thrust          [N]
numtries = 0;       %number of tries player has had
no_ship = 0;        %Represents whether explosion has happened
choose_exp = 0;     %Shows whether expedition has been chosen

%% TITLE SCREEN
%Create figure window using function
[fig1] = create_fig;

%Create ACTUAL Title screen with functions
%Create axes limits
xmax = 2000;
ymax = 2000;
%Create axes for Title screen (same as for gameplay)
axis([0,xmax,0,ymax]);
axis manual
axis off
tscale = 50; %Initialize tscale, scale variable used for title screen

pause(0.5)

%Create LUNAR
[L1_p,U_p,N1_p,St1_p,R1_p,P] = LUNAR(tscale);

%Create LANDER!
[L2_p,St2_p,N2_p,D_p,E_p,R2_p] = LANDER(tscale,P); 
 
%Create titleship
[tship_p,OFFx,OFFy,Q,r,ax,bx,cx,dx,tshipy,s_scale,tscale] =...
    create_tship(tscale);

%Turn tship
[tship_p,tship,tshiprotate,Q,OFFx,OFFy] =...
    turn_tship(s_scale,OFFx,OFFy,r,Q,ax,bx,cx,dx,tshipy,tship_p);

%Movement of tship and gradual scale of exclamation mark (fire trail)
[tship_p,expm_p] = zoom(tship_p,tship,tshiprotate,Q,OFFx,OFFy,tscale);

pause(1) %briefly to appreciate beautiful title screen animation

%Delete ACTUAL title screen
deleteLUNAR(L1_p,U_p,N1_p,St1_p,R1_p)
deleteLANDER(L2_p,St2_p,N2_p,D_p,E_p,R2_p)
delete_tship_expm(tship_p,expm_p)

%Reset Q value because I was lazy and used it for the ACTUAL title screen
Q = 0; %[deg]

%% SHIP SELECTION SCREEN
%Create infinite delay loop to allow user to remain on Title screen until
%pressing the "Ready for Takeoff?" button

%Register KeyPressFcn to allow figure window to respond to keypresses
set(fig1, 'KeyPressFcn', @KeyPress);

%Create "Gear Up!" button
Gearbutton = uicontrol('Style','pushbutton','String','Gear Up!',...
    'Position',[325,740,150,30],'Callback',@Start,'FontSize',14);

%Create variable representing condition of being on ship selection screen
s_sel = 1;
%Create variable to represent whether a ship has been selected
choose = 0;
%Created variable to represent the selected ship
selected = 1;

%Create buttons for choosing a ship
[goleft_butt] = goleft_button;
[goright_butt] = goright_button;

%Plot first ship first
[ship,ship_patch,hship,wship] = create_ship(x,y,s_sel,xmax,ymax);
ship_name = uicontrol('Style','text','String','Classic Lander',...
    'Position',[300,180,200,50],'FontSize',18);

%Create text box for ship selection screen
ChooseBox = uicontrol('Style','text','String','CHOOSE YOUR SHIP:',...
    'Position',[200,600,400,50],'FontSize',28);

%Initialize
P = 0.00575; %pause variable
nmax_shift = 15; %number of iterations
button_on = 1;
dshift = (xmax-500)/nmax_shift; %shift size - this is a dependent
%variable, DON'T CHANGE IT
%Scrolls ships left and right
while choose ~= 1
    if selected == 1
        if    right == 1
            %Display image and name of current ship
            [ship,ship_patch,hship,wship] =...
                create_ship(x,y,s_sel,xmax,ymax,left,right);
            ship_name = uicontrol('Style','text','String',...
            'Classic Lander','Position',[300,180,200,50],'FontSize',18);
            %"Scroll"
            button_on = 0; %turn button off
            Rshift = 0;
            for n = 1:nmax_shift
                Rshift = Rshift-dshift;
                set(s3_patch,'XData',s3(1,:)+xmax/2+Rshift,...
                             'YData',s3(2,:)+ymax/2);
                set(ship_patch,'XData',ship(1,:)+xmax+500+Rshift,...
                               'YData',ship(2,:)+ymax/2);
                pause(P)
            end
            button_on = 1; %turn button on
            %Delete image and name of previous ship
            delete(s3_patch)
            delete(s3_name)
            %Reset right value
            right = 0;
        elseif left == 1
            %Display image and name of current ship
            [ship,ship_patch,hship,wship] =...
                create_ship(x,y,s_sel,xmax,ymax,left,right);
            ship_name = uicontrol('Style','text','String',...
            'Classic Lander','Position',[300,180,200,50],'FontSize',18);
            %"Scroll" away previous ship
            button_on = 0; %turn button off
            Lshift = 0; %reset Lshift
            for n = 1:nmax_shift
                Lshift = Lshift+dshift;
                set(s2_patch,'XData',s2(1,:)+xmax/2+Lshift,...
                             'YData',s2(2,:)+ymax/2);
                set(ship_patch,'XData',ship(1,:)-500+Lshift,...
                             'YData',ship(2,:)+ymax/2);
                pause(P)
            end
            button_on = 1; %turn button on
            %Delete image and name of previous ship
            delete(s2_patch)
            delete(s2_name)
            %Reset left value
            left = 0;
        else
            %Do not change anything if player has not pressed button
        end
    elseif selected == 2
        if    right == 1
            %Display image and name of current ship
            [s2,s2_patch,hs2,ws2] = create_s2(x,y,s_sel,xmax,ymax);
            s2_name = uicontrol('Style','text','String',...
            'Wrecking Ball','Position',[300,180,200,50],'FontSize',18);
            %"Scroll" away previous ship
            button_on = 0; %turn button off
            Rshift = 0; %reset Rshift
            for n = 1:nmax_shift
                Rshift = Rshift-dshift;
                set(ship_patch,'XData',ship(1,:)+xmax/2+Rshift,...
                               'YData',ship(2,:)+ymax/2);
                set(s2_patch,'XData',s2(1,:)+xmax+500+Rshift,...
                             'YData',s2(2,:)+ymax/2);
                pause(P)
            end
            button_on = 1; %turn button on
            %Delete image and name of previous ship
            delete(ship_patch)
            delete(ship_name)
            %Reset right value
            right = 0;
        elseif left == 1
            %Display image and name of current ship
            [s2,s2_patch,hs2,ws2] = create_s2(x,y,s_sel,xmax,ymax);
            s2_name = uicontrol('Style','text','String',...
            'Wrecking Ball','Position',[300,180,200,50],'FontSize',18);
            %"Scroll"
            button_on = 0; %reset button_on
            Lshift = 0; %reset Lshift
            for n = 1:nmax_shift
                Lshift = Lshift+dshift;
                set(s3_patch,'XData',s3(1,:)+xmax/2+Lshift,...
                             'YData',s3(2,:)+ymax/2);
                set(s2_patch,'XData',s2(1,:)-500+Lshift,...
                             'YData',s2(2,:)+ymax/2);
                pause(P)
            end
            button_on = 1; %turn button on
            %Delete image and name of previous ship
            delete(s3_patch)
            delete(s3_name)
            %Reset left value
            left = 0;
        else
            %Do not change anything if player has not pressed button
        end
    elseif selected == 3
        if    right == 1
            %Display image and name of current ship
            [s3,s3_patch,hs3,ws3] = create_s3(x,y,s_sel,xmax,ymax);
            s3_name = uicontrol('Style','text','String',...
            'Space Missile','Position',[300,180,200,50],'FontSize',18);
            %"Scroll" away previous ship
            button_on = 0; %turn button off
            Rshift = 0; %reset Rshift
            for n = 1:nmax_shift
                Rshift = Rshift-dshift;
                set(s2_patch,'XData',s2(1,:)+xmax/2+Rshift,...
                             'YData',s2(2,:)+ymax/2);
                set(s3_patch,'XData',s3(1,:)+xmax+500+Rshift,...
                             'YData',s3(2,:)+ymax/2);
                pause(P)
            end
            button_on = 1; %turn button on
            %Delete image and name of previous ship
            delete(s2_patch)
            delete(s2_name)
            %Reset right value
            right = 0;
        elseif left == 1
            %Display image and name of current ship
            [s3,s3_patch,hs3,ws3] = create_s3(x,y,s_sel,xmax,ymax);
            s3_name = uicontrol('Style','text','String',...
            'Space Missile','Position',[300,180,200,50],'FontSize',18);
            %"Scroll"
            button_on = 0; %turn button off
            Lshift = 0; %reset Lshift
            for n = 1:nmax_shift
                Lshift = Lshift+dshift;
                set(ship_patch,'XData',ship(1,:)+xmax/2+Lshift,...
                               'YData',ship(2,:)+ymax/2);
                set(s3_patch,'XData',s3(1,:)-500+Lshift,...
                             'YData',s3(2,:)+ymax/2);
                pause(P)
            end
            button_on = 1; %turn button on
            %Delete image and name of previous ship
            delete(ship_patch)
            delete(ship_name)
            %Reset left value
            left = 0;
        else
            %Do not change anything if player has not pressed button
        end
    end
pause(.1)
end

%Delete patch and name textbox for most recently selected ship
if selected == 1
    delete(ship_patch)
    delete(ship_name)
elseif selected == 2
    delete(s2_patch)
    delete(s2_name)
elseif selected == 3
    delete(s3_patch)
    delete(s3_name)
end

%Delete text box and buttons for ship selection screen
delete(ChooseBox)
delete(Gearbutton)
delete(goleft_butt);
delete(goright_butt);
%Reset start variable value so function handle @start can be reused
start = 0;

%% EXPEDITIONS SCREEN
%Create text boxes and buttons for expd screen using function
[ChooseExpBox,MoonButt,MarsButt,Mooninfo1,Mooninfo2a,...
 Mooninfo2b,Mooninfo3a,Mooninfo3b,Marsinfo1a,Marsinfo1b,...
 Marsinfo1c,Marsinfo2a,Marsinfo2b,Marsinfo3a,Marsinfo3b,...
 Marsinfo3c] = expd_boxes;

%Create patches for expd screen using function
[expd_bg1_p,expd_bg2_p,tmoon_p,tmars_p] = planetpatches;

%delay loop until player chooses expedition
while choose_exp ~= 1
    pause(0.1)
end

%Delete everything from expd screen using function
delete_expd(ChooseExpBox,MoonButt,MarsButt,Mooninfo1,Mooninfo2a,...
            Mooninfo2b,Mooninfo3a,Mooninfo3b,Marsinfo1a,Marsinfo1b,...
            Marsinfo1c,Marsinfo2a,Marsinfo2b,Marsinfo3a,Marsinfo3b,...
            Marsinfo3c,expd_bg1_p,expd_bg2_p,tmoon_p,tmars_p);

%Adjust physics and fuel values depending on chosen expedition
[gm,Rm,mship,mfuel,mtot,mfuel0] = setplanet(mars);

%% SETUP

%Change s_sel so create-ship functions function normally
s_sel = 0;

%Create Title screen with explanation of game controls and mission
[Readybutton,Mission,CONTROLS,ctr1,ctr2,ctr3,ctr4,GL] = tscreen;

%Create infinite delay loop to allow user to remain on Title screen until
%pressing the "Ready for Takeoff?" button
while start ~= 1
    pause(0.1);
end
start = 0;
%Reset start variable value so function handle @start can be reused

%Delete Title screen text boxes
deletetscreen(Readybutton,Mission,CONTROLS,ctr1,ctr2,ctr3,ctr4,GL);

%Create plot using function
[bg] = create_axis;

%Create "setting" using function
[bg_patch,xmoon,ymoon,moon_patch] = create_setting(bg,mars);

%Allow coordinates for selected ship to be used in game
if     selected == 1
    [ship,ship_patch,hship,wship] = create_ship(x,y,s_sel,xmax,ymax);
elseif selected == 2
    [s2,s2_patch,hs2,ws2] = create_s2(x,y,s_sel,xmax,ymax);
    ship = s2;
    ship_patch = s2_patch;
    hship = hs2;
    wship = ws2;
elseif selected == 3
    [s3,s3_patch,hs3,ws3] = create_s3(x,y,s_sel,xmax,ymax);
    ship = s3;
    ship_patch = s3_patch;
    hship = hs3;
    wship = ws3;
end

%Create fire patch coords using function
    [f,fp1,fp2,flength,rf1,rf2,tf1,tf2,fx]...
    = create_firepatch(x,y,selected);

%Create coords for explosion animation using functions
[expl_org] = create_expl;
[expl2_org] = create_expl2;
[expl3a_org,expl3b_org,expl3c_org] = create_expl3;

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

if numtries > 0
    %Respawn fire patches
    fy = [0,flength,0];
    f = 10*[fx;fy];
    fp1 = patch(f(1,:)+x+rf1*cosd(tf1),...
                f(2,:)+y+rf1*sind(tf1),'y');
    fp2 = patch(f(1,:)+x+rf2*cosd(tf2),...
                f(2,:)+y+rf2*sind(tf2),'y');
end

while y-hship > h_moon %game plays while base of ship is above moon surface
                       %and above 0 (if off screen)
    
    %Calculate acceleration, velocity, and position using function
    [Ax,Ay,x,y,Vx,Vy] = LanderAccel(x,y,Vx,Vy,Q,mtot,T,gm,Rm,dt);
    
    %Delete "old" fire patch so updated one can take its place
    delete(fp1)
    delete(fp2)
    
    %THIS NEEDS TO BECOME A FUNCTION
    %Update size of fire using function (not using function yet)
    if mfuel > 0
        %Update y-cords with new flength value
        fy = [0,flength,  0];
    else
        %No fuel means no fire
        flength = 0;
        fy = [0,flength,0];
    end
    %Concatenate x and y coordinates
    f = 10*[fx;fy];
    %Use patch function to create patch
    fp1 = patch(f(1,:)+x+rf1*cosd(tf1),...
            f(2,:)+y+rf1*sind(tf1),'y');
    fp2 = patch(f(1,:)+x+rf2*cosd(tf2),...
            f(2,:)+y+rf2*sind(tf2),'y'); 
    
    %Update position of ship and fire
    set(ship_patch,'XData',ship(1,:)+x,...
                   'YData',ship(2,:)+y);
    set(fp1,'XData',f(1,:)+x+rf1*cosd(tf1+Q),...
            'YData',f(2,:)+y+rf1*sind(tf1+Q));
    set(fp2,'XData',f(1,:)+x+rf2*cosd(tf2+Q),...
            'YData',f(2,:)+y+rf2*sind(tf2+Q));
    
    %Create a rotation matrix to update ship's orientation
    RtnMtrx = [cosd(-Q), sind(-Q); -sind(-Q), cosd(-Q)];
    %Multiply the ship coordinates by the rotation matrix to reorient it
    shiprotate = RtnMtrx*ship;
    frotate = RtnMtrx*f;
    
    %Update orientation of ship and fire
    set(ship_patch,'XData',shiprotate(1,:)+x,'YData',shiprotate(2,:)+y);
    set(fp1,'XData',frotate(1,:)+x+rf1*cosd(tf1+Q),...
            'YData',frotate(2,:)+y+rf1*sind(tf1+Q))
    set(fp2,'Xdata',frotate(1,:)+x+rf2*cosd(tf2+Q),...
            'YData',frotate(2,:)+y+rf2*sind(tf2+Q))
    
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
[point] = evaluate(Vy,mfuel,x,Q,ymoon,mars);

%Set velocity and thrust values equal to zero now that program is finished
%(program doesn't update them outside of loop)
Vx = 0;
Vy = 0;
T = 0;
flength = 0;
%Update text boxes
%Velocity in x direction
VxnumBox.String = num2str(Vx,'%4.2f');
%Velocity in y direction
VynumBox.String = num2str(Vy,'%4.2f');
%Thrust magnitude
TnumBox.String = num2str(T,'%4.2f');
%Delete fire patch if it hasn't been already
delete(fp1)
delete(fp2)

%Create animation depending on whether player has won
if point ~= 1
%Player has lost, ship explodes
animate_expl(expl_org,expl2_org,expl3a_org,expl3b_org,expl3c_org,...
    x,y,ship_patch);
no_ship = 1; %change no_ship
else
%Player has won, american flag is planted  
pause(0.5)
[pole_p,flag_p,fbox_p,stripe1_p,stripe2_p,stripe3_p] =...
    flag(x,y,hship,wship);
end

%Add to numtries variable
numtries = numtries+1;
 
%Spawn try again button if conditions are met
if point ~= 1 && tabutton ~= 1
    %Function to spawn tryagain button
    tryagain_button
    %Change value of tabutton
    tabutton = 1;
end

while tryagain ~= 1 && point ~= 1 && main ~= 1
    %Reset variables
    x = x0;
    y = y0;
    mfuel = mfuel0;
    Q = 0;
    flength = 0;
   %Pause until user presses "Try Again!" button
    pause(0.1);
end

%Re-plot ship_patch if player has lost
if no_ship == 1
    %Use ship selection function to plot whichever ship was selected
    if     selected == 1
        [ship,ship_patch,hship] = create_ship(x,y,s_sel,xmax,ymax);
    elseif selected == 2
        [s2,s2_patch,hs2] = create_s2(x,y,s_sel,xmax,ymax);
        ship = s2;
        ship_patch = s2_patch;
        hship = hs2;
    elseif selected == 3
        [s3,s3_patch,hs3] = create_s3(x,y,s_sel,xmax,ymax);
        ship = s3;
        ship_patch = s3_patch;
        hship = hs3;
    end
    no_ship = 0; %change no_ship
end

%Reorient ship and place at starting point
RtnMtrx = [cosd(-Q), sind(-Q); -sind(-Q), cosd(-Q)];
shiprotate = RtnMtrx*ship;
set(ship_patch,'XData',shiprotate(1,:)+x,'YData',shiprotate(2,:)+y);
end

%Reset tryagain variable
tryagain = 0;

%% TITLE SCREEN FUNCTIONS

function [L1_p,U_p,N1_p,St1_p,R1_p,P] = LUNAR(tscale)
P = 0.125; %Pause variable
L1 = tscale*[ 1, 2, 6, 6, 3, 3;
            28,20,23,25,22,29];   
L1_p = patch(L1(1,:),L1(2,:),'y');
pause(P)
U = tscale*[ 6,8.5,10,12,12,11, 9,11,11,10, 9,7.5;
           30, 25,24,26,28,31,31,27,26,25,26, 30];
U_p = patch(U(1,:),U(2,:),'y');
pause(P)
N1 = tscale*[14,15,15,19,21,19,18.5,15;
            24,24,28,23,31,32,  25,31];
N1_p = patch(N1(1,:),N1(2,:),'y');
pause(P)
St1 = tscale*[21,24,25,26,28,  26,26,24,22,23;
             23,25,22,26,27,27.5,31,28,29,26];
St1_p = patch(St1(1,:),St1(2,:),'y');
pause(P)
R1 = tscale*[27,29,29,31,32,30,29,30,31,31,29,30,32,32,30;
            22,21,25,21,22,25,26,28,28,26,26,25,25,29,30];
R1_p = patch(R1(1,:),R1(2,:),'y');
pause(P)
end

function [L2_p,St2_p,N2_p,D_p,E_p,R2_p] = LANDER(tscale,P)
L2 = tscale*[ 6, 9,8,6.5, 5, 4;
             12,14,16, 13,19,18];
L2_p = patch(L2(1,:),L2(2,:),'y');
pause(P)
St2 = tscale*[10,12,14,13,15,13,12,11,09,11;
              14,16,14,17,19,19,21,19,19,17];
St2_p = patch(St2(1,:),St2(2,:),'y');   
pause(P)
N2 = tscale*[15,16,17,18,20,18,18,17;
             14,14,17,14,19,19,16,19];
N2_p = patch(N2(1,:),N2(2,:),'y');
pause(P)
D = tscale*[22,20,19,21,20,22,23,21,19,19,21,24,22;
            20,15,13,14,14,19,18,14,13,13,13,18,20];
D_p = patch(D(1,:),D(2,:),'y');    
pause(P)
E = tscale*[22,25,25,23,24,26,27,24.5,25,27,28,25;
            12,12,13,13,15,15,16,  16,17,17,18,18];
E_p = patch(E(1,:),E(2,:),'y');
pause(P)
R2 = tscale*[26,27,28,30,31,29,28,29,30,30,28,29,31,31,29;
             11,10,13,10,11,13,14,16,16,14,14,13,13,16,17];
R2_p = patch(R2(1,:),R2(2,:),'y');
pause(P)
end

function [tship_p,OFFx,OFFy,Q,r,ax,bx,cx,dx,tshipy,s_scale,tscale] =...
    create_tship(tscale)
%Offset values must be used AFTER coordinates are made (around (0,0) so the
%varying scaling works)
OFFx = (37.5)*tscale;
OFFy = 24*tscale;

%Create coords of ship
[tship,tshipx,tshipy,s_scale] = tshipcords(tscale,OFFx,OFFy);
%Create tship_p
tship_p = patch(tshipx+OFFx,tshipy+OFFy,'b');

%Make ship "move"
Q = 0; %Initialize starting angle [deg]
%"save" original coords for tshipx and tshipy
tshipx_org = tshipx;
tshipy_org = tshipy;
 tship_org = tship;
for n = 1:15
    Q = n*255/15; %Change Q by increment
    ratio = n/15; %Change ratio by increment
    for n = 1:length(tshipx)
        tshipx(n) = ratio*tshipx_org(n);
        tshipy(n) = ratio*tshipy_org(n);
    end
    %Concatenate x and y coords
    tship = [tshipx;tshipy];
    
    %Update tship size
    set(tship_p,'XData',tship(1,:)+OFFx,'YData',tship(2,:)+OFFy);
    %Create a rotation matrix to update ship's orientation
    
    %Create rotation matrix to rotate tship
    RtnMtrx = [cosd(-Q), sind(-Q); -sind(-Q), cosd(-Q)];
    %Multiply the ship coordinates by the rotation matrix to reorient it
    tshiprotate = RtnMtrx*tship;
    
    %Update orientation of tship
    set(tship_p,'XData',tshiprotate(1,:)+OFFx,...
                'YData',tshiprotate(2,:)+OFFy);

    pause(0.05)
end

    %Creates coordinates for tship
    function [tship,tshipx,tshipy,s_scale] = tshipcords(tscale,OFFx,OFFy)
            r = sqrt(29);
    %Create coordinates for circle
    %A- curve
    thmin = -atand(5/2);
    thmax = -thmin;
    dth = (thmax-thmin)/11.2; %this 11.2 was originally 10... it works for
                              %now, the issue was the curve C not making it
                              %all the way around.  Not sure why, but this
                              %seemed to fix it. May cause future issues...
                              %we'll see.
    nmax = (thmax-thmin)/dth+1;
    for n = 1:nmax
        th(n) = thmin+(n-1)*dth;
        ax(n) = r*cosd(th(n));
        ay(n) = r*sind(th(n));
    end
    %B- curve
    bx = [ 2, 2,-2,-2];
    by = [ 5, 7, 7, 5];
    %C- curve
    thmin = thmin+180;
    thmax = thmax+180;
    dth = (thmax-thmin)/dth+1;
    nmax - (thmax-thmin)/dth+1;
    for n = 1:nmax
        th(n) = thmin+(n-1)*dth;
        cx(n) = r*cosd(th(n));
        cy(n) = r*sind(th(n));
    end    
    %D- curve
    dx = -bx;
    dy = -by;

    %Concatenate curves
    tshipx = [ax,bx,cx,dx];
    tshipy = [ay,by,cy,dy];

    %Scale
    s_scale = tscale/sqrt(3.5^2+2.5^2);
    for n = 1:length(tshipx)
        tshipx(n) = s_scale*tshipx(n);
        tshipy(n) = s_scale*tshipy(n);
    end

    %Concatenate x and y coords
    tship = [tshipx;tshipy];
    end
end

function [tship_p,tship,tshiprotate,Q,OFFx,OFFy] =...
    turn_tship(s_scale,OFFx,OFFy,r,Q,ax,bx,cx,dx,tshipy,tship_p)
%Create original value vectors (setup for 1st part)
ax_org = ax*s_scale;
cx_org = cx*s_scale;

%Create newtship coords (setup for 2nd part)
[newtshipx,newtshipy,newtship] = newcoords(s_scale);

%Create limits for turnR, ratio for the turning-of-the-ship process
for n = 1:16
    if n <= 8
        %Ship is compressing
    turnR = 11-n;
    %"Squish" ax and cx
    for n = 1:length(ax)
        ax(n) = ax_org(n)/10*(turnR+0);
        cx(n) = cx_org(n)/10*(turnR+0);
    end
    %Adjust bx and dx accordingly
    bx(1) = ax(length(ax)); %line up with last value in ax vector array
    bx(2) = bx(1); %horizontal shift to maintain shape
    bx(3) = bx(1)-4*s_scale; %horizontal shift to maintain shape
    bx(4) = bx(3);
    %dx values are same, but order is switched
    dx(1) = bx(4);
    dx(2) = bx(4);
    dx(3) = bx(1);
    dx(4) = bx(1);
    
    %Re-concatenate using function
    [tshipx,tship] = concatenate(ax,bx,cx,dx,tshipy);
    else
        turnR = n-5;
        %Ship is expanding
        tshipx = newtshipx*(turnR/11);
        tship = [tshipx;newtshipy];
    end
    
    %Reorient ship using Rotation Matrix
    RtnMtrx = [cosd(-Q), sind(-Q); -sind(-Q), cosd(-Q)];
    %Multiply the ship coordinates by the rotation matrix to reorient it
    tshiprotate = RtnMtrx*tship;

    %Update orientation of tship
    set(tship_p,'XData',tshiprotate(1,:)+OFFx,...
                'YData',tshiprotate(2,:)+OFFy);
    pause(.005)
end

    function [tshipx,tship] = concatenate(ax,bx,cx,dx,tshipy)
    %Concatenate curves
    tshipx = [ax,bx,cx,dx];

    %Concatenate x and y coords
    tship = [tshipx;tshipy];
    end
    function [newtshipx,newtshipy,newtship] = newcoords(s_scale)
        newtshipx = ...
        s_scale*[ 6, 4,-3,-4,-6,-7,-7,-6,-4,  -4,-6,-7,-7,-6,-4,-3, 4, 6];
        newtshipy = ...
        s_scale*[ 3, 5, 5, 6, 6, 7, 3, 4, 4,  -4,-4,-3,-7,-6,-6,-5,-5,-3];
        newtship = [newtshipx;newtshipy];
    end

end

function [tship_p,expm_p] =...
    zoom(tship_p,tship,tshiprotate,Q,OFFx,OFFy,tscale)
%Create expm coords
expmx_org = tscale*[0, -4.75, -1.75];
expmy_org = tscale*[0,-11,-11.75];
%expm = [expmx;expmy];

%Initialize slope
%slope_path = 7/2; %(1200-500)/(1875-1675)
dx_org = -2;
dy_org = -7;
for n = 1:100
    dx = dx_org*n;
    dy = dy_org*n;
    
    %Reorient ship using Rotation Matrix
    %(this must be done because the new ship coords are oriented
    %differently from the start; they must be re-oriented everytime you use
    %them
    RtnMtrx = [cosd(-Q), sind(-Q); -sind(-Q), cosd(-Q)];
    %Multiply the ship coordinates by the rotation matrix to reorient it
    tshiprotate = RtnMtrx*tship;
    
    %Update position of tship
    set(tship_p,'XData',tshiprotate(1,:)+dx+OFFx,...
                'YData',tshiprotate(2,:)+dy+OFFy);
    if n >= 10
    %Scale expm coords
    expmx = expmx_org*(n-9)/(100-9);
    expmy = expmy_org*(n-9)/(100-9);
    %Offset values must be used AFTER coordinates are made (around (0,0) so the
    %varying scaling works)
        if n == 10
        %Create expm patch
        expm_p = patch(expmx+OFFx,expmy+OFFy,'y');
        else
        %Update expm patch
        set(expm_p,'XData',expmx+OFFx,...
                   'YData',expmy+OFFy);
        end
    end        
    pause(.00575)
end

end

function deleteLUNAR(L1_p,U_p,N1_p,St1_p,R1_p)
delete(L1_p)
delete(U_p)
delete(N1_p)
delete(St1_p)
delete(R1_p)
end
function deleteLANDER(L2_p,St2_p,N2_p,D_p,E_p,R2_p)
delete(L2_p)
delete(St2_p)
delete(N2_p)
delete(D_p)
delete(E_p)
delete(R2_p)
end
function delete_tship_expm(tship_p,expm_p)
delete(tship_p)
delete(expm_p)
end

%% SHIP SELECTION SCREEN FUNCTIONS

function [goleft_butt] = goleft_button
%Create "<" button
goleft_butt = uicontrol('Style','pushbutton','String','<',...
          'Position',[220,200,50,30],'Callback',@go_left,'FontSize',24);
end
function [goright_butt] = goright_button
%Create ">" button
goright_butt = uicontrol('Style','pushbutton','String','>',...
          'Position',[530,200,50,30],'Callback',@go_right,'FontSize',24);
end
function go_right(pushbutton,~)
%Give access to globals
global selected right left button_on

if button_on == 1 %THIS IF STATEMENT DOESN'T SEEM TO HELP... :/
%Center ship to right
    if     selected == 1 %1 up to 2
        selected = 2;
    elseif selected == 2 %2 up to 3
        selected = 3;
    elseif selected == 3 %3 up to 1
        selected = 1;
    end
end
%Update right and left values
right = 1;
left = 0;
end
function go_left(pushbutton,~)
%Give access to globals
global selected right left button_on

%if button_on == 1
%Center ship to left
        if     selected == 1 %1 down to 3
            selected = 3;
        elseif selected == 2 %2 down to 1
            selected = 1;
        elseif selected == 3 %3 down to 2
            selected = 2;
        end
%end
%Update right and left values
left = 1;
right = 0;
end

function [ship,ship_patch,hship,wship] =...
    create_ship(x,y,s_sel,xmax,ymax,left,right)
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
xship = [0,-1,-2,-2,-3,-3,-4,-1,-2,-2,   2, 2, 1, 4, 3, 3, 2, 2, 1];
%ycords:
yship = [4, 4, 3, 1, 0,-3,-4,-4,-3,-2,  -2,-3,-4,-4,-3, 0, 1, 3, 4];
%Gap signifies middle of ship

hship = 40; %distance between center and bottom of ship [m]
wship = hship; %distance between center and side of ship [m]

if s_sel ~=1
    ship = 10*[xship;yship]; %concatenate x and y cords and scale it
    %use patch function to create ship patch
    ship_patch = patch(ship(1,:)+x,ship(2,:)+y,'b'); %ship "spawns" at
    %top-center of the plot and is blue
else
    ship = 100*[xship;yship];
    ship_patch = patch(ship(1,:)+xmax/2,ship(2,:)+ymax/2,'b'); %ship 
    %"spawns" at center of screen, is blue, and is enlarged
end

end
function [s2,s2_patch,hs2,ws2] =...
    create_s2(x,y,s_sel,xmax,ymax,left,right)
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

%Initialize
r = sqrt(13);
th1 = atand(2/3);
th2 = 90-2*th1;
thmin = -th1;
thmax = 180+1*th1;
dth = 5;
nmax = (thmax-thmin)/dth+1;
%Create coordinates that make up shape
%A- curve
for n = 1:nmax
    th(n) = thmin+dth*(n-1);
    xs2a(n) = r*cosd(th(n));
    ys2a(n) = r*sind(th(n));
end
%B- manual
xs2b = [-3,-3,-3.75,-3.75,-1.25,-1.25,-2,-2];
ys2b = [-2,-4,-4,-5,-5,-4,-4,-3];
%Reinitialize
thmin = 180+th1+th2;
thmax = thmin+2*th1;
nmax = (thmax-thmin)/dth+1;
%C- curve
for n = 1:nmax
    th(n) = thmin+dth*(n-1);
    xs2c(n) = r*cosd(th(n));
    ys2c(n) = r*sind(th(n));
end
%D- manual
xs2d = [ 2, 2, 1.25, 1.25, 3.75, 3.75, 3, 3];
ys2d = [-3,-4,-4,-5,-5,-4,-4,-2];
       
hs2 = 50; %distance between center and bottom of ship [m]
ws2 = 40; %distance between center and side of ship [m]

if s_sel ~= 1
    %Concatenate all four parts: A, B, C, D
    s2 = 10*[xs2a,xs2b,xs2c,xs2d;
           ys2a,ys2b,ys2c,ys2d];
    s2_patch = patch(s2(1,:)+x,s2(2,:)+y,'c'); %ship "spawns" at
    %top-center of the plot and is cyan
else
    s2 = 100*[xs2a,xs2b,xs2c,xs2d;...
              ys2a,ys2b,ys2c,ys2d];
    s2_patch = patch(s2(1,:)+xmax/2,s2(2,:)+ymax/2,'c'); %ship "spawns" at
    %middle of the plot, is cyan, and is enlarged
end
end
function [s3,s3_patch,hs3,ws3] =...
    create_s3(x,y,s_sel,xmax,ymax,left,right)
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
xs3 = [ 2, 2, 1, 1, 2, 2, 1, 1,  0,  -1,-1,-2,-2,-1,-1,-2,-2];
%ycords:
ys3 = [-5,-4,-3,-1,-1, 0, 1, 3,  5,   3, 1, 0,-1,-1,-3,-4,-5]; 
%Gap signifies middle of ship

hs3 = 50; %distance between center and bottom of ship [m]
ws3 = 25; %distance between center and side of ship [m]
%(it's actually 20, but the flag looks better when ws3=25)

if s_sel ~= 1
    s3 = 10*[xs3;ys3]; %concatenate x and y cords and scale it
    %use patch function to create ship patch
    s3_patch = patch(s3(1,:)+x,s3(2,:)+y,'r'); %ship "spawns" at
    %top-center of the plot and is red
else
    s3 = 100*[xs3;ys3]; %concatenate x and y cords and scale it
    %use patch function to create ship patch
    s3_patch = patch(s3(1,:)+xmax/2,s3(2,:)+ymax/2,'r'); %ship "spawns" at
    %center of the plot, is red, and is enlarged
end
end

%% EXPEDITIONS SCREEN FUNCTIONS

function chooseMoon(pushbutton,~)
%Give access to globals
global mars choose_exp
mars = 0; %shows that mars was not chosen
choose_exp = 1; %shows that expedition was chosen
end
function chooseMars(pushbutton,~)
%Give access to globals
global mars choose_exp
mars = 1; %shows that mars was chosen
choose_exp = 1; %shows that expedition was chosen
end

function [ChooseExpBox,MoonButt,MarsButt,Mooninfo1,Mooninfo2a,...
          Mooninfo2b,Mooninfo3a,Mooninfo3b,Marsinfo1a,Marsinfo1b,...
          Marsinfo1c,Marsinfo2a,Marsinfo2b,Marsinfo3a,Marsinfo3b,...
          Marsinfo3c] = expd_boxes
%Create text box
ChooseExpBox = uicontrol('Style','text','String','EXPEDITIONS:',...
    'Position',[200,700,400,50],'FontSize',28);
%Create moon option button
MoonButt = uicontrol('Style','pushbutton','String','The Moon',...
    'Position',[150,350,200,50],'Callback',@chooseMoon,'FontSize',20);
%Create mars option button
MarsButt = uicontrol('Style','pushbutton','String','Mars',...
    'Position',[450,350,200,50],'Callback',@chooseMars,'FontSize',20);
%Create expedition info textboxes
Mooninfo1 = uicontrol('Style','text','String',...
    'For the Plutos of astronauts.','Position',...
    [100,280,300,50],'FontSize',12);
Mooninfo2a = uicontrol('Style','text','String',...
    'Intense gravity? What the',...
    'Position',[100,200,300,50],'FontSize',12);
Mooninfo2b = uicontrol('Style','text','String',...
    'heck is that?',...
    'Position',[100,180,300,50],'FontSize',12);
Mooninfo3a = uicontrol('Style','text','String',...
    'You have enough fuel to feed the',...
    'Position',[100,120,300,50],'FontSize',12);
Mooninfo3b = uicontrol('Style','text','String',...
    ' entire cast of Pixar''s Cars.',...
    'Position',[100,100,300,50],'FontSize',12);
Marsinfo1a = uicontrol('Style','text','String',...
    'For astronauts whose middle names are:',...
    'Position',[400,280,300,50],'FontSize',12);
Marsinfo1b = uicontrol('Style','text','String',...
    '"SPACE EXPLORATION ''TIL DEATH"',...
    'Position',[400,260,300,50],'FontSize',12);
Marsinfo1c = uicontrol('Style','text','String',...
    '(in all caps).','Position',[400,240,300,50],'FontSize',12);
Marsinfo2a = uicontrol('Style','text','String',...
    'Magnitude of gravity that''s higher',...
    'Position',[400,200,300,50],'FontSize',12);
Marsinfo2b = uicontrol('Style','text','String',...
    'than Snoop Dogg.',...
    'Position',[400,180,300,50],'FontSize',12);
Marsinfo3a = uicontrol('Style','text','String',...
    'if destination == far away        ',...
    'Position',[400,120,300,50],'FontSize',12);
Marsinfo3b = uicontrol('Style','text','String',...
    '      fuel_at_arrival = not a lot;',...
    'Position',[400,100,300,50],'FontSize',12);
Marsinfo3c = uicontrol('Style','text','String',...
    'end                                        ',...
    'Position',[400,80,300,50],'FontSize',12);
end

function [expd_bg1_p,expd_bg2_p,tmoon_p,tmars_p] = planetpatches
%Create background patch for expedition screen
r_bg = 60*sqrt(5^2+5^2);
thmin = 0;
thmax = 360;
dth = 6;
for n = 1:61
    th(n) = thmin+dth*(n-1);
    x_bg(n) = r_bg*cosd(th(n));
    y_bg(n) = r_bg*sind(th(n));
end
%Use patch function
expd_bg1_p = patch(x_bg+500,y_bg+1450,'k');
expd_bg2_p = patch(x_bg+1510,y_bg+1450,'k');
%Create patch for moon
rtmoon = 30*sqrt(5^2+5^2);
thmin = 45;
thmax = -135;
dth = (thmax-thmin)/30;
for n = 1:31
    th(n) = thmin+dth*(n-1);
    xtmoon(n) = rtmoon*cosd(th(n));
    ytmoon(n) = rtmoon*sind(th(n));
end
rtmoon = 30*sqrt(11^2+1^2);
thmin = -90+atand(1/11);
thmax = -atand(1/11);
dth = (thmax-thmin)/30;
offx = -600*3/10;
offy =  600*3/10;
for n = 32:62
    th(n) = thmin+dth*(n-32);
    xtmoon(n) = rtmoon*cosd(th(n))+offx;
    ytmoon(n) = rtmoon*sind(th(n))+offy;
end
%Use patch function
tmoon_p = patch(xtmoon+500,ytmoon+1450,'w');
%Create patch for mars
rtmars = 30*sqrt(5^2+5^2);
thmin = 0;
thmax = 360;
dth = 6;
for n = 1:61
    th(n) = thmin+dth*(n-1);
    xtmars(n) = rtmars*cosd(th(n));
    ytmars(n) = rtmars*sind(th(n));
end
%Use patch function
tmars_p = patch(xtmars+1510,ytmars+1450,'r');
end

function delete_expd(ChooseExpBox,MoonButt,MarsButt,Mooninfo1,Mooninfo2a,...
          Mooninfo2b,Mooninfo3a,Mooninfo3b,Marsinfo1a,Marsinfo1b,...
          Marsinfo1c,Marsinfo2a,Marsinfo2b,Marsinfo3a,Marsinfo3b,...
          Marsinfo3c,expd_bg1_p,expd_bg2_p,tmoon_p,tmars_p)
%Delete everything on expeditions screen
delete(ChooseExpBox)
delete(MoonButt)
delete(MarsButt)
delete(Mooninfo1)
delete(Mooninfo2a)
delete(Mooninfo2b)
delete(Mooninfo3a)
delete(Mooninfo3b)
delete(Marsinfo1a)
delete(Marsinfo1b)
delete(Marsinfo1c)
delete(Marsinfo2a)
delete(Marsinfo2b)
delete(Marsinfo3a)
delete(Marsinfo3b)
delete(Marsinfo3c)
delete(expd_bg1_p)
delete(expd_bg2_p)
delete(tmoon_p)
delete(tmars_p)
end

function [gm,Rm,mship,mfuel,mtot,mfuel0] = setplanet(mars)
%Change values for mars
if mars == 1
    %Use Mars' gravity and radius
    gm = 3.71;    %[m/s^2]
    Rm = 3.389E6; %[m]
    %Initialize mass fuel and other mass values
    mship = 4900;       %mass of the ship     [kg]
    mfuel = 10300;      %mass of the fuel     [kg]
    mfuel0 = 10300;     %initial mass of fuel [kg]
    mtot = mship+mfuel; %total mass of the ship (sum of ship and fuel) [kg]
else
    %Use Moon's gravity and radius
    gm = 1.62;    %[m/s^2]
    Rm = 1.738E6; %[m]
    %Initialize mass fuel and other mass values
    mship = 4900;       %mass of the ship         [kg]
    mfuel = 15300;      %mass of the fuel         [kg]
    mfuel0 = 15300;     %initial mass of the fuel [kg]
    mtot = mship+mfuel; %total mass of the ship (sum of ship and fuel) [kg]
end
end

%% SETUP FUNCTIONS

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

%Title (was scrapped, keeping it in code in case)
%uicontrol('Style','text','String','***LUNAR***LANDER***',...
%   'Position',[100,750,600,50],'FontSize',30);

%"Ready for Takeoff?" button
Readybutton = uicontrol('Style','pushbutton','String','Ready for Takeoff?',...
    'Position',[325,740,150,30],'Callback',@Start,'FontSize',10);
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
%Define background, bg
bg = [0, xmax, xmax,    0;...
      0,    0, xmax, xmax]; %[m]
end

function [bg_patch,xmoon,ymoon,moon_patch] = create_setting(bg,mars)
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

%Create stars for background after background is plotted and before moon
%surface is plotted
createstars;

%Create moon patch
%Generate x-coordinates of lunar surface
xmoon = linspace(0,2000,2001); %[m]
%Generate y-coordinates of lunar surface and create patch
if mars == 1
    [ymoon] = genymars; %[m]
    moon_patch = patch(xmoon,ymoon,'r');
else
    [ymoon] = genymoon; %[m]
    moon_patch = patch(xmoon,ymoon,'w');
end

%Use nested function to create ymoon
function [ymoon] = genymoon
    %NOTE: randomized = psuedo-randomized
%Generate randomized coords
%Create "empty" array to be filled with randomized y coords
y = zeros(1,19);
%Create first array (101 elements, not 100)
y(1) = randi(13);
ymoon = 100*linspace(0,y(1),101);
%Create randomized x coords for flat surfaces
flat1 =       2+randi(2);
flat2 = flat1+1+randi(4);
flat3 = flat2+1+randi(4);
flat4 = flat3+1+randi(4);
%Create the "middle" coords
for a = 2:19
    %Make sure that there are 3 flat surfaces
    if a == flat1 || a == flat2 || a == flat3 || a == flat4
        y(a) = y(a-1); %make y = previous y value so surface is flat
    elseif y(a-1) >= 13 %max height, y must go down
        y(a) = y(a-1)-randi(3); %decrease
    elseif y(a-1) <= 5 %min height, y must go up
        y(a) = y(a-1)+randi(3); %increase 
    elseif a < flat1 %before first flat
        y(a) = y(a-1)+randi(3); %increase
    elseif a > flat1 && a < flat2 %between first and second flats
        y(a) = y(a-1)+randi(3); %increase
    elseif a > flat2 && a < flat3 %between second and third flats
        y(a) = y(a-1)-randi(3); %increase
    elseif a > flat3 && a < flat4 %between third and fourth flats
        y(a) = y(a-1)+randi(3); %increase
    else %after fourth flat
        y(a) = y(a-1)-randi(3); %decrease
    end
    newcoords = 100*linspace(y(a-1),y(a),100); %create 100 element
    %array connecting previous y value to current y value
    ymoon = [ymoon,newcoords];
end
%Create last array
ymoon = [ymoon,100*linspace(y(19),0,100)];
end

%Use nested function to create ymars
function [ymoon] = genymars
    %NOTE: randomized = psuedo-randomized
%Generate randomized coords
%Create "empty" array to be filled with randomized y coords
y = zeros(1,19);
%Create first array (101 elements, not 100)
y(1) = randi(13);
ymoon = 100*linspace(0,y(1),101);
%Create randomized x coords for flat surfaces
flat1 =       2+randi(2);
flat2 = flat1+2+randi(3);
flat3 = flat2+2+randi(2);
flat4 = flat3+2+randi(3);
%Create the "middle" coords
for a = 2:19
    %Make sure that there are 3 flat surfaces
    if a == flat1 || a == flat2 || a == flat3 || a == flat4
        y(a) = y(a-1); %make y = previous y value so surface is flat
    elseif y(a-1) >= 13 %max height, y must go down
        y(a) = y(a-1)-randi(3); %decrease
    elseif y(a-1) <= 5 %min height, y must go up
        y(a) = y(a-1)+randi(3); %increase 
    elseif a < flat1 %before first flat
        y(a) = y(a-1)+randi(3); %increase
    elseif a > flat1 && a < flat2 %between first and second flats
        y(a) = y(a-1)+randi(3); %increase
    elseif a > flat2 && a < flat3 %between second and third flats
        y(a) = y(a-1)-randi(3); %increase
    elseif a > flat3 && a < flat4 %between third and fourth flats
        y(a) = y(a-1)+randi(3); %increase
    else %after fourth flat
        y(a) = y(a-1)-randi(3); %decrease
    end
    %Analyze type of segment
    if a == flat1 || a == flat2 || a == flat3 || a == flat4
        newcoords = 100*linspace(y(a-1),y(a),80); %flats are less wide
    elseif abs(flat1-a) == 1 || abs(flat2-a) == 1 ||...
           abs(flat3-a) == 1 || abs(flat4-a) == 1
       %if segment is adjacent to flat surface
       newcoords = 100*linspace(y(a-1),y(a),110); %add 10 extra elements on
       %each adjacent segment to make up for shorter flats
    else
    newcoords = 100*linspace(y(a-1),y(a),100); %create 100 element
    %array connecting previous y value to current y value
    end
    %Add new coords to ymoon array
    ymoon = [ymoon,newcoords];
end
%Create last array
ymoon = [ymoon,100*linspace(y(19),0,100)];
end

%Use nested function to create the OG ymoon
function [ymoon] = genymoon_OG
%Generate psuedo-random coords

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
    function createstars
    %Create x and y coords for star patch
    star_x = [ 4, 1, 0,-1,-4,-1, 0, 1, 4];
    star_y = [ 0, 1, 4, 1, 0,-1,-4,-1, 0];
    %Concatenate x and y coords
    star_org = 2/3*[star_x;star_y];
    %Create array of 100 randomized star positions and plot 100 stars of
    %varying size
        for n = 1:100
            star_pos_x = randi(2000);
            star_pos_y = randi(2000);
            %Shift star upwards if too low
            if star_pos_y < 500
                star_pos_y = star_pos_y+500+randi(1000);
            end
            star = randi(3)*star_org;
            patch(star(1,:)+star_pos_x,...
                  star(2,:)+star_pos_y,'w');
        end
    end
end

function [f,fp1,fp2,flength,rf1,rf2,tf1,tf2,fx]...
    = create_firepatch(x,y,selected)
%Creates fire patch by creates x and y coordinate arrays, concatenating
%them, then using the patch function to create the patch

%Inputs:  none
%Outputs: fp   - the fire patch that appears when there is thrust
%         fpl  - fire patch length, variable used to vary fire patch length

flength = 0;  %length of fire starts at 0
%Change coords depending on ship selected
if selected == 1
    %Use coords for ship (1st ship)
    %Create x and y coordinates
    fx = [0,      1, 2];
    fy = [0,flength, 0]; %Fire patch has area = 0 to start because at the
                         %start, thrust is off
    %Create offset values
    offx1 = -35;
    offy1 = -40;
    offx2 =  15;
    offy2 = -40;
elseif selected == 2
    %s2 coordinates (same as 1st ship)
    fx = [0,      1, 2];
    fy = [0,flength, 0];
    %s2 offset values
    offx1 = -35;
    offy1 = -50;
    offx2 =  15;
    offy2 = -50;
elseif selected == 3
    %s3 coordinates
    fx = [0,    1.5, 3];
    fy = [0,flength, 0];
    %s3 offset values
    offx1 = -15;
    offy1 = -50;
    offx2 = -15;
    offy2 = -50;
end

%Concatenate x and y coordinates
f = 10*[fx;fy];
    
%Radii and angles of fire patches origin to ship origin
rf1 = sqrt(offx1^2+offy1^2);
rf2 = sqrt(offx2^2+offy2^2);
tf1 = 180+abs(atand(offy1/offx1));
tf2 = 270+abs(atand(offx2/offy2));

%Make fire patches overlap for s3
if selected == 3
    %s3 radius and theta values are same for BOTH patches because they overlap
    %to appear as a single patch
    rf2 = rf1;
    tf2 = tf1;
end

%Use patch function to create patch
fp1 = patch(f(1,:)+x+rf1*cosd(tf1),...
            f(2,:)+y+rf1*sind(tf1),'y');

fp2 = patch(f(1,:)+x+rf2*cosd(tf2),...
            f(2,:)+y+rf2*sind(tf2),'y');
end

function [FuelBox,FuelunitsBox,AltBox,AltunitsBox,...
           VxBox,VxunitsBox,VyBox,VyunitsBox,TBox,TunitsBox] = text_boxes
%Creates text boxes and displays them in the figure window, next to the
%plot using the UI Controls

%Inputs: none
%Outputs: none

%Create Start button using UIControls
uicontrol('Style','pushbutton','String','Start Mission!',...
          'Position',[325,740,150,30],'Callback',@Start,'FontSize',14)

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
VyunitsBox = uicontrol('Style','text','String','m/s',...
          'Position',[0,225,100,50],'FontSize',12);
%Thrust magnitude
TBox = uicontrol('Style','text','String','Thrust:',...
        'Position',[0,150,100,50],'FontSize',14);
%Units of thrust magnitude
TunitsBox = uicontrol('Style','text','String','N',...
          'Position',[0,100,100,50],'FontSize',12);
end

function [expl_org] = create_expl
%Initial coords of curves (relative to center)
x1 = [ 3, 0,-2,-4,-3, 2, 2, 4]; %initial x
y1 = [-3,-4,-4, 1, 3, 5, 2, 0]; %initial y

%Offset values of curve centers
offx_arr = [ 0,-3,-5,-3, 1, 2, 5, 3];
offy_arr = [-4,-3, 1, 4, 5, 3,-1,-3];

%min and max values of theta
thmin_arr = [ -45,270,180+atan2d(4,2),180-atan2d(1,4),...
     135,90-atan2d(2,5), 45,  0];
thmax_arr = [-135,180,180-atan2d(4,2), 90-atan2d(1,4),...
      45,  -atan2d(2,5),-45,-90];

%Initialize nmin and nmax
nmin = -10;
nmax = 0;
for a = 1:8
%Calculate radius, angles, and offset values
r = sqrt(y1(a)^2+x1(a)^2);
thmin = thmin_arr(a);
thmax = thmax_arr(a);
dth = (thmax-thmin)/10;
offx = offx_arr(a);
offy = offy_arr(a);
nmin = nmin+11;
nmax = nmax+11;
    for n = nmin:nmax
        th(n) = thmin+dth*(n-(11*(a-1)+1));
        expl_x(n) = r*cosd(th(n))+offx;
        expl_y(n) = r*sind(th(n))+offy;
    end
end
%Concatenate x and y coords
expl_org = 2*[expl_x;expl_y];
end

function [expl2_org] = create_expl2
%Create coords of expl2
%Initial coords of curves (relative to center)
x1 = [ 3,2, 1/2,-1,-3,-1,-1/2]; %initial x
y1 = [ 1,0,-5/2,-3, 0, 2, 5/2]; %initial y

%Offset values of curve centers
offx_arr = [-2,-5,-11/2,-2,4, 5, 7/2];
offy_arr = [-6,-3,  3/2, 5,4,-1,-9/2];

%min and max values of theta
thmin_arr = [   atand(1/3), 0,-90+atand(1/5),-90-atand(1/3),...
    -180, 90+atand(1/2), 90+atand(1/5)];
thmax_arr = [90+atand(1/3),90,    atand(1/5),   -atand(1/3),...
     -90,180+atand(1/2),180+atand(1/5)];

%Initialize nmin and nmax
nmin = -10;
nmax = 0;
for a = 1:7
%Calculate radius, angles, and offset values
r = sqrt(y1(a)^2+x1(a)^2);
thmin = thmin_arr(a);
thmax = thmax_arr(a);
dth = (thmax-thmin)/10;
offx = offx_arr(a);
offy = offy_arr(a);
nmin = nmin+11;
nmax = nmax+11;
    for n = nmin:nmax
        th(n) = thmin+dth*(n-(11*(a-1)+1));
        expl_x2(n) = r*cosd(th(n))+offx;
        expl_y2(n) = r*sind(th(n))+offy;
    end
end
%Concatenate x and y coords
expl2_org = 20*[expl_x2;expl_y2];
end

function [expl3a_org,expl3b_org,expl3c_org] = create_expl3
[expl3a_org] = create_expl3a;
[expl3b_org] = create_expl3b;
[expl3c_org] = create_expl3c;

    function [expl3a_org] = create_expl3a %done
    %Initial coords of curves (relative to center)
    x1 = [-1, 1,-1, 1]; %initial x
    y1 = [ 2,-2,-2, 2]; %initial y
    
    %Offset values of curve centers
    offx_arr = [ 3, 0, 0, 0]; %[-2,-5,-5,-5];
    offy_arr = [-3, 0, 4, 0]; %[-6,-3, 1,-3];
    
    %min and max values of theta
    thmin_arr = [  90+atand(1/2), -90+atand(1/2),270-atand(1/2),...
        90-atand(1/2),90-atand(1/2)];
    thmax_arr = [-180-atand(1/2),-270+atand(1/2),-90+atand(1/2),...
          -atand(1/2),  -atand(1/2)];
    
    %Initialize nmin and nmax
    nmin = -20;
    nmax = 0;
    for a = 1:4
    %Calculate radius, angles, and offset values
    r = sqrt(y1(a)^2+x1(a)^2);
    thmin = thmin_arr(a);%atand(abs(y1(a)/x1(a)))+thmin_off;
    thmax = thmax_arr(a);%atand(abs(y2(a)/x2(a)))+thmax_off;
    dth = (thmax-thmin)/20;
    offx = offx_arr(a);
    offy = offy_arr(a);
    nmin = nmin+21;
    nmax = nmax+21;
        for n = nmin:nmax
            th(n) = thmin+dth*(n-(21*(a-1)+1));
            expl3a_x(n) = r*cosd(th(n))+offx;
            expl3a_y(n) = r*sind(th(n))+offy;
        end
    end
    %Concatenate x and y coords
    expl3a_org = 20*[expl3a_x;expl3a_y];
    end
    function [expl3b_org] = create_expl3b
    %Initial coords of curves (relative to center)
    x1 = [ 1,-1,-1,-2]; %initial x
    y1 = [-2, 2, 1,-1]; %initial y
    
    %Offset values of curve centers
    offx_arr = [-3, 0, 3, 0];
    offy_arr = [ 1, 0, 0, 0];
    
    %min and max values of theta
    thmin_arr = [270+atand(1/2),90+atand(1/2), 90+atand(1/2),...
        -atand(1/2)];
    thmax_arr = [    atand(1/2),   atand(1/2),-90-atand(1/2),...
        -180+atand(1/2)];
    
    %Initialize nmin and nmax
    nmin = -20;
    nmax = 0;
    for a = 1:4
    %Calculate radius, angles, and offset values
    r = sqrt(y1(a)^2+x1(a)^2);
    thmin = thmin_arr(a);%atand(abs(y1(a)/x1(a)))+thmin_off;
    thmax = thmax_arr(a);%atand(abs(y2(a)/x2(a)))+thmax_off;
    dth = (thmax-thmin)/20;
    offx = offx_arr(a);
    offy = offy_arr(a);
    nmin = nmin+21;
    nmax = nmax+21;
        for n = nmin:nmax
            th(n) = thmin+dth*(n-(21*(a-1)+1));
            expl3b_x(n) = r*cosd(th(n))+offx;
            expl3b_y(n) = r*sind(th(n))+offy;
        end
    end
    %Concatenate x and y coords
    expl3b_org = 20*[expl3b_x;expl3b_y];
    end
    function [expl3c_org] = create_expl3c
    %Initial coords of curves (relative to center)
    x1 = [-2, 2]; %initial x
    y1 = [-1, 1]; %initial y
    
    %Offset values of curve centers
    offx_arr = [1/2,-1/2];
    offy_arr = [3/2,-3/2];
    
    %min and max values of theta
    thmin_arr = [180+atand(1/2),     atand(1/2)];
    thmax_arr = [-90+atand(1/2),-270+atand(1/2)];
    
    %Initialize nmin and nmax
    nmin = -20;
    nmax = 0;
    for a = 1:2
    %Calculate radius, angles, and offset values
    r = sqrt(y1(a)^2+x1(a)^2);
    thmin = thmin_arr(a);%atand(abs(y1(a)/x1(a)))+thmin_off;
    thmax = thmax_arr(a);%atand(abs(y2(a)/x2(a)))+thmax_off;
    dth = (thmax-thmin)/20;
    offx = offx_arr(a);
    offy = offy_arr(a);
    nmin = nmin+21;
    nmax = nmax+21;
        for n = nmin:nmax
            th(n) = thmin+dth*(n-(21*(a-1)+1));
            expl3c_x(n) = r*cosd(th(n))+offx;
            expl3c_y(n) = r*sind(th(n))+offy;
        end
    end
    %Concatenate x and y coords
    expl3c_org = 20*[expl3c_x;expl3c_y];
    end
end

%% GAMEPLAY FUNCTIONS

function Start(pushbutton, ~)
%Allows user to change start variable value by pressing the "Start Mission!"
%button

%Inputs: pushbutton - source of push button event
%        ~          - I'm not sure what this means.  Maybe it is the push
%                     button event?
%Outputs: none

%Give access to global variables
global start choose

%Change start value
start = 1;

%Change choose value (this can be added to this function because choose
%only changes once)
choose = 1;
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

function [h_moon] = calc_h_moon(x,ymoon)
%Calculates the height of the moon directly below the ship

%Inputs:  x      - x-position [m]
%         ymoon  - array of y-coordinates for moon patch      [m]
%Outputs: h_moon - y coordinate of moon directly beneath ship [m]

    if x >= 0 && x <= 2000 
        %Center of ship is within axes limits, calculate normally
        %Calculate counter, n, by rounding x value and using scale variable
        n = round(x);
        %Set h_moon equal to the y-value of moon corresponding to ship's x-value
        h_moon = ymoon(n); %[m]
    else
        %Center of ship is not within axes limits, set h_moon = 0
        h_moon = 0;
    end
end

function [Ax,Ay,x,y,Vx,Vy] =...
    LanderAccel(x,y,Vx,Vy,Q,mtot,T,gm,Rm,dt)
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
    %gm = 1.62;      % {m/s^2} I am using external gm value instead
    %Rm = 1.738E6;   % {m}

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
    
    %Calc position vectors using kinematics formulas and Euler's method
    x = x+Vx*dt+1/2*Ax*dt^2; %[m]
    y = y+Vy*dt+1/2*Ay*dt^2; %[m]
    
    %Calc velocity vectors using kinematics formulas and Euler's method
    Vx = Vx*dt+1/2*Ax*dt^2; %[m/s]
    Vy = Vy*dt+1/2*Ay*dt^2; %[m/s]
    %The above formulas are completely wrong...nice.  The following ones
    %are correct, but make the gameplay terrible because everything else in
    %the code regarding physics was programmed around my incorrect
    %formulas.
    %Vx = Vx+Ax*dt; %[m/s]
    %Vy = Vy+Ay*dt; %[m/s]
end

function ...
    animate_expl(expl_org,expl2_org,expl3a_org,expl3b_org,expl3c_org,...
    x,y,ship_patch)
%Create initial patch for expl and in_expl (inner explosion)
expl_p = patch(1/10*expl_org(1,:)+x,1/10*expl_org(2,:)+y,'y');
in_expl_p = patch(3/40*expl_org(1,:)+x,3/40*expl_org(2,:)+y,'r');
%Make it expand
for n = 1:10
    expl = expl_org*n;
    in_expl = expl_org*n*3/4;
    set(expl_p,'XData',expl(1,:)+x,'YData',expl(2,:)+y)
    set(in_expl_p,'XData',in_expl(1,:)+x,'YData',in_expl(2,:)+y)
    pause(0.01)
end
%Create initial patch for expl2
expl2_p = patch(1/10*expl2_org(1,:),1/10*expl2_org(2,:),'k');
for n = 1:10
    expl2 = expl2_org*n/10;
    set(expl2_p,'XData',expl2(1,:)+x,'YData',expl2(2,:)+y)
    pause(0.01)
end
%Delete patches of expl, expl_in, expl2
delete(expl_p)
delete(in_expl_p)
delete(expl2_p)
%Delete ship and fire patches
delete(ship_patch)
%Create initial patch for expl3
expl3a_p = patch(expl3a_org(1,:)-50,expl3a_org(2,:)-30,'y');
expl3b_p = patch(expl3b_org(1,:)+10,expl3b_org(2,:)+40,'y');
expl3c_p = patch(expl3c_org(1,:)+45,expl3c_org(2,:)-25,'y');
%Shrink expl3 patches
for n = 1:10
    if n<6
        expl3a = expl3a_org/n;
        expl3b = expl3b_org/n;
        expl3c = expl3c_org/n;
    elseif n == 6
        expl3a_org = expl3a;
        expl3b_org = expl3b;
        expl3c_org = expl3c;
        expl3a = expl3a_org/(n-4);
        expl3b = expl3b_org/(n-4);
        expl3c = expl3c_org/(n-4);
    else
        expl3a = expl3a_org/(n-4);
        expl3b = expl3b_org/(n-4);
        expl3c = expl3c_org/(n-4);
    end
    set(expl3a_p,'XData',expl3a(1,:)-50+x,'YData',expl3a(2,:)-30+y);
    set(expl3b_p,'XData',expl3b(1,:)+10+x,'YData',expl3b(2,:)+40+y);
    set(expl3c_p,'XData',expl3c(1,:)+45+x,'YData',expl3c(2,:)-25+y);
    pause(0.05)
end
%Delete patches of expl3
delete(expl3a_p)
delete(expl3b_p)
delete(expl3c_p)
end

function [pole_p,flag_p,fbox_p,stripe1_p,stripe2_p,stripe3_p] =...
    flag(x,y,hship,wship)
%Initialize
offx = x+wship;
offy = y-hship;
%Create pole coords and patch
polex = [  0,  5,  5,  0];
poley_org = 1/10*[  0,  0,110,110];
pole_p = patch(polex+offx,poley_org+offy,'w');

%Raise pole
for n = 2:10
    %Change y coords
    poley = n*poley_org;
    %Set new patches
    set(pole_p,'XData',polex+offx,...
               'YData',poley+offy);
    pause(0.01)
end
    pause(0.25)
    
%Create flag coords and patch
%Red background
    flagx_org = 1/10*[  5, 80, 80,  5];
    flagy = [ 60, 60,110,110];
    flag_p = patch(flagx_org+offx,flagy+offy,'r');
%Blue corner box
    fboxx_org = 1/10*[  5, 40, 40,  5];
    fboxy = [ 90, 90,110,110];
    fbox_p = patch(fboxx_org+offx,fboxy+offy,'b');
%White stripes
    stripe12x_org = 1/10*[  5, 80, 80, 5];
    stripe3x_org = 1/10*[ 40, 80, 80, 40];
    stripe1y = [ 65, 65, 75, 75];
    stripe2y = [ 80, 80, 90, 90];
    stripe3y = [ 95, 95,105,105];
    stripe1_p = patch(stripe12x_org+offx,stripe1y+offy,'w');
    stripe2_p = patch(stripe12x_org+offx,stripe2y+offy,'w');
    stripe3_p = patch(stripe3x_org+40+offx,stripe3y+offy,'w');
    
%Extend flag
for n = 2:10
    %Change x coords
    flagx = flagx_org*n;
    fboxx = fboxx_org*n;
    stripe12x = stripe12x_org*n;
    stripe3x = stripe3x_org*n;
    %Set new patches
    set(flag_p,'XData',flagx+offx,...
               'YData',flagy+offy);
    set(fbox_p,'XData',fboxx+offx,...
               'YData',fboxy+offy);
    set(stripe1_p,'XData',stripe12x+offx,...
               'YData',stripe1y+offy);
    set(stripe2_p,'XData',stripe12x+offx,...
               'YData',stripe2y+offy);
    set(stripe3_p,'XData',stripe3x+offx,...
               'YData',stripe3y+offy);
    pause(0.005)
end
end

function [point] = evaluate(Vy,mfuel,x,Q,ymoon,mars)
%Uses conditional operators to evaluate how the player has done.
%Depeneding on the conditions of the ship's landing, a message is displayed
%in the figure window using UI Controls

%Inputs:  Vy    - velocity in y-direction of ship [m/s]
%         mfuel - mass of fuel                    [kg]
%         x     - x-position of ship              [m]
%         Q     - angle from vertical of ship to vertical of plot [deg]
%Outputs: none

%Evaluate player's conditions and output message in text box
%Evaluate angle, Q using nested function
[Qeval] = evalQ(Q);
if x < 0 || x > 2000
    %the player has managed to finish the game off of the screen
    uicontrol('Style','text','String','WOW!  You managed to land OFF SCREEN!',...
        'Position',[100,35,600,50],'FontSize',12);
    uicontrol('Style','text','String','Maybe you should try a different game...',...
        'Position',[100,15,600,50],'FontSize',12);
    point = 0;
elseif ymoon(round(x)) == ymoon(round(x)+1)...
        || ymoon(round(x)) == ymoon(round(x)-1)
    %Player has landed on a flat surface
    if Qeval == 1 && Vy>-6*10^(-1)
        %Player has landed ship oriented upright, with low velocity
        %Player has met all requirements to win, display message
        if mars ~= 1
            uicontrol('Style','text','String','Congratulations! You have successfully landed on the Moon!',...
            'Position',[100,35,600,50],'FontSize',12);
        else
            uicontrol('Style','text','String','Congratulations! You have successfully landed on Mars!',...
            'Position',[100,35,600,50],'FontSize',12);
        end
        uicontrol('Style','text','String','Good luck getting back home with only',...
        'Position',[95,15,400,50],'FontSize',12);
        uicontrol('Style','text','String',num2str(mfuel,'%7.2f'),...
        'Position',[435,15,60,50],'FontSize',12);
        uicontrol('Style','text','String','kg of fuel left! XD',...
        'Position',[495,15,150,50],'FontSize',12);
        %Add to point counter, point
        point = 1;
    else
        %Player has not met requirements for a smooth landing, display
        %message
        uicontrol('Style','text','String','That landing was too harsh, your ship got destroyed! Try again!',...
        'Position',[100,35,600,50],'FontSize',12);
        %Keep the same point value
        point = 0;
    end
else
    %Player has not landed on one of four flat surfaces, display message
    uicontrol('Style','text','String','Your ship cannot land on a slope that steep! Try again!',...
    'Position',[100,35,600,50],'FontSize',12);
    %Keep the same point value
    point = 0;
end
    function [Qeval] = evalQ(Q)
        if  Q   /360 == round( Q   /360) ||...
           (Q-2)/360 == round((Q-2)/360) ||...
           (Q+2)/360 == round((Q+2)/360)
           %Respective descriptions of situations:
           %Q is a multiple of 360
           %Q is >= 2 less than a multiple of 360
           %Q is <= 2 more than a multiple of 360
            Qeval = 1;
        else
            %Q is not within 2 degrees of a multiple of 360
            Qeval = 0;
        end
    end
end

function tryagain_button
%Create Try Again button using UIControls
uicontrol('Style','pushbutton','String','Try Again!',...
          'Position',[325,740,150,30],'Callback',@Tryagain,'FontSize',14)
end

function Tryagain(pushbutton, ~)
%Allows user to change start variable value by pressing the "Start Mission!"
%button

%Inputs: pushbutton - source of push button event
%        ~          - I'm not sure what this means.  Maybe it is the push
%                     button event?
%             point - I don't remember why this is an input
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
global x y Q T mfuel dT Tmax flength selected

%Reads the name of the pressed key
keyname = event.Key;

%Check which key was pressed and respond
%Player is on title screen    
    switch keyname
%    case 'd'
%        %Center ship to right
%        if     selected == 1 %1 up to 2
%            selected = 2;
%        elseif selected == 2 %2 up to 3
%            selected = 3;
%        elseif selected == 3 %3 up to 1
%            selected = 1;
%        end
%    case 'a'
%        %Center ship to left
%        if     selected == 1 %1 down to 3
%            selected = 3;
%        elseif selected == 2 %2 down to 1
%            selected = 1;
%        elseif selected == 3 %3 down to 2
%            selected = 2;
%        end
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
            %Change fire patch using nested function
            flength = flength-0.5;
        end
    case 'downarrow'
        %Decrease Thrust if Thrust is not already 0
        if T > 0
            T = T-dT; %[N]
            %Change fire patch using nested function
            flength = flength+0.5;
        end
    end
end
