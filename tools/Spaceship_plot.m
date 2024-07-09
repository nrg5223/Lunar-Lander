%Clear Workspace
clear
clc

%Determine screen size
scrsize = get(0, 'ScreenSize');
%[1,1,1920,1080]

%Create Figure Window
    %Set bottom-left corner location
    %Set width and height
    %Create Figure window
    fig1 = figure('Position', [50,10,600,600]);
    
%Create x and y axis centered in figure
    %Define axis limits
    xmax = 150;
    xmin = -xmax;
    ymax = xmax;
    ymin = -ymax;
    %Create axis with equal aspect ratio, turn off display of axis and
    %prevent MATLAB from auto-rescaling the axis
    axis([xmin,xmax,ymin,ymax])
    axis equal
    axis off
    axis manual

    Xs = 5;
    Ys = 5;
    
%Create 'patch' for spaceship
    %Coordinates of polygon defining spaceship shape
    SHIP = [0.00, 1.00, 4.00,1.00,0.00,-1.00,-4.00,-1.00;... 
           -2.50,-0.75,-1.50,1.50,5.00, 1.50,-1.50,-0.75];
    %Scale ship to make it more visible
    SHIPscale = 2;
    SHIP = SHIP*SHIPscale; %shortcut for array multiplic
    %Create patch, move it to point (Xs,Ys) & color it blue
    Shp_Patch = patch(SHIP(1,:)+Xs,SHIP(2,:)+Ys,'b');
    %Variable Shp_Patch contains the 'handle' to the patch
    %The colon in the column position of the SHIP array means "all"