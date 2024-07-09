%Grosskopf_N_S08_E2
%Animate Spaceship Instructions

%Clear Workspace
clear
clc

%Determine screen size
scrsize = get(0, 'ScreenSize')
%[1,1,1920,1080]

%Create Figure Window
    %Set bottom-left corner location
    %Set width and height
    %Create Figure window
    swidth = scrsize(3);
    sheight = scrsize(4);
    fig1 = figure('Position', [0,0,swidth/2,sheight]);
    
%Create x and y axis centered in figure
    %Define axis limits
    xmax = 100;
    xmin = -xmax;
    ymax = xmax;
    ymin = -ymax;

    axis([xmin,xmax,ymin,ymax])
    %--axis equal %sets axis values to inputted values
    %axis off %turns off display so it isnt in the image
    %--axis manual %allows inputted axis values to be used

%Create 'patch' for spaceship
    %Coordinates of polygon defining spaceship shape
    %COCK = [2,2,2.2,2.2,1, -1,-2.2,-2.2,-2,-2,-3,-4,-5,  -4,-3,  -1,-1,   3, 4, 5, 4, 3,2;...
    %        1,5,5.5,  6,6,  6,   6, 5.5, 5, 1, 0,-1,-2,  -3,-4,-4.5,-4,-4.5,-4,-3,-2,-1,0];
    %--SHIP = [0.00, 1.00, 4.00,1.00,0.00,-1.00,-4.00,-1.00;... 
    %--       -2.50,-0.75,-1.50,1.50,5.00, 1.50,-1.50,-0.75];
    %SHIP = COCK; %if you want it to
    
    %Scale ship to make it more visible
    %--SHIPscale = 4;
    %--SHIP = SHIP*SHIPscale; %shortcut for array multiplic
    
    %Initialize data
    %--r_orbit = 100; %radius of orbit
    %--theta_min = 0; %minimum angle, theta
    %--theta_max = 8*pi; %maximum angle, theta
    %--dtheta = 0.01; %change in theta
    %--n_max = (theta_max-theta_min)/dtheta+1; %highest n value in for loop
    %--Xs = 0; %Initial change in x position of ship
    %--Ys = 0; %Initial change in y position of ship

 %Create patch, move it to point (Xs,Ys) & color it blue
    %Define Shp_Patch in term of SHIP, Xs, and Ys
    %--Shp_Patch = patch(SHIP(1,:)+Xs,SHIP(2,:)+Ys,'b');
    %Variable Shp_Patch contains the 'handle' to the patch
    %The colon in the column position of the SHIP array means "all"

%Create for loop to calculate Xs & Ys as functions of theta & r_orbit   
%for b = 1:n_maxr
    %r(b) = r_min+(b-1)*dr;
%for n = 1:n_max %n_max
    %--theta(n) = theta_min+(n-1)*dtheta; %calculates array of theta values
    %--Xs(n) = r_orbit*cos(theta(n)); %calculates array of x position values
    %--Ys(n) = r_orbit*sin(theta(n)); %calculates array of y position values
    
    %Within the for loop, use the set function to update the coordinates of
    %Shp_Patch, relocating it for each value of theta
    %--set(Shp_Patch,'Xdata',SHIP(1,:)+Xs(n),'Ydata',SHIP(2,:)+Ys(n));
    %Add a pause after the set so you can see the movement
    
    %Add rotation of the ship as it orbits
    %Create a rotation matrix
    %--RtnMtrx = [cos(-theta(n)) sin(-theta(n)); -sin(-theta(n)) cos(-theta(n))];
    %Multiply the SHIP coordinates by the rotation matrix before moving it
    %to (Xs,Ys)
    %--SHIP2 = RtnMtrx*SHIP; %Matrix multiplication
    %--set(Shp_Patch,'Xdata',SHIP2(1,:)+Xs(n),'Ydata',SHIP2(2,:)+Ys(n))
    pause(0.05) %pauses between each loop so the motion can be seen
%end
%end
