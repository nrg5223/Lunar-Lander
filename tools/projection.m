% Projection Mini-Project
% Nicholas Grosskopf
% Started 8/15/22
clear
clc
close all

%% INITIALIZE

% GLOBAL
global x_pos y_pos quit

% COLORS
darkred = [0.9,0,0];
brightred = [1,0,0];
gray_wbluetint = [0.8,0.8,0.85];
gold = [1,0.85,0];
black = [0,0,0];
white = [1,1,1];
brightblue = [0,0,1];
brightgreen = [0,1,0];
darkgreen = [0,0.5,0];
pink = [1,0.4,0.6];
orange = [1,0.5,0];
purple = [0.7,0,0.7];
blue = [0,0,1];

% FIG AND AXES
fig_color = white;
fig_scale = 100;
axes_color = gray_wbluetint;

% SQUARE
sq_scale = 10;
sq_color = white;

% DOTS
markersize = 10;
max_x = 20;
max_y = max_x;
min_x = -max_x;
min_y = min_x;
num_points = 10; % num points in 1d arrays for x and y coords
inc = (max_x - min_x) / num_points;


% KEYBOARD
x_pos = 0;
y_pos = 0;
quit = 0;

% PHYSICS
dt = 0.1;


%% MAIN PROGRAM

% Creates figure window on left half of screen
[fig] = create_fig(fig_color);

% Create test axes
x_min = -1 * fig_scale;
x_max =  1 * fig_scale;
y_min = x_min;
y_max = x_max;
% Creates axes with... 
axes('XLim', [x_min, x_max], 'YLim', [y_min, y_max],... % given limits
    'color', axes_color);%,...                          % color
%    'XTickLabels',[],'YTickLabels',[]);                 % blank graph marks
grid on

% Create square patch
%x_coords = sq_scale * [-1, -1,  1,  1];
%y_coords = sq_scale * [-1,  1,  1, -1];
sq_pos = [0, 0];
%square = patch(x_coords + sq_pos(1),...
%               y_coords + sq_pos(2),... 
%               sq_color);

% Create array of dots in the shape of a square
    % Create 1d array from 0 to sq_scale with num_points points
    xs_dots = linspace(0, sq_scale, num_points);
    ys_dots = linspace(0, sq_scale, num_points);
    % Create 2d array of dots that make the square shape
    dots = zeros(num_points, num_points);
    dots_x = zeros(num_points, num_points);
    dots_y = zeros(num_points, num_points);
    x_count = 1;
    y_count = 1;
    for x = min_x:inc:max_x
        for y = min_y:inc:max_y
%            fprintf("(%3.1f, %3.1f)\n", x, y); % to confirm the points are
%            %correct
            dots(x_count, y_count) = line(x, y, 'marker','.',...
                                    'markersize',markersize,'color', blue);
            dots_x(x_count, y_count) = x;
            dots_y(x_count, y_count) = y;
            y_count = y_count + 1;
        end
        x_count = x_count + 1;
    end

% Main loop
while quit ~= 1
    % Set the square's position and update the patch
%    sq_pos = sq_scale * [x_pos, y_pos];
%    set(square, 'XData', x_coords + sq_pos(1),...
%                'YData', y_coords + sq_pos(2));

    % Set the dots' positions and update the dots (which are lines)
    x_count = 1;
    y_count = 1;
    for x = min_x:inc:max_x
        for y = min_y:inc:max_y
            set(dots(x_count, y_count),...
                'XData', dots_x(x_count, y_count) + x_pos,...
                'YData', dots_y(x_count, y_count) + y_pos);
            y_count = y_count + 1;
        end
        x_count = x_count + 1;
    end
    
    % Pause so you can see changes in the UI
    pause(dt);
end

%%%% End of program %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FUNCTIONS
function [fig1] = create_fig(color)
%Creates a figure window on the left half of the screen using screen
%dimensions

%Get screen dimensions
scrsize = get(0, 'ScreenSize'); % = [1,1,1920,1080]
swidth = scrsize(3);
sheight = scrsize(4);

% Create figure window
fig1 = figure('Position', [0,0,swidth/2,sheight], 'color', color,...
    'KeyPressFcn', @keyboard);
end

function keyboard(figure,event)
global x_pos y_pos quit
switch event.Key
    case 'leftarrow'
        x_pos = x_pos - 1;
        disp('x');
    case 'rightarrow'
        x_pos = x_pos + 1;
    case 'downarrow'
        y_pos = y_pos - 1;
    case 'uparrow'
        y_pos = y_pos + 1;
    case 'x'
        quit = 1;
end
end