
%compute mean flow field with respect to the moving frame of reference of
%the ball

video_name = "2024-07-10 - 14.19.31.avi";
obj = VideoReader(video_name);
tmp = readFrame(obj);                           %here I take the size of the frames

[ylengy,xlengy,~] = size(tmp);                  %I store the size of the frames, first Y (# rows) then X (# columns)
disp(ylengy);
disp(xlengy);

%length of video
%%lengy = obj.NumFrames;                     %I store how many images
start_frame = 500;
end_frame = 2100;
freq_frame = 1;

lengy = (end_frame - start_frame) / freq_frame;


%% first track ball

circleRadius = 100; % Change the radius as needed

%imageSize = [2000, 4000]; % Change the image dimensions as needed       %also here first I store Y then X
%I think I can do like this
imageSize = [500, 500];

% Create a meshgrid representing the image coordinates
[xx, yy] = meshgrid(1:imageSize(2), 1:imageSize(1));

% Generate a binary image with a white circle on a black background
binaryImage = (xx - imageSize(2)/2).^2 + (yy - imageSize(1)/2).^2 <= circleRadius^2;

% Create a 2D Gaussian filter to smooth the reference circle
%gaussianFilter = fspecial('gaussian', [100, 100], 50); % Adjust parameters as needed

% Smooth the binary circle image using the Gaussian filter
smoothedReference = imgaussfilt(double(binaryImage), 10);
%smoothedReference = conv2(double(binaryImage), gaussianFilter, 'same');

% Normalize the smoothed reference image to [0, 1]
smoothedReference = smoothedReference / max(smoothedReference(:));

% Display the binary image
%figure;
%imshow(smoothedReference);
%title('Binary Image of Circle');

centerpos = NaN(lengy,2); %create array for position of ball

figure;
offs = 1;

for i=offs:(lengy)
    frame_index = start_frame + freq_frame * i;
    tmp = read(obj,frame_index);

    % Convert the image to double precision for processing
    img_double = im2double(tmp);
    
    % Extract the channels (assuming the image is in RGB format)
    green_channel = img_double(:, :, 2);
    red_channel = img_double(:, :, 1);
    blue_channel = img_double(:,:,3);
    
    % Threshold to get a binary image
    %binary_green = green_channel > 0.5;
    %binary_red = red_channel > 0.25;
    binary_selection = (blue_channel < .35).*(red_channel);

    %Normalize the reference image to [0, 1]
    referenceImage = binary_selection / max(binary_selection(:)) ;
    
    % Compute the autocorrelation using normxcorr2
    correlationResult = normxcorr2(smoothedReference,referenceImage);
    
    % Display the autocorrelation result
    imshow(correlationResult, 'colormap', jet);
    title('Autocorrelation Result');
    
    % Find the peak correlation location
    [maxValue, maxIndex] = max(correlationResult(:));
    [centerpos(i,1), centerpos(i,2)] = ind2sub(size(correlationResult), maxIndex);
    
    % Display the location of the peak correlation
    hold on;
    plot(centerpos(i,2), centerpos(i,1), 'r+', 'MarkerSize', 10);
    hold off;

    pause(.2)

end
clear obj

save('center_position')

%% The interesting region is only from frame 600 to 1000 in our case.
% here we fix the center position, which is slowly drifting of 5 pixels in
% this region. Note again that the x position is stored in the second entry

load('center_position')

for i=600:1000
    centerpos(i,2) = centerpos(i,2) - round((i-600)*5/400);
end


%% Get radial velocity at the center

interstep = 5;
xmin = 475;     %controls the left edge
xmax = 25;     %controls the right edge of the grid

ymin = 475;    %controls the upper edge
ymax = 25;     %controls the lower edge

offs = 1;
slengy = lengy-offs;

Uq_equat_arr = NaN(length(xmax:interstep:xmin),2*slengy);
Vq_equat_arr = NaN(length(xmax:interstep:xmin),2*slengy);

for i=offs:(lengy-1)
    frame_index = start_frame + freq_frame * i;
    
    % plot(centerpos(i,2)-imageSize(2)/2, centerpos(i,1)-imageSize(1)/2, 'r+', 'MarkerSize', 100);  %centering fix; somehow there is an offset; presumably due to size of convolution image
    
    % process only subset of grid points close to ball
    xs = (centerpos(i,2)-xmax):-interstep:(centerpos(i,2)-xmin);
    ys = (centerpos(i,1)-ymax):-interstep:(centerpos(i,1)-ymin);

    [Xs,Ys] = meshgrid(xs,ys);

    % PIV data extraction from .mat export    
    PIVx = cell2mat(x(i));
    PIVy = cell2mat(y(i));
    PIVu = cell2mat(u_original(i));
    PIVv = cell2mat(v_original(i));
    
    %make interpolations in new grid subset, over finer mesh that is
    %co-moving with the ball
    Uq = interp2(PIVx,PIVy,PIVu,Xs,Ys,'cubic');
    Vq = interp2(PIVx,PIVy,PIVv,Xs,Ys,'cubic');

    %I save only two lines to an array for later averaging
    %the one in the middle
    Uq_equat = Uq(int16((xmin-xmax)/(2*interstep)),:);
    Vq_equat = Vq(int16((xmin-xmax)/(2*interstep)),:);
   
    Uq_equat_arr(:,2*i) = Uq_equat';
    Vq_equat_arr(:,2*i) = Vq_equat';

    %the one just above
    Uq_equat = Uq(int16((xmin-xmax)/(2*interstep))+1,:);
    Vq_equat = Vq(int16((xmin-xmax)/(2*interstep))+1,:);
   
    Uq_equat_arr(:,2*i+1) = Uq_equat';
    Vq_equat_arr(:,2*i+1) = Vq_equat';

end

frame_to_sec = 1/71.16;
pix_to_cm = 4/175;
pix_per_frame_to_cm_per_sec = pix_to_cm / frame_to_sec;

%either this
Uqarrmean = nanmean(Uq_equat_arr(:,315:1000),2)*pix_per_frame_to_cm_per_sec;
Vqarrmean = nanmean(Vq_equat_arr(:,315:1000),2)*pix_per_frame_to_cm_per_sec;

Radial_velocity = (Uqarrmean.^2 + Vqarrmean.^2).^(0.5)

%or this
%Radial_velocity_arr = (Uq_equat_arr.^2 + Vq_equat_arr.^2).^(0.5);

%Radial_velocity_mean = nanmean(Radial_velocity_arr(:,315:1000),2)*pix_per_frame_to_cm_per_sec



