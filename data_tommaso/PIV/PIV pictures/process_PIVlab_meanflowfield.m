
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

%%


centervelocities = NaN(lengy-1,2);
%note that centerpos(i,1) is the y component, while centerpos(i,2) is the x
%component, so the velocities are switched (1 is y velocity and 2 is x
%velocity)
for i=offs:(lengy-1)
    centervelocities(i,1) = centerpos(i+1,1) - centerpos(i,1);
    centervelocities(i,2) = centerpos(i+1,2) - centerpos(i,2);
end

%%
obj = VideoReader(video_name);
% find settings of interpolation grid and make array space for it.
interstep = 15;
xmin = 500;     %controls the left edge
xmax = 0;     %controls the right edge of the grid

ymin = 500;    %controls the upper edge
ymax = 000;     %controls the lower edge

offs = 1;
slengy = lengy-offs;

Uqarr = NaN(length(xmax:interstep:xmin),length(ymax:interstep:ymin),slengy);
Vqarr = NaN(length(xmax:interstep:xmin),length(ymax:interstep:ymin),slengy);

%run through images, plot image, PIV data and the interpolated data. Make a
%video out of it.

% Specify the file name and video settings
outputFileName = 'output_movie.avi';  % Change the file name as needed
frameRate = 10;  % Adjust the frame rate as needed

% Create a VideoWriter object
videoObj = VideoWriter(outputFileName, 'Motion JPEG AVI');
videoObj.FrameRate = frameRate;

% Open the video file for writing
open(videoObj);

for i=offs:(lengy-1)
    frame_index = start_frame + freq_frame * i;
    tmp = read(obj,frame_index);
    grey_tmp = rgb2gray(tmp);
    
    imshow(grey_tmp, 'colormap', jet);
    
    hold on;
    plot(centerpos(i,2)-imageSize(2)/2, centerpos(i,1)-imageSize(1)/2, 'r+', 'MarkerSize', 100);  %centering fix; somehow there is an offset; presumably due to size of convolution image
    
    % process only subset of grid points close to ball
    xs = (centerpos(i,2)-xmax):-interstep:(centerpos(i,2)-xmin);
    ys = (centerpos(i,1)-ymax):-interstep:(centerpos(i,1)-ymin);

    [Xs,Ys] = meshgrid(xs,ys);

    scatter(Xs,Ys,'^m')
    
    % PIV data extraction from .mat export    
    PIVx = cell2mat(x(i));
    PIVy = cell2mat(y(i));
    PIVu = cell2mat(u_original(i));
    PIVv = cell2mat(v_original(i));
    
    %plot the gridpoints where data was obtained
    scatter(PIVx,PIVy,'ok')

    %make a velcity field plot
    quiver(PIVx,PIVy,PIVu,PIVv,10)
    
    %make interpolations in new grid subset, over finer mesh that is
    %co-moving with the ball
    Uq = interp2(PIVx,PIVy,PIVu,Xs,Ys,'cubic');
    Vq = interp2(PIVx,PIVy,PIVv,Xs,Ys,'cubic');

    %show the interpolated velocity field; check consistency
    quiver(Xs,Ys,Uq,Vq,10)
    
    %save it to an array for later averaging
    Uqarr(:,:,i) = Uq';
    Vqarr(:,:,i) = Vq';
    hold off;

    % Capture the current figure as a frame
    currentFrame = getframe(gcf);
    
    % Write the frame to the video file
    writeVideo(videoObj, currentFrame);
    
    % Close the current figure to prepare for the next one
    clf;
end

close(gcf);

% Close the video file
close(videoObj);
clear obj

%%
% conversion factors
frame_to_sec = 1/71.16;
pix_to_cm = 4/175;
pix_per_frame_to_cm_per_sec = pix_to_cm / frame_to_sec;

% compute mean of field in comoving frame

Ball_mean_velocity = nanmean(centervelocities(315:1000,:),1)*pix_per_frame_to_cm_per_sec ;
ball_vertical_velocity = (centerpos(1600,1) - centerpos(1,1))/1600*pix_per_frame_to_cm_per_sec ;

Uqarrmean = nanmean(Uqarr(:,:,315:1000),3)*pix_per_frame_to_cm_per_sec;
Vqarrmean = nanmean(Vqarr(:,:,315:1000),3)*pix_per_frame_to_cm_per_sec;




%Uqarrmean = Uqarrmean*pix_per_frame_to_cm_per_sec;
%Vqarrmean = Vqarrmean*pix_per_frame_to_cm_per_sec;

for i=1:34
    for j=1:34
        if (i-17.5)^2 + (j-17.5)^2 < 38
            Uqarrmean(i,j) = Ball_mean_velocity(2) ;
            Vqarrmean(i,j) = Ball_mean_velocity(1) ;
        end
        %if ((i-17.5)^2 + (j-17.5)^2 > 36) * ((i-17.5)^2 + (j-17.5)^2 < 44)
        %    Uqarrmean(i,j) = 10;
        %    Vqarrmean(i,j) = 10;
        %end
    end
end

size_x = [0 (xmin-xmax)*pix_to_cm];
size_y = [(ymin-ymax)*pix_to_cm 0];

subplot 121
imagesc(size_x, size_y,transpose(squeeze(Uqarrmean)),[-0.025, 0.025])
colorbar('eastoutside')
colormap turbo
title('Velocity left/right [cm/s]')

subplot 122
imagesc(size_x, size_y,transpose(squeeze(Vqarrmean)),[-0.04, 0.08])
colorbar('eastoutside')
title('Velocity with/against gravity [cm/s]')

