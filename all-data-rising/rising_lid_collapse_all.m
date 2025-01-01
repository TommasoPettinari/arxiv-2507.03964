

%% read data on stress and size dependence of rising behavior

fnames = dir('data_stress_size\*.dat');
numfiles = length(fnames);

lidweightarr = NaN(numfiles,1);
sizearr = NaN(numfiles,1);
addedweightarr = NaN(numfiles,1);
lengarr = NaN(numfiles,1);

timearr = NaN(numfiles, 50); % make enough space
displarr = NaN(numfiles, 50);
errarr = NaN(numfiles, 50);

% for creep curves
slopearr = NaN(numfiles,1);
offsetarr = NaN(numfiles,1);
errorarr = NaN(numfiles,2);

% for lin curves
slopearr_lin = NaN(numfiles,1);
offsetarr_lin = NaN(numfiles,1);
errorarr_lin = NaN(numfiles,2);
conf_lin = NaN(numfiles,2,2);

%settings to extract best fit exponent
bexpy = 0.2; % lower exp
expyr = 0.025; % step exponent
eexpy = 1.4; % max exp
expyarr = NaN(numfiles,1+round((eexpy-bexpy)/expyr)); % save all exp fit rsquares here
stdarr = NaN(numfiles,1+round((eexpy-bexpy)/expyr)); % stats of residuals
kurtarr = NaN(numfiles,1+round((eexpy-bexpy)/expyr)); % kurtosis stats of residuals
skewarr = NaN(numfiles,1+round((eexpy-bexpy)/expyr)); % skewness stats of residuals
distopyarr = NaN(numfiles,1+round((eexpy-bexpy)/expyr)); % distance to optimum

% extract prefactors for all nonlin fits (new for these data)
fitexpyarr = bexpy:expyr:eexpy;
slopearr_nonlin = NaN(numfiles,length(fitexpyarr));
offsetarr_nonlin = NaN(numfiles,length(fitexpyarr));
exparr_nonlin = NaN(numfiles,length(fitexpyarr));
conf_nonlin = NaN(numfiles,2,2);

bestexparr = NaN(numfiles,1); 
optexparr = NaN(numfiles,1);

typearr = []; % save type: cyl or sphere
lidonarr = []; % extract whether lid is on or not

for i=1:numfiles
    
    if ~mod(i,2), fprintf('numfile # %d\n', i); end
    
    % units all in centimeters and minutes!!!!!!!!!!
    readfile = char(fnames(i).name);
    readfiletotal = char(horzcat(fnames(i).folder,'\',fnames(i).name));
    
    
    typearr{i} = readfile(2:5);
    sizearr(i) = str2num(readfile(9:12))/100;
    
    tmp = importdata(readfiletotal);
    tmpdata = tmp.data;
    tmpstring = tmp.textdata;
    tmpstring = tmpstring{1};
    datl = size(tmpdata);
    
    
    %extract lid status
    if strfind(readfile,'l'),
        lidonarr(i) = 1;
    else
        lidonarr(i) = 0;
    end

    
    lengarr(i) = datl(1);
    timearr(i,1:datl) = tmpdata(:,1);
    % subtract the first data point; the beginning is arbitrary. This procedure is different from previous works as Tom is now video-analyzing data.
    displarr(i,1:datl) = tmpdata(:,2)-tmpdata(1,2); 

    if contains(tmpstring, 'err') % if there is an error estimate, it will be in the third column.
        errarr(i,1:datl(1)) = tmpdata(:,3);
    end
    
    %%%%%%%%%%%% we first assume square root behavior:
    % delta = sqrt(Dx)
    % delta^2 = Dx
    % x = time
    % delta = displacement

    xtmp = [timearr(i,1:lengarr(i))];
    deltatmp = [displarr(i,1:lengarr(i))-displarr(i,1)]; %do linear fit. Prefactor is D


    % fit poly1, offset should be zero hence added point at origin
    [fitobj,gof2] = fit(xtmp',deltatmp','poly1');

    tmp = coeffvalues(fitobj);    
    %prefactor is slope on y^2(t) plot
    slopearr(i) = tmp(1);
    offsetarr(i) = tmp(2);
    errorarr(i,1) = gof2.rmse;%/datl(1);
    errorarr(i,2) = gof2.rsquare;
    
    %%%%%%%%%%%%%%%% now fit the linear parts. may not exist for lid
    %%%%%%%%%%%%%%%% experiments
    deltatmp = [displarr(i,1:lengarr(i))];
    
    % fit poly1, offset should be zero hence added point at origin
    [fitobj,gof2] = fit(xtmp',deltatmp','poly1');
    
    %figure
    %plot(fitobj,xtmp',deltatmp')
    
    tmp = coeffvalues(fitobj);    
    slopearr_lin(i) = tmp(1);
    offsetarr_lin(i) = tmp(2);
    errorarr_lin(i,1) = gof2.rmse;%/datl(1);
    errorarr_lin(i,2) = gof2.rsquare;
    conf_lin(i,:,:) = confint(fitobj,0.95); 
    
    %%%%%%%%%%%%%%%% fit a range of different exponents just to see which
    %%%%%%%%%%%%%%%% might work best
    xtmp = timearr(i,2:lengarr(i));
    deltatmp = displarr(i,2:lengarr(i))-displarr(i,1);

    ctr = 0;
    for expy = bexpy:expyr:eexpy
        ctr = ctr + 1;
        opty = fitoptions('power2','Lower', [0,expy,-.0010], 'Upper', [Inf,expy,.0010]);
        %opty = fitoptions('power1','Lower', [0,expy], 'Upper', [Inf,expy]);
       %[fitobj,gof2] = fit(xtmp',deltatmp','power2',opty);
       [fitobj,gof2] = fit(xtmp',deltatmp','power2',opty);
       expyarr(i,ctr) = gof2.rsquare;
       %thought: check stats of residuals
       resi = fitobj(xtmp)-deltatmp';
       stdarr(i,ctr) = nanstd(resi);
       kurtarr(i,ctr) = kurtosis(resi);
       skewarr(i,ctr) = skewness(resi);

       coiffure = coeffvalues(fitobj);
       exparr_nonlin(i,ctr) = expy; 
       slopearr_nonlin(i,ctr) = coiffure(1);
       offsetarr_nonlin(i,ctr) = coiffure(3);
    end

    [~,indexy] = max(expyarr(i,:));
    bestexparr(i) = bexpy+(indexy-1)*expyr;
    
    distopyarr(i,:) = sqrt(skewarr(i,:).^2 + (kurtarr(i,:)-3).^2 + (gradient(expyarr(i,:)).^2));
    [~,indexy] = min(distopyarr(i,:));
    optexparr(i) = bexpy+(indexy-1)*expyr;
end

% fix zeros in timearr; some files contain spurious zeros
timearr(find(timearr == 0)) = NaN;
timearr(:,1) = 0;

fprintf('done\n')

close all

%% some admin stuff: extract buoyancy data

% buoyancy force comes from the weight of displaced volume minus the weight
% from intruder

% weight of displaced volume: rho * 4/3pir^3
% assume water density = 1000
% convert to grams for easier calcs
% sizearr has diameters in mm

wdisplvol = 1000*1000*0.333*pi()*4*(sizearr/2000).^3; % in grams, weight of displaced volume

% read the names and properties of the experimental data files
tmpdata = readmatrix(['data_stress_size\files_numbers.txt']);
tmpnames = importdata(['data_stress_size\files.txt']);

weightarr = NaN(1,numfiles);

% add all experimental info to the database
indy = [];
for i=1:numfiles
    
    if ~mod(i,2), fprintf('numfile # %d\n', i); end
    readfile = char(fnames(i).name);

    %walk through file list in tmpnames
    for j=1:numfiles+1
        if strcmp(readfile,char(tmpnames(j))),
            indy = [indy j]; 
            weightarr(i) = tmpdata(j,1); % weight of ball comes from filename, in grams
            surfy = pi()*(sizearr(i)/2000)^2; % surface area compute from diameter
            buoyarr(i) = 9.81*(wdisplvol(i)-weightarr(i))/1000/surfy; % buoyancy stress is g*grams/surface area
            if ~buoyarr(i); fprintf('error\n'); end % if something was wrong in the prep arrays, it will perhaps be printed
        end
    end

end

% Some more admin for buoyancy versus sinking rate analyses

% weight of displaced volume: rho * 4/3 pi r^3
% assume water density = 1000 kg/m3. Not technically correct due to
% presence of hydrogels, but good enough

% compute volume of all spheres
volsy = 0.333*pi()*4*(sizearr/2000).^3;
% weight of displaced volume
wdisplvol = 1000*volsy;
% density of spheres; convert weight to kg as well
denspart = weightarr./(1000*volsy');
% the delta rho calculation, using 1kg/L as reference for density of liquid
deltarho = 1000-denspart;
% computing the effective viscosity from the slope
etaconv = 2.*deltarho'.*(sizearr/2000).^2*9.81./9;

%% FIGURE 1 is a picture / description of the experiments

%% FIGURE 2 Single linear and non-linear to introduce results.

figure(1)

% select the data sets that are pairs of lid & no-lid:
% 2 & 24
% 3 & 25
% 4 & 26 <-- cherry picked for paper because it shows the clearest difference
% 5 & 27
% 6 % 28
% 7 & 29
% 8 & 30

maskex = 26; % the filename=hist12 example as picked

for i=1:1:length(maskex)
    
    scatter(timearr(maskex(i),:),displarr(maskex(i),:), '^k')
    hold on
    % a 5% was considered reasonable for these measurements, based on pixel
    % inaccuracy estimated in positioning
    errorbar(timearr(maskex(i),:),displarr(maskex(i),:),0.05*displarr(maskex(i),:),'k')
    
    %for optical improvement, also include a line between the data points
    plot(timearr(maskex(i),:),displarr(maskex(i),:), '-.k')

    % plot reference line from nonlin fit. #33 is just the linear option
    % in this set
    plot(timearr(maskex(i),:),slopearr_nonlin(maskex(i),33)*timearr(maskex(i),:),'-k')

    fprintf('lidstress = %d \n', buoyarr(maskex(i)))

end

%  then plot data from the same stress runs but then with a lid
maskex = maskex - 22; % select the one with lid

for i=1:1:length(maskex)
    
    scatter(timearr(maskex(i),:),displarr(maskex(i),:), 'ob')
    hold on
    errorbar(timearr(maskex(i),:),displarr(maskex(i),:),0.05*displarr(maskex(i),:),'b')

    plot(timearr(maskex(i),:),displarr(maskex(i),:), '-.b')
    fprintf('lidstress = %d \n', buoyarr(maskex(i)))

end


box on
xlabel('t [sec]')
ylabel('\delta [mm]')
text(2500,10,'(a)','FontSize',18)
xlim([0 300])
ylim([0, 110])
ax=gca;
ax.FontSize = 14;

%save
print('2-excomp', '-dpdf')


%% FIGURE 3 (a) t-series for pingpong with phi 6.35 3(b) slopes vs buoyancy for this data.
 

% first option: find all data with rising sphere data, and no lid present
% on the surface of the packing
lidchoice = 0;
masklid = find(lidonarr == lidchoice);
% excluding the non-PPball
masklid = find(lidonarr == lidchoice & sizearr' == 40);
% excluding all the runs after the lid/no-lid varianst with single phi for
% hydrogel density
masklid = [20,21,24:36];

% plot reference line from nonlin fit. #33 is just the linear option
% in this set
indexchoice = 33;
alphy = 1/fitexpyarr(indexchoice);


figure(31)

% show an overall family of pp ball data. color is stress
subplot 211
cmap = colormap(parula(length(unique(buoyarr(masklid)))));

for i=1:1:length(masklid)
    
    indexer = find(unique(buoyarr(masklid))==buoyarr(masklid(i)));
    
    scatter(timearr(masklid(i),:),displarr(masklid(i),:), 'k','MarkerEdgeColor',cmap(indexer,:))
    hold on
    plot(timearr(masklid(i),:),displarr(masklid(i),:), 'k','Color',cmap(indexer,:))

    if ~isnan(errarr(masklid(i))) 
  %      errorbar(timearr(masklid(i),:),displarr(masklid(i),:),errarr(masklid(i),:),'k')%,'Color',cmap(indexer,:))
        fprintf('errorrrrrrrrrr %d\n', expy)
    end
        
        

    
end

box on
xlabel('t [sec]')
ylabel('\delta [mm]')
text(1300,110,'(a)','FontSize',18)
xlim([0 1500])
ylim([0, 125])
ax=gca;
ax.FontSize = 14;
hold off
c = colorbar('EastOutside');
c.Label.String = '\sigma_S [Pa]';
caxis([0,250]);

% display the stress dependence of the effective viscosity
subplot 212

% masklid has the selection as chosen above
minerr = etaconv(masklid)./slopearr(masklid)-etaconv(masklid)./conf_lin(masklid,2,1);
maxerr = etaconv(masklid)./conf_lin(masklid,1,1)-etaconv(masklid)./slopearr(masklid);

errorbar(buoyarr(masklid),etaconv(masklid)./slopearr(masklid),minerr,maxerr,'ok')
hold on
scatter(buoyarr(masklid),etaconv(masklid)./slopearr(masklid), 90, buoyarr(masklid), 'o','filled')
errorbar(buoyarr(masklid),etaconv(masklid)./slopearr(masklid),minerr,maxerr,'.k')

box on
xlabel('\sigma_S [N/m^2]')
ylabel('\eta_{\rm eff} [Pa\cdot s]')
text(235,60,'(b)','FontSize',18)
set(gca,'Yscale','log')
ylim([.5,100]);
xlim([120,250]);
ax=gca;
ax.FontSize = 14;

hold on
xtmp = 10:1:250;
plot(xtmp,7E3*exp(-xtmp./27),'-.k')
hold off 

print(horzcat(('3-dense_ex'), num2str(lidchoice,2)), '-dpdf', '-bestfit')

%% Fig 4(a) Collapse for all pingpong ball data (6.35,6.25,6) 4(b) Slopes for all ping pong ball data

% excluding the non-PPball
masklid = find(lidonarr == lidchoice & sizearr' == 40);

cmap = colormap(parula(length(unique(buoyarr(masklid)))));

figure(33)
% show the rescaling of the sped with eta effective
subplot 211

for i=1:1:length(masklid) 
    
    indexer = find(unique(buoyarr(masklid))==buoyarr(masklid(i)));
    

    %scatter(timearr(masklid(i),:),(1./slopearr(masklid(i)))*displarr(masklid(i),:).^alphy, 'k','MarkerEdgeColor',cmap(indexer,:))
    %hold on
    %plot(timearr(masklid(i),:),(1./slopearr(masklid(i)))*displarr(masklid(i),:).^alphy, 'k','Color',cmap(indexer,:))

    scatter(timearr(masklid(i),:),(displarr(masklid(i),:)./slopearr_nonlin(masklid(i),indexchoice)).^alphy, 'k','MarkerEdgeColor',cmap(indexer,:))
    hold on

    if ~isnan(errarr(masklid(i))) 
        fprintf('errorrrrrrrrrr %d\n', expy)
    end
        
        
end

plot([1,10000],[1,10000],'-.k')

box on
xlabel('t [sec]')
ylabel(horzcat('(\delta\eta_{eff})^{',num2str(alphy),'} [A.U]'))
ax=gca;
ax.FontSize = 14;
hold off
xlim([1 10000])
ylim([1 10000])
set(gca,'Yscale','log')
set(gca,'Xscale','log')

text(3250,10,'(a)','FontSize',18)

c = colorbar('EastOutside');
c.Label.String = '\sigma_S [Pa]';
caxis([0,250]);

subplot 212

%scatter(buoyarr(13:57),etaconv(13:57)./slopearr(13:57), 90, buoyarr(13:57), 'o','filled')
% or exclude some
scatter(buoyarr(masklid),etaconv(masklid)./slopearr(masklid), 90, buoyarr(masklid), 'o','filled')

hold on

box on
xlabel('\sigma_S [N/m^2]')
ylabel('\eta_{\rm eff} [Pa\cdot s]')
text(225,30,'(b)','FontSize',18)
set(gca,'Yscale','log')
ylim([.01,100]);
xlim([30,250]);
ax=gca;
ax.FontSize = 14;

% plot three reference lines, slopes and offsets evaluated manually
% this concerns three different samples (6.35, 6.25, 6.0 g/L) so different
% slopes are expected
xtmp = 10:1:250;
plot(xtmp,7E3*exp(-xtmp./27),'-.k')
plot(xtmp,3E3*exp(-xtmp./17),'-.k')
plot(xtmp,3E1*exp(-xtmp./12),'-.k')


hold off

legend('lid','no lid','\propto e^{-\sigma/\sigma_0}','Location','northwest')

%save
print(horzcat(('4-risingex_rescale'), num2str(lidchoice,2)), '-dpdf', '-bestfit')

%% Figure 5 PIV data, Stokes solution etc

% 
% look elsewhere in repository

%% FIGURE 6 Temperature dependent sinking rate; 
% 
% look elsewhere in repository

%% FIGURE 7 Visualization of the surface flow
% 
% just sample images from typical experiments

%%  FIGURE 8 High speed video processing showing no stick slip i.e. smooth motion  (long timescale)  

figure(20) 

readfile = char("data_tommaso/no_stick_slip/smooth_rising.dat");
centerpos = readmatrix(readfile);

hs_pix2m = 40/730; %730 pixels correspond to 40 mm
displacement_z = centerpos(:,1);
displacement_x = centerpos(:,2);

pos0 = [.12,.14,.85,.82];
pos1 = [.17,.6,.3,.3];
pos2 = [.63,.20,.3,.3];
ax0 = axes('Position',pos0,'Box','on')
ax1 = axes('Position',pos1,'Box','on')
ax2 = axes('Position',pos2,'Box','on')

axes(ax0)

offset = displacement_z(1);
n_frames = length(displacement_z);
time = linspace(1, n_frames, n_frames);
displacement_z(1168) = NaN;      %removing one weird tracking outlier
scatter(time,abs(displacement_z-offset)*hs_pix2m)
hold on
plot(time,abs(displacement_z-offset)*hs_pix2m)


%plot(hs_time,gradient((hs_displacement_z-hs_displacement_z(1))*hs_pix2m,hs_time))
box on
xlim([-.5, 1600])
ylim([0, 130])
xlabel('t [sec]')
ylabel('\delta [mm]')
text(1000,120,'(a)','FontSize',25)
ax=gca;
ax.FontSize = 14;

axes(ax1)

%hist((hs_displacement_z(2:65)-mean(hs_displacement_z(2:65)))*hs_pix2m)

scatter(time,abs(displacement_z-offset)*hs_pix2m)
hold on
plot(time,abs(displacement_z-offset)*hs_pix2m)
plot([0,500],[10,50]+0.75,'-k')
text(35,45,'(b)','FontSize',14)



%plot(hs_time,gradient((hs_displacement_z-hs_displacement_z(1))*hs_pix2m,hs_time))
box on
xlim([-.5 500])
ylim([0, 55])
%xlabel('t [sec]')
%ylabel('\delta [mm]')

axes(ax2)

rangezoom = 450:650;

scatter(time(rangezoom),abs(displacement_z(rangezoom)-offset)*hs_pix2m)
hold on
plot(time(rangezoom),abs(displacement_z(rangezoom)-offset)*hs_pix2m)

box on
xlim([480, 650])
ylim([47,60])
%xlabel('t [sec]')
%ylabel('\delta [m]')
text(500,58,'(c)','FontSize',14)


%save
print('8-nostickslip_tommaso', '-dpdf')

%% FIGURE 9 Sinking towards a liquid layer - Fluid boundary at bottom  

figure(3)

pix_to_mm = 40/733;

fnames = dir('data_tommaso/sinking_with_liquid/*.dat');
numfiles = length(fnames);
interval = [5.0, 5.0, 10.0];
etas = [0.3, 0.3, 0.055];
colors = ["#0072BD","#EDB120","#77AC30"];

subplot 211

for i=1:numfiles
    
    readfile = char(strcat(fnames(i).folder,'/',fnames(i).name));
    tmp = readmatrix(readfile);
    
    displacement = tmp(:,1);
    offset = displacement(1);
    n_datapoints = length(displacement);
    timearry = linspace(0, (n_datapoints - 1)*(interval(i)), n_datapoints);

    scatter(timearry,(displacement-offset)*pix_to_mm, "MarkerEdgeColor",colors(i))
    hold on
    plot(timearry,displacement)
   
end

box on
xlabel('t [sec]')
ylabel('\delta [mm]')
%text(3250,.015,'(b)','FontSize',18)
%set(gca,'Yscale','log')
%set(gca,'Xscale','log')
xlim([0 2000])
ylim([0, 150])
ax=gca;
ax.FontSize = 14;
hold off

text(1800,30,'(a)','FontSize',14)

subplot 212

for i=1:numfiles
    
    readfile = char(strcat(fnames(numfiles - i + 1).folder,'/',fnames(numfiles - i + 1).name));
    tmp = readmatrix(readfile);
    
    displacement = tmp(:,1);
    offset = displacement(1);
    n_datapoints = length(displacement);
    timearry = linspace(0, (n_datapoints - 1)*(interval(numfiles - i + 1)), n_datapoints);

    scatter(timearry,(displacement-offset)*pix_to_mm/etas(numfiles - i + 1), "MarkerEdgeColor", colors(numfiles - i + 1))
    hold on
        
end

plot([1,2000],[1,2000],'-.k')

box on
%plot([10,3000],[10,3000],'-.k')
xlabel('t [sec]')
ylabel(horzcat('\delta/U [sec]'))
%text(3250,500,'(c)','FontSize',18)
ax=gca;
ax.FontSize = 14;
hold off
xlim([0 2000])
ylim([0, 3000])
%set(gca,'Yscale','log')
%set(gca,'Xscale','log')

text(1800,500,'(b)','FontSize',14)


%save
print('9-sinking with layer', '-dpdf', '-bestfit')

%% FIGURE 10 verify nonlinearity with lid 
% (a) Time series for cases with lid plus (b) collapse 

indexchoice = 23;

lidchoice = 1;
alphy = 1/fitexpyarr(indexchoice);

masklid = find(lidonarr == lidchoice);
% excluding the non-PPball
masklid = find(lidonarr == lidchoice & sizearr' == 40);

figure(10)

subplot 211
cmap = colormap(parula(length(unique(buoyarr(masklid)))));

for i=1:1:length(masklid)
    
    indexer = find(unique(buoyarr(masklid))==buoyarr(masklid(i)));
    

    scatter(timearr(masklid(i),:),displarr(masklid(i),:), 'k','MarkerEdgeColor',cmap(indexer,:))
    hold on
    plot(timearr(masklid(i),:),displarr(masklid(i),:), 'k','Color',cmap(indexer,:))

    if ~isnan(errarr(masklid(i))) 
  %      errorbar(timearr(masklid(i),:),displarr(masklid(i),:),errarr(masklid(i),:),'k')%,'Color',cmap(indexer,:))
        fprintf('errorrrrrrrrrr %d\n', expy)
    end
        
        

    
end

%legend('','exp 1', 'exp 0.5','bla','Location','southeast')

box on
xlabel('t [sec]')
ylabel('\delta [mm]')
%text(3250,.015,'(b)','FontSize',18)
%set(gca,'Yscale','log')
%set(gca,'Xscale','log')
xlim([0 4000])
ylim([.01, 110])
ax=gca;
ax.FontSize = 14;
hold off
c = colorbar('EastOutside');
c.Label.String = '\sigma_S [Pa]';
%colorbar('Ticks',[0,100,200,300],...
%         'TickLabels',{'0','100','200','300'})
caxis([0,250]);

text(3500,90,'(a)','FontSize',18)


subplot 212

for i=1:1:length(masklid) 
    
        indexer = find(unique(buoyarr(masklid))==buoyarr(masklid(i)));
    

    scatter(timearr(masklid(i),:),(displarr(masklid(i),:)./slopearr_nonlin(masklid(i),indexchoice)).^alphy, 'k','MarkerEdgeColor',cmap(indexer,:))
    hold on
    
    %scatter(timearr(masklid(i),:),(1./slopearr(masklid(i)))*displarr(masklid(i),:).^alphy, 'k','MarkerEdgeColor',cmap(indexer,:))
    %hold on
    %plot(timearr(masklid(i),:),(1./slopearr(masklid(i)))*displarr(masklid(i),:).^alphy, 'k','Color',cmap(indexer,:))

    if ~isnan(errarr(masklid(i))) 
  %      errorbar(timearr(masklid(i),:),(1./slopearr(masklid(i)))*displarr(masklid(i),:).^alphy,errarr(masklid(i),:),'k')%,'Color',cmap(indexer,:))
        fprintf('errorrrrrrrrrr %d\n', expy)
    end
        
        
end

plot([1,20000],0.8*[1,20000],'-.k')

box on
%plot([10,3000],[10,3000],'-.k')
xlabel('t [sec]')
ylabel(horzcat('(\delta/\eta_{eff})^{',num2str(alphy),'} [A.U]'))
%text(3250,500,'(c)','FontSize',18)
ax=gca;
ax.FontSize = 14;
hold off

caxis([0,250])
xlim([5 20000])
%ylim([1, 2E6])
set(gca,'Yscale','log')
set(gca,'Xscale','log')

text(10000,10,'(b)','FontSize',18)


%save
print(horzcat(('10-checknl'), num2str(lidchoice,2)), '-dpdf', '-bestfit')

%% FIGURE 11 plot of all stress dependence of buoyancy speed
% PP ball, smaller balls, lid, no-lid cases.

figure(11)

colormap parula

% plot data with lid, with linear fit (which doesn't work very well), 
% to show that the new lid data still shows the
% exponential stress dependence

plot(buoyarr(1:12),etaconv(1:12)./slopearr(1:12), '^',...    
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[1,1,0])
hold on

% plot(buoyarr(13:57),etaconv.*slopearr(13:57), 'o',...    
%     'LineWidth',2,...
%     'MarkerSize',10,...
%     'MarkerEdgeColor',[0,0,0],...
%     'MarkerFaceColor',[0.5,0.,0.2])
scatter(buoyarr(13:57),etaconv(13:57)./slopearr(13:57), 90, sizearr(13:57),'o','filled')

box on
xlabel('\sigma_S [N/m^2]')
ylabel('\eta_{\rm eff} [Pa\cdot s]')
%text(3250,.015,'(b)','FontSize',18)
set(gca,'Yscale','log')
ylim([.002,2E2])
xlim([30,250])
ax=gca;
ax.FontSize = 14;

c = colorbar('eastoutside');
c.Label.String = 'Intruder diameter [mm]';
c.Ticks = [20,25,40];

hold on

% plot three reference lines, slopes and offsets evaluated manually
xtmp = 10:1:250;
plot(xtmp,1E4*exp(-xtmp./23),'-.k')
plot(xtmp,1E3*exp(-xtmp./17),'-.k')
plot(xtmp,1E2*exp(-xtmp./10),'-.k')
plot(xtmp,5E2*exp(-xtmp./10),'-.k')

hold off

legend('lid','no lid','\propto e^{-\sigma/\sigma_0}','Location','northwest')

% save
print('11-risingspeed-both', '-dpdf', '-bestfit')



%%%%%%%%%%%%%%%%%%%%%% END OF PAPER PICTURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% START OF SOME USEFUL SUPPLEMENTARY PICTURES %%%%%%%%%%%%%%%%%%%%%

%% hydrogel level vs density

file_name = 'data_tommaso/rising_hydrogel/data.dat';

masses = [];
height_hydro = [];
height_water = [];
normalized_height = [];
volume = [];
density = [];

% Open the file and read the lines
fid = fopen(file_name, 'r');
lines = textscan(fid, '%s', 'Delimiter', '\n');
lines = lines{1};
fclose(fid);

% Remove the header line
lines(1) = [];

% Process each line
for count = 1:length(lines)
    obs = str2double(strsplit(lines{count}));
    
    masses(count) = obs(1);
    height_hydro(count) = obs(2);
    height_water(count) = obs(3);
    volume(count) = obs(4);
    
    normalized_height(count) = height_hydro(count) / height_water(count);
    density(count) = masses(count) / volume(count);
end

% Create the figure and plot
figure;
hold on;

box on

ylabel('h_h / h_w');
xlabel('\rho [g/L]');
xlim([0.5 3.5])
ylim([-0.1, 1.1])
ax=gca;
ax.FontSize = 14;

% Generate and plot the theoretical line
x1 = linspace(0.5, 3, 100);
y1 = 0.34 * x1;
plot(x1, y1, 'g--');

% Plot the data
plot(density, normalized_height, 'k+');

% Add horizontal and vertical lines
yline(1, 'r--');
xline(2.95, 'b');

hold off;

% save
print('level-hydrogel', '-depsc','-r600')
print('level-hydrogel', '-dpdf', '-bestfit')

%%%%% END OF SUPPL. PICTURES %%%%

