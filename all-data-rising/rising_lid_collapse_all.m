

%% load previous data from other papers and extract from it:
%- exponent arrays
%- sample data for comparison with current work



load("PRLdata.mat")

%check displarr
%line 11, 39, 86
displarr_ex = displarr([11,39,86],:);
timearr_ex = timearr([11,39,86],:);

save("examplecreep.mat", "displarr_ex", "timearr_ex", "-mat")

bestexparr_surftens = bestexparr;

clear bestexparr

%% then the SM paper
load("SMdata.mat")

bestexparr_lid = bestexparr;

bestexparr_combined = [bestexparr_lid ; bestexparr_surftens];

save("bestexparr_combined.mat", "bestexparr_combined","-mat")

clear all

%% some details regarding conversions etc for current data
D2si = 0.01^2/60; % convert creep slopes fitted to SI units (measurements are in cm/min)
Eta2si = 0.01/60; % convert linear slopes fitted to SI units (measurements are in cm/min)

w2pa_23 = 9.81*0.001/(3.1415*0.0115^2); % surface area of sphere 23mm

%% read data

fnames = dir('*.dat');
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

% extract prefactors for all non fits (new for these data)
fitexpyarr = bexpy:expyr:eexpy;
slopearr_nonlin = NaN(numfiles,length(fitexpyarr));
offsetarr_nonlin = NaN(numfiles,length(fitexpyarr));
exparr_nonlin = NaN(numfiles,length(fitexpyarr));

bestexparr = NaN(numfiles,1); 
optexparr = NaN(numfiles,1);

typearr = []; % save type: cyl or sphere
lidonarr = []; % extract whether lid is on or not

for i=1:numfiles
    
    if ~mod(i,2), fprintf('numfile # %d\n', i); end
    
    % units all in centimeters and minutes!!!!!!!!!!
    readfile = char(fnames(i).name);
    
    %lidweightarr(i) = str2num(readfile(2:4));
    typearr{i} = readfile(2:5);
    %lidonarr = [lidonarr, readfile(1)];
    sizearr(i) = str2num(readfile(9:12))/100;
    %addedweightarr(i) = str2num(readfile(8 :9));
    
    tmp = importdata(readfile);
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

    %figure;
    %plot(fitobj,xtmp,deltatmp)
    
    
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

%buoyancy force comes from the weight of displaced volume minus the weight
%from intruder

%weight of displaced volume: rho * 4/3pir^3
% assume water density = 1000
% convert to grams for easier calcs
% sizearr has diameters in mm

wdisplvol = 1000*1000*0.333*pi()*4*(sizearr/2000).^3; % in grams, weight of displaced volume

tmpdata = readmatrix(['files_numbers.txt']);
tmpnames = importdata(['files.txt']);

weightarr = NaN(1,numfiles);

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
            if ~buoyarr(i); fprintf('RRRRRR\n'); end % if something was wrong in the prep arrays, it will perhaps be printed
        end
    end

end

%% high speed data: extract sample from Tommaso's measurements July '24


tmp = importdata('2024-07-11 - 15.22.mqa');

hs_time = tmp.data(:,2); % already in seconds
hs_displacement_x = tmp.data(:,3);
hs_displacement_z = tmp.data(:,end);

hs_pix2m = 100/520; % the calibration is 500 pixels for 100mm.


%% FIGURE 1 make comparison of previous data and current data.

% relies on data extracted in first step of code. Extract delta(t)
% experiments at similar time/displacement levels

figure(1)

%subplot(2,2,[1 4])
%pos2 = [0.05 0.08 0.9 0.85];
%subplot('Position',pos2)

% first plot data from the current rising sphere data in black triangles
% select the data sets that are pairs of lid & no-lid:
% 2 & 24
% 3 & 25
% 4 & 26
% 5 & 27
% 6 % 28
% 7 & 29
% 8 & 30

maskex = 26; %the hist12 example Tom selected
%maskex = find(lidonarr == 0); % if all data needs to be displayed

for i=1:1:length(maskex)
    
    scatter(timearr(maskex(i),:),displarr(maskex(i),:), '^k')
    hold on
    plot(timearr(maskex(i),:),displarr(maskex(i),:), '-.k')

    plot(timearr(maskex(i),:),slopearr_nonlin(maskex(i),33)*timearr(maskex(i),:),'-k')

    fprintf('lidstress = %d \n', buoyarr(maskex(i)))

end

%  then plot data from the same stress runs but then with a lid
maskex = maskex - 22; % select the ones with lid
%maskex = find(lidonarr == 1);

for i=1:1:length(maskex)
    
    scatter(timearr(maskex(i),:),displarr(maskex(i),:), 'ob')
    hold on
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

%tmp = importdata('appar.jpg');

%pos1 = [0.05 0.08 0.4 0.85];
%subplot('Position',pos1)

%ax1 = axes('Position',[0.15 0.6 0.2 0.30]);

%timelapsy = imagesc(tmp);
%xticks([])
%yticks([])


%save
print('excomp', '-depsc','-r600')
print('excomp', '-dpdf')

%% demonstrate no stick slip

figure(10) 

pos0 = [.12,.14,.85,.82];
pos1 = [.17,.6,.3,.3];
pos2 = [.63,.20,.3,.3];
ax0 = axes('Position',pos0,'Box','on')
ax1 = axes('Position',pos1,'Box','on')
ax2 = axes('Position',pos2,'Box','on')

axes(ax0)

offsety = mean(hs_displacement_z(2:65));

scatter(hs_time-0.901409,(hs_displacement_z-offsety)*hs_pix2m)
hold on
plot(hs_time-0.901409,(hs_displacement_z-offsety)*hs_pix2m)


%plot(hs_time,gradient((hs_displacement_z-hs_displacement_z(1))*hs_pix2m,hs_time))
box on
xlim([-.5, 32])
ylim([0, 160])
xlabel('t [sec]')
ylabel('\delta [mm]')
text(27,130,'(a)','FontSize',25)
ax=gca;
ax.FontSize = 14;

axes(ax1)

%hist((hs_displacement_z(2:65)-mean(hs_displacement_z(2:65)))*hs_pix2m)

scatter(hs_time-0.901409,(hs_displacement_z-offsety)*hs_pix2m)
hold on
plot(hs_time-0.901409,(hs_displacement_z-offsety)*hs_pix2m)
plot([0,3],[4,21]+0.75,'-k')
text(-.3,12,'(b)','FontSize',14)



%plot(hs_time,gradient((hs_displacement_z-hs_displacement_z(1))*hs_pix2m,hs_time))
box on
xlim([-.5 2])
ylim([0, 16])
%xlabel('t [sec]')
%ylabel('\delta [mm]')

axes(ax2)

rangezoom = 500:600;

scatter(hs_time(rangezoom),(hs_displacement_z(rangezoom)-hs_displacement_z(1))*hs_pix2m)
hold on
plot(hs_time(rangezoom),(hs_displacement_z(rangezoom)-hs_displacement_z(1))*hs_pix2m)

box on
xlim([7.1, 8.4])
ylim([34,41])
%xlabel('t [sec]')
%ylabel('\delta [m]')
text(8.1,36,'(c)','FontSize',14)


%save
print('nostickslip', '-depsc','-r600')
print('nostickslip', '-dpdf')

%% FIGURE 3 make overview figure without lid 
% 

lidchoice = 0;


masklid = find(lidonarr == lidchoice);

if lidchoice % 1 is with lid 
    alphy = 2;
    prefy = 100;
else % alphy = 1.5 gives a better collapse, but the asymptote really seems to be a slope of 1.
    alphy = 1;
    prefy = 1;
end

indexchoice = 33;
alphy = 1/fitexpyarr(indexchoice);


figure(3)

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
%text(3250,.015,'(b)','FontSize',18)
%set(gca,'Yscale','log')
%set(gca,'Xscale','log')
xlim([0 4000])
ylim([0, 150])
ax=gca;
ax.FontSize = 14;
hold off
c = colorbar('EastOutside')
c.Label.String = '\sigma_S [Pa]';
%colorbar('Ticks',[0,100,200,300],...
%         'TickLabels',{'0','100','200','300'})
caxis([0,250])


subplot 212

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

plot([1,10000],prefy*[1,10000],'-.k')

box on
%plot([10,3000],[10,3000],'-.k')
xlabel('t [sec]')
ylabel(horzcat('(\delta/\eta_{eff})^{',num2str(alphy),'} [A.U]'))
%text(3250,500,'(c)','FontSize',18)
ax=gca;
ax.FontSize = 14;
hold off
%xlim([5 200])
%ylim([5, 150])
set(gca,'Yscale','log')
set(gca,'Xscale','log')

%save
print(horzcat(('risingex'), num2str(lidchoice,2)), '-depsc','-r600')
print(horzcat(('risingex'), num2str(lidchoice,2)), '-dpdf', '-bestfit')

%% FIGURE 2 plot stress dependence of buoyancy speed

figure(2)

colormap parula

% weight of displaced volume: rho * 4/3pir^3
% assume water density = 1000 kg/m3. Not technically correct due to
% presence of hydrogels, but good enough


volsy = 0.333*pi()*4*(sizearr/2000).^3;
wdisplvol = 1000*volsy;
denspart = weightarr./(1000*volsy');
deltarho = 1000-denspart;
etaconv = 2.*deltarho*(sizearr/2000).^2*9.81./9;

plot(buoyarr(1:12),etaconv./slopearr(1:12), '^',...    
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[1,1,0])
hold on
%scatter(buoyarr(1:12),etaconv.*slopearr(1:12), 90, sizearr(1:12),'^','filled')


% plot(buoyarr(13:57),etaconv.*slopearr(13:57), 'o',...    
%     'LineWidth',2,...
%     'MarkerSize',10,...
%     'MarkerEdgeColor',[0,0,0],...
%     'MarkerFaceColor',[0.5,0.,0.2])
scatter(buoyarr(13:57),etaconv./slopearr(13:57), 90, sizearr(13:57),'o','filled')



box on
xlabel('\sigma_S [N/m^2]')
ylabel('\eta_{\rm eff} [Pa\cdot s]')
%text(3250,.015,'(b)','FontSize',18)
set(gca,'Yscale','log')
%set(gca,'Xscale','log')
%xlim([0 2500])%xlim([20,6000])
ylim([1,1E4])%ylim([3E-3,.15])
ax=gca;
ax.FontSize = 14;

c = colorbar('eastoutside');
c.Label.String = 'Intruder diameter [mm]';
c.Ticks = [20,25,40];

hold on
xtmp = 10:1:250;
plot(xtmp,1E6*exp(-xtmp./23),'-.k')
plot(xtmp01E5*exp(-xtmp./17),'-.k')
plot(xtmp,1E4*exp(-xtmp./10),'-.k')
plot(xtmp,5E4*exp(-xtmp./10),'-.k')

hold off

legend('lid','no lid','\propto e^{-\sigma/\sigma_0}','Location','northwest')




% save
%print('risingspeed', '-depsc','-r600')
print('risingspeed-both', '-dpdf', '-bestfit')

    
%% FIGURE 4 verify sqrt with lid 
% 

indexchoice = 23;

lidchoice = 1;
alphy = 1/fitexpyarr(indexchoice);

masklid = find(lidonarr == lidchoice);

figure(4)

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

if lidchoice
    title('with lid')
else
    title('without lid')
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
colorbar('EastOutside')
%colorbar('Ticks',[0,100,200,300],...
%         'TickLabels',{'0','100','200','300'})
%caxis([0,350])


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
colorbar('EastOutside')
colorbar('Ticks',[0,100,200,300],...
         'TickLabels',{'0','100','200','300'})
caxis([0,350])
xlim([5 20000])
%ylim([1, 2E6])
set(gca,'Yscale','log')
set(gca,'Xscale','log')

title(horzcat(('exp ='), num2str(1/alphy,2)))

%save
print(horzcat(('checksq'), num2str(lidchoice,2)), '-depsc','-r600')
print(horzcat(('checksq'), num2str(lidchoice,2)), '-dpdf', '-bestfit')

%% make figure for speeds 
% 

lidchoice = 1;
exymple = 4;

masklid = find(lidonarr == lidchoice);

if lidchoice % 1 is with lid 
    alphy = 2;
    prefy = 100;
else % no lid
    alphy = 1.47;
    prefy = 10;
end


figure(8)

subplot 221
cmap = colormap(parula(length(unique(buoyarr(masklid)))));

for i=1:1:length(masklid)
    
    indexer = find(unique(buoyarr(masklid))==buoyarr(masklid(i)));
    

    scatter(timearr(masklid(i),:),gradient(displarr(masklid(i),:),timearr(masklid(i),:)), 'k','MarkerEdgeColor',cmap(indexer,:))
    hold on
    plot(timearr(masklid(i),:),gradient(displarr(masklid(i),:),timearr(masklid(i),:)), 'k','Color',cmap(indexer,:))

    if ~isnan(errarr(masklid(i))) 
        %errorbar(timearr(masklid(i),:),displarr(masklid(i),:),errarr(masklid(i),:),'k')%,'Color',cmap(indexer,:))
        fprintf('errorrrrrrrrrr %d\n', expy)
    end
          
end

%plot([1,3000],0.1*[1,3000],'-.k')
%plot([1,3000],4*[1,3000].^0.5,'-k')

if lidchoice
    title('with lid')
else
    title('without lid')
end

%legend('','exp 1', 'exp 0.5','bla','Location','southeast')

box on
xlabel('t [sec]')
ylabel('d\delta/dt [mm/s]')
%text(3250,.015,'(b)','FontSize',18)
set(gca,'Yscale','log')
set(gca,'Xscale','log')
xlim([5 10000])
ylim([2E-3, 2E0])
ax=gca;
ax.FontSize = 14;
hold off
colorbar('EastOutside')
colorbar('Ticks',[0,100,200,300],...
         'TickLabels',{'0','100','200','300'})
%caxis([0,350])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% single example


subplot 222
cmap = colormap(parula(length(unique(buoyarr(masklid)))));

for i=exymple:exymple
    
    indexer = find(unique(buoyarr(masklid))==buoyarr(masklid(i)));
    

    scatter(timearr(masklid(i),:),gradient(displarr(masklid(i),:),timearr(masklid(i),:)), 'k','MarkerEdgeColor',cmap(indexer,:))
    hold on
    plot(timearr(masklid(i),:),gradient(displarr(masklid(i),:),timearr(masklid(i),:)), 'k','Color',cmap(indexer,:))

    if ~isnan(errarr(masklid(i))) 
        errorbar(timearr(masklid(i),:),displarr(masklid(i),:),errarr(masklid(i),:),'k')%,'Color',cmap(indexer,:))
        fprintf('errorrrrrrrrrr %d\n', expy)
    end
        
        

    
end

%plot([1,3000],0.1*[1,3000],'-.k')
%plot([1,3000],4*[1,3000].^0.5,'-k')

if lidchoice
    title('with lid')
else
    title('without lid')
end

%legend('','exp 1', 'exp 0.5','bla','Location','southeast')

box on
xlabel('t [sec]')
ylabel('d\delta/dt [mm/s]')
%text(3250,.015,'(b)','FontSize',18)
%xlim([0 3000])
ylim([0,1])
ax=gca;
ax.FontSize = 14;
hold off
colorbar('EastOutside')
colorbar('Ticks',[0,100,200,300],...
         'TickLabels',{'0','100','200','300'})


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% as function of delta

subplot(2,2,[3 4])
cmap = colormap(parula(length(unique(buoyarr(masklid)))));

for i=1:1:length(masklid)
    
    indexer = find(unique(buoyarr(masklid))==buoyarr(masklid(i)));
    

    scatter(displarr(masklid(i),:),gradient(displarr(masklid(i),:),timearr(masklid(i),:)), 'k','MarkerEdgeColor',cmap(indexer,:))
    hold on
    plot(displarr(masklid(i),:),gradient(displarr(masklid(i),:),timearr(masklid(i),:)), 'k','Color',cmap(indexer,:))

    if ~isnan(errarr(masklid(i))) 
        %errorbar(timearr(masklid(i),:),displarr(masklid(i),:),errarr(masklid(i),:),'k')%,'Color',cmap(indexer,:))
        fprintf('errorrrrrrrrrr %d\n', expy)
    end
        
        

    
end

%plot([1,3000],0.1*[1,3000],'-.k')
%plot([1,3000],4*[1,3000].^0.5,'-k')

if lidchoice
    title('with lid')
else
    title('without lid')
end

%legend('','exp 1', 'exp 0.5','bla','Location','southeast')

box on
xlabel('\delta [mm]')
ylabel('d\delta/dt [mm/s]')
%text(3250,.015,'(b)','FontSize',18)
set(gca,'Yscale','log')
%set(gca,'Xscale','log')
%xlim([5 10000])
ylim([1E-3, 2E0])
ax=gca;
ax.FontSize = 14;
hold off
colorbar('EastOutside')
colorbar('Ticks',[0,100,200,300],...
         'TickLabels',{'0','100','200','300'})
%caxis([0,350])


%save
print(horzcat(('fig-vel-lid ='), num2str(lidchoice,2)), '-depsc','-r600')
print(horzcat(('fig-vel-lid ='), num2str(lidchoice,2)), '-dpdf', '-bestfit')

%% make figure for single example speeds combine lid/no lid 

lidchoice = 0;
exymple = 7;
masklid = find(lidonarr == lidchoice);

if lidchoice % 1 is with lid 
    alphy = 2;
    prefy = 100;
else
    alphy = 1.47;
    prefy = 10;
end

fprintf('lid buoyarr = %d\n', buoyarr(masklid(exymple)))


figure(9)
cmap = colormap(parula(length(unique(buoyarr(masklid)))));

for i=exymple:exymple
    
    indexer = find(unique(buoyarr(masklid))==buoyarr(masklid(i)));
    

    scatter(displarr(masklid(i),:),gradient(displarr(masklid(i),:),timearr(masklid(i),:)), 'o','MarkerEdgeColor',cmap(indexer,:))
    hold on
    %plot(timearr(masklid(i),:),gradient(displarr(masklid(i),:),timearr(masklid(i),:)), 'o','Color',cmap(indexer,:))

    if ~isnan(errarr(masklid(i))) 
        errorbar(timearr(masklid(i),:),displarr(masklid(i),:),errarr(masklid(i),:),'o')%,'Color',cmap(indexer,:))
        fprintf('errorrrrrrrrrr %d\n', expy)
    end
        
        

    
end

lidchoice = 0;
masklid = find(lidonarr == lidchoice);
map = colormap(parula(length(unique(buoyarr(masklid)))));

for i=exymple:exymple
    
    indexer = find(unique(buoyarr(masklid))==buoyarr(masklid(i)));
    

    scatter(displarr(masklid(i),:),gradient(displarr(masklid(i),:),timearr(masklid(i),:)), '>',...
                                                        'MarkerEdgeColor',cmap(indexer,:),'MarkerFaceColor',cmap(indexer,:))
    hold on
    %plot(timearr(masklid(i),:),gradient(displarr(masklid(i),:),timearr(masklid(i),:)), '>','Color',cmap(indexer,:))

    if ~isnan(errarr(masklid(i))) 
        errorbar(timearr(masklid(i),:),displarr(masklid(i),:),errarr(masklid(i),:),'>')%,'Color',cmap(indexer,:))
        fprintf('errorrrrrrrrrr %d\n', expy)
    end  

    
end

fprintf('no lid buoyarr = %d\n', buoyarr(masklid(exymple)))


box on
xlabel('\delta [mm]')
ylabel('d\delta/dt [mm/s]')
%text(3250,.015,'(b)','FontSize',18)
%xlim([0 3000])
ylim([0,.2])
ax=gca;
ax.FontSize = 14;
hold off
%colorbar('EastOutside')
%colorbar('Ticks',[0,100,200,300],...
%         'TickLabels',{'0','100','200','300'})

legend('lid','no lid','Location','northeast')



%save
print(horzcat(('fig-comb-lid ='), num2str(lidchoice,2)), '-depsc','-r600')
print(horzcat(('fig-comb-lid ='), num2str(lidchoice,2)), '-dpdf', '-bestfit')


%% ADDITIONAL FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% at the moment these are not included in the article, but show useful data
% nonetheless.

%% plot distribution of best fitted exponents

% for just the rising data:
% trend is clearly with peak around 1, but sigma is large. 
% what can be included is the data from the first paper and the second
% paper: those have clear peaks around 0.5 (but also broad).


% get data from previous experiments
load("bestexparr_combined.mat")

% the historical sinking data has many 0.2 values as exponent. These skew
% the statistics. Otherwise they are peaked around 0.5
bestexparrall = [bestexparr_combined; bestexparr;];

figure(1)

subplot 121

plot(1:256,bestexparrall(1:256), 'o')
hold on
plot(257:length(bestexparrall),bestexparrall(257:end), 'r>')
hold off

xlabel('experiment #')
ylabel('best exponent [-]')
legend('lid','no lid')



subplot 122

binny = bexpy:3*expyr:eexpy;

histogram(bestexparrall(1:244),binny)
hold on
histogram(bestexparrall(245:256),binny)
histogram(bestexparrall(257:end)',binny)

hold off

% save
print('exponent-lid', '-depsc','-r600')
print('exponent-lid', '-dpdf', '-bestfit')

%% check exponent dependence on stress

figure

surfpp = (3.1415*0.0215^2);
plot(buoyarr(1:12)./surfpp,bestexparr(1:12), 'o')
hold on
plot(buoyarr(13:end)./surfpp,bestexparr(13:end), 'r>')

xlabel('Buoyancy Stress [N/m^2]')
ylabel('best exponent [-]')
legend('lid','no lid')

% save
%print('risingspeed', '-depsc','-r600')
print('stress-vs-exponent', '-dpdf', '-bestfit')

%% END OF USEFUL FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)

% first plot the D dynamics for both cases

%subplot(2,2,[1 2])

tmpspherestress = unique(spherestressarr(masksp));
cmap = colormap(jet(length(tmpspherestress)));


for fixmass  = 0:600
 
    for i=1:length(masksp)
        
        indexer = find(unique(spherestressarr(masksp))==spherestressarr(masksp(i)));


        if addedweightarr(masksp(i)) == fixmass,
            psh1 = scatter([lidstressarr(masksp(i))],[slopearr(masksp(i))*D2si],[30],[spherestressarr(masksp(i))],'Fill',...
            'MarkerEdgeColor',cmap(indexer,:),...
            'MarkerFaceColor',cmap(indexer,:),...
            'LineWidth',1.5);
        
            if lidonarr(i) == 'l'
                psh1.Marker = 'd';
            end

            hold on
            %plot(timearr(maskcyl30(i),:)*60,0.01*(timearr(maskcyl30(i),:)*slopearr(maskcyl30(i))).^0.5,'-.','Color',cmap(indexer,:))
              
            sigL = 50:1:600; % draw a line
            sigS = spherestressarr(masksp(i));
            Dfit = D0*exp(prefy*(sigS./sigS0 - (sigL - sigLb)./sigL0));
            plot(sigL,Dfit,'Color',cmap(indexer,:))
        end

    end
    
end


xlim([0,600])
ylim([1e-9,1e-4])
set(gca,'Yscale','log')
xlabel('surface stress [Pa]')
ylabel('D [m^2/s]')
box on
hold off
c = colorbar('EastOutside','Ticks',[0,300,600],...
         'TickLabels',{'0','300','600'});
c.Label.String = 'Stress on sphere [Pa]';
caxis([0,600]);



%text(7.0,1E-5,'\sigma = 355 Pa','Color','blue','FontSize',14)
%text(6.6,7E-7,'\sigma = 247 Pa','Color','red','FontSize',14)


% save
print('fig2', '-depsc','-r600')
print('fig2', '-dpdf', '-bestfit')


%% how does D depend on added weights on sphere, for different 

figure(3)

% first sphere data

subplot(2,1,1);

tmplidstress = unique(lidstressarr(masksp));
cmap = colormap(jet(length(tmplidstress)));

%xtmp = 300:1:600;
%plot(xtmp,0.65E-10*exp(xtmp./35),'-.k')
%hold on


for i=1:length(masksp)
        
    indexer = find(unique(lidstressarr(masksp))==lidstressarr(masksp(i)));
    
    if lidweightarr(masksp(i)) < lintreshsp, %== 6.96,
        %scatter([spherestressarr(masksp(i))],[slopearr(masksp(i))*D2si],[90],[lidstressarr(masksp(i))],'Fill');
        psh1 = scatter([spherestressarr(masksp(i))],[slopearr(masksp(i))*D2si],[90],'Fill',...
            'MarkerEdgeColor',cmap(indexer,:),...
            'MarkerFaceColor',cmap(indexer,:),...
            'LineWidth',1.5');
            if lidonarr(i) == 'l'
                psh1.Marker = 'd';
            end
            hold on
    end
    
    sigL = lidstressarr(masksp(i));
    if sigL >= 200 % above threshold, lid stress dep exists
        D0 = 5E-11;
        sigS0 = 45;
        sigL0 = 30;
        sigLb = 275;
        sigS = 50:1:1000;
        Dfit = D0*exp(prefy*(sigS./sigS0 - (sigL - sigLb)./sigL0)); 
        plot(sigS,Dfit,'Color',cmap(indexer,:)); 
    else % below threshold, no lid stress dependence
        sigS0 = 22;
        sigS = 0:1:400;
        Dfit = D0*exp(prefy*(sigS./sigS0)); 
        plot(sigS,Dfit,'Color',cmap(indexer,:)); 
    end

end

plot([0,1000],[1E-7,1E-7],'-.k')

xlim([0,1000])
ylim([1e-9,1E-4])
set(gca,'Yscale','log')
xlabel('sphere stress [Pa]')
ylabel('D [m^2/s]')
%text(485,1E-6,'(a)','FontSize',18)
box on
colorbar('EastOutside')
colorbar('Ticks',[0,100,200,300,400,500],...
         'TickLabels',{'0','100','200','300','400','500'})
caxis([min(tmplidstress),max(tmplidstress)]);
ax=gca;
ax.FontSize = 14;

hold off

subplot(2,1,2);

tmplidstress = unique(lidstressarr(masksp));
cmap = colormap(jet(length(tmplidstress)));

%xtmp = 300:1:600;
%plot(xtmp,0.65E-10*exp(xtmp./35),'-.k')
%hold on


for i=1:length(masksp)
        
    indexer = find(unique(lidstressarr(masksp))==lidstressarr(masksp(i)));
    
    %if lidweightarr(masksp(i)) > 100,
        %scatter([spherestressarr(masksp(i))],[slopearr(masksp(i))*D2si],[90],[lidstressarr(masksp(i))],'Fill');
        
    %end
    
    sigL = lidstressarr(masksp(i));
    if sigL >= 0 % above threshold, lid stress dep exists
        D0 = 5E-11;
        sigS0 = 45;
        sigL0 = 30;
        sigLb = 250;
        sigS = 50:1:1000;
        Dfit = D0*exp(prefy*(sigS./sigS0 - (sigL - sigLb)./sigL0)); 
        %plot(sigS,Dfit,'Color',cmap(indexer,:)); 

        % rescalecoeff = D0*exp(prefy*(- (sigL - sigLb)./sigL0)); % works
        % nicely, but doesn't capture change in slope
        rescalecoeff = D0*exp(prefy*(- (sigL)./sigL0));
        
        %or try more subtle rescaling, in which sigL0 depends on sigL
        %D0 = 5E-11;
        %sigS0 = 45;
        %sigL0 = 24;
        %sigLb = 275;
        %sigS = 50:1:1000;
        %Dfit = D0*exp(prefy*(sigS./sigS0 - (sigL - sigLb)./sigL0));
        % rescalecoeff = D0*exp(prefy*(-sigL./(sigL0+0.011*sigL)));

        
        psh1 = scatter([spherestressarr(masksp(i))],[slopearr(masksp(i))*D2si]./rescalecoeff,[90],'Fill',...
        'MarkerEdgeColor',cmap(indexer,:),...
        'MarkerFaceColor',cmap(indexer,:),...
        'LineWidth',1.5');
        if lidonarr(i) == 'l'
            psh1.Marker = 'd';
        end
        hold on
    end

end

%plot([0,1000],[1E-7,1E-7],'-.k')

%xlim([0,1000])
ylim([1e2,1E15])
set(gca,'Yscale','log')
xlabel('sphere stress [Pa]')
ylabel('D/(D_0exp((\sigma_L-\sigma_{Lb})/\sigma_{L0})) [-]')
%text(485,1E-6,'(a)','FontSize',18)
box on
colorbar('EastOutside')
colorbar('Ticks',[0,100,200,300,400,500],...
         'TickLabels',{'0','100','200','300','400','500'})
caxis([min(tmplidstress),max(tmplidstress)]);
ax=gca;
ax.FontSize = 14;

hold off


% save
print('fig3', '-depsc','-r600')
print('fig3', '-dpdf', '-bestfit')

%% how does D depend on added weights on sphere, for different 

figure(4)

% first sphere data

%bluey = subplot(2,1,1);

%xtmp = 300:1:600;
%plot(xtmp,0.65E-10*exp(xtmp./35),'-.k')
%hold on


for i=1:length(masksp)

    if lidweightarr(masksp(i)) < lintreshsp, %== 6.96,
        scatter([lidstressarr(masksp(i))],[slopearr(masksp(i))*D2si],[30],[spherestressarr(masksp(i))],'Fill')

        hold on
    end

    sigS = spherestressarr(masksp(i));
    if sigS >= 0 % above threshold, lid stress dep exists
        D0 = 5E-11;
        sigS0 = 45;
        sigL0 = 30;
        sigLb = 275;
        sigL = 0:1:500;
        Dfit = D0*exp(prefy*(sigS./sigS0 - (sigL - sigLb)./sigL0)); 
        plot(sigL,Dfit,'Color',cmap(indexer,:)); 
    else % below threshold, no lid stress dependence
        %sigS0 = 22;
        %sigS = 0:1:400;
        %Dfit = D0*exp(prefy*(sigS./sigS0)); 
        %plot(sigS,Dfit,'Color',cmap(indexer,:)); 
    end

end

%xlim([340,510])
%ylim([5e-7,1E-4])
set(gca,'Yscale','log')
xlabel('lid stress [Pa]')
ylabel('D [m^2/s]')
%text(485,1E-6,'(a)','FontSize',18)
box on
colorbar('EastOutside')
colorbar('Ticks',[0,100,500,1000],...
         'TickLabels',{'0','100','500','1000'})
clim([0,1000])
ax=gca;
ax.FontSize = 14;



% save
%print('fig4', '-depsc','-r600')
%print('fig4', '-dpdf', '-bestfit')

%% other stuff
bluey = subplot(2,1,2);


for fixmass  = 0:1:60
 
    for i=1:length(masksp)

        if addedweightarr(masksp(i)) == fixmass,
            scatter([lidstressarr(masksp(i))],[slopearr(masksp(i))*D2si]/(platysp*exp(spherestressarr(masksp(i))/35)),[30],[spherestressarr(masksp(i))],'Fill')

            hold on
        end

    end

end
%xlim([6.7,7.17])
%ylim([1e-7,1E-4])
set(gca,'Yscale','log')
xlabel('lid stress [Pa]')
ylabel('D/D_0 exp(\sigma/\sigma_s) [-]')
%text(7.1,4E-5,'(c)','FontSize',18)
box on
colorbar('EastOutside')
%caxis([350,450])
ax=gca;
ax.FontSize = 14;

MapRed = [linspace(0,1,256).^1',zeros(256,2)];
MapBlue = [zeros(256,2),linspace(0,1,256)'];
colormap(bluey, MapBlue);

% save
print('fig3', '-depsc','-r600')
print('fig3', '-dpdf', '-bestfit')


%% SI figure log-log scale version Fig 1


figure(5)

subplot 221
cmap = colormap(parula(length(unique(lidweightarr(masksp)))));

for i=1:4:length(masksp)
    
    %indexer = find(lidarruniq==lidarr(masksp(i)));
    indexer = find(unique(lidweightarr(masksp))==lidweightarr(masksp(i)));
    
    if lidweightarr(masksp(i)) < lintreshsp
       
        %scatter(timearr(masksp(i),:)*60,displarr(masksp(i),:)*0.01,'MarkerEdgeColor',cmap(indexer,:))
        scatter(timearr(masksp(i),:)*60,displarr(masksp(i),:)*0.01, 'k'); %,'MarkerEdgeColor',cmap(indexer,:))
        hold on
        plot(timearr(masksp(i),:)*60,displarr(masksp(i),:)*0.01, 'k')%,'Color',cmap(indexer,:))

        if ~isnan(errarr(masksp(i))) 
            %errorbar(timearr(masksp(i),:)*60,displarr(masksp(i),:)*0.01,errarr(masksp(i),:)*0.01,'Color',cmap(indexer,:))
            errorbar(timearr(masksp(i),:)*60,displarr(masksp(i),:)*0.01,errarr(masksp(i),:)*0.01,'k')%,'Color',cmap(indexer,:))
            fprintf('errorrrrrrrrrr %d\n', expy)
        end
        
        
    end
    
end

set(gca,'Xscale','log')
set(gca,'Yscale','log')

plot([10,100000],0.0015*[10,100000].^0.5,'-.k')
plot([10,100000],0.0015*[10,100000].^0.4,':k')
plot([10,100000],0.0015*[10,100000].^0.6,':k')


box on
xlabel('t [sec]')
ylabel('\delta [m]')
text(2000,.0035,'(a)','FontSize',18)
%set(gca,'Yscale','log')
%set(gca,'Xscale','log')
xlim([20, 5000])
ylim([0.002,0.2])
ax=gca;
ax.FontSize = 14;
hold off



subplot 222

for i=length(masksp):-1:1
    
    %indexer = find(lidarruniq==lidarr(masksp(i)));
    indexer = find(unique(lidweightarr(masksp))==lidweightarr(masksp(i)));
    
    if lidweightarr(masksp(i)) < lintreshsp
       
        %plot(timearr(masksp(i),:)*60,displarr(masksp(i),:).^2/(10000*slopearr(masksp(i))*D2si),'-o','Color',cmap(indexer,:),'MarkerSize', 6)
        plot(timearr(masksp(i),:)*60,displarr(masksp(i),:).^2/(10000*slopearr(masksp(i))*D2si),'-ok')
        %loglog(timearr(masksp(i),:)*60,displarr(masksp(i),:).^2/(10000*slopearr(masksp(i))*D2si),'-ok')

        %convert D to proper units, but also the war time and displacement data
        hold on
        
    end
    
end
plot([10,10000],[10,10000],'-.k')
xlabel('t [sec]')
ylabel('\delta^2/D [sec]')
text(2000,8.5,'(b)','FontSize',18)
xlim([20, 5000])
ylim([4,5000])
ax=gca;
ax.FontSize = 14;
hold off

set(gca,'Xscale','log')
set(gca,'Yscale','log')

% save
print('figSI', '-depsc','-r600')
print('figSI', '-dpdf', '-bestfit')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Sieve: make the first overview figure 


figure(1)

subplot 221
%cmap = colormap(parula(length(unique(lidarr(masksv)))));
cmap = colormap(parula(10));

for i=1:1:length(masksv)
    
    %indexer = find(lidarruniq==lidarr(masksv(i)));
    indexer = find(lidarruniq==lidweightarr(masksv(i)));
    
    if lidweightarr(masksv(i)) < 5.5 %< lintreshsp
      
        scatter(timearr(masksv(i),:)*60,displarr(masksv(i),:)*0.01,'MarkerEdgeColor',cmap(i,:),'MarkerFaceColor',cmap(i,:))

        hold on
        plot(timearr(masksv(i),:)*60,displarr(masksv(i),:)*0.01,'Color',cmap(indexer,:))

        if ~isnan(errarr(masksv(i))) 
            errorbar(timearr(masksv(i),:)*60,displarr(masksv(i),:)*0.01,errarr(masksv(i),:)*0.01,'Color',cmap(indexer,:))
            fprintf('errorrrrrrrrrr %d\n', expy)
        end
        
        
    end
    
end

%plot(6*[10,1000],0.001*[10,1000].^0.6, '-.k')
%plot(6*[10,1000],0.001*[10,1000].^0.5, '-.k')

box on
xlabel('t [sec]')
ylabel('\delta [m]')
%text(3250,.015,'(b)','FontSize',18)
%set(gca,'Yscale','log')
%set(gca,'Xscale','log')
xlim([0 4000])%xlim([20,6000])
ylim([0,0.1])%ylim([3E-3,.15])
hold off



subplot 222

for i=length(masksv):-1:1
    
    %indexer = find(lidarruniq==lidarr(masksv(i)));
    indexer = find(unique(lidweightarr(masksv))==lidweightarr(masksv(i)));
    
    if lidweightarr(masksv(i)) < 5.5; %lintreshsp
       
        plot(timearr(masksv(i),:)*60,displarr(masksv(i),:).^2/(10000*slopearr(masksv(i))*D2si),'Color',cmap(i,:))
        %convert D to proper units, but also the war time and displacement data
        hold on
        
    end
    
end
plot([0,4000],[0,4000],'-.k')
xlabel('t [sec]')
ylabel('\delta^2/D [sec]')
%text(3250,500,'(c)','FontSize',18)
hold off

subplot 223


xtmp = 300:1:600;
plot(xtmp,0.5E-10*exp(xtmp./50),'-.k')
hold on

for i=1:length(masksv)

    if lidweightarr(masksv(i)) < 5.5 
        scatter([spherestressarr(masksv(i))],[slopearr(masksv(i))*D2si],[30],[lidweightarr(masksv(i))],'Fill')

        hold on
    end

end

%xlim([340,510])
%ylim([5e-7,1E-4])
set(gca,'Yscale','log')
xlabel('\sigma [Pa]')
ylabel('D [m^2/s]')
%text(485,1E-6,'(a)','FontSize',18)
box on
%colorbar('EastOutside')
%colorbar('Ticks',[6.8,6.9,7.0],...
%         'TickLabels',{'6.8','6.9','7.0'})
%caxis([350,450])
ax=gca;
ax.FontSize = 14;

%% %%%%%%%%%%%%%%%%%%%%%TALK FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make the first overview figure 


figure(1)

subplot 222
cmap = colormap(parula(length(unique(lidweightarr(masksp)))));

for i=1:4:length(masksp)
    
    %indexer = find(lidarruniq==lidarr(masksp(i)));
    indexer = find(unique(lidweightarr(masksp))==lidweightarr(masksp(i)));
    
    if lidweightarr(masksp(i)) < lintreshsp
       
        %scatter(timearr(masksp(i),:)*60,displarr(masksp(i),:)*0.01,'MarkerEdgeColor',cmap(indexer,:))
        scatter(timearr(masksp(i),:)*60,displarr(masksp(i),:)*0.01, 'k'); %,'MarkerEdgeColor',cmap(indexer,:))
        hold on
        plot(timearr(masksp(i),:)*60,displarr(masksp(i),:)*0.01, 'k')%,'Color',cmap(indexer,:))

        if ~isnan(errarr(masksp(i))) 
            %errorbar(timearr(masksp(i),:)*60,displarr(masksp(i),:)*0.01,errarr(masksp(i),:)*0.01,'Color',cmap(indexer,:))
            errorbar(timearr(masksp(i),:)*60,displarr(masksp(i),:)*0.01,errarr(masksp(i),:)*0.01,'k')%,'Color',cmap(indexer,:))
            fprintf('errorrrrrrrrrr %d\n', expy)
        end
        
        
    end
    
end

%plot(6*[10,1000],0.001*[10,1000].^0.6, '-.k')
%plot(6*[10,1000],0.001*[10,1000].^0.5, '-.k')

box on
xlabel('t [sec]')
ylabel('\delta [m]')
text(3250,.015,'(b)','FontSize',18)
%set(gca,'Yscale','log')
%set(gca,'Xscale','log')
xlim([0 4000])%xlim([20,6000])
ylim([0,0.1])%ylim([3E-3,.15])
ax=gca;
ax.FontSize = 14;
hold off



subplot 224

for i=length(masksp):-1:1
    
    %indexer = find(lidarruniq==lidarr(masksp(i)));
    indexer = find(unique(lidweightarr(masksp))==lidweightarr(masksp(i)));
    
    if lidweightarr(masksp(i)) < lintreshsp
       
        %plot(timearr(masksp(i),:)*60,displarr(masksp(i),:).^2/(10000*slopearr(masksp(i))*D2si),'-o','Color',cmap(indexer,:),'MarkerSize', 6)
        plot(timearr(masksp(i),:)*60,displarr(masksp(i),:).^2/(10000*slopearr(masksp(i))*D2si),'-ok')
        %convert D to proper units, but also the war time and displacement data
        hold on
        
    end
    
end
plot([0,4000],[0,4000],'-.k')
xlabel('t [sec]')
ylabel('\delta^2/D [sec]')
text(3250,500,'(c)','FontSize',18)
ax=gca;
ax.FontSize = 14;
hold off

% save
print('figtalk1', '-dpng','-r600')

%% FIG TALK 2

figure(2)

% first plot the D dynamics for both cases

subplot(2,2,[1 2])


plot([6.3,crtreshsp],[platysp,platysp], '-.k') % for spheres
hold on
plot([crtreshsp,crtreshsp],[1E-3,1E3], ':k')


for fixmass  = 0:0 % 1:60
 
    for i=1:length(masksp)

        if addedweightarr(masksp(i)) == fixmass,
            psh1 = scatter([lidweightarr(masksp(i))],[slopearr(masksp(i))*D2si],[30],'b','Fill')
            %plot(timearr(maskcyl30(i),:)*60,displarr(maskcyl30(i),:)*0.01,'Color',cmap(indexer,:))
            %convert D to proper units, but also the war time and displacement data
            hold on
            %plot(timearr(maskcyl30(i),:)*60,0.01*(timearr(maskcyl30(i),:)*slopearr(maskcyl30(i))).^0.5,'-.','Color',cmap(indexer,:))
        end

    end
    
    xtmp = crtreshsp:0.01:7.2;
    plot(xtmp,platysp*exp(-sloperhosp*(xtmp-crtreshsp)),'-k')

end
text(7.0,1E-5,'\sigma = 355 Pa','Color','blue','FontSize',14)

xlim([6.51,7.17])
ylim([1e-8,5E-5])
set(gca,'Yscale','log')
xlabel('\rho [g/L]')
ylabel('D [m^2/s]')
box on

subplot(2,2,[3 4])

plot([6.3,crtreshsp],[platysp,platysp], '-.k') % for spheres
hold on
plot([crtreshsp,crtreshsp],[1E-3,1E3], ':k')


for fixmass  = 0:0 % 1:60
 
    for i=1:length(masksp)

        if addedweightarr(masksp(i)) == fixmass,
            psh1 = scatter([lidweightarr(masksp(i))],[slopearr(masksp(i))*D2si],[30],'b','Fill')
            %plot(timearr(maskcyl30(i),:)*60,displarr(maskcyl30(i),:)*0.01,'Color',cmap(indexer,:))
            %convert D to proper units, but also the war time and displacement data
            hold on
            %plot(timearr(maskcyl30(i),:)*60,0.01*(timearr(maskcyl30(i),:)*slopearr(maskcyl30(i))).^0.5,'-.','Color',cmap(indexer,:))
        end

    end
    
    xtmp = crtreshsp:0.01:7.2;
    plot(xtmp,platysp*exp(-sloperhosp*(xtmp-crtreshsp)),'-k')

end

% also for the cylinder
plot([crtreshsp,crtreshsp],[1E-9,1E-3], ':k')


plot([6.3,crtreshsp],[platycyl,platycyl], '-.k') % for cyl

for fixmass  = 0:0
 
    for i=1:length(maskcyl30)

        if addedweightarr(maskcyl30(i)) == fixmass,
            pcyl1 = scatter([lidweightarr(maskcyl30(i))],[slopearr(maskcyl30(i))*D2si],[30],'rs', 'Fill')

            %plot(timearr(maskcyl30(i),:)*60,displarr(maskcyl30(i),:)*0.01,'Color',cmap(indexer,:))
            %convert D to proper units, but also the war time and displacement data
            hold on
            %plot(timearr(maskcyl30(i),:)*60,0.01*(timearr(maskcyl30(i),:)*slopearr(maskcyl30(i))).^0.5,'-.','Color',cmap(indexer,:))
        end

    end
    
    xtmp = crtreshsp:0.01:8.1;
    plot(xtmp,platycyl*exp(-sloperhocyl*(xtmp-crtreshsp)),'-k')


 

end
plot([crtreshsp,crtreshsp],[1E-9,1E-3], ':k')

xlim([6.51,7.17])
ylim([1e-8,5E-5])
set(gca,'Yscale','log')
xlabel('\rho [g/L]')
ylabel('D [m^2/s]')
box on
legend([psh1,pcyl1], 'sphere','cylinder','Location','southwest')


text(7.0,1E-5,'\sigma = 355 Pa','Color','blue','FontSize',14)
text(6.6,7E-7,'\sigma = 247 Pa','Color','red','FontSize',14)

% save
print('figtalk2', '-dpng','-r600')

%% FIG TALK 3

figure(3)

% first sphere data

bluey = subplot(2,2,1);

xtmp = 300:1:600;
plot(xtmp,0.65E-10*exp(xtmp./35),'-.k')
hold on

for i=1:length(masksp)

    if lidweightarr(masksp(i)) < lintreshsp, %== 6.96,
        scatter([spherestressarr(masksp(i))],[slopearr(masksp(i))*D2si],[30],[lidweightarr(masksp(i))],'Fill')

        hold on
    end

end

xlim([340,510])
ylim([5e-7,1E-4])
set(gca,'Yscale','log')
xlabel('\sigma [Pa]')
ylabel('D [m^2/s]')
box on
colorbar('EastOutside')
colorbar('Ticks',[6.8,6.9,7.0],...
         'TickLabels',{'6.8','6.9','7.0'})
%caxis([350,450])
ax=gca;
ax.FontSize = 14;

bluey = subplot(2,2,3);

plot([6.3,crtreshsp],[Dsigmsp,Dsigmsp], '-.k'); 
hold on
%plot([6.3,crtreshsp],[6.5E-5,6.5E-5], '-.k') 
plot([crtreshsp,crtreshsp],[1E-9,1E-3], ':k')
xtmp = crtreshsp:0.01:7.2;
%plot(xtmp,1.5E-5*exp(-12*(xtmp-crtreshsp)),'k')
plot(xtmp,Dsigmsp*exp(-sloperhosp*(xtmp-crtreshsp)),'-k')
%plot(xtmp,2E-4*exp(-12*(xtmp-crtreshsp)),'k')


for fixmass  = 0:1:60
 
    for i=1:length(masksp)

        if addedweightarr(masksp(i)) == fixmass,
            scatter([lidweightarr(masksp(i))],[slopearr(masksp(i))*D2si]/(platysp*exp(spherestressarr(masksp(i))/35)),[30],[spherestressarr(masksp(i))],'Fill')

            hold on
        end

    end

end
xlim([6.7,7.17])
ylim([1e-7,1E-4])
set(gca,'Yscale','log')
xlabel('\rho [g/L]')
ylabel('D/D_0 exp(\sigma/\sigma_s) [-]')
box on
colorbar('EastOutside')
caxis([350,450])
ax=gca;
ax.FontSize = 14;

MapRed = [linspace(0,1,256).^1',zeros(256,2)];
MapBlue = [zeros(256,2),linspace(0,1,256)'];
colormap(bluey, MapBlue);

% save
print('figtalk3', '-dpng','-r600')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)

% first plot the D dynamics for both cases

subplot(2,2,[1 2])


plot([6.3,crtreshsp],[platysp,platysp], '-.k') % for spheres
hold on
%plot([crtreshsp,crtreshsp],[1E-3,1E3], ':k')


for fixmass  = 0:0 % 1:60
 
    for i=1:length(masksv)

        if addedweightarr(masksv(i)) == fixmass,
            scatter([lidweightarr(masksv(i))],[slopearr(masksv(i))*D2si],[30],'g','Fill')
            %plot(timearr(maskcyl30(i),:)*60,displarr(maskcyl30(i),:)*0.01,'Color',cmap(indexer,:))
            %convert D to proper units, but also the war time and displacement data
            hold on
            %plot(timearr(maskcyl30(i),:)*60,0.01*(timearr(maskcyl30(i),:)*slopearr(maskcyl30(i))).^0.5,'-.','Color',cmap(indexer,:))
        end

    end
    
    xtmp = crtreshsp:0.01:7.2;
    plot(xtmp,platysp*exp(-sloperhosp*(xtmp-crtreshsp)),'-k')

end

% also for the cylinder
plot([crtreshsp,crtreshsp],[1E-9,1E-3], ':k')


plot([6.3,crtreshsp],[platycyl,platycyl], '-.k') % for cyl

for fixmass  = 0:0
 
    for i=1:length(maskcyl30)

        if addedweightarr(maskcyl30(i)) == fixmass,
            scatter([lidweightarr(maskcyl30(i))],[slopearr(maskcyl30(i))*D2si],[30],'rs', 'Fill')

            %plot(timearr(maskcyl30(i),:)*60,displarr(maskcyl30(i),:)*0.01,'Color',cmap(indexer,:))
            %convert D to proper units, but also the war time and displacement data
            hold on
            %plot(timearr(maskcyl30(i),:)*60,0.01*(timearr(maskcyl30(i),:)*slopearr(maskcyl30(i))).^0.5,'-.','Color',cmap(indexer,:))
        end

    end
    
    xtmp = crtreshsp:0.01:8.1;
    plot(xtmp,platycyl*exp(-sloperhocyl*(xtmp-crtreshsp)),'-k')


 

end
plot([crtreshsp,crtreshsp],[1E-9,1E-3], ':k')

xlim([6.31,7.17])
ylim([1e-9,5E-3])
set(gca,'Yscale','log')
xlabel('\rho [g/L]')
ylabel('D [m^2/s]')
box on

%text(7.0,1E-5,'\sigma = 355 Pa','Color','blue','FontSize',14)
%text(6.6,7E-7,'\sigma = 247 Pa','Color','red','FontSize',14)























%% D rescaling BOTH
figure(5);

cmap = colormap(parula(length(lidarruniq)));

for i=1:length(maskcyl30)
    
    if lidweightarr(maskcyl30(i)) >= crtreshsp && lidweightarr(maskcyl30(i)) < lintreshsp
    
        indexer = find(lidarruniq==lidweightarr(maskcyl30(i)));
        
        %scatter(addedweightarr(maskcyl30(i)),vertscale_manual(i).*slopearr(maskcyl30(i)),30,cmap(indexer,:), 'Fill','s')
        scatter(spherestressarr(maskcyl30(i)),...
                slopearr(maskcyl30(i))*D2si,70,cmap(indexer,:), 'Fill','s')
        hold on
    
    end

end

% now do the sphere data

for i=1:length(masksv)
    
    if lidweightarr(masksv(i)) >= crtreshsp && lidweightarr(masksv(i)) < lintreshsp
    
        indexer = find(lidarruniq==lidweightarr(masksv(i)));
        
        %scatter(addedweightarr(masksv(i)),vertscale_manual(i).*slopearr(masksp(i)),30,cmap(indexer,:), 'Fill','s')
        scatter(spherestressarr(masksp(i)),...
                slopearr(masksp(i))*D2si,50,cmap(indexer,:), 'Fill')
        hold on
    
    end

end
caxis([crtreshsp,lintreshsp])
colorbar('EastOutside')

set(gca,'Yscale','log')
xlabel('\sigma [Pa]')
ylabel('D [m^2/s]')
%ylim([.005,10])
%xlim([0,60])
box on

%% D rescaling CYL30
figure(12);

cmap = colormap(parula(length(lidarruniq)));

% low dens data

% subplot 221             
%              
% for i=1:length(maskcyl30)
%     
%     if lidarr(maskcyl30(i)) < crtreshsp
%     
%         indexer = find(lidarruniq==lidarr(maskcyl30(i)));
%         %scatter(addedweightarr(maskcyl30(i)),vertscale_manual(i).*slopearr(maskcyl30(i)),30,cmap(indexer,:), 'Fill','s')
%         scatter(stressarr(maskcyl30(i)),...
%                 slopearr(maskcyl30(i))*D2si,30,cmap(indexer,:), 'Fill','s')
%         hold on
%     
%     end
% 
% end
% 
% set(gca,'Yscale','log')
% %set(gca,'Xscale','log')
% title('<6.82')
% xlabel('\sigma [Pa]')
% ylabel('D [m^2/s]')
% %ylim([.5,10])
% %xlim([0,5])
% box on

subplot 231

for i=1:length(maskcyl30)
    
    if lidweightarr(maskcyl30(i)) >= crtreshsp % && lidarr(maskcyl30(i)) < lintreshsp
    
        indexer = find(lidarruniq==lidweightarr(maskcyl30(i)));
        
        %scatter(addedweightarr(maskcyl30(i)),vertscale_manual(i).*slopearr(maskcyl30(i)),30,cmap(indexer,:), 'Fill','s')
        scatter(spherestressarr(maskcyl30(i)),...
                slopearr(maskcyl30(i))*D2si,30,cmap(indexer,:), 'Fill','s')
        hold on
    
    end

end


% colorbar('EastOutside',...
%     'Ticks',[0,0.25,0.5,0.75,1],...
%     'TickLabels',{num2str(mindens,2),...
%                     num2str(mindens+0.25*densrange,2),...
%                     num2str(mindens+0.5*densrange,2),...
%                     num2str(mindens+0.75*densrange,2),...
%                     num2str(maxdens,2)})


set(gca,'Yscale','log')
%set(gca,'Xscale','log')
title('Cylinder 30')
xlabel('\sigma [Pa]')
ylabel('D [m^2/s]')
%ylim([.005,10])
%xlim([0,60])
box on


%rescale cylinder data
subplot 232  

vertscale_manual = [ones(4,1)*1; ... %6.82
                 ones(8,1)*4; ... % 6.87
                 ones(4,1)*16; ... % 6.92
                 ones(8,1)*50; ... % 6.97
                 ones(5,1)*80; ... % 7.02
                 ones(4,1)*400; ... % 7.09
                 ones(4,1)*8000; ... % 7.04
                 ones(4,1)*150000; ... % 7.80
                 ones(9,1)*300000]; % 8.27

prefacarr(10) = 1;
prefacarr(12) = 4;
prefacarr(14) = 16;
prefacarr(16) = 50;
prefacarr(17) = 80;
prefacarr(20) = 400;
%prefacarr(x) = 8000;
%prefacarr(x) = 150000;
%prefacarr(x) = 300000;
             
for i=1:length(maskcyl30)
    
    if lidweightarr(maskcyl30(i)) >= crtreshsp && lidweightarr(maskcyl30(i)) < lintreshsp
    
        indexer = find(lidarruniq==lidweightarr(maskcyl30(i)));
             
        %prefact = 1E5*(lidarr(maskcyl30(i))-crtreshsp)^0.5;
        %scatter(addedweightarr(maskcyl30(i)),vertscale_manual(i).*slopearr(maskcyl30(i)),30,cmap(indexer,:), 'Fill','s')
        scatter(spherestressarr(maskcyl30(i)),...
                prefacarr(indexer)*slopearr(maskcyl30(i))*D2si,30,cmap(indexer,:), 'Fill','s')
        hold on
    
    end

end


set(gca,'Yscale','log')
%set(gca,'Xscale','log')
xlabel('\sigma [Pa]')
ylabel('a(\rho)D [m^2/s]')
%ylim([.005,10])
%xlim([0,60])
box on

subplot 233

for i=1:length(maskcyl30)
    
    if lidweightarr(maskcyl30(i)) >= crtreshsp && lidweightarr(maskcyl30(i)) < lintreshsp
    
        indexer = find(lidarruniq==lidweightarr(maskcyl30(i)));
             
        scatter(lidweightarr(maskcyl30(i)),prefacarr(indexer),30,cmap(indexer,:),'Fill','s')
        
        hold on
    
    end

end
caxis([crtreshsp,lintreshsp])
colorbar('EastOutside')

set(gca,'Yscale','log')
%set(gca,'Xscale','log')
xlabel('\rho [gr/L]')
ylabel('a(\rho)')
%ylim([.005,10])
%xlim([0,60])
box on

% now do the sphere data

subplot 234

for i=1:length(masksp)
    
    if lidweightarr(masksp(i)) >= crtreshsp && lidweightarr(masksp(i)) < lintreshsp
    
        indexer = find(lidarruniq==lidweightarr(masksp(i)));
        
        %scatter(addedweightarr(masksp(i)),vertscale_manual(i).*slopearr(masksp(i)),30,cmap(indexer,:), 'Fill','s')
        scatter(spherestressarr(masksp(i)),...
                slopearr(masksp(i))*D2si,30,cmap(indexer,:), 'Fill')
        hold on
    
    end

end


set(gca,'Yscale','log')
title('Sphere')
xlabel('\sigma [Pa]')
ylabel('D [m^2/s]')
%ylim([.005,10])
%xlim([0,60])
box on


subplot 235  

prefacarr(11) = 1;
prefacarr(13) = 2;
prefacarr(15) = 4;
prefacarr(17) = 7;
prefacarr(18) = 9;
prefacarr(19) = 14;
prefacarr(20) = 18;


             
for i=1:length(masksp)
    
    if lidweightarr(masksp(i)) >= crtreshsp && lidweightarr(masksp(i)) < lintreshsp
    
        indexer = find(lidarruniq==lidweightarr(masksp(i)));

             
        %prefact = 1E5*(lidarr(masksp(i))-crtreshsp)^0.5;
        %scatter(addedweightarr(masksp(i)),vertscale_manual(i).*slopearr(masksp(i)),30,cmap(indexer,:), 'Fill','s')
        scatter(spherestressarr(masksp(i)),...
                prefacarr(indexer)*slopearr(masksp(i))*D2si,30,cmap(indexer,:), 'Fill')
        hold on
    
    end

end


set(gca,'Yscale','log')
%set(gca,'Xscale','log')
xlabel('\sigma [Pa]')
ylabel('a(\rho)D [m^2/s]')
%ylim([.005,10])
%xlim([0,60])
box on

subplot 236

for i=1:length(masksp)
    
    if lidweightarr(masksp(i)) >= crtreshsp && lidweightarr(masksp(i)) < lintreshsp
    
        indexer = find(lidarruniq==lidweightarr(masksp(i)));
             
        scatter(lidweightarr(masksp(i)),prefacarr(indexer),30,cmap(indexer,:),'Fill')
        
        hold on
    
    end

end

caxis([crtreshsp,lintreshsp])
colorbar('EastOutside')

set(gca,'Yscale','log')
%set(gca,'Xscale','log')
xlabel('\rho [gr/L]')
ylabel('a(\rho)')
%ylim([.005,10])
%xlim([0,60])
box on

% save
print('fig3', '-depsc','-r600')
print('fig3', '-dpdf')


%% OLD FIG 1


figure(2)

% plot examples from entire database
[sortie,sortindex] = sort(slopearr);

subplot 221
cmap = colormap(parula(length(unique(lidweightarr))));

for expy = 1:14:length(sortie)
    
    indexer = find(unique(lidweightarr)==lidweightarr(expy));
    indexer = indexer(1);

    plot(timearr(expy,:)*60,displarr(expy,:)*0.01,'Color',cmap(indexer,:))
    
    if ~isnan(errarr(expy))
        errorbar(timearr(expy,:)*60,displarr(expy,:)*0.01,errarr(expy,:)*0.01,'Color',cmap(indexer,:))
        fprintf('errorrrrrrrrrr %d\n', expy)
    end
    %convert D to proper units, but also the war time and displacement data
    hold on
    plot(timearr(expy,:)*60,0.01*(timearr(expy,:)*slopearr(expy)).^0.5,'-.','Color',cmap(indexer,:))

end

xlabel('t [sec]')
ylabel('\delta [m]')
title('Examples')
hold off
box on

%plot colorbar

[sortie,sortindex] = sort(slopearr);

subplot 223
cmap = colormap(parula(length(unique(lidweightarr))));

for expy = 2:14:length(sortie)
    
    indexer = find(unique(lidweightarr)==lidweightarr(expy));
    indexer = indexer(1);

 
    plot(timearr(expy,:)*60,displarr(expy,:)*0.01,'Color',cmap(indexer,:))
    
    if ~isnan(errarr(expy))
        errorbar(timearr(expy,:)*60,displarr(expy,:)*0.01,errarr(expy,:)*0.01,'Color',cmap(indexer,:))
        fprintf('errorrrrrrrrrr %d\n', expy)
    end
    %convert D to proper units, but also the war time and displacement data
    hold on
    plot(timearr(expy,:)*60,0.01*(timearr(expy,:)*slopearr(expy)).^0.5,'-.','Color',cmap(indexer,:))


end

xlabel('t [sec]')
ylabel('\delta [m]')
hold off
box on
%set(gca,'Yscale','log')
%set(gca,'Xscale','log')


colorbar('EastOutside',...
    'Ticks',[0,0.25,0.5,0.75,1],...
    'TickLabels',{num2str(mindens,2),...
                    num2str(mindens+0.25*densrange,2),...
                    num2str(mindens+0.5*densrange,2),...
                    num2str(mindens+0.75*densrange,2),...
                    num2str(maxdens,2)})

                
% rescaling of delta^2/D

subplot 222
cmap = colormap(parula(length(unique(lidweightarr(masksp)))));

for i=1:length(masksp)
    
    %indexer = find(lidarruniq==lidarr(masksp(i)));
    indexer = find(unique(lidweightarr(masksp))==lidweightarr(masksp(i)));
       
    plot(timearr(masksp(i),:)*60,displarr(masksp(i),:).^2/(10000*slopearr(masksp(i))*D2si),'Color',cmap(indexer,:))
    %convert D to proper units, but also the war time and displacement data
    hold on
    
end
plot([0,6000],[0,6000],'-.k')
xlabel('t [sec]')
ylabel('\delta^2/D [sec]')
hold off
title('23mm sphere')
                
                
% then plot the D dynamics

subplot 224

platy = 1.6E-5;
plot([6.3,6.82],[platy,platy], '-.k') 
hold on
plot([6.82,6.82],[1E-3,1E3], ':k')


for fixmass  = 0:0 % 1:60
 
    for i=1:length(masksp)

        if addedweightarr(masksp(i)) == fixmass,
            scatter([lidweightarr(masksp(i))],[slopearr(masksp(i))*D2si],[30],[addedweightarr(masksp(i))],'Fill')

            %plot(timearr(maskcyl30(i),:)*60,displarr(maskcyl30(i),:)*0.01,'Color',cmap(indexer,:))
            %convert D to proper units, but also the war time and displacement data
            hold on
            %plot(timearr(maskcyl30(i),:)*60,0.01*(timearr(maskcyl30(i),:)*slopearr(maskcyl30(i))).^0.5,'-.','Color',cmap(indexer,:))
        end

    end
    
    xtmp = 6.82:0.01:7.2;
    plot(xtmp,platy*exp(-14*(xtmp-6.82)),'k')

    xlim([6.7,7.2])
    ylim([1e-8,1E-4])
    set(gca,'Yscale','log')
    xlabel('\rho [gr/L]')
    ylabel('D [m^2/s]')
    box on
    %colorbar('EastOutside')


end

% save
%print('example-sp-sq', '-depsc','-r600')


%%

% figure D for sphere and cyl 35

%figure(2)


subplot 223


%plot([6.82,6.82],[1E-3,1E3], '-.k')

for fixmass  = 0:1:60
 
    for i=1:length(maskcyl35)

        if addedweightarr(maskcyl35(i)) == fixmass,
            scatter([lidweightarr(maskcyl35(i))],[slopearr(maskcyl35(i))],[30],[addedweightarr(maskcyl35(i))],'Fill')

            hold on
        end

    end

    xlim([6.4,6.9])
    ylim([1e-3,100])
    set(gca,'Yscale','log')
    xlabel('\rho [gr/L]')
    ylabel('D [m^2/s]')
    title('Cyl 35mm')
    box on
    colorbar('EastOutside')


end

subplot 224


%plot([6.82,6.82],[1E-3,1E3], '-.k')

for fixmass  = 0:1:60
 
    for i=1:length(masksp)

        if addedweightarr(masksp(i)) == fixmass,
            scatter([lidweightarr(masksp(i))],[slopearr(masksp(i))],[30],[addedweightarr(masksp(i))],'Fill')

            hold on
        end

    end

    xlim([6.7,7.2])
    ylim([1e-1,100])
    set(gca,'Yscale','log')
    xlabel('\rho [gr/L]')
    ylabel('D [m^2/s]')
    title('Sphere 23mm')
    box on
    colorbar('EastOutside')


end

% save
print('example-sq-35-sp-rescale', '-dpng','-r600')



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% plotting fits exponents

figure(1)

subplot 131
imagesc(bexpy:expyr:eexpy,maskcyl30,distopyarr(maskcyl30,:))
colorbar('EastOutside')
xlabel('exponent')
ylabel('measurement #')
title('cylinder 30')

subplot 132
imagesc(bexpy:expyr:eexpy,maskcyl35,distopyarr(maskcyl35,:))
colorbar('EastOutside')
xlabel('exponent')
ylabel('measurement #')
title('cylinder 35')


subplot 133
imagesc(bexpy:expyr:eexpy,masksp,distopyarr(masksp,:))
colorbar('EastOutside')
xlabel('exponent')
ylabel('measurement #')
title('sphere 23')

%% plotting fits exponents

figure(1)

subplot 131
imagesc(bexpy:expyr:eexpy,maskcyl30,expyarr(maskcyl30,:))
colorbar('EastOutside')
xlabel('exponent')
ylabel('measurement #')
title('cylinder 30')

subplot 132
imagesc(bexpy:expyr:eexpy,maskcyl35,expyarr(maskcyl35,:))
colorbar('EastOutside')
xlabel('exponent')
ylabel('measurement #')
title('cylinder 35')


subplot 133
imagesc(bexpy:expyr:eexpy,masksp,expyarr(masksp,:))
colorbar('EastOutside')
xlabel('exponent')
ylabel('measurement #')
title('sphere 23')


hold off
% save
print('exponent_qual_all', '-dpng','-r600')

%% plot fits exponent as a function of rho and mass

figure(2)

subplot 131
%[X,Y] = meshgrid(lidarr(maskcyl30),addedweightarr(maskcyl30));
scatter(lidweightarr(maskcyl30),addedweightarr(maskcyl30),30,optexparr(maskcyl30), 'Fill')
xlabel('\rho [gr/L]')
ylabel('added mass [gr]')
title('Cyl 30mm')
box on
colorbar('Southoutside')
caxis([bexpy,0.8])

subplot 132
%[X,Y] = meshgrid(lidarr(maskcyl30),addedweightarr(maskcyl30));
scatter(lidweightarr(maskcyl35),addedweightarr(maskcyl35),30,optexparr(maskcyl35), 'Fill')
xlabel('\rho [gr/L]')
ylabel('added mass [gr]')
title('Cyl 35mm')
box on
colorbar('Southoutside')
caxis([bexpy,eexpy])


subplot 133
%[X,Y] = meshgrid(lidarr(maskcyl30),addedweightarr(maskcyl30));
scatter(lidweightarr(masksp),addedweightarr(masksp),30,optexparr(masksp), 'Fill')
xlabel('\rho [gr/L]')
ylabel('added mass [gr]')
title('Sphere 23mm')
box on
colorbar('Southoutside')
caxis([bexpy,eexpy])

hold off
% save
print('exponent_qual_all_perprobe', '-dpng','-r600')

%% plot fits exponent as a function of rho FOR GIVEN MASS

figure(24)

subplot 131

maskX = find(addedweightarr(maskcyl30) == 0);
%[X,Y] = meshgrid(lidarr(maskcyl30),addedweightarr(maskcyl30));
scatter(lidweightarr(maskcyl30(maskX)),addedweightarr(maskcyl30(maskX)),30,bestexparr(maskcyl30(maskX)), 'Fill')
xlabel('\rho [gr/L]')
ylabel('added mass [gr]')
title('Cyl 30mm')
box on
colorbar('Southoutside')
caxis([bexpy,eexpy])

subplot 132
maskX = find(addedweightarr(maskcyl35) == 0);
%[X,Y] = meshgrid(lidarr(maskcyl30),addedweightarr(maskcyl30));
scatter(lidweightarr(maskcyl35(maskX)),addedweightarr(maskcyl35(maskX)),30,bestexparr(maskcyl35(maskX)), 'Fill')
xlabel('\rho [gr/L]')
ylabel('added mass [gr]')
title('Cyl 35mm')
box on
colorbar('Southoutside')
caxis([bexpy,eexpy])


subplot 133
maskX = find(addedweightarr(masksp) == 0);
%[X,Y] = meshgrid(lidarr(maskcyl30),addedweightarr(maskcyl30));
scatter(lidweightarr(masksp(maskX)),addedweightarr(masksp(maskX)),30,bestexparr(masksp(maskX)), 'Fill')
xlabel('\rho [gr/L]')
ylabel('added mass [gr]')
title('Sphere 23mm')
box on
colorbar('Southoutside')
caxis([bexpy,eexpy])

hold off
% save
print('exponent_qual_all_perprobe_fM', '-dpng','-r600')

%% 

figure(24)

maskX = find(addedweightarr(maskcyl30) == 0);
%[X,Y] = meshgrid(lidarr(maskcyl30),addedweightarr(maskcyl30));
scatter(lidweightarr(maskcyl30(maskX)),slopearr(maskcyl30(maskX)),30,bestexparr(maskcyl30(maskX)), 'Fill')
set(gca,'Yscale','log')
xlabel('\rho [gr/L]')
ylabel('D')
title('Cyl 30mm')
box on


%% histogram of all exponents

figure(34);
cdfplot(bestexparr);
xlabel('exponent')


% save
print('exponent_all_cdf', '-dpng','-r600')

%% collapse all displacement data

% plot collapse of square root behavior
% call displacement \delta

figure(20);

subplot 131
cmap = colormap(parula(length(lidarruniq)));
for i=1:length(maskcyl30)
    
    indexer = find(lidarruniq==lidweightarr(maskcyl30(i)));
       
    plot(timearr(maskcyl30(i),:)*60,displarr(maskcyl30(i),:).^2/(10000*slopearr(maskcyl30(i))*D2si),'Color',cmap(indexer,:))
    %convert D to proper units, but also the war time and displacement data
    hold on
    
end
plot([0,6000],[0,6000],'-.k')
xlabel('t [sec]')
ylabel('\delta^2/D [sec]')
hold off
title('30mm cyl')

subplot 132
cmap = colormap(parula(length(lidarruniq)));
for i=1:length(maskcyl35)
    
    indexer = find(lidarruniq==lidweightarr(maskcyl35(i)));
       
    plot(timearr(maskcyl35(i),:)*60,displarr(maskcyl35(i),:).^2/(10000*slopearr(maskcyl35(i))*D2si),'Color',cmap(indexer,:))
    %convert D to proper units, but also the war time and displacement data
    hold on
    
end
plot([0,6000],[0,6000],'-.k')
xlabel('t [sec]')
ylabel('\delta^2/D [sec]')
hold off
title('35mm cyl')

subplot 133
cmap = colormap(parula(length(lidarruniq)));
for i=1:length(masksp)
    
    indexer = find(lidarruniq==lidweightarr(masksp(i)));
       
    plot(timearr(masksp(i),:)*60,displarr(masksp(i),:).^2/(10000*slopearr(masksp(i))*D2si),'Color',cmap(indexer,:))
    %convert D to proper units, but also the war time and displacement data
    hold on
    
end
plot([0,6000],[0,6000],'-.k')
xlabel('t [sec]')
ylabel('\delta^2/D [sec]')
hold off
title('23mm sphere')

% save
print('allcollapse', '-dpng','-r600')

%% EXAMPLE collapse all displacement data FIXED MASS

% plot collapse of square root behavior
% call displacement \delta

fixmass = 0;

figure(200);

subplot 121
cmap = colormap(parula(length(lidarruniq)));
for i=1:length(maskcyl30)
    
    indexer = find(lidarruniq==lidweightarr(maskcyl30(i)));
       
    if addedweightarr(maskcyl30(i)) == fixmass,
        plot(timearr(maskcyl30(i),:)*60,displarr(maskcyl30(i),:)*0.01,'Color',cmap(indexer,:))
        %convert D to proper units, but also the war time and displacement data
        hold on
        plot(timearr(maskcyl30(i),:)*60,0.01*(timearr(maskcyl30(i),:)*slopearr(maskcyl30(i))).^0.5,'-.','Color',cmap(indexer,:))
    end
    
end
%plot([0,6000],[0,6000],'-.k')
xlabel('t [sec]')
ylabel('\delta [m]')
hold off

subplot 122
cmap = colormap(parula(length(lidarruniq)));
for i=1:length(maskcyl30)
    
    indexer = find(lidarruniq==lidweightarr(maskcyl30(i)));
       
    if addedweightarr(maskcyl30(i)) == fixmass,
        plot(timearr(maskcyl30(i),:)*60,displarr(maskcyl30(i),:).^2/(10000*slopearr(maskcyl30(i))*D2si),'Color',cmap(indexer,:))
        %convert D to proper units, but also the war time and displacement data
        hold on
    end
    
end
plot([0,6000],[0,6000],'-.k')
xlabel('t [sec]')
ylabel('\delta^2/D [sec]')
hold off

% save
print('allcollapse_example', '-dpng','-r600')


%% make plots of all displacement data per density
% get around collision of same dens for 30 / 35 by multiplying

[dataID,ia,ic] = unique(lidweightarr.*sizearr);

for j=1:length(dataID)
   
    figure
    masky = find(ic == j);
    cmap = colormap(parula(max(addedweightarr(masky))+1));
    for i=1:length(masky)
    
        indexer = find(lidarruniq==lidweightarr(masky(i)));      
        plot(timearr(masky(i),:)*60,displarr(masky(i),:).^2/(10000*slopearr(masky(i))*D2si),'Color',cmap(1+addedweightarr(masky(i)),:))
        %convert D to proper units, but also the war time and displacement data
        hold on
    
    end
    caxis([0, max(addedweightarr(masky))+1])
    c = colorbar('EastOutside');
    c.Label.String = 'Added weight';
    plot([0,6000],[0,6000],'-.k')
    xlabel('t [sec]')
    ylabel('\delta^2/D [sec]')
    titlestr = horzcat('rho-',num2str(lidweightarr(masky(1))),'-size-',num2str(sizearr(masky(1))),'.png');
    title(titlestr)
    print(titlestr, '-dpng','-r600')
end

%% first plot different D
hFig = figure(1);

cmap = colormap(parula(length(lidarruniq)));


for i=1:numfiles
    
    indexer = find(lidarruniq==lidweightarr(i));
    
    if sizearr(i) == 23,
        scatter(addedweightarr(i),slopearr(i),30,indexer,'Fill')
    elseif sizearr(i) == 30,
        scatter(addedweightarr(i),slopearr(i),30,indexer,'Fill','s')
    else
        scatter(addedweightarr(i),slopearr(i),30,indexer,'Fill','^')
    end
    
    hold on
    
end


hold on

set(gca,'Yscale','log')

xlabel('Added weight [grams]')
ylabel('D [cm^2/min]')
%xlim([25,100])

box on

axis('square')
set(gcf,'PaperPositionMode','manual')
set(gcf,'PaperPosition',[0 0 10 10])


hold off
% save
print('D-overview', '-dpng','-r600')



%% then rescaling of vertical axis

figure(5);
cmap = colormap(parula(length(lidarruniq)));


vertscale_manual = [ones(4,1)*1; ... %6.82
                 ones(8,1)*4; ... % 6.87
                 ones(4,1)*16; ... % 6.92
                 ones(8,1)*50; ... % 6.97
                 ones(5,1)*80; ... % 7.02
                 ones(4,1)*400; ... % 7.09
                 ones(4,1)*8000; ... % 7.04
                 ones(4,1)*150000; ... % 7.80
                 ones(9,1)*300000];

subplot 131             
cmap = colormap(parula(length(lidarruniq)));
             
for i=1:length(maskcyl30)
    
    indexer = find(lidarruniq==lidweightarr(maskcyl30(i)));
    scatter(addedweightarr(maskcyl30(i)),vertscale_manual(i).*slopearr(maskcyl30(i)),30,cmap(indexer,:), 'Fill','s')
    hold on

end


set(gca,'Yscale','log')
%set(gca,'Xscale','log')

xlabel('Added weight [grams]')
ylabel('a(\rho)*D [cm^2/min]')
ylim([.5,5E6])
xlim([0,60])
box on

% then rescaling of vertical axis for cyl 35

vertscale_manual_35 = [ones(5,1)*1; ... %6.58
                 ones(7,1)*4; ... % 6.63
                 ones(10,1)*10; ... % 6.70
                 ];

subplot 132
cmap = colormap(parula(length(lidarruniq)));

for i=1:length(maskcyl35)
    
    indexer = find(lidarruniq==lidweightarr(maskcyl35(i)));
    scatter(addedweightarr(maskcyl35(i)),vertscale_manual_35(i).*slopearr(maskcyl35(i)),30,cmap(indexer,:), 'Fill','^')
    hold on

end


set(gca,'Yscale','log')
%set(gca,'Xscale','log')

xlabel('Added weight [grams]')
ylabel('a(\rho)*D [cm^2/min]')
ylim([.5,5E6])
xlim([0,60])
box on

% then rescaling of vertical axis for sphere 23

vertscale_manual_23 = [ones(5,1)*1; ... %6.58
                 ones(5,1)*3; ... % 6.63
                 ];

subplot 133
cmap = colormap(parula(length(lidarruniq)));

for i=1:length(masksp)
    
    indexer = find(lidarruniq==lidweightarr(masksp(i)));
    scatter(addedweightarr(masksp(i)),vertscale_manual_23(i).*slopearr(masksp(i)),30,cmap(indexer,:), 'Fill')
    hold on

end


set(gca,'Yscale','log')
%set(gca,'Xscale','log')

xlabel('Added weight [grams]')
ylabel('a(\rho)*D [cm^2/min]')
ylim([.5,5E6])
xlim([0,60])
box on
print('D-rescale-all', '-dpng','-r600')

%% plot shift factors

figure(6)

cmap = colormap(parula(length(lidweightarr)));

scatter(lidweightarr(maskcyl30),vertscale_manual,30,lidweightarr(maskcyl30),'Fill','s')
hold on
scatter(lidweightarr(maskcyl35),vertscale_manual_35,30,lidweightarr(maskcyl35), 'Fill','^')
scatter(lidweightarr(masksp),vertscale_manual_23,30,lidweightarr(masksp), 'Fill')


set(gca,'Yscale','log')
%set(gca,'Xscale','log')

xlabel('\rho [grams/L]')
ylabel('a(\rho)')
%xlim([25,100])

box on

axis('square')
set(gcf,'PaperPositionMode','manual')
set(gcf,'PaperPosition',[0 0 10 10])


hold off
% save
print('rescale-fact', '-dpng','-r600')
