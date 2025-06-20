
%% some details regarding conversions etc 

Eta2si = 0.01/60; % convert linear slopes fitted to SI units (measurements are in cm/min)

ballrad = 0.02; % ball radius in m

% weight of displaced volume: rho * 4/3pir^3
% assume water density = 1000 kg/m3. Not technically correct due to
% presence of hydrogels, but good enough

wdisplvol = 1000*0.333*pi()*4*(ballrad).^3;

% some other parameters
ppweight = 0.00270; % in kg, weight of pp ball
denspp = ppweight/(0.333*pi()*4*(ballrad).^3); % dens of pp ball
deltarho = 1000-denspp;

%% read data

fnames = dir('data_tempdep\*.dat');
numfiles = length(fnames);

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

bestexparr = NaN(numfiles,1); 
optexparr = NaN(numfiles,1);

temparr = NaN(numfiles,1);

for i=1:numfiles
    
    if ~mod(i,2), fprintf('numfile # %d\n', i); end
    
    % units all in centimeters and minutes!
    readfile = char(fnames(i).name);

    readfiletotal = char(horzcat(fnames(i).folder,'\',fnames(i).name));
    
  
    tmp = importdata(readfiletotal);
    tmpdata = tmp.data;
    tmpstring = tmp.textdata;
    tmpstring = tmpstring{1};
    datl = size(tmpdata);

    temparr(i) = str2num(readfile(7:8));
    
   
    lengarr(i) = datl(1);
    timearr(i,1:datl) = tmpdata(:,1);
    displarr(i,1:datl) = tmpdata(:,2);

    if contains(tmpstring, 'err') % if there is an error estimate, it will be in the third column.
        errarr(i,1:datl(1)) = tmpdata(:,3);
    end
    
    %%%%%%%%%%%% we first assume square root behavior:

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
        opty = fitoptions('power2','Lower', [0,expy,-10], 'Upper', [Inf,expy,10]);
       [fitobj,gof2] = fit(xtmp',deltatmp','power2',opty);
       expyarr(i,ctr) = gof2.rsquare;
       resi = fitobj(xtmp)-deltatmp';
       stdarr(i,ctr) = nanstd(resi);
       kurtarr(i,ctr) = kurtosis(resi);
       skewarr(i,ctr) = skewness(resi);  
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

%% read the temperature water data

% contains data from an engineering toolbox website, and an Anton Paar
% provided reference, both of which are in accordance with the reference
% from Journal of Physical and Chemical Reference Data

T = readtable('data_tempdep\watervisc.xlsx'); 


%% FIGURE WITH TEMP DEPENDENT VISCOSITY DATA

% buoyancy force comes from the weight of displaced volume minus the weight
% from intruder

weightarr = NaN(1,numfiles);

indy = [];
for i=1:numfiles
    
    if ~mod(i,2), fprintf('numfile # %d\n', i); end


    %walk through file list in tmpnames
    for j=1:numfiles+1

            indy = [indy j]; 
            weightarr(i) = ppweight; 
            surfy = pi()*(ballrad)^2 % - all pingpong ball
            buoyarr(i) = 9.81*(wdisplvol-weightarr(i))/surfy; %buoyancy stress in Pa (about 240 Pa for PP)

        
    end

end


%% plot stress dependence of buoyancy speed

figure(5)

% slope to visc conversion with equation from main manuscript.

% eta = (2/9)*deltarho*r^2*g * (1/v)

etaconv = (2/9)*deltarho*ballrad^2*9.81;

normy_hydro = etaconv./(Eta2si*slopearr); 

% slopearr is now in cm/min; conversion via Eta2si
% normalize on the viscosity value at 20 degrees, which turns out to be 0.8244E5
plot(temparr,etaconv./(Eta2si*slopearr)/normy_hydro(7), 'o',...    
    'LineWidth',1,...
    'MarkerSize',7,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0,0,1])

hold on

% now plot water viscosity data. Automatically normalized at value 
% for 20 degrees, as the values were in mPas.

celly = T.('C_1');
velly = T.('dyn_1');

plot(celly,velly/velly(16),'>',...    
    'LineWidth',1,...
    'MarkerSize',7,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0,0,0])

%manually insert Arrhenius trends with manually determined coefficients
xtmp = 15:0.1:30;
plot(xtmp,4.6*exp(-xtmp./11),'-b')
plot(celly,1.65*exp(-celly./41)/velly(16),'-k')

box on
xlabel('Temperature [^oC]')
ylabel('\eta/\eta_{17} [-]')
xlim([15 30])
ylim([0.3,1.1])
ax=gca;
ax.FontSize = 14;
hold on
xtmp = 16:.1:29;

hold off

% save
print('risingspeed-temp', '-dpdf', '-bestfit')

