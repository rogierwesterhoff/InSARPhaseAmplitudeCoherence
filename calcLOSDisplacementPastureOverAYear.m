%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
cc

baresoil=false;
xlsfilename = 'i:\GroundWater\Research\Smart\Satellite RS\Papers_Reports_Presentations\Journal_Papers\GRL_InSAR_phase_disturbances\YearInAPaddock.xlsx';
% xlsfilename = 'i:\GroundWater\Research\Smart\Satellite RS\Papers_Reports_Presentations\Journal_Papers\GRL_InSAR_phase_disturbances\YearOnALevee.xlsx';
xlstabname = 'Sheet1';
% freq=5.405e9; % C-band
freq=10e9; % X-band
lambda = 3e8/freq; % wavelength (needed later)
numscenarios = 50;
% imode_inc=1; % incident wave
imode_refl=2; % reflection
imode_trans=3; % transmitted
Ei0=1;
% eps_soil=3; % the reference soil. A dry sand
% eps_soil=24-1i*5;% the references soil: silty clay Hallikainen
% eps_soil=75-1i*30;% fresh-brackish water-saturated
% eps_soil=5.55-1i*1.25;% vegetation
% eps_soil=2-1i*1;% light vegetation

[num,~,raw]=xlsread(xlsfilename,xlstabname);
dateStrings = datetime(raw(2:size(num,1)+1,1),'ConvertFrom','excel');
dateNums = datenum(dateStrings);
dateVecs = datevec(dateNums);

% THIS XTICK IS DONE MANUALLY
xtick = ones(15,1);
for imonth=1:15
    xtick(imonth) = datenum(2012,imonth+11,1);
end
xticklabel={'D','J','F','M','A','M','J','J','A','S','O','N','D','J','F'};

soilMoisture = num(:,1);

percLOSDeformationE = zeros(size(soilMoisture,1),numscenarios); % first one should be zero
percLOSDeformationH = zeros(size(soilMoisture,1),numscenarios); % first one should be zero
Hcoh = nan(size(soilMoisture));
Ecoh = nan(size(soilMoisture));
angleE02 = nan(size(soilMoisture));
angleH02 = nan(size(soilMoisture));

reps=num(1,2); ieps = num(1,3); % tmp
eps_soil=reps(1)-1i*ieps(1);
[E1,H1,E01,H01] = em_propagation_2lyrs(1,eps_soil,Ei0,false,freq);
    
for iscenario = 1:numscenarios
    if iscenario == 1
        noisefactor = 0;
    else
        noisefactor = 0.05;
    end
    
    reps = num(:,2)+noisefactor*randn(size(num,1),1).*num(:,2);
    ieps = num(:,3)+noisefactor*randn(size(num,1),1).*num(:,3);
   
    if baresoil
        for i=1:length(reps)
            
            if i==1
                Eu1=ones(10,1)*E01(imode_refl);
                Hu1=ones(10,1)*H01(imode_refl);
            else
                Eu1=ones(10,1)*E02a(imode_refl);
                Hu1=ones(10,1)*H02a(imode_refl);
            end
            
            eps_soil=reps(i)-1i*ieps(i);
            [E2a,H2a,E02a,H02a] = em_propagation_2lyrs(1,eps_soil,1,false,freq);
            
            angleE02(i)=round(180*angle(E02a(imode_refl))/pi);
            angleH02(i)=round(180*angle(H02a(imode_refl))/pi);
            
            Eu2=ones(10,1)*E02a(imode_refl);
            Hu2=ones(10,1)*H02a(imode_refl);
            
            Hcomplconjprod=nan(size(Hu1));
            Ecomplconjprod=nan(size(Eu1));
            
            for iCoherence=1:length(Hu1)
                Hcomplconjprod(iCoherence)=Hu1(iCoherence)*conj(Hu2(iCoherence));
                Ecomplconjprod(iCoherence)=Eu1(iCoherence)*conj(Eu2(iCoherence));
            end
            
            Hcoh(i)=sum(Hcomplconjprod)/(sqrt(sum((abs(Hu1)).^2))*sqrt(sum((abs(Hu2)).^2)));
            Ecoh(i)=sum(Ecomplconjprod)/(sqrt(sum((abs(Eu1)).^2))*sqrt(sum((abs(Eu2)).^2)));
            
        end
    else
        % vegetation properties
        grassHeight = num(:,4);
        repsveg = num(:,6)+noisefactor*randn(size(num,1),1).*num(:,6);
        iepsveg = num(:,7)+noisefactor*randn(size(num,1),1).*num(:,7);        
          
        
        for i=1:length(reps)
            
            if i~=1
                Eu1=ones(10,1)*E02c(imode_trans);
                Hu1=ones(10,1)*H02c(imode_trans);
            end
            
            numWaves = grassHeight(i)/lambda;
            eps_veg=repsveg(i)-1i*iepsveg(i);
            
            % a: travel through vegetation.
            [E2a,H2a,E02a,H02a] = em_propagation_2lyrs(1,eps_veg,Ei0,false,freq,numWaves);
            damp1=abs(E2a(3,end))./abs(E2a(3,1)); % 1= not damped, 0=totally damped
            % b: reflection on soil. Input = row 3, plus damping factor.
            [E2b,H2b,E02b,H02b] = em_propagation_2lyrs(eps_veg,eps_soil,damp1*E02a(3),false,freq); % transmitting signal from vegetation as input (E02a(3))
            % travel back through vegetation. Input = row 2, plus damping factor
            damp2=abs(E2b(2,end))./abs(E2b(2,1)); % 1= not damped, 0=totally damped
            %% SO CHECK THIS: if Z_vegetation is higher than Z_soil, the reflection is different. That is rarely the case though (thick wet vegetation on completely dry ground is rare).
            %             if E02b(2)/E02b(1)<0
            %                 externalrefl=true; % E field is reversed. H field not. External reflection. Z2 > Z1.
            %                 internalrefl=true; % H field is reversed. E field not. Internal reflection. Z1 > Z2.
            %             end
            [E2c,H2c,E02c,H02c] = em_propagation_2lyrs(eps_veg,1,damp2*E02b(2),false,freq,numWaves); %
            
            if i==1
                Eu1=ones(10,1)*E02c(imode_trans);
                Hu1=ones(10,1)*H02c(imode_trans);
            end
            
            angleE02(i)=round(180*angle(E02c(imode_trans))/pi);
            angleH02(i)=round(180*angle(H02c(imode_trans))/pi);
            
            Eu2=ones(10,1)*E02c(imode_trans);
            Hu2=ones(10,1)*H02c(imode_trans);
            
            Hcomplconjprod=nan(size(Hu1));
            Ecomplconjprod=nan(size(Eu1));
            
            for iCoherence=1:length(Hu1)
                Hcomplconjprod(iCoherence)=Hu1(iCoherence)*conj(Hu2(iCoherence));
                Ecomplconjprod(iCoherence)=Eu1(iCoherence)*conj(Eu2(iCoherence));
            end
            
            Hcoh(i)=sum(Hcomplconjprod)/(sqrt(sum((abs(Hu1)).^2))*sqrt(sum((abs(Hu2)).^2)));
            Ecoh(i)=sum(Ecomplconjprod)/(sqrt(sum((abs(Eu1)).^2))*sqrt(sum((abs(Eu2)).^2)));
            
        end
    end
    absHcoh=abs(Hcoh);
    absEcoh=abs(Ecoh);
    
    percLOSDeformationE(:,iscenario) = lambda*(angleE02-angleE02(1))/360;
    percLOSDeformationH(:,iscenario) = lambda*(angleH02-angleH02(1))/360;
    
    if iscenario == 1
        reps2plot = reps;
        ieps2plot = ieps;
        repsveg2plot = repsveg;
        iepsveg2plot = iepsveg;
    end
end

%% PLOT

hfig=figure;
set(hfig, 'Position', [150, 150, 1000, 1100]);
subplot(5,1,1)
yyaxis left
plot(dateNums,soilMoisture,'LineWidth',1.5); hold on;
ax=gca;
ax.XTick=xtick;
ax.XTickLabel=xticklabel;
xlim([xtick(1)-5 xtick(end)])
ylim([0 0.5])
ylabel(['M_v (m^3/m^3)'])
ax.FontSize = 14;
if ~baresoil
    yyaxis right
    plot(dateNums,grassHeight,'LineWidth',1.5)
    ylim([0 0.5])
    ylabel(['grass height (m)'])
end
grid on

subplot(5,1,2)
yyaxis left
plot(dateNums,reps2plot,'LineWidth',1.5); hold on;
ylim([0 30])
ax=gca;
ax.XTick=xtick;
ax.XTickLabel=xticklabel;
xlim([xtick(1)-5 xtick(end)])
ylabel(['soil \epsilon'''])
yyaxis right
plot(dateNums,ieps2plot,'LineWidth',1.5)
ylim([0 10])
ax.FontSize = 14;
ylabel(['soil \epsilon'''''])
grid on

subplot(5,1,3)
yyaxis left
plot(dateNums,repsveg2plot,'LineWidth',1.5); hold on;
ylabel(['vegetatiion \epsilon'''])
ylim([0 6])
xlim([xtick(1)-5 xtick(end)])
yyaxis right
plot(dateNums,iepsveg2plot,'LineWidth',1.5)
ylabel(['vegetation \epsilon'''''])
ylim([0 6])
ax=gca;
ax.XTick=xtick;
ax.XTickLabel=xticklabel;
ax.FontSize = 14;
grid on

%% X-band repeat
percLOSDeformationE2xband = zeros(size(soilMoisture,1),numscenarios); % first one should be zero
percLOSDeformationH2xband = zeros(size(soilMoisture,1),numscenarios); % first one should be zero
Hcoh = nan(size(soilMoisture));
Ecoh = nan(size(soilMoisture));
angleE02 = nan(size(soilMoisture));
angleH02 = nan(size(soilMoisture));

reps=num(1,2); ieps = num(1,3); % tmp
eps_soil=reps(1)-1i*ieps(1);
[E1,H1,E01,H01] = em_propagation_2lyrs(1,eps_soil,Ei0,false,freq);

freq=5.405e9; % C-band
% freq=10e9; % X-band
lambda = 3e8/freq; % wavelength (needed later)

for iscenario = 1:numscenarios
    if iscenario == 1
        noisefactor = 0;
    else
        noisefactor = 0.05;
    end
    
    reps = num(:,2)+noisefactor*randn(size(num,1),1).*num(:,2);
    ieps = num(:,3)+noisefactor*randn(size(num,1),1).*num(:,3);
   
    if baresoil
        for i=1:length(reps)
            
            if i==1
                Eu1=ones(10,1)*E01(imode_refl);
                Hu1=ones(10,1)*H01(imode_refl);
            else
                Eu1=ones(10,1)*E02a(imode_refl);
                Hu1=ones(10,1)*H02a(imode_refl);
            end
            
            eps_soil=reps(i)-1i*ieps(i);
            [E2a,H2a,E02a,H02a] = em_propagation_2lyrs(1,eps_soil,1,false,freq);
            
            angleE02(i)=round(180*angle(E02a(imode_refl))/pi);
            angleH02(i)=round(180*angle(H02a(imode_refl))/pi);
            
            Eu2=ones(10,1)*E02a(imode_refl);
            Hu2=ones(10,1)*H02a(imode_refl);
            
            Hcomplconjprod=nan(size(Hu1));
            Ecomplconjprod=nan(size(Eu1));
            
            for iCoherence=1:length(Hu1)
                Hcomplconjprod(iCoherence)=Hu1(iCoherence)*conj(Hu2(iCoherence));
                Ecomplconjprod(iCoherence)=Eu1(iCoherence)*conj(Eu2(iCoherence));
            end
            
            Hcoh(i)=sum(Hcomplconjprod)/(sqrt(sum((abs(Hu1)).^2))*sqrt(sum((abs(Hu2)).^2)));
            Ecoh(i)=sum(Ecomplconjprod)/(sqrt(sum((abs(Eu1)).^2))*sqrt(sum((abs(Eu2)).^2)));
            
        end
    else
        % vegetation properties
        grassHeight = num(:,4);
        repsveg = num(:,6)+noisefactor*randn(size(num,1),1).*num(:,6);
        iepsveg = num(:,7)+noisefactor*randn(size(num,1),1).*num(:,7);        
          
        
        for i=1:length(reps)
            
            if i~=1
                Eu1=ones(10,1)*E02c(imode_trans);
                Hu1=ones(10,1)*H02c(imode_trans);
            end
            
            numWaves = grassHeight(i)/lambda;
            eps_veg=repsveg(i)-1i*iepsveg(i);
            
            % a: travel through vegetation.
            [E2a,H2a,E02a,H02a] = em_propagation_2lyrs(1,eps_veg,Ei0,false,freq,numWaves);
            damp1=abs(E2a(3,end))./abs(E2a(3,1)); % 1= not damped, 0=totally damped
            % b: reflection on soil. Input = row 3, plus damping factor.
            [E2b,H2b,E02b,H02b] = em_propagation_2lyrs(eps_veg,eps_soil,damp1*E02a(3),false,freq); % transmitting signal from vegetation as input (E02a(3))
            % travel back through vegetation. Input = row 2, plus damping factor
            damp2=abs(E2b(2,end))./abs(E2b(2,1)); % 1= not damped, 0=totally damped
            %% SO CHECK THIS: if Z_vegetation is higher than Z_soil, the reflection is different. That is rarely the case though (thick wet vegetation on completely dry ground is rare).
            %             if E02b(2)/E02b(1)<0
            %                 externalrefl=true; % E field is reversed. H field not. External reflection. Z2 > Z1.
            %                 internalrefl=true; % H field is reversed. E field not. Internal reflection. Z1 > Z2.
            %             end
            [E2c,H2c,E02c,H02c] = em_propagation_2lyrs(eps_veg,1,damp2*E02b(2),false,freq,numWaves); %
            
            if i==1
                Eu1=ones(10,1)*E02c(imode_trans);
                Hu1=ones(10,1)*H02c(imode_trans);
            end
            
            angleE02(i)=round(180*angle(E02c(imode_trans))/pi);
            angleH02(i)=round(180*angle(H02c(imode_trans))/pi);
            
            Eu2=ones(10,1)*E02c(imode_trans);
            Hu2=ones(10,1)*H02c(imode_trans);
            
            Hcomplconjprod=nan(size(Hu1));
            Ecomplconjprod=nan(size(Eu1));
            
            for iCoherence=1:length(Hu1)
                Hcomplconjprod(iCoherence)=Hu1(iCoherence)*conj(Hu2(iCoherence));
                Ecomplconjprod(iCoherence)=Eu1(iCoherence)*conj(Eu2(iCoherence));
            end
            
            Hcoh(i)=sum(Hcomplconjprod)/(sqrt(sum((abs(Hu1)).^2))*sqrt(sum((abs(Hu2)).^2)));
            Ecoh(i)=sum(Ecomplconjprod)/(sqrt(sum((abs(Eu1)).^2))*sqrt(sum((abs(Eu2)).^2)));
            
        end
    end
    absHcoh=abs(Hcoh);
    absEcoh=abs(Ecoh);
    
    percLOSDeformationE2xband(:,iscenario) = lambda*(angleE02-angleE02(1))/360;
    percLOSDeformationH2xband(:,iscenario) = lambda*(angleH02-angleH02(1))/360;
    
    if iscenario == 1
        reps2plot = reps;
        ieps2plot = ieps;
        repsveg2plot = repsveg;
        iepsveg2plot = iepsveg;
    end
end

subplot(5,1,4:5)
% plot(percLOSDeformationE*1000,'LineWidth',1.5); hold on;
plot(dateNums,percLOSDeformationE*1000,'-','LineWidth',0.5,'Color',[0.9 0.9 0.9]); hold on
plot(dateNums,percLOSDeformationE2xband*1000,'-.','LineWidth',0.5,'Color',[0.8 0.8 0.8]);
% plot(percLOSDeformationH*1000,'LineWidth',0.5);
ylim([-10 10])
xlim([xtick(1)-5 xtick(end)])
%% TODO: plot X-band and C-band
ylabel(['Perceived LOS deformation (mm)'])
ax=gca;
ax.XTick=xtick;
ax.YTick=-10:2:10;
ax.XTickLabel=xticklabel;
ax.FontSize = 14;
grid on

leg2=plot(dateNums,percLOSDeformationE2xband(:,1)*1000,'-','LineWidth',2,'Color',[0.2 0.5 0.2]); 
leg1=plot(dateNums,percLOSDeformationE(:,1)*1000,'-','LineWidth',2,'Color',[0.1 0.1 0.1]); 

% plot(percLOSDeformationH*1000,'LineWidth',0.5);
str={'X-band','C-band'};
hleg = legend([leg1,leg2],str,'Location','Northwest');