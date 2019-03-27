function [E,H,E0,H0] = em_propagation_2lyrs(eps_1,eps_2,varargin)

% plotornot=true;
% eps_1=1;
% eps_2=3-1i*1;
% varargin(1)=Ei0=1 by defeault;
%               Other options are giving this as a 2col vector, 
%               with Hi0 the second input.
% varargin(2)=plotornot=false by default

%OUTPUT
% E and H field in 3 cases, each case one row:
%   1 - Original incoming field in medium 1 (e.g. air)
%   2 - Reflected field from wave travelling from medium 1 to medium 2
%   3 - Transmitted field in medium 2
%
%   including the z-vector, indicating distance.
% nargin

if nargin < 6
    numberoflambdas=2;  
else
    numberoflambdas=varargin{4};
end

if nargin < 5
     freq=5.405e9;  % freq in Hz (C-band)
else
    freq=varargin{3};
end

if nargin < 4
     plotornot=false;  
else
    plotornot=logical(varargin{2});
end

if nargin < 3
    Ei0=1;
else
    EandHi0=varargin{1};
end
Ei0=EandHi0(1);
    
% START user input
% freq=5e9;%Hz
% freq=5.405e9; % Hz, C-band Sentinel-1
% definitions of MEDIUM 1 and 2 % column vectors for layer 1 (col 1) and layer 2(col 2)
% eps_soil=3;
% eps_soil=75-1i*30;%saturated sand (fresh-brackish)
% eps_soil=24-1i*5;%silty clay Hallikainen
% eps_veg=5.55-1i*1.25; % vegetation layer
% END user input

epsilon_0=8.854187817e-12;
mu_0=4e-7*pi;
% mu=mu_0*[1.00000037,1];%col vector for layer1 and layer2

mu=mu_0*[1,1];
epsilon=epsilon_0*[eps_1,eps_2];%col vector for layer1 and layer2 %dry sand
omega=2*pi*freq; %angular frequency
eta=sqrt(mu./epsilon);
k=omega*sqrt(mu.*epsilon);
lambda=1./(freq.*sqrt(mu.*epsilon));% m

% Calculate E0 and H0s
% calculate wavefield in three media

% PROPAGATION IN MEDIUM 1
numzsteps=200;
% Ei0=1;
t=0;
% calculate travel through the first medium

%% incoming field
z=zeros(length(lambda),numzsteps); % calculate z for both media
% stick to figure 32-2 of Lorrain and Corson and its x-y-z conventions
% Eq. 32-1 to 32-3
for iz=1:numzsteps
    z(:,iz)=abs(numberoflambdas*lambda(:))*(iz-1)/numzsteps;
end

E_I=Ei0.*exp(1i*(omega*t-k(1).*z(1,:)));

if length(EandHi0)<2
    Hi0=Ei0/eta(1);
    H_I=E_I/eta(1);    
else
    Hi0=EandHi0(2);
    H_I=Hi0.*exp(1i*(omega*t-k(1).*z(1,:)));
end

%% transmission and reflection from medium 1 to 2
% total fields in two media:
Er0=Ei0*(k(1)/k(2)-1)/(k(1)/k(2)+1); % 32-14 for E normal to incidence plane
Hr0=-Hi0*(k(1)/k(2)-1)/(k(1)/k(2)+1);

Et0=Ei0*(2*k(1)/k(2))/(k(1)/k(2)+1); % 32-15 for E normal to incidence plane
Ht0=(eta(1)/eta(2))*Hi0*(2*k(1)/k(2))/(k(1)/k(2)+1);

E0(1)=Ei0;
E0(2)=Er0;
E0(3)=Et0;

H0(1)=Hi0;
H0(2)=Hr0;
H0(3)=Ht0;

E_R=Er0*exp(1i*(omega*t+k(1).*z(1,:)));
E_T=Et0*exp(1i*(omega*t-k(2)*z(2,:)));
H_R=Hr0*exp(1i*(omega*t+k(1).*z(1,:)));
H_T=Ht0*exp(1i*(omega*t-k(2)*z(2,:)));

% check if boundary conditions are met
Emed1=Ei0+Er0;
Emed2=Et0;
if round(Emed1*1e6) ~= round(Emed2*1e6)
    warning(['Boundary condition not met: E in medium 1 = ',num2str(Emed1),' E in medium 2 = ',num2str(Emed2)])
end

Hmed1=Hi0+Hr0;
Hmed2=Ht0;
if round(Hmed1*1e5) ~= round(Hmed2*1e5)
    warning(['Boundary condition not met: H in medium 1 = ',num2str(Hmed1),' H in medium 2 = ',num2str(Hmed2)])
end

E=zeros(3,numzsteps);H=E; % 1=incoming, 2=reflected, 3=transmitted
E(1,:)=E_I;
H(1,:)=H_I;
E(2,:)=E_R;
H(2,:)=H_R;
E(3,:)=E_T;
H(3,:)=H_T;

zplot=z;
zplot(2,:)=z(1,:);
zplot(3,:)=z(2,:);

if plotornot
    for istep=1:1:3
        PAR1=E(istep,:);
        PAR2=H(istep,:);
        PAR1_0=E0(istep);
        PAR2_0=H0(istep);
        zz=zplot(istep,:);
        zlim=[0 max(zplot(:))];
        
        % % % % back to 2D array
        % % % if length(size(PAR1))>2
        % % % PAR1=reshape(PAR1,[numzsteps numtsteps]);
        % % % PAR2=reshape(PAR2,[numzsteps numtsteps]);
        % % % end
        % % %
        % calculate the phase angles in degrees
        if istep==1
            angle1=round(180*angle(PAR1_0)/pi);
            angle2=round(180*angle(PAR2_0)/pi);
            Htitle=['Phase H_{i} = ',num2str(angle2),'^{\circ}'];
            Etitle=['Phase E_{i} = ',num2str(angle1),'^{\circ}'];
        end
        if istep==2
            angle1=round(180*angle(PAR1_0)/pi);
            angle2=round(180*angle(PAR2_0)/pi);
            Htitle=['Phase H_{r} = ',num2str(angle2),'^{\circ}'];
            Etitle=['Phase E_{r} = ',num2str(angle1),'^{\circ}'];
        end
        if istep==3
            angle1=round(180*angle(PAR1_0)/pi);
            angle2=round(180*angle(PAR2_0)/pi);
            Htitle=['Phase H_{t} = ',num2str(angle2),'^{\circ}'];
            Etitle=['Phase E_{t} = ',num2str(angle1),'^{\circ}'];
        end
        
        %% START PLOTTING
        ssize=get(0,'ScreenSize');
        hfig = figure;
        set(hfig, 'Position', [ssize(1)+10, ssize(2)+50+10, ssize(3)-100, ssize(4)/2]);
        
        subplot(1,7,1);
        h=compass(PAR1(1,1));
        title(Etitle,'Fontsize',16)
        set(h,'Linewidth',2,'Color',[0 0 1]);
        
        rPAR1=real(PAR1);
        rPAR2=real(PAR2);
        subplot(1,7,(3:5));
        % [AX,H1,H2] = plotyy(t,rPAR1(1,:),t,rPAR2(1,:),'plot');
        [AX,H1,E2] = plotyy(zz,rPAR1,zz,rPAR2,'plot');
        set(AX,'XLim',zlim);
        set(get(AX(1),'Ylabel'),'String','E_{x} (E_{0})','Fontsize',16,'Color',[0 0 1]);
        set(get(AX(2),'Ylabel'),'String','H_{y} (E_{0})','Fontsize',16,'Color',[1 0 0]);
        set(AX(1),'YColor',[0 0 1],'Fontsize',16,'ylim',[-1 1],'ytick',-1:0.2:1);
        set(AX(2),'YColor',[1 0 0],'Fontsize',16,'ylim',[-4e-3 4e-3],'ytick',-4e-3:1e-3:4e-3);
        set(H1,'LineStyle',':','LineWidth',2.5,'Color',[0 0 1]);
        set(E2,'LineStyle',':','LineWidth',2.5,'Color',[1 0 0]);
        % title(['step ',num2str(ilayer),' ',titletext(ilayer)],'Fontsize',18);
        % xlabel('time (s)','Fontsize',16);
        xlabel('z (m)','Fontsize',16);
        
        subplot(1,7,7);
        h2=compass(PAR2(1,1));
        title(Htitle,'Fontsize',16)
        set(h2,'Linewidth',2,'Color',[1 0 0]);
        
        % filename=['i:\GroundWater\Research\Smart\Satellite RS\Latex\images\EandH_step',num2str(ilayer)];
        % saveas(hfig,[filename,'.eps'],'eps')%use eps later
    end
end