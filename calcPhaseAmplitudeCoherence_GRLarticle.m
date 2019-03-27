% Code belonging to manucscript 'Explanation of InSAR phase disturbances by seasonal characteristics of soil and vegetation'.
% Rogier Westerhoff and Moira Steyn-Ross
%   This script is used to calculate phases, amplitudes and coherence of microwave radar

cc

baresoil=true;
plotcoherence=true;

freq=5.405e9;
numlambdas = 3;

% imode_inc=1; % incident wave
imode_refl=2; % reflection
imode_trans=3; % transmitted
Ei0=1;
% eps_soil=3; % the reference soil. A dry sand
% eps_soil=24-1i*5;% the references soil: silty clay Hallikainen
eps_soil=12-1i*2; minmax = [-180 -130];% the references soil: sandly loam (Mv=0.2) Hallikainen
% eps_soil=75-1i*30;% fresh-brackish water-saturated
% eps_soil=5.55-1i*1.25;% vegetation
% eps_soil=2-1i*1;% light vegetation

[E1,H1,E01,H01] = em_propagation_2lyrs(1,eps_soil,Ei0,false);

if baresoil
    % bare soil
    % calculate the E and H field and coherence
    if plotcoherence
        reps=2:0.1:30;
        ieps=0:0.1:10;
        plotornot=false;
    else
        % % testing,
        plotornot=true;
        reps=3;
        ieps=0;
    end
    for ireal=1:length(reps)
        for iimag=1:length(ieps)
            eps_soil=reps(ireal)-1i*ieps(iimag);
            
            [E2a,H2a,E02a,H02a] = em_propagation_2lyrs(1,eps_soil,1,plotornot);
            E02(iimag,ireal)=E02a(imode_refl);
            H02(iimag,ireal)=H02a(imode_refl);
            angleE02(iimag,ireal)=round(180*angle(E02a(imode_refl))/pi);
            angleH02(iimag,ireal)=round(180*angle(H02a(imode_refl))/pi);
            
            Eu1=ones(10,1)*E01(imode_refl);
            Eu2=ones(10,1)*E02(iimag,ireal);
            Hu1=ones(10,1)*H01(imode_refl);
            Hu2=ones(10,1)*H02(iimag,ireal);
            
            for i=1:length(Hu1)
                Hcomplconjprod(i)=Hu1(i)*conj(Hu2(i));
                Ecomplconjprod(i)=Eu1(i)*conj(Eu2(i));
            end
            Hcoh(iimag,ireal)=sum(Hcomplconjprod)/(sqrt(sum((abs(Hu1)).^2))*sqrt(sum((abs(Hu2)).^2)));
            Ecoh(iimag,ireal)=sum(Ecomplconjprod)/(sqrt(sum((abs(Eu1)).^2))*sqrt(sum((abs(Eu2)).^2)));
        end
    end
    xtick=[2,reps(reps==10):10:reps(end)];
    ytick=ieps(1):5:ieps(end);
    xplot=[3,24,12,5.55,2]; % different soils, incl vegetation
    yplot=[0,5,2,1.25,1];  % different soils, incl vegetation
else
    
    if plotcoherence
        reps=1:0.1:10;
        ieps=0.02:0.02:3;
        plotornot=false;

    else
        % % testing,
        plotornot=true;
        reps=5;
        ieps=1;
    end
    for ireal=1:length(reps)
        for iimag=1:length(ieps)
            externalrefl=false;
%                          eps_veg = 5.55-1i*1.25;
            eps_veg=reps(ireal)-1i*ieps(iimag);
            % first, travel through vegetation.
            [E2a,H2a,E02a,H02a] = em_propagation_2lyrs(1,eps_veg,Ei0,false,freq,numlambdas);
            % reflection on soil. Input = row 3, plus damping factor.
            damp1=abs(E2a(3,end))./abs(E2a(3,1)); % 1= not damped, 0=totally damped
            [E2b,H2b,E02b,H02b] = em_propagation_2lyrs(eps_veg,eps_soil,damp1*E02a(3),false,freq,numlambdas);
            % travel back through vegetation. Input = row 2, plus damping factor
            damp2=abs(E2b(2,end))./abs(E2b(2,1)); % 1= not damped, 0=totally damped
%              if E02b(2)/E02b(1)>0
%                  warning('E field reversed, since Z1>Z2')
%              end
%                 externalrefl=true; % E field is reversed. H field not. External reflection. Z2 > Z1.
%                 internalrefl=true; % H field is reversed. E field not. Internal reflection. Z1 > Z2.
%             end
            [E2c,H2c,E02c,H02c] = em_propagation_2lyrs(eps_veg,1,damp2*E02b(2),false,freq,numlambdas); 
            
            E02(iimag,ireal)=E02c(imode_trans);
            H02(iimag,ireal)=H02c(imode_trans);
            angleE02(iimag,ireal)=round(180*angle(E02c(imode_trans))/pi);
            angleH02(iimag,ireal)=round(180*angle(H02c(imode_trans))/pi);
            
            Eu1=ones(10,1)*E01(imode_refl);
            Eu2=ones(10,1)*E02(iimag,ireal);
            Hu1=ones(10,1)*H01(imode_refl);
            Hu2=ones(10,1)*H02(iimag,ireal);
            
            for i=1:length(Hu1)
                Hcomplconjprod(i)=Hu1(i)*conj(Hu2(i));
                Ecomplconjprod(i)=Eu1(i)*conj(Eu2(i));
            end
            Hcoh(iimag,ireal)=sum(Hcomplconjprod)/(sqrt(sum((abs(Hu1)).^2))*sqrt(sum((abs(Hu2)).^2)));
            Ecoh(iimag,ireal)=sum(Ecomplconjprod)/(sqrt(sum((abs(Eu1)).^2))*sqrt(sum((abs(Eu2)).^2)));
            
        end
    end
    xtick=[1,reps(reps==2):1:reps(end)];
    ytick=0:0.5:ieps(end);
    xplot=[5.55,2]; % different vegetation
    yplot=[1.25,1]; % different vegetation
end
realHcoh=real(Hcoh);
realEcoh=real(Ecoh);
angleHcoh=180*imag(Hcoh)/pi;
angleEcoh=180*imag(Ecoh)/pi;

if plotcoherence
    % Optional TODO: build mask from real(Ecoh < 0.3). 
    % Not needed now, since no values < 0.3
    
    ssize=get(0,'ScreenSize');
    hfig = figure;
    %     cm=colormap('gray');
    cm=flipud(colormap(jet(128)));
    cm=colormap(cm);
    
    angleE02(abs(realEcoh)<0.3)=nan;
    E02(abs(realEcoh)<0.3)=nan;
    
    set(hfig, 'Position', [ssize(1)+100, ssize(2)+150, 1000, 700]);
    subplot(2,1,1)
    if baresoil
        plot_imagesc_rogier_angle(angleE02,reps,ieps,'\Phi E_r (Degrees)',cm); hold on
    else
        plot_imagesc_rogier_angle(angleE02,reps,ieps,'\Phi E_r (Degrees)',cm,minmax); hold on
    end
    plot(xplot,yplot,'wo','MarkerSize',6,'MarkerFaceColor',[0.8 0.8 0.8])
    ax=gca;
    ax.XTick=xtick;ax.XTickLabel=xtick; 
    ax.YTick=ytick;ax.YTickLabel=ytick; 
    ax.FontSize=14;
    subplot(2,1,2)
    plot_imagesc_rogier_angle(abs(E02),reps,ieps,'|E_r| / E_{0,air}',cm)
    ax=gca;
    ax.XTick=xtick;ax.XTickLabel=xtick; 
    ax.YTick=ytick;ax.YTickLabel=ytick; 
    ax.FontSize=14;
    hold off
    
 
% % plot of amplitude and phase of E02 and H02 of the reflected wave at the boundary in medium 1
%     ssize=get(0,'ScreenSize');
%     hfig = figure;
%     %     cm=colormap('gray');
%     cm=flipud(colormap(jet(128)));
%     cm=colormap(cm);
%     set(hfig, 'Position', [ssize(1)+10, ssize(2)+50, ssize(3)*.7, ssize(4)*.42]);
%     subplot(2,2,1)
%     plot_imagesc_rogier_angle(angleE02,reps,ieps,'\Phi E_0 (Degrees)',cm); hold on
%     plot(xplot,yplot,'wo','MarkerSize',6,'MarkerFaceColor','w')
%     ax=gca;
%     ax.XTick=xtick;ax.XTickLabel=xtick; 
%     ax.YTick=ytick;ax.YTickLabel=ytick; 
%     ax.FontSize=14;
%     subplot(2,2,2)
%     plot_imagesc_rogier_angle(angleH02,reps,ieps,'\Phi H_0(Degrees)',cm)
%     ax=gca;
%     ax.XTick=xtick;ax.XTickLabel=xtick; 
%     ax.YTick=ytick;ax.YTickLabel=ytick; 
%     ax.FontSize=14;
%     subplot(2,2,3)
%     plot_imagesc_rogier_angle(abs(E02),reps,ieps,'|E_0| / E_{0,air}',cm)
%     ax=gca;
%     ax.XTick=xtick;ax.XTickLabel=xtick; 
%     ax.YTick=ytick;ax.YTickLabel=ytick; 
%     ax.FontSize=14;
%     subplot(2,2,4)
%     plot_imagesc_rogier_angle(abs(H02)/H02a(1),reps,ieps,'|H_0| / H_{0,air}',cm)
%     ax=gca;
%     ax.XTick=xtick;ax.XTickLabel=xtick; 
%     ax.YTick=ytick;ax.YTickLabel=ytick; 
%     ax.FontSize=14;
%     hold off
% %     tightfig;
    
%     ssize=get(0,'ScreenSize');
%     hfig = figure;
%     %     cm=colormap('gray');
%     cm=flipud(colormap(jet(128)));
%     cm=colormap(cm);
%     set(hfig, 'Position', [ssize(1)+10, ssize(2)+50, ssize(3)*.7, ssize(4)*.42]);
%     subplot(2,2,1)
%     plot_imagesc_rogier_angle((realEcoh),reps,ieps,'Re(\gamma)_E',cm) ; hold on
%     plot(xplot,yplot,'wo','MarkerSize',6,'MarkerFaceColor','w')
%     ax=gca;
%     ax.XTick=xtick;ax.XTickLabel=xtick; 
%     ax.YTick=ytick;ax.YTickLabel=ytick; 
%     ax.FontSize=14;
%     subplot(2,2,2)
%     plot_imagesc_rogier_angle((realHcoh),reps,ieps,'Re(\gamma)_H',cm)
%     ax=gca;
%     ax.XTick=xtick;ax.XTickLabel=xtick; 
%     ax.YTick=ytick;ax.YTickLabel=ytick; 
%     ax.FontSize=14;
%     subplot(2,2,3)
%     plot_imagesc_rogier_angle((angleEcoh),reps,ieps,'\angle\gamma_E',cm)
%     ax=gca;
%     ax.XTick=xtick;ax.XTickLabel=xtick; 
%     ax.YTick=ytick;ax.YTickLabel=ytick; 
%     ax.FontSize=14;
%     subplot(2,2,4)
%     plot_imagesc_rogier_angle((angleHcoh),reps,ieps,'\angle\gamma_H',cm)
% %     set(gca,'Xtick',xtick,'XTickLabel',xtick)
% %     set(gca,'Ytick',ytick,'YTickLabel',ytick)
%     ax=gca;
%     ax.XTick=xtick;ax.XTickLabel=xtick; 
%     ax.YTick=ytick;ax.YTickLabel=ytick; 
%     ax.FontSize=14;
%     hold off
% %     tightfig;
end