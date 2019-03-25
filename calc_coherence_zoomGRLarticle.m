%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
cc

baresoil=false;
plotcoherence=true;

% imode_inc=1; % incident wave
imode_refl=2; % reflection
imode_trans=3; % transmitted
Ei0=1;
eps_soil=12-1i*2; % the reference soil. A dry sand
% eps_soil=24-1i*5;% the references soil: silty clay Hallikainen

[E1,H1,E01,H01] = em_propagation_2lyrs(1,eps_soil,Ei0,false);

if baresoil
    % bare soil
    % calculate the E and H field and coherence
    if plotcoherence
        reps=2:0.5:30;
        ieps=0:0.5:20;
        plotornot=false;
    else
        % % testing,
%         plotornot=true;
%         reps=3;
%         ieps=0;
    end
    for ireal=1:length(reps)
        for iimag=1:length(ieps)
            eps_soil=reps(ireal)-1i*ieps(iimag);
            [E2a,H2a,E02a,H02a] = em_propagation_2lyrs(1,eps_soil,1,false);
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
    ytick=ieps(1):10:ieps(end);
else
    if plotcoherence
        reps=1:0.1:10;
        ieps=0:0.02:2;
        plotornot=false;
    else
        % % testing,
        plotornot=true;
        reps=5;
        ieps=0;
    end
    for ireal=1:length(reps)
        for iimag=1:length(ieps)
            externalrefl=false;
            %             eps_veg = 5.55-1i*1.25;
            eps_veg=reps(ireal)-1i*ieps(iimag);
            % first, travel through vegetation.
            [E2a,H2a,E02a,H02a] = em_propagation_2lyrs(1,eps_veg,Ei0,plotornot);
            % reflection on soil. Input = row 3, plus damping factor.
            damp1=abs(E2a(3,end))./abs(E2a(3,1)); % 1= not damped, 0=totally damped
            [E2b,H2b,E02b,H02b] = em_propagation_2lyrs(eps_veg,eps_soil,damp1*E02a(3),plotornot);
            % travel back through vegetation. Input = row 2, plus damping factor
            damp2=abs(E2b(2,end))./abs(E2b(2,1)); % 1= not damped, 0=totally damped
            if E02b(2)/E02b(1)<0
                externalrefl=true; % E field is reversed. H field not. External reflection. Z2 > Z1.
                internalrefl=true; % H field is reversed. E field not. Internal reflection. Z1 > Z2.
            end
            [E2c,H2c,E02c,H02c] = em_propagation_2lyrs(eps_veg,1,damp2*E02b(2),plotornot,externalrefl); % reverseE
            
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
    ytick=ieps(1):0.5:ieps(end);
end
realHcoh=real(Hcoh);
realEcoh=real(Ecoh);
angleHcoh=180*imag(Hcoh)/pi;
angleEcoh=180*imag(Ecoh)/pi;

if plotcoherence
    % plot of amplitude and phase of E02 and H02 of the reflected wave at the boundary in medium 1
    ssize=get(0,'ScreenSize');
    hfig = figure;
    %     cm=colormap('gray');
    cm=flipud(colormap(jet(128)));
    cm=colormap(cm);
    set(hfig, 'Position', [ssize(1)+10, ssize(2)+50, ssize(3)*.7, ssize(4)*.42]);
    subplot(2,2,1)
    plot_imagesc_rogier_angle(angleE02,reps,ieps,'\Phi E_0 (Degrees)',cm)
    set(gca,'Xtick',xtick,'XTickLabel',xtick)
    set(gca,'Ytick',ytick,'YTickLabel',ytick)
    subplot(2,2,2)
    plot_imagesc_rogier_angle(angleH02,reps,ieps,'\Phi H_0(Degrees)',cm)
    set(gca,'Xtick',xtick,'XTickLabel',xtick)
    set(gca,'Ytick',ytick,'YTickLabel',ytick)
    subplot(2,2,3)
    plot_imagesc_rogier_angle(abs(E02),reps,ieps,'|E_0| / E_{0,air}',cm)
    set(gca,'Xtick',xtick,'XTickLabel',xtick)
    set(gca,'Ytick',ytick,'YTickLabel',ytick)
    subplot(2,2,4)
    plot_imagesc_rogier_angle(abs(H02)/H02a(1),reps,ieps,'|H_0| / H_{0,air}',cm)
    set(gca,'Xtick',xtick,'XTickLabel',xtick)
    set(gca,'Ytick',ytick,'YTickLabel',ytick)

    ssize=get(0,'ScreenSize');
    hfig = figure;
    %     cm=colormap('gray');
    cm=flipud(colormap(jet(128)));
    cm=colormap(cm);
    set(hfig, 'Position', [ssize(1)+10, ssize(2)+50, ssize(3)*.7, ssize(4)*.42]);
    subplot(2,2,1)
    plot_imagesc_rogier_angle((realEcoh),reps,ieps,'Re(\gamma)_E',cm) 
    set(gca,'Xtick',xtick,'XTickLabel',xtick)
    set(gca,'Ytick',ytick,'YTickLabel',ytick)
    subplot(2,2,2)
    plot_imagesc_rogier_angle((realHcoh),reps,ieps,'Re(\gamma)_H',cm)
    set(gca,'Xtick',xtick,'XTickLabel',xtick)
    set(gca,'Ytick',ytick,'YTickLabel',ytick)
    subplot(2,2,3)
    plot_imagesc_rogier_angle((angleEcoh),reps,ieps,'\angle\gamma_E',cm)
    set(gca,'Xtick',xtick,'XTickLabel',xtick)
    set(gca,'Ytick',ytick,'YTickLabel',ytick)
    subplot(2,2,4)
    plot_imagesc_rogier_angle((angleHcoh),reps,ieps,'\angle\gamma_H',cm)
    set(gca,'Xtick',xtick,'XTickLabel',xtick)
    set(gca,'Ytick',ytick,'YTickLabel',ytick)
    tightfig;
end