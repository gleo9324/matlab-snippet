%%%%%%%%%%% energy harvester parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% equations parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = 1;                                                % modal mass [m]
Cp = e33*b*Lp/hp;                                     % equivalent piezo capacity [F]
% Theta = -(Ep*d31*b)/(1-nip^2)*(hp/2+hs-z0);         % electro-mechanical coefficient [As/m]
Theta = -(Ep*d31*b)/(2*hp)*(hc^2-hb^2);               % electro-mechanical coefficient [As/m]

dC = 0;
A = 1;                                                % forcing amplitude [Kg*m/s^2] 

for i=1:3
    omr(i) = 2*pi*freq(i);                            % resonance pulsation [rad/s]
end
R = [100,10^3,10^4,10^5,10^6];                        % resistance [Ohm] 100,1000,10000,100000,10^6

df = 0.1;               % frequency step size

%%%%%%%%%%%%%% Analytical solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
nf = 15000;

 vfft2 = zeros([length(R),nf]);
 v2 = zeros([length(R),nf]);
 vfft = zeros([1,nf]);
 v = zeros([1,nf]);
 om = zeros([1,nf]);
%  G = zeros([1,10000]);
%  H = zeros([1,10000]);
 Tau = 0;

 
for p=1:length(R)
    
    Tau = Cp*R(p);
    
 for j=1:nf
     
    omj(j) = j*df*2*pi;
    
    G = zeros([1,nf]);
    H = zeros([1,nf]);
    
    for k=1:nmodes
        
        Chi(k) = Theta*chi(k);
        Phi(k) = Theta/Cp*chi(k);%-(d31*Yp*hpc*hp)/(e33*L)*chi(k);
        Zitar(k) = cs*omr(k)/2+ca/(2*omr(k));%C(k);
%        G(j) = G(j)+(-1i*m(0)*omj(j)*Phi(k)*phi3(k))/(omr(k)^2-omj(j)^2+1i*2*Zitar(k)*omr(k)*omj(j));
        G(j) = G(j)+(-1i*m*omj(j)*Phi(k)*phi_tre(k))/(omr(k)^2-omj(j)^2+1i*2*Zitar(k)*omr(k)*omj(j));
        H(j) = H(j)+(1i*omj(j)*Chi(k)*Phi(k))/(omr(k)^2-omj(j)^2+1i*2*Zitar(k)*omr(k)*omj(j));
        
    end
    
    v(j) = G(j)/(H(j)+(1+1i*omj(j)*Tau)/Tau);
    vfft(j) = abs(v(j));
    v2(p,j) = v(j);
    vfft2(p,j) = vfft(j);
        
 end
 
    figure(p)
    plot(omj/(2*pi),vfft,'r');
    hold on
        
end

% Simulink session %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% punti vicini solo alle frequenze di interesse
% dom1 = [0.9,0.95,0.96,0.97,0.98,0.99,1,1.01,1.02,1.03,1.04,1.05,1.1];
% % dom2 = [0.96,0.97,0.98,0.99,1,1.01,1.02,1.03,1.04,1.05];
% dom2 = [0.9,0.95,0.96,0.97,0.98,0.99,1,1.01,1.02,1.03,1.04,1.05,1.1];
% dom3 = [0.9,0.95,0.96,0.97,0.98,0.99,1,1.01,1.02,1.03,1.04,1.05,1.1];


% maggior numero di punti
 dom1 = [0,0.6,0.7,0.8,0.9,0.95,0.96,0.97,0.98,0.99,1,1.01,1.02,1.03,1.04,1.05,1.1,1.2,1.3,1.4,1.7,2,2.3,2.6,3,3.5,4,4.5];
dom2 = [0.7,0.8,0.9,0.95,0.96,0.97,0.98,0.99,1,1.01,1.02,1.03,1.04,1.05,1.1,1.2,1.3,1.4,1.7,2,2.3,2.6];
dom3 = [0.9,0.95,0.96,0.97,0.98,0.99,1,1.01,1.02,1.03,1.04,1.05,1.1,1.2];

for k=1:3
    if (k==1)
        for j=1:length(dom1)
        omj2(j) = dom1(j)*omr(k);
        th(j) = 10^(-7);
        nfft(j) = 10;
        end
    elseif(k==2)
        for j=1:length(dom2)
            omj2(length(dom1)+j) = dom2(j)*omr(k);
            th(length(dom1)+j) = 10^(-6);
            nfft(length(dom1)+j) = 50;
        end
%       if(k==2)
%         for j=1:length(dom2)
%             omj2(j) = dom2(j)*omr(k);
%             th(j) = 10^(-6);
%             nfft(j) = 50;
%         end
    elseif(k==3)
        for j=1:length(dom3)
            omj2(length(dom1)+length(dom2)+j) = dom3(j)*omr(k);
            th(length(dom1)+length(dom2)+j) = 10^(-7);
            nfft(length(dom1)+length(dom2)+j) = 50;
        end
    end
end

% th = 10^(-7);
% nfft = 10;

for p=1:length(R)
 for j=1:length(omj2)
  
%     omj2(j) = j*df*2*pi;
     
    options = odeset('MaxStep',10^(-4));
    tspan = [0 2];
    y0 = [0 0 0 0 0 0 0];
    
    [t,y] = ode45(@(t,y) odefcn(t,y,cmod,kmod,Theta,phi_due,phi_uno,A,fmod,omj2(j),chi,Cp,R(p)),tspan,y0,options);
    
   v = y(:,7);
       
    % Tagliare segnale

    if (j>1)
        [pks,locs] = findpeaks(v,t);
        prova = zeros(size(pks));
        
        for i=2:length(locs)
            prova(i-1) = abs(pks(i)-pks(i-1));
            if prova(i-1) < th(j)
                itrans = find(t==locs(i));
                break
            end
        end
    else
        itrans = 1;
    end
    
       
    Y = v(itrans:end);
    x = t(itrans:end);

%     figure(j)
%     plot(x,Y,'x',t,v);
%     plot(t,v);

%     figure(1000)
%     plot(t,y(:,1))
    %xlim([0 0.5]);
%%%%% fft  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    n1 = 2^(nextpow2(length(x)*nfft(j)));           % fft funziona meglio con potenze di 2. Più si aumenta, più è precisa, ma c'è comunque il rumore.   
    
    Z1 = fft(Y,n1);                              % creo trasformata calcolata in un numero di punti pari alla lunghezza del vettore tempo
    Z2 = abs(Z1/length(Y));                      % prendo il modulo del numero immaginario e lo moltiplico per il sampling rate (che non è considerato diverso da 1 in MATLAB). 
                                                 % Se dt=1 allora divido per length(x)
    Z3 = 2*Z2(1:n1/2+1);                         % taglio metà del segnale, altrimenti ho tutto speculare. Ottengo la semiampiezza
    f1 = (1/(x(end)-x(1)))*(0:n1/2);             % creo il vettore frequenza che è l'inverso del periodo
    

    % Dal segnale tagliato il transitorio
    Vsimfft2(j) = max(Z3);
    fft_V2(p,j) = Vsimfft2(j);
    
   
 end   
end


for p=1:length(R)
      
    [pks10,locs10] = findpeaks(fft_V2(p,:));
    [pks30,locs30] = findpeaks(vfft2(p,:));
    
    figure(p*100)

    plot(omj/(2*pi),vfft2(p,:),'r',omj2/(2*pi),fft_V2(p,:)/(2*pi),'g',omj(locs30)/(2*pi),pks30,'r*',omj2(locs10)/(2*pi),pks10/(2*pi),'go');
    
    hold on


    legend('Erturk','V/Acc','fftCutSign')
    
end


if length(R)==5

[pks60,locs60] = findpeaks(fft_V2(1,:));
[pks70,locs70] = findpeaks(fft_V2(2,:));
[pks80,locs80] = findpeaks(fft_V2(3,:));
[pks90,locs90] = findpeaks(fft_V2(4,:));
[pks100,locs100] = findpeaks(fft_V2(5,:));


figure(1000)

hold on

p2(1) = plot(omj2(1)/(2*pi),fft_V2(1,1),'ks','MarkerSize',25);
p2(2) = plot(omj(1)/(2*pi),vfft2(1,1),'k','MarkerSize',25);

% analitica
p1 = semilogy(omj/(2*pi),vfft2(1,:),'b',omj/(2*pi),vfft2(2,:),'r',omj/(2*pi),vfft2(3,:),'g',omj/(2*pi),vfft2(4,:),'m',omj/(2*pi),vfft2(5,:),'c','MarkerSize',25);

% numerica
semilogy(omj2/(2*pi),fft_V2(1,:),'bs',omj2/(2*pi),fft_V2(2,:),'rs',omj2/(2*pi),fft_V2(3,:),'gs',omj2/(2*pi),fft_V2(4,:),'ms',omj2/(2*pi),fft_V2(5,:),'ks','MarkerSize',5);

ylim([10^(-8) 10^2]);
xlim([0 1000]);
xlabel('Frequency [Hz]','FontSize',30)
ylabel('Voltage/Acceleration [V/(m/s^2)]','FontSize',30)
set(gca,'YScale', 'log','FontSize',30);
set(gcf, 'Position', get(0, 'Screensize'));
pbaspect([9 5 1])

% semilogy(omj2(locs60)/(2*pi),pks60,'bx',omj2(locs70)/(2*pi),pks70,'rx',omj2(locs80)/(2*pi),pks80,'gx',omj2(locs90)/(2*pi),pks90,'mx',omj2(locs100)/(2*pi),pks100,'kx');

ah1 = axes('position',get(gca,'position'),'visible','off');

legend1 = legend(p1,'R = 10^2 Ohm','R = 10^3 Ohm','R = 10^4 Ohm','R = 10^5 Ohm','R = 10^6 Ohm','FontSize',25,'Box','off');
legend2 = legend(ah1,p2(1:2),'Present theory','Erturk','Position',[0.825 0.75 0.05 0.07],'FontSize',25,'Box','off');

% legend('R = 10^2 Ohm','R = 10^3 Ohm','R = 10^4 Ohm','R = 10^5 Ohm','R = 10^6 Ohm');
end


