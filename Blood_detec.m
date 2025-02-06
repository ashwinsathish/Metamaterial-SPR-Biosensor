%%%%%%%%%% for tm binary grating glass metal  metalgrating dielectric
%%%%%%%%%% wavelength modulation
clear;
close all;
clc;
clf;
hold all;

c=3e8;

count=1;
w=66:0.008:75;

iii_range = w;
total_iterations = length(iii_range);
td = 2500e-9; % Array of different td values
IR12_nd_133 = [];
nd_range = [1.368, 1.376, 1.381, 1.385, 1.387, 1.392, 1.390, 1.395, 1.399, 1.401];

figure
title('SPR curves');
hold all;

for nd=nd_range
    for iii=iii_range
        lambda0=1.550e-6;
        lamdac = 2.4511e-5;
        lamdap = 1.0657e-7;
        epsilonreal = 1-((lambda0.^2.*(lamdac)^2)./(lamdap^2.*(lambda0.^2+lamdac^2)));
        epsilonim = ((lambda0.^3.*(lamdac))./(lamdap^2.*(lambda0.^2+lamdac^2)));

        n2 = sqrt((sqrt(epsilonreal.^2+epsilonim.^2)+epsilonreal)./2);
        k2 = sqrt((sqrt(epsilonreal.^2+epsilonim.^2)-epsilonreal)./2);

        Numords = 101;  % number of diffractive orders maintained
        nc = 1.426;   % CaF2 prism
        ns = 1;  % region 3 substrate refractive index
        Ngrat = 8;   % number of grating slices
        period = 100e-9;  % grating period in microns

        if (nd == 1.368)
            detec_text = 'Healthy Hela';
        elseif (nd == 1.376)
            detec_text = 'Healthy Jurkat';
        elseif (nd == 1.381)
            detec_text = 'Healthy PC12';
        elseif (nd == 1.385)
            detec_text = 'Healthy MDA-MB-231';
        elseif (nd == 1.387)
            detec_text = 'Healthy MCF-7';  
        elseif (nd == 1.392)
            detec_text = 'Cervical cancer (Hela)';
        elseif (nd == 1.390)
            detec_text = 'Blood cancer (Jurkat)';
        elseif (nd == 1.395)
            detec_text = 'Adrenal cancer (PC12)';
        elseif (nd == 1.399)
            detec_text = 'Breast cancer (MDA-MB-231)';
        elseif (nd == 1.401)
            detec_text = 'Breast cancer (MCF-7)';
        end

        %nd = 1.33;
        nm = n2-1i*k2; % Aluminium refractive index
        A = 3.44904;
        A1 = 2271.88813;
        A2 = 3.39538;
        t1 = 0.058304;
        t2 = 0.30384;

        neth = 1.36;
        nALO = 1.746;
        teth = 0e-9;
        %tALO = kkk*5e-9;
        mos2_layers = 0;
        j = sqrt(-1);
        c1 = 5.446e6;

        nsi = A + A1*exp(-1.55/t1) + A2*exp(-1.55/t2);
        nx = 1; % dummy variable
        nmos2 = 3.647;
        ngr =  3.0 - 1.1491i;
        nTi = 2.7043 - 3.7657i;
        nAg = 0.0562 - 4.2776i;
        nBFO = 2.968;
        nBTO = 2.30;
        %nmeta = -sqrt((-4.0 + 1i*0.001)*(-2.4 + 1i*0.001));
        nmeta = -3.0984 - 0.001i;
        ntio2 = 2.4328;

        nr = [nm, nBTO, nm, nmeta, ntio2, nmos2, nx, nd];  % Ridge 
        depth = [30e-9, 5e-9, 10e-9, 160e-9, 2e-9, 0.65e-9, 0e-9,td];  % structure 7

        ng = nr; % Groove
        Filfac = [.5 .5 .5 .5 .5 .5 .5 .5 ];  % fill factor for ridges
        Disp = [0 0 0 0 0 0 0 0];   % ridge displacement in a frac                                                                                                                                                                                                                                              tion of period

        theta0=iii;  %   angle of incidence
        phi0=0;    % azimuthal angle of incidence
        deg=pi/180; 
        Nmax=(Numords-1)/2; % highest order number retained
        I=(-Nmax:Nmax)';  % I is the order index
        p=(Numords+1)/2;  % index of zeroth order
        theta0=theta0*deg;
        phi0=phi0*deg;  % converting in radians
        epsc=nc^2;    % relative permittivity
        epss=ns^2;

        k0=2*pi/lambda0;  % free space vector
        K=2*pi/period;    % grating vector

        kc=k0*nc;
        kx=kc*sin(theta0)*cos(phi0)-I*K;  % region1 wave vector components

        ky=kc*sin(theta0)*sin(phi0)*ones(size(kx));
        kzc=sqrt(kc^2-kx.^2-ky.^2);
        bad_indices=find((real(kzc)-imag(kzc))<0);  
        kzc(bad_indices)=-kzc(bad_indices);              
        ks=k0*ns;
        kzs=sqrt(ks^2-kx.^2-ky.^2);
        bad_indices=find((real(kzs)-imag(kzs))<0);  % region3 wavevector
        kzs(bad_indices)=-kzs(bad_indices); 

        %%%%%% define some  matrices and vectors %%%%%%%
        Zv=zeros(Numords,1);
        Dv=Zv;
        Dv(p)=1;
        Zv2=[Zv;Zv];
        Eye=eye(Numords);  % identity matrix
        Kx=diag(kx/k0);
        Kxsq=Kx.^2;
        Kzc=diag(kzc/k0);
        Kzcsq= Kzc.^2;
        Kzs=diag(kzs/k0);
        Kzssq= Kzs.^2;
        M=Numords-1;
        temp1=Eye/ns;
        fmat=Eye;
        gmat=j*Kzs/ns^2;

        for ii=Ngrat:-1:1
            epsg=ng(ii).^2;  % groove permittivity
            epsr=nr(ii).^2;   % ridge permittivity
            epsG=(1-Filfac(ii))*epsg+Filfac(ii)*epsr;  % average grating
            iepsG=(1-Filfac(ii))/epsg+Filfac(ii)/epsr; 
            Sinc=sin(pi*Filfac(ii)*(1:M))./(pi*(1:M));
            vm=(epsr-epsg)*fliplr(Sinc);
            v0=epsG;
            vp=(epsr-epsg)*Sinc;
            v=[vm v0 vp].*exp(+j*2*pi*Disp(ii)*(-M:M));
            ivm=(1/epsr-1/epsg)*fliplr(Sinc);
            iv0=iepsG;
            ivp=(1/epsr-1/epsg)*Sinc;
            iv=[ivm iv0 ivp].*exp(+j*2*pi*Disp(ii)*(-M:M));
            Epsilon=toeplitz(fliplr(v(1:Numords)),v(Numords:2*Numords-1));
            Alpha=toeplitz(fliplr(iv(1:Numords)),iv(Numords:2*Numords-1));

            clear Sinc v  vm  v0  vp 
            B=Kx*(Epsilon\Kx)-Eye;  % cofficient matrix
            [W,V]=eig(Alpha\B);    % W is the eigen vector and V are the eigen values
            Q=sqrt(V);
            M0=Alpha*W*Q;
            E=expm(-k0*Q*depth(ii));
            v=[W W;M0,-M0]\[fmat;gmat];
            temp2=v(1:Numords,:)\E;
            temp3=E*v(Numords+1:2*Numords,:)*temp2;
            temp1=temp1*temp2;
            fmat=W+W*temp3;
            gmat=M0-M0*temp3;

        end

        gfi=gmat/fmat;
        RHS=-gfi(:,p);
        RHS(p)=RHS(p)+j*kzc(p)/k0/epsc;
        Rs=(gfi+j*Kzc/nc^2)\RHS;
        Ts=(temp1/fmat)*(Rs+Dv)*nc;
     

       IR1=(abs(Rs).^2).*real(kzc./kzc(p));
       IT1=(abs(Ts).^2).*real(kzs./kzc(p));

       e=sum(IT1);
       f=sum(IR1);
       g=1-e-f;
       IT12(count)=e;
       IR12(count)=f; 
       loss(count)=g;

       count=count+1;
 
       progress = (count/total_iterations)*100;
       clc;
       fprintf('nd = %.3f\t', nd);
       fprintf(detec_text);
       fprintf('\t\t');
       fprintf('\nProgress: [%s%s] %.2f%%\r', repmat('=',1,floor(progress/2)), repmat(' ',1,50-floor(progress/2)), floor(progress));

    end

    %%%%% Plotting results %%%%%
  
    
    % Plotting the curve for current nd value
    plot(w, IR12, 'DisplayName', ['nd = ', num2str(nd), ' ', detec_text]);
    
    % Resetting count for next iteration
    count=1;
    
end    

% Setting plot title and labels
title('SPR curves');
xlabel('Angle');
ylabel('Reflectivity');
legend('Location', 'southwest');

% Resetting minVals and minWs for next td value
minVals = [];
minWs = [];
    
