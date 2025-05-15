% Periodically driven-dissipative SPINFULL free fermion
Lspin = 16; %system size (total number of spinful fermion sites)
L=2*Lspin; %system size in c_{j,up}=c_{2j-1}, c_{j,down} = c_{2j} representation.
eta = randi([0 1],L,1); %L*1 vector storing the initial occupation numbers.
if mod(sum(eta),2) ~= 0
   eta(1) = mod(eta(1)+1,2); % make sure that there are even number of particles
end
%eta = ones(L,1);
%eta(L/2+1:L) = 0;
phi = eye(L); % A set of basis wave functions stored in the COLUMNS.
Corr = zeros(L); % Correlation matrix at t=0
for j = 1:L
    for k = 1:L
        Corr(j,k) = (conj(phi(j,:)).*phi(k,:))*eta;
    end
end
C0 = zeros(2*L); %Covariance matrix at t=0;
for j = 1:L
    for k = 1:L
        C0(2*j,2*k) = 2i*imag(Corr(j,k));
        C0(2*j-1,2*k-1) = 2i*imag(Corr(j,k));
        C0(2*j-1,2*k) = 2i*real(Corr(j,k)) - 1i*(j-k==0);
        C0(2*j,2*k-1) = -2i*real(Corr(j,k)) + 1i*(j-k==0);
    end
end

% Free fermion Hamiltonian
J=1;
h=0; 

Hc = zeros(L);%Hamiltonian in canonical fermions
% alpha = 0.85;
for j = 1:L-2
    Hc(j,j+2) = J;
    Hc(j+2,j) = J';
    Hc(j,j) = h;
end
Hc(L-1,L-1) = h;
Hc(L,L) = h;

H0 = zeros(2*L);%Hamiltonian in majorana fermions
for j = 1:L
    for k = 1:L
        H0(2*j-1,2*k-1) = (Hc(j,k) - Hc(k,j))/8;
        H0(2*j,2*k-1) = 1i*(Hc(j,k) + Hc(k,j))/8;
        H0(2*j-1,2*k) = -1i*(Hc(j,k) + Hc(k,j))/8;
        H0(2*j,2*k) = (Hc(j,k) - Hc(k,j))/8;
    end
end

% add some pairing in the bulk
pair = 0.3+0.4i;
Hp = zeros(2*L);
for j = 1:Lspin
    Hp(4*j-3,4*j-1) = (pair-pair')/8;
    Hp(4*j-2,4*j-1) = 1i*(pair + pair')/8;
    Hp(4*j-3,4*j) = 1i*(pair + pair')/8;
    Hp(4*j-2,4*j) = -(pair - pair')/8;
end
Hp = Hp+Hp';
H0 = H0+Hp;

% Driving Hamiltonian 
Delta = 1; %drive amplitude
Delta1 = 0;
omega = 2; %frequency
Hd = zeros(2*L);
Hd(1,3) = (Delta-Delta')/8 + (Delta1-Delta1')/8;
Hd(2,3) = 1i*(Delta + Delta')/8 + 1i*(Delta1 + Delta1')/8;
Hd(1,4) = 1i*(Delta + Delta')/8 - 1i*(Delta1 + Delta1')/8;
Hd(2,4) = -(Delta - Delta')/8 + (Delta1 - Delta1')/8;
Hd = Hd + Hd';



%%
% Dissipator
l = zeros(8,2*L);

Gamma1 = 0.5; l(1,1) = Gamma1*1/2; l(1,2) = Gamma1*1i/2;%L1 = c1^\dag
Gamma2 = 1.5; l(2,1) = Gamma2*1/2; l(2,2) = -Gamma2*1i/2;%L2 = c1
Gamma3 = 0.5; l(3,3) = Gamma3*1/2; l(3,4) = Gamma3*1i/2;%L3 = c2^\dag
Gamma4 = 1.5; l(4,3) = Gamma4*1/2; l(4,4) = -Gamma4*1i/2;%L4 = c2
Gamma5 = 0.8; l(5,2*L-3) = Gamma5*1/2;  l(5,2*L-2) = Gamma5*1i/2;%L5 = c{L-1}^\dag
Gamma6 = 1.2; l(6,2*L-3) = Gamma6*1/2;  l(6,2*L-2) = -Gamma6*1i/2;%L6 = c{L-1}
Gamma7 = 0.8; l(7,2*L-1) = Gamma7*1/2;  l(7,2*L) = Gamma7*1i/2;%L7 = cL^\dag
Gamma8 = 1.2; l(8,2*L-1) = Gamma8*1/2;  l(8,2*L) = -Gamma8*1i/2;%L8 = cL


% Set up the dynamic equation of C
M = zeros(2*L);
for j = 1:8
    M = M+l(j,:).'*conj(l(j,:));
end

Mr = real(M);
Mi = imag(M);

% Set up the dephasing terms
gamma = 0.1;



tmax = 1000; % a long enough time for the system to relax
tst = 5*2*pi/omega;
tspan  = [0 tmax/2 linspace(tmax,tst+tmax,1000)];
C0 = reshape(C0,[],1);
[t, C] = ode45(@(t,C) myODE(t,C,H0,Mr,Mi,Hd,omega,L,gamma), tspan, C0);
C = reshape(C,[],2*L,2*L);


%%

%occupation number on a site vs. t
site = 7;
occ = zeros(length(t)-2,1);
for j = 3:length(t)
    occ(j-2) = 1/2 + 1i*C(j,2*site,2*site-1)/2;
end


Cavg = squeeze(mean(C(3:length(t),:,:),1));

Ccorravg = zeros(L);
for j = 1:L
    for k = 1:L
        Ccorravg(j,k) = (Cavg(2*j-1,2*k-1) + Cavg(2*j,2*k) + 2*(j-k==0) + 1i*Cavg(2*j,2*k-1) - 1i*Cavg(2*j-1,2*k))/4;
    end
end

% occlist1 = [];
% for site = 1:L
%     occlist1 = [occlist1;1/2 + 1i*Cavg(2*site,2*site-1)/2];
% end

%We have two occupations
occlist1 = zeros(L/2,1); %spin up occupations on sites
occlist2 = zeros(L/2,1); %spin down occupations on sites
for site = 1:L/2
    occlist1(site) = Ccorravg(2*site-1,2*site-1);
    occlist2(site) = Ccorravg(2*site,2*site);
end

%We have two separated currents
curr1 = zeros(L/2+1,1); %spin up currents
curr2 = zeros(L/2+1,1); %spin down currents
curr1(1) = 2*Gamma1^2*(1-Ccorravg(1,1)) - 2*Gamma2^2*Ccorravg(1,1);
curr2(1) = 2*Gamma3^2*(1-Ccorravg(2,2)) - 2*Gamma4^2*Ccorravg(2,2);
curr1(L/2+1) = -2*Gamma5^2*(1-Ccorravg(L-1,L-1)) + 2*Gamma6^2*Ccorravg(L-1,L-1);
curr2(L/2+1) = -2*Gamma7^2*(1-Ccorravg(L,L)) + 2*Gamma8^2*Ccorravg(L,L);
for site = 1:L/2-1
    curr1(site+1) = 1i*(J*Ccorravg(2*site-1,2*site+1) - J'*Ccorravg(2*site+1,2*site-1));
    curr2(site+1) = 1i*(J*Ccorravg(2*site,2*site+2) - J'*Ccorravg(2*site+2,2*site));
end
currs = curr1 - curr2;

% currentlist = zeros(L+1,2); 
% currentlist(:,1) = 0:L;
% currentlist(1,2) = 2*Gamma1^2*(1-Ccorravg(1,1)); %current goes into site 1 from left bath, seems to off by a factor of 2
% for site = 1:L-1
%     currentlist(site+1,2) = 1i*J*(Ccorravg(site,site+1) - Ccorravg(site+1,site));
% end
% currentlist(L+1,2) = 2*Gamma4^2*Ccorravg(L,L);  %current from site L to right bath, seems to off by a factor of 2


figure
plot(t(3:length(t)),occ)
title("<n_j>(t)")
xlabel("t")
ylabel("<n_j>")

figure
tiledlayout(2,1)
title('Occupation number profile (time averaged)')

nexttile
plot(1:L/2,occlist1,'-o')
%plot(2:L/2-1,occlist1(2:L/2-1),'-o')
%ylim([0.4 0.5])
xlabel("j")
ylabel("<n_j>")
title('Occupation number profile (time averaged) - spin up')

nexttile
plot(1:L/2,occlist2,'-o')
%plot(2:L/2-1,occlist2(2:L/2-1),'-o')
%ylim([0.4 0.5])
xlabel("j")
ylabel("<n_j>")
title('Occupation number profile (time averaged) - spin down')

% slope_calc = gamma*Gamma1^2/(Gamma1^4+1+(L/2-1)*gamma*Gamma1^2);
% title(['Occupation number profile (time averaged). Calculated slope = ',num2str(slope_calc)])


figure
tiledlayout(3,1)

nexttile
plot(0:L/2,curr1,'-o')
xlabel("j")
ylabel("<J_j>")
title('Particle current profile (time averaged) - spin up')


nexttile
plot(0:L/2,curr2,'-o')
xlabel("j")
ylabel("<J_j>")
title('Particle current profile (time averaged) - spin down')

nexttile
plot(0:L/2,currs,'-o')
xlabel("j")
ylabel("<J_j>")
title('Spin current profile (time averaged)')



% current_calc = Gamma1^2/(Gamma1^4+1+(L/2-1)*gamma*Gamma1^2);
% title(['Particle current profile (time averaged). Calculated current = ',num2str(current_calc)])


% figure
% plot(currentlist(:,1),currentlist(:,2),'-o')
% current_calc = Gamma1^2/(Gamma1^4+1+(L-1)*gamma*Gamma1^2);
% title(['Particle current profile (time averaged). Calculated current = ',num2str(current_calc)])
% xlabel("j")
% ylabel("<J_j>")



%plot the time-averaged correlation matrix
% xx=kron(0:L-1,[1,0])+kron(1:L,[0,1]);
% [X,Y]=meshgrid(xx,xx);
% Z=kron(log(abs(Ccorravg)),[1,1;1,1]);% log scale
% % Zm=max(max(abs(Z))); Z=Z/Zm;
% figure;
% surf(X,Y,Z);
% view(90,90);
% xlim([0,L]);ylim([0,L]);
% xticks([0.5,L-0.5]);yticks([0.5,L-0.5]);xticklabels({'1',L});yticklabels({'1',L});
% caxis([-30,0]);colormap('jet');colorbar;

%%

function dC = myODE(t,Ct,H0,Mr,Mi,Hd,omega,L,gamma)
    Ct = reshape(Ct,2*L,2*L);
    Ht = H0+Hd*cos(omega*t);
    %Ht = H0+Hd*(0*cos(omega*t)+1); %with a dc offset
    X = 4*(1i*Ht+Mr);
    Y = 8*Mi;
    alt = zeros(2*L-1,1);
    alt(1:2:2*L-1) = 1;
    dephasing = diag(alt,-1) + diag(alt,1) - ones(2*L);
    dC = -X*Ct - Ct*X.' - 1i*Y + 2*Ct.*dephasing*gamma;
    dC = reshape(dC,[],1);
end
