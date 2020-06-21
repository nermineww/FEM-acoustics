% -------------------------------------------------------------------------
% Universidade Federal de Santa Catarina
% Programa de Pós-Graduação em Engenharia Mecânica
% Laboratório de Vibrações e Acústica - LVA
% Disciplina: Vibroacústica Computacional - EMC410104
% Professor: Dr. Júlio A. Cordioli
% Mestrando: Augusto R. C. de Sousa, Eng
% Matrícula: 201604230
% -------------------------------------------------------------------------
% Este código utiliza o Método dos Elementos Finitos para a avaliação do
% sistema proposto no Trabalho 3 da disciplina.
% -------------------------------------------------------------------------
% Início do código
% -------------------------------------------------------------------------

%% Cavidade
%%
clear all; clc;

Lx = 0.24;
Ly = 0.24;
Lz = 0.54;
d = 0.02;

rhoac = 1.1; %% Densidade aproximada do ar T = 20°C
co = 343;   %% Velocidade do som no ar T = 20°C, p_o = 101,325 kPa

a = d/2;

Nx = round(Lx/d); %Número de elementos em x
Ny = round(Ly/d); %Número de elementos em y
Nz = round(Lz/d); %Número de elementos em z
Nnx = Nx+1; %Número de nós em x
Nny = Ny+1; %Número de nós em y
Nnz = Nz+1; %Número de nós em z
Ntotal = Nx*Ny*Nz; %Número total de elementos
Nn = (Nx+1)*(Ny+1)*(Nz+1); %Número total de nós
N = 8; %Número de nós por elemento

%% Matrizes de inércia e rigidez acústica do elemento
syms xi1 xi2 xi3

Be = [1 -1 -1 -1 1 1 1 -1;
      1 1 -1 -1 -1 -1 1 1;
      1 1 1 -1 1 -1 -1 -1;
      1 -1 1 -1 -1 1 -1 1;
      1 -1 -1 1 1 -1 -1 1;
      1 1 -1 1 -1 1 -1 1;
      1 1 1 1 1 1 1 1;
      1 -1 1 1 -1 -1 1 -1;];
  
g = [1; xi1; xi2; xi3; xi1*xi2; xi1*xi3; xi2*xi3; xi1*xi2*xi3];
g1 = diff(g,xi1);
g2 = diff(g,xi2);
g3 = diff(g,xi3);

He = eval(inv(Be')*int(int(int((((1/a^2)*g1*g1') + ((1/a^2)*g2*g2') + ((1/a^2)*g3*g3'))*(a^3),xi1,-1,1),xi2,-1,1),xi3,-1,1)*inv(Be));
Qe = eval((1/co^2)*inv(Be')*int(int(int(g*g'*(a^3),xi1,-1,1),xi2,-1,1),xi3,-1,1)*inv(Be));

%% Matriz de conectividade da cavidade para o problema proposto
MconAc = zeros(Nx*Ny,N);
for k = 1:Nz
    for j = 1:Ny
    for i = 1:Nx
        for n = 1:N
            aux1 = (k-1)*(Nx*Ny);
            aux2 = (k-1)*(Nx+1)*(Ny+1);
            if n==1
                MconAc (i + Nx*(j-1)+aux1,n) = ((i+Ny*(j-1)) + (j-1)) + aux2;
            else
                if n==2
                    MconAc (i + Nx*(j-1)+aux1,n) = ((i+Ny*(j-1)) + (j-1)+1) + aux2;
                else
                    if n==3
                        MconAc (i + Nx*(j-1)+aux1,n) = ((i+Ny*(j-1)) + (j-1)+(Nx+1)) + aux2;
                    else
                        if n==4
                            MconAc (i + Nx*(j-1)+aux1,n) = ((i+Ny*(j-1)) + (j-1)+(Nx+2)) + aux2;
                        else
                            if n==5
                                MconAc (i + Nx*(j-1)+aux1,n) = ((i+Ny*(j-1)) + Nx*Nx+Nx+Ny+1+(j-1)) + aux2;
                            else
                                if n==6
                                    MconAc (i + Nx*(j-1)+aux1,n) = ((i+Ny*(j-1)) + Nx*Nx+Nx+Ny+1+(j-1)+1) + aux2;
                                else
                                    if n==7
                                        MconAc (i + Nx*(j-1)+aux1,n) = ((i+Ny*(j-1)) +  Nx*Nx+Nx+Ny+1+(j-1)+(Nx+1)) + aux2;
                                    else
                                        if n==8
                                            MconAc (i + Nx*(j-1)+aux1,n) = ((i+Ny*(j-1)) + Nx*Nx+Nx+Ny+1+(j-1)+(Nx+2)) + aux2;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    end
end

Ngl = 1; %Número de graus de liberdade por nó = 1
for i=1:N
   for j=1:N  
       Qe_s(j+(i-1)*N).sub = Qe((1+(j-1)*Ngl):(Ngl+(j-1)*Ngl),(1+(i-1)*Ngl):(Ngl+(i-1)*Ngl));
       He_s(j+(i-1)*N).sub = He((1+(j-1)*Ngl):(Ngl+(j-1)*Ngl),(1+(i-1)*Ngl):(Ngl+(i-1)*Ngl));  
   end
end

%% Pre-alocação das matrizes de massa e rigidez do sistema
Q = zeros(Nn);
H = zeros(Nn); 
%% Montagem das Matrizes Globais
for i=1:Ntotal
    for k=1:N
        for j=1:N
            Q(((1+(MconAc(i,k)-1)):(Ngl+(MconAc(i,k)-1))),((1+(MconAc(i,j)-1)):(Ngl+(MconAc(i,j)-1)))) = Qe_s((k-1)+1+(j-1)*N).sub + Q(((1+(MconAc(i,k)-1)*Ngl):(Ngl+(MconAc(i,k)-1)*Ngl)),((1+(MconAc(i,j)-1)*Ngl):(Ngl+(MconAc(i,j)-1)*Ngl)));
            H(((1+(MconAc(i,k)-1)):(Ngl+(MconAc(i,k)-1))),((1+(MconAc(i,j)-1)):(Ngl+(MconAc(i,j)-1)))) = He_s((k-1)+1+(j-1)*N).sub + H(((1+(MconAc(i,k)-1)*Ngl):(Ngl+(MconAc(i,k)-1)*Ngl)),((1+(MconAc(i,j)-1)*Ngl):(Ngl+(MconAc(i,j)-1)*Ngl)));
        end
    end
end
%% Extraindo autovalores e autovetores
Qsparse = sparse(Q);
Hsparse = sparse(H);

Modos = 40;
nn = 1;
    [A,Wn_Cavidade]=eigs(Hsparse,Qsparse,Modos,'SM'); %Extraindo os Autovetores (A) e Autovalores (Wn^2) de M e K
    resultados.Wn_Cavidade = diag(Wn_Cavidade).^(0.5)/(2*pi); % freq. naturais;
    [Fn_Cavidade.fn, ind.fn] = sort(resultados.Wn_Cavidade);

%% Criando a malha

xp1 = 0:d:Lx;
yp1 = 0:d:Ly;
zp1 = 0:d:Lz;

[X1,Y1,Z1] = meshgrid(xp1,yp1,zp1);

% Cálculo da pressão
Pref = zeros(1,Modos);
for m = 1:Modos
    N=ind.fn(m,1);
    Pressure = zeros(Nny,Nnx,Nnz);
    Pref(m) = max(max(max(abs(A(:,N)))));
    for k=1:Nnz
        for j=1:Nny   
            for i=1:Nnx
                    Pressure(i,j,k) =  A(i+(j-1)*Nnx+(k-1)*(Nnx*Nny),N)/Pref(m);
            end
        end
    end
    P(m).p = Pressure;
end

n=1;
Nm = 5;

%% Figuras da cavidade

% a) Primeiras 5 frequências naturais e modos da cavidade
figure(1)
for i = 1:1:(Nm-1)
subplot(2,2,i)
slice(X1, Y1, Z1, P(i).p, 0, 0, 0.54);
shading interp
set(gca, 'FontSize', 20)
xlabel('Lx (m)')
ylabel('Ly (m)')
zlabel('Lz (m)')
title(['Forma modal ' num2str(i) ' da cavidade, f_n = ' num2str(Fn_Cavidade.fn(i)) ' Hz']);
axis equal
colorbar
end

figure(2)
slice(X1, Y1, Z1, P(5).p, 0, 0, 0.54);
shading interp
set(gca, 'FontSize', 20)
xlabel('Lx (m)')
ylabel('Ly (m)')
zlabel('Lz (m)')
title(['Forma modal ' num2str(5) ' da cavidade, f_n = ' num2str(Fn_Cavidade.fn(5)) ' Hz']);
axis equal
colorbar

%% Placa
%% Propriedades da porta
rhostr = 7850; %Densidade do material
E = 205E9; %Módulo de elasticidade do material
poisson = 0.28; %Coeficiente de Poisson
EP = E/(1-poisson^2); 
G = E/(2*(1+poisson)); %Módulo de cisalhamento (ou de rigidez) do material
h = 0.002; %Espessura da porta
I = (h^3)/12; %Momento de inércia de área (segundo)
Lx = 0.2; %Largura em x
Ly = 0.5; %Largura em y

%% Propriedades da malha e dos elemento
dx = 0.02; %Tamanho do elemento na direção x
dy = 0.02; %Tamanho do elemento na direção y
a = dx/2;
b = dy/2;
Nx = round((Lx)/dx); %Elementos em x
Ny = round((Ly)/dy); %Elementos em y
Ne = Nx*Ny; %Número total de elementos
Nn_x = Nx+1; %Nós em x
Nn_y = Ny+1; %Nós em y
Nnt = Nn_x*Nn_y; %Número total de nós
Nn_e = 4;%Nós por elemento
Ngl_n = 3;%Graus de Liberdade por nó
Ntgls = Nnt*Ngl_n;%Total de Graus de Liberdade do sistema
%% Matriz Ae do elemento, de acordo com Petyt (2010)
[Ae] = [1 -1 -1 1 1 1 -1 -1 -1 -1 1 1;
        0 0 1/b 0 -1/b -2/b 0 1/b 2/b 3/b -1/b -3/b;
        0 -1/a 0 2/a 1/a 0 -3/a -2/a -1/a 0 3/a 1/a;
        1 1 -1 1 -1 1 1 -1 1 -1 -1 -1;
        0 0 1/b 0 1/b -2/b 0 1/b -2/b 3/b 1/b 3/b;
        0 -1/a 0 -2/a 1/a 0 -3/a 2/a -1/a 0 3/a 1/a;
        1 1 1 1 1 1 1 1 1 1 1 1;
        0 0 1/b 0 1/b 2/b 0 1/b 2/b 3/b 1/b 3/b;
        0 -1/a 0 -2/a -1/a 0 -3/a -2/a -1/a 0 -3/a -1/a;
        1 -1 1 1 -1 1 -1 1 -1 1 -1 -1;
        0 0 1/b 0 -1/b 2/b 0 1/b -2/b 3/b -1/b -3/b;
        0 -1/a 0 2/a -1/a 0 -3/a 2/a -1/a 0 -3/a -1/a];

%% Matriz constitutiva do material , [D]
D = [EP EP*poisson 0; EP*poisson EP 0; 0 0 G];
%% Polinômio p(csi,eta) e vetor {X}
syms csi eta
p = [1; csi; eta; csi^2; csi*eta; eta^2; csi^3; csi^2*eta; csi*eta^2; eta^3; csi^3*eta; csi*eta^3];
X(1,1:12) = (1/(a^2))*[0 0 0 2 0 0 6*csi 2*eta 0 0 6*csi*eta 0];
X(2,1:12) = (1/(b^2))*[0 0 0 0 0 2 0 0 2*csi 6*eta 0 6*csi*eta];
X(3,1:12) = (2/(a*b))*[0 0 0 0 1 0 0 2*csi 2*eta 0 3*(csi^2) 3*(eta^2)];

%% Matrizes elementares de Massa e Rigidez e vetor de forças unitário
Me = eval(inv(Ae')*(int(int(rhostr*h*p*p'*a*b,csi,-1,1),eta,-1,1))*inv(Ae));
Ke = eval(inv(Ae')*(int(int(I*X'*D*X*a*b,csi,-1,1),eta,-1,1))*inv(Ae));
fe  = eval(inv(Ae')*(int(int(p*a*b,csi,-1,1),eta,-1,1)));
%% Matriz de conectividade
MConS = zeros(Ne,Nn_e+1);
aux = 0;
for i=1:Ne
    if i<=Ne
    if i>Nx && rem(i-1,Nx)==0 && i<=Ne
    aux = aux + 1;
    end
        for j=1:Nn_e
            MConS(i,1) = i;
                if j<=2
                    MConS(i,j+1) = j+(i-1)+aux;
                else
                    if j==3
                        MConS(i,j+1) = j+(i-1)+aux+Nx;
                    else
                        MConS(i,j+1) = (j-1) + (Nx-1) + (i-1) +  aux;
                    end
                end 
        end    
    end
end
%% Determinação das matrizes de massa e rigidez para cada nó

for i=1:Nn_e
   for j=1:Nn_e
    
      Me_s(j+(i-1)*Nn_e).sub = Me((1+(j-1)*Ngl_n):(3+(j-1)*Ngl_n),(1+(i-1)*Ngl_n):(3+(i-1)*Ngl_n));
      Ke_s(j+(i-1)*Nn_e).sub = Ke((1+(j-1)*Ngl_n):(3+(j-1)*Ngl_n),(1+(i-1)*Ngl_n):(3+(i-1)*Ngl_n));
   
   end
end

%% Matrizes globais
Kglobal = zeros(Ntgls);
Mglobal = zeros(Ntgls);
Mconec = MConS(:,2:(Nn_e+1));

for i=1:Ne
    for k=1:Nn_e
        for j=1:Nn_e
            Mglobal(((1+(Mconec(i,k)-1)*3):(3+(Mconec(i,k)-1)*3)),((1+(Mconec(i,j)-1)*3):(3+(Mconec(i,j)-1)*3))) = Me_s((k-1)+1+(j-1)*Nn_e).sub + Mglobal(((1+(Mconec(i,k)-1)*3):(3+(Mconec(i,k)-1)*3)),((1+(Mconec(i,j)-1)*3):(3+(Mconec(i,j)-1)*3))); 
            Kglobal(((1+(Mconec(i,k)-1)*3):(3+(Mconec(i,k)-1)*3)),((1+(Mconec(i,j)-1)*3):(3+(Mconec(i,j)-1)*3))) = Ke_s((k-1)+1+(j-1)*Nn_e).sub +  Kglobal(((1+(Mconec(i,k)-1)*3):(3+(Mconec(i,k)-1)*3)),((1+(Mconec(i,j)-1)*3):(3+(Mconec(i,j)-1)*3)));    
        end
    end
end

%% Graus de liberdade referente às dobradiças e fechadura

for j = 1:1:Ny 
        FirstEy(j) = 1 + Nx*(j-1);
        SecondEy(j) = FirstEy(j)+1;
        LastEy(j) = Nx*(j);
        BeforeLastEy(j) = LastEy(j)-1;
end
% Calcula quais linhas e elementos está aplicado o primiero engaste -
% dobradiça inferior
L1 = round(0.1/dy + 1); E1 = FirstEy(L1); E11 = SecondEy(L1);
L2 = round((0.1+dy)/dy); E2 = FirstEy(L2); E22 = SecondEy(L2);

% Calcula quais linhas e elementos está aplicado o segundo engaste -
% dobradiça superior
L3 = round((0.5-0.1)/dy); E3 = FirstEy(L3); E33 = SecondEy(L3);
L4 = round((0.5-(0.1+dy))/dy); E4 = FirstEy(L4); E44 = SecondEy(L4);

%Calcula quais as linhas e elementos está aplicado o engaste da fechadura
L5 = round((0.5-0.22)/dy); E5 = LastEy(L5); E55 = BeforeLastEy(L5);
L6 = round((0.5-(0.22+dy))/dy); E6 = LastEy(L6); E66 = BeforeLastEy(L6);
L7 = round((0.5-(0.22+(2*dy)))/dy); E7 = LastEy(L7); E77 = BeforeLastEy(L7);
L8 = round((0.5-(0.22+(3*dy)))/dy); E8 = LastEy(L8); E88 = BeforeLastEy(L8);

% Nós que devem ser retirados devido às dobradiças
N1 = Mconec(E1,:); N11 = Mconec(E11,:);            %Nós do elemento 1...
N2 = Mconec(E2,:); N22 = Mconec(E22,:);
N3 = Mconec(E3,:); N33 = Mconec(E33,:);
N4 = Mconec(E4,:); N44 = Mconec(E44,:);
N5 = Mconec(E5,:); N55 = Mconec(E55,:);
N6 = Mconec(E6,:); N66 = Mconec(E66,:);
N7 = Mconec(E7,:); N77 = Mconec(E77,:);
N8 = Mconec(E8,:); N88 = Mconec(E88,:);

Nos = [N1 N3 N5 N6];

Nos = unique(Nos);
% Os GL a serem excluidos
GL = zeros(1,length(Nos)*Ngl_n); % Criando o vetor de tamanhao necessário para alocar GL's
for j=1:1:length(Nos)
   GL(1,j*3) = Nos(j)*3;
   GL(1,(j*3)-1) = Nos(j)*3 - 1;
   GL(1,(j*3)-2) = Nos(j)*3 - 2;
end

    for j=1:length(GL)
        Mglobal(:,GL(j)-j+1) = []; %Tirando a coluna do nó restrito da matriz M
        Mglobal(GL(j)-j+1,:) = []; %Tirando a linha do nó restrito da Matriz M
        Kglobal(:,GL(j)-j+1) = []; %Tirando a coluna do nó restrito da matriz K
        Kglobal(GL(j)-j+1,:) = []; %Tirando a linha do nó restrito da Matriz K
    end
    
%% Autovalor e Autovetor
Modos = [10 20 30 40];
Nfn = 20;
M = sparse(Mglobal);
K = sparse(Kglobal);
[Avet, Wn_Placa] = eigs(K,M,Modos(4),'SM');
resultados.Wn_Placa = diag(Wn_Placa).^(0.5)/(2*pi); % freq. naturais;
[Fn_Placa.fn, ind.fn] = sort(resultados.Wn_Placa);
%% Incluindo os GL das condicoes de contorno
   for j=1:length(GL)
        [s1,s2]=size(Avet);
        nline=zeros(1,s2);
        V1=Avet(1:(GL(j)-1),:);
        V2=Avet(GL(j):s1,:);
        Avet=cat(1,V1,nline,V2);
   end
    
%% Excluindo Rotações
  Rot=Avet(1:3:end,:);
  
%% Criando a malha
Lxp = 0:dx:Lx;  
Lyp = 0:dy:Ly;   
       
[X,Y] = meshgrid(Lxp,Lyp);

%% Selecionando as formas modais a serem plotadas

for k=1:Nfn
modo=Rot/max(max(Rot));
        if sum(modo(1:10,k))>=0
            modo(:,k)=modo(:,k);
        else
            modo(:,k)=-modo(:,k);
        end
eval(['f_m',num2str(k),'=vec2mat(modo(:,',num2str(ind.fn(k,1)),'),Nn_x);'])
end
%% Cálculo e representação das 5 primeiras frequências naturais da placa

%Modos 1 a 4
figure(3)
for n = 1:1:length(Modos);
subplot(2,2,n)
set(gcf,'papersize',[20 10])
surf(X,Y,eval(['f_m',num2str(n)]))
grid on;
box off; 
colorbar;
set(gca,'fontsize',20)
xlabel('Lx (m)');
ylabel('Ly (m)');
zlabel('Amplitude');
title(['Forma modal ' num2str(n) ' da placa, f_n =' num2str(Fn_Placa.fn(n)) ' Hz']);
end

%Modo 5
figure(4)
set(gcf,'papersize',[20 10])
surf(X,Y,eval(['f_m',num2str(5)]))
grid on;
box off; 
colorbar;
set(gca,'fontsize',20)
xlabel('Lx (m)');
ylabel('Ly (m)');
zlabel('Amplitude');
title(['Forma modal ' num2str(5) ' da placa, f_n = ' num2str(Fn_Placa.fn(5)) ' Hz']);

%% Acoplamento
syms xi1 xi2 xi3 csi eta
xi2 = -1;
xi1 = csi;
xi3 = eta;
g = subs(g);
Se = eval(inv(Ae')*int(int(p*g'*a*b,csi,-1,1),eta,-1,1)*inv(Be));
Re = Se';
%% Determinação das submatrizes para cada nó

NAc = 1; % Número de graus de liberdade acústicos
NSt = 3; % Número de graus de liberdade da estrutura
NnAc = 8; % Número de nós por elemento da cavidade acústica
NnS = 4; % Número de nós por elemento da placa
NnSt = 858; % Numero de graus de liberdade estrutural
Nn = Nnx*Nny*Nnz; % Numero de graus de liberdade da cavidade

for i=1:NnAc
   for j=1:NnS
       Se_s(j+(i-1)*NnS).sub = Se((1+(j-1)*NSt):(NSt+(j-1)*NSt),(1+(i-1)*NAc):(NAc+(i-1)*NAc));
       Re_s(j+(i-1)*NnS).sub = Re((1+(i-1)*NAc):(NAc+(i-1)*NAc),(1+(j-1)*NSt):(NSt+(j-1)*NSt));  
   end
end

%% Montagem das matrizes globais
% Determinando os elementos da cavidade novamente
LxAc = 0.24;
LyAc = 0.24;
d = 0.02;
NxAc = LxAc/d; %Número de elementos em x
NyAc = LyAc/d; %Número de elementos em y
Id_lS = MConS;
Id_lAc = MconAc;
S_parc = zeros(NnSt,Nn);

for i=1:NxAc*NyAc
    for k=1:NnS
        for j=1:NnAc
            S_parc(((1+(Id_lS(i,k)-1)*NSt):(NSt+(Id_lS(i,k)-1)*NSt)),((1+(Id_lAc(i,j)-1)*NAc):(NAc+(Id_lAc(i,j)-1)*NAc))) = Se_s((k-1)+1+(j-1)*NnS).sub + S_parc(((1+(Id_lS(i,k)-1)*NSt):(NSt+(Id_lS(i,k)-1)*NSt)),((1+(Id_lAc(i,j)-1)*NAc):(NAc+(Id_lAc(i,j)-1)*NAc)));       
        end
    end
end

 R_parc = S_parc'; % Matrizes R e S
 
for j=1:length(GL)
        R_parc(:,GL(j)-j+1) = [];
end 

R = R_parc;
S = R';

O = zeros(size(R'));

M1 = [M O; -R Q];
M2 = [K S; O' H];

Nfn = 100; %Numero de fn a serem calculadas
M1sparse=sparse(M1);
M2sparse=sparse(M2);
    [A_Acoplado,Wn_Acoplado]=eigs(M2sparse,M1sparse,Nfn,'sm'); %Extraindo os Autovetores (A) e Autovalores (Wn^2) de M e K
    resultados.Wn_Acoplado = diag(Wn_Acoplado).^(0.5)/(2*pi); % freq. naturais;
    [Fn_Acoplado.fn, ind.fn_Acoplado] = sort(resultados.Wn_Acoplado);

%% Calculo do NPS dentro da cavidade e definindo vetor de força e posição do ponto para calculo de NPS
px=0.16; py=0; pz=0.42; Nnx_st = 11; Ngln=3; d=0.02; %coordenadas da aplicação da força
no_st=(Nnx_st*(pz/0.02) + (px/0.02)+1)*Ngln-4; %Calculo do nó onde foi aplicada a força
forca=zeros(length(M),1);
forca(no_st,1)= 1; %Força unitária no ponto A

px_ac=0.1; py_ac=0.16; pz_ac=0.34; %coordenadas da coleta de pressão
no_ac = (px_ac/d + 1) + (NxAc+1)*(py_ac/d) + (pz_ac/d)*((NxAc+1)*(NyAc+1)); %Nó da coleta de pressão

freq = 100:2:700;
omega = 2*pi*freq;
eta = 0.03;

Lp = zeros(length(Q),length(freq));
for i=1:length(freq)
    w = ((-omega(i)^2*M) + (1+1i*eta)*K)\forca;
    Lp(:,i) = (-omega(i)^2*w'*S)/(-omega(i)^2*Q + H);
end

figure(5)
semilogx(freq,dlmread('COMSOL2.txt'));
hold on
semilogx(freq,20*log10(abs(Lp(no_ac,:))/(20*10^(-6))))
set(gca,'fontsize',25)
legend('NPS - MATLAB', 'NPS - COMSOL');
ylabel('NPS (dB)');
xlabel('Frequência (Hz)');
xlim([100,700])
grid on

figure(6)
semilogx(freq,dlmread('COMSOL.txt'));
hold on
semilogx(freq,dlmread('COMSOL2.txt'));
hold on
semilogx(freq,20*log10(abs(Lp(no_ac,:))/(20*10^(-6))))
set(gca,'fontsize',25)
legend('NPS - COMSOL (sem abertura)', 'NPS - COMSOL (com abertura)','NPS - MATLAB (sem abertura)');
ylabel('NPS (dB)');
xlabel('Frequência (Hz)');
xlim([100,700])
grid on

clearvars -except Q H M K Fn_Cavidade Fn_Placa Fn_Acoplado