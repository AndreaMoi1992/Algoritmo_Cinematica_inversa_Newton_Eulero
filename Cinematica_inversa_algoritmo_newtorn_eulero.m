clear all
% Algoritmo cinematica inversa

% Si richiede di implementare l'algoritmo
% di Newton che risolva la cinematica inversa del punto quando l'EE
% percorre la traiettoria C := (x+xc)^2+(z+zc)^2-r^2 = 0;
% xc = a1; zc = (a2 + r); r = 1/4 [u:m].
% a1=0.5; a2=1

% Traiettoria richiesta
phi=0:0.04:2*pi;
x=0.5+0.25*cos(phi);
z=5/4+0.25*sin(phi);

% Valori di partenza End Effector
d_0=0.3;
alfa_0=pi/3;

% Matrice di partenza per il calcolo
q_0=[d_0; alfa_0];
andamento=q_0;

% Errore di approssimazione
errore_minimo=10^(-12);

for i=1:length(z)
    % Inizializzo l'errore iniziale
    deltx=[1;1];
    % Calcolo per iterazioni
    while norm(deltx)>=norm(errore_minimo)
        
        % Valuto la posizione del punto dell'end effector rispetto a quella
        % voluta
        f1=x(i)-(0.5+cos(q_0(2)));
        f2=z(i)-(q_0(1)+sin(q_0(2)));
        F=[f1;f2];
        
        J11=0;
        J12=sin(q_0(2));
        J21=-1;
        J22=-cos(q_0(2));
        J=[J11 J12;
            J21 J22];
        % Calcolo il delta
        deltx=-(J\F);
        % Aggiorno il valore dei punti di partenza per l'iterazione
        % successiva
        q_0=q_0+deltx;    
    end
    % Memorizzo la posizione d e l'angolo tra il braccio 2 e 1
    andamento(:,i+1)=q_0;

end

% Calcolo movimento braccio
% Dimenzione bracci
a1=0.5;
a2=1;


num_punti=length(andamento(1,:));
traiettoria=length(x);

% Stampo la lunghezza del braccio d e angolo tra braccio 2 e 3 
figure(4)
subplot(2,1,1)
plot(1:num_punti,andamento(1,:)),grid on,title('Andamento d')
subplot(2,1,2)
plot(1:num_punti,andamento(2,:)),grid on,title('Variazione angolo [rad]')

% Verifica che la posizione richiesta sia la stessa seguita dall'end
% effector
rx=a1+a2*cos(andamento(2,:));
rz=andamento(1,:)+a2*sin(andamento(2,:));

for i=2:traiettoria
    figure(10)
    plot([x(1:i-1),x(i)],[z(1:i-1),z(i)],'k-'),grid on,title('Traiettoria end effector con traiettoria richiesta');
    hold
    
end

for i=1:3:num_punti
figure(10)
    X0 = [0 0];
    Z0 = [0 0];
    x0 = 0; z0=andamento(1,i);
    X1 = [0 x0];
    Z1 = [0 z0];
    x2 = a1; z2=andamento(1,i);
    X2 = [x0 x2];
    Z2 = [z0 z2]; 
    x3 = a1+a2*cos(andamento(2,i)); z3 =andamento(1,i)+a2*sin(andamento(2,i));
    X3 = [x2 x3];
    Z3 = [z2 z3];
Coords= line(X1, Z1,'LineWidth',2,'Color','blue','Marker','o','MarkerSize',2,'MarkerFaceColor','blue');
Coords= line(X2, Z2,'LineWidth',2,'Color','blue','Marker','o','MarkerSize',2,'MarkerFaceColor','blue');
Coords= line(X3, Z3,'LineWidth',2,'Color','k','Marker','o','MarkerSize',2,'MarkerFaceColor','k');
 pause(.3)
  pause on
end
grid on
  axis('equal');
  drawnow

% Calcolo Matrici rototraslazione
syms d theta a1 a2

DH=[
    pi/2, a1, d 0;
    0, a2 0 theta
    ];

A=Calcolo_denavit(DH);

A_02=eye(4,4)*A(:,:,1)*A(:,:,2);
disp('Matrice A_01')
disp(A(:,:,1))
disp('Matrice A_12')
disp(A(:,:,2))
disp('Matrice A_02')
disp(A_02)