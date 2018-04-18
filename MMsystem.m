%Gillespie SSA for Michaelis-Menten System

%%Hits an error when species reaches 0, needs a while loop to quit once subtrate has ran out. Could
% do with rewriting as a proper function

maxN=10^4;
X=[100;100;0;0];
k=[10,5,5]; %Rate Constants
v1=[-1,-1,1,0]';
v2=[1,1,-1,0]';
v3=[1,0,-1,1]';

N=0; %Reaction counter
time=0;
timeplot=zeros(1,maxN);
Xplot=zeros(4,maxN);
timeplot(1)=0;
TimeSpent=zeros(1,601);

while N<maxN
    
    alpha=[k(1)*X(1)*X(2),k(2)*X(3),k(3)*X(3)];
    alpha0=sum(alpha);
    u=rand(); %rand for reaction
    r=rand(); %rand for timestep
    tau=(1/alpha0)*log(1/u);
    time=time+tau;
    
    Xplot(:,N+1)=X';
    
     
    if(r<alpha(1)/alpha0)       %R1 has occurred
        X=X+v1;%+S(:,1);
    elseif(r<sum(alpha(1:2))/alpha0) %R2 has occurred
        X=X+v2;%S(:,2);
    else                             %R3 has occurred
        X=X+v3;%S(:,3);
    end
    N=N+1;

end

Xplot;
%TimeSpent=TimeSpent/(sum(TimeSpent));
%bar(0:600,TimeSpent);
plot(Xplot(1,:))
hold on
plot(Xplot(2,:),'r')
plot(Xplot(3,:),'g')
plot(Xplot(4,:),'y')

