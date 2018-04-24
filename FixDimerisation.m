%Gillespie SSA for Reversible Dimerisation System
maxN=10^6;
X=[100;100];
k=[142,1,880,92.8,10,500,6]; %Rate Constants
v1=[1,0]';
v2=[-1,0]';
v3=[1,0]';
v4=[-1,0]';
v5=[-2,1]';
v6=[2,-1]';
v7=[0,-1]';   %Stoichiometric vectors

N=0; %Reaction counter
time=0;
TimeSpentS=zeros(1,900);

while N<maxN
    
    alpha=[k(1)*X(2),k(2)*X(1)*X(2),k(3),k(4)*X(1),k(5)*X(1)*(X(1)-1),k(6)*X(2),k(7)*X(2)];
    alpha0=sum(alpha);
    u=rand(); %rand for reaction
    r=rand(); %rand for timestep
    tau=(1/alpha0)*log(1/u);
    time=time+tau;
    S=X(1);+2*X(2);
    TimeSpentS(S+1)=TimeSpentS(S+1)+tau;
    
    if(r<alpha(1)/alpha0)            %R1 has occurred
        X=X+v1;
    elseif(r<sum(alpha(1:2))/alpha0) %R2 has occurred
        X=X+v2;
    elseif(r<sum(alpha(1:3))/alpha0) %R3 has occurred
        X=X+v3;
    elseif(r<sum(alpha(1:4))/alpha0) %R4 has occurred
        X=X+v4;
    elseif(r<sum(alpha(1:5))/alpha0) %R5 has occurred
        X=X+v5;
    elseif(r<sum(alpha(1:6))/alpha0) %R6 has occurred
        X=X+v6;
    else                        %R7 has occurred
        X=X+v7;
    end
    N=N+1;
end

TimeSpentS=TimeSpentS/(sum(TimeSpentS));
plot(TimeSpentS);
