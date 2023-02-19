clc
clear
n=1001//input("enter number of datapoints : ")
h=6.626D-34
h_bar=1.054*(1e-34)
l=input("enter value of l :") //s wave=0;p wave=1;d wave=2
b=14
rmin=0//input("enter r_min : ")
rmax=b*1D-15//input("enter r_max : ")
counter =4
e=1.6D-19 
m=1.67*10^(-27)
E0 = 9*1D9
r=linspace(rmin,rmax,n)//range
rnew =r(2:n-1)'
dr=r(2)-r(1)
C=-((h/(2*%pi))^2)/(2*m*(dr^2)*e)
//rtp calculation
w=5.34*(1D+21)
for i=1:4
    rtp(i)=sqrt((2*i-1)*(h/(2*%pi))/(m*w))
end
V=zeros(1,n)
for i=1:n
    V(i)=(0.5*m*w*w*(r(i)^2)+((l*(l+1)*(h_bar)^2)/(2*m*(r(i))^2)))/e  
     //in eV
end
A=eye(n-2,n-2)
v=diag(V(2:n-1))
//disp(A)
//----------------------By inbuilt command-------------------------
//tic()
D=(-2*C)*ones(n-2,1)
A1=diag(D)
//disp(A1)
C1=C*ones(n-3,1)
A2=diag(C1,1)
//disp(A2)
A3=diag(C1,-1)
//disp(A3)
//Hermitian matrix (tridiagonal)
H=A1+A2+A3+v
//disp("By inbuilt : ",H)
[a1 a2]=spec(H)
Z=spec(H)
for i=1:n
    Eth(i)=(((2*(i-1))+l+(3/2))*(h/(2*%pi))*w)/e
end

//--------------------------Normalizations & Expectation values-------------------------
vnew = V(2:n-1)
for i = 1:4
    eigenvec = a1(:,i).*a1(:,i)
    Nrm(i)=inttrap(rnew,eigenvec)
    u_r(:,i)=(1/sqrt(Nrm(i))).*a1(:,i)
    neigen(:,i) = u_r(:,i).*u_r(:,i)
    neigen1(:,i) = rnew.*neigen(:,i)
    neigen2(:,i) = rnew.*neigen1(:,i)
    neigen3(:,i) = vnew'.*neigen(:,i)
end
for i=1:counter
    ex_r(i)=inttrap(rnew,neigen1(:,i))
    ex_r2(i)= inttrap(rnew,neigen2(:,i))
    ex_V(i) = inttrap(rnew,neigen3(:,i))
end
for i = 1:n-3
    for j = 1:counter
        r2(i) = (rnew(i)+rnew(i+1))/2
        mid_u(i,j) = (u_r(i,j)+u_r(i+1,j))/2
        diff_u(i,j) = (u_r(i+1,j)-u_r(i,j))/dr
    end
end

for i = 1:n-4
    for j = 1:counter
        r3(i) = (r2(i) + r2(i+1))/2
        mid2_u(i,j) = (u_r(i+2,j) + 2*u_r(i+1,j)+ u_r(i,j))/4
        diff2_u(i,j) = (u_r(i+2,j) - 2*u_r(i+1,j) + u_r(i,j))/(dr*dr)
    end
end
//----------------Uncertainty check--------------------------------------
Un = (4*%pi)/h

for i=1:counter
    y2(:,i) = mid_u(:,i).*diff_u(:,i)
    y3(:,i) = mid2_u(:,i).*diff2_u(:,i)
    ex_p(i) = -1*%i*(h/(2*%pi))*inttrap(r2,y2(:,i))
    ex_p2(i) = -1*(h/(2*%pi))**2*inttrap(r3,y3(:,i))
    sig_r(i) = sqrt(ex_r2(i) - (ex_r(i)*ex_r(i)))
    sig_p(i) = sqrt(ex_p2(i) - (ex_p(i)*ex_p(i)))
    ex_K(i) = ex_p2(i)/(2*m*e)
    un(i) = Un*(sig_r(i).* sig_p(i))
end
//-----------------total energy------------------------
E = ex_V + ex_K
disp("Expectation value <r> =",ex_r)
disp("Expectation value <r2> =",ex_r2)
disp("Expectation value <p> =",ex_p)
disp("Expectation value <p2> =",ex_p2)
disp("Expectation value <V> =",ex_V)
disp("Expectation value <KE> =",ex_K)
disp("Standard deviation of r =",sig_r)
disp("Standard deviation of p =",sig_p)
disp("Uncertanity Product (hbar/2) = ",un)
disp("E_Th, E_comp(<KE>+<V>) - ")
disp(" E_th(eV)         E_comp(eV) (<KE>+<V>) Energy")
disp(string(Eth(1:3))+"        "+string(E(1:3)))
//Comparison between theo & Calculated value
//Z=Z(1:counter)
for i=1:3
    ratio(i)=(Z(i))/Eth(i)
    ratio2(i)=(Z(i)/((h_bar*w)/e))
end

show_window(1)
for i=1:3
    subplot(1,3,i)
    plot(rnew*(1D+15),u_r(:,i),'r')
    title("U(r) vs r plot ",'color','brown','Fontsize','5','Fontname','2')
    xlabel('r (fm) ------> ','color','brown','Fontsize','4','Fontname','4')
    ylabel('U(r) ------> ','color','brown','Fontsize','4','Fontname','4')
end
show_window(2)
for i=1:3
    subplot(1,3,i)
    plot(rnew*(1D+15),neigen(:,i))
    title("|U(r)|^2 vs r plot ",'color','brown','Fontsize','5','Fontname','2')
    xlabel('r (fm) ------> ','color','brown','Fontsize','4','Fontname','4')
    ylabel('|U(r)|^2 ------> ','color','brown','Fontsize','4','Fontname','4')
end
//----------Radial R(r)------------------
for i=1:3
    R(:,i)=u_r(:,i)./rnew
end

show_window(4)
for i=1:3
    subplot(1,3,i)
    plot(rnew*(1D+15),R(:,i),'g')
    title("R(r) vs r plot ",'color','brown','Fontsize','5','Fontname','4')
    xlabel('r (fm) ------> ','color','brown','Fontsize','4','Fontname','4')
    ylabel('R(r) ------> ','color','brown','Fontsize','4','Fontname','4')
end
