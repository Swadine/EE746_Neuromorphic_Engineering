function [V_,I_,avg_]=LIF(n,m)

h=0.1e-3;
a=0.1;
I_c=2.7e-9;
V_T = 20e-3;
E_L = -70e-3;

V=ones(n,m)*E_L;
I=ones(n,m);
for i = 1:n
  I(i,:)=I_c*(1+i*a);
end

V(:,1)=RK2(V(:,1), h,I(:,1));
c = zeros(n);
avg = zeros(n);
prev = zeros(n);
for a = 2:m 
   V(:,a)=RK2(V(:,a-1), 0.1e-3,I(:,a));
   for j=1:n
       if(V(j,a) > V_T)
           %V(j,a-1) = V(j,a);
           V(j,a) = E_L;
           if(c(j)>0)   
               avg(j) = avg(j) + (a-prev(j)-avg(j))/(c(j)+1);
           end 
           c(j) = c(j) + 1;
           prev(j) = a;
       end
   end
   
end
V_=V;
I_=I;
avg_=avg;
end

function k2 = RK2(V_vec, h, I_app)
k1 = V_vec + f_(V_vec, I_app) * (h/2);
k2 = k1 + f_(k1, I_app) * (h/2);
end

function deriv = f_(V, I_app)
C = 300e-12;
E_L = -70e-3;
G = 30e-9;
deriv = (1/C) * (I_app - G*(V - E_L));
end



