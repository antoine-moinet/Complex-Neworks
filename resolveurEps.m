clear
eps = 0.001;

#R = input("sum top edge = ");

for i = 1:30
	
	#beta = 0.1*i;
	#A = (beta)/(eps^-(beta)-1);
	gam = 0.1*i+0.001; 

	#s(i) = fsolve(@(t) quad(@(c) A*c^(-beta-1)*sum([0:R].^2 .*(erf(([0:R]+1)/sqrt(c*t))-erf([0:R]/sqrt(c*t)))),eps,1) -quad(@(c) A*c^(-	beta-1)*sum([0:R].*(erf(([0:R]+1)/sqrt(c*t))-erf([0:R]/sqrt(c*t)))),eps,1)^2 + quad(@(c) A*c^(-beta-1)*sum([0:R].*(erf(([0:R]+1)/sqrt(c*t))-erf([0:R]/sqrt(c*t)))),eps,1)-1, 1.);
	s(i) = fsolve(@(x) (1-gam)/(1-eps^(1-gam))*quad(@(a) a^(-gam)*(((1-a)/(1-eps))^x+eps),eps,1) - (1-gam)/(1-eps^(1-gam))*quad(@(a) a^(1-gam),eps,1),1.);  
	
endfor

op=fopen('multiperc.dat','w');

for j=1:30
	fprintf(op,"exponent[%d]  = %f; \n", j, s(j));
endfor

fclose(op);



 
