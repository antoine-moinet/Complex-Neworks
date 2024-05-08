clear
cmax = 10^5;
cmin = 1;

R = input("sum top edge = ");
first_guess = input("first guess = ")

for i = 1:80
	
	beta = 0.025*i;
	A = (beta)/(cmin^-beta-cmax^-(beta));
#	sum_of_r = 0;
#	for k=1:R
#		sum_of_r = sum_of_r + k*

# molloy reed criterion

	s(i) = fsolve(@(t) quad(@(c) A*c^(-beta-1)*sum([0:R].^2 .*(erf(([0:R]+1)/sqrt(c*t))-erf([0:R]/sqrt(c*t)))),cmin,cmax) + 3*quad(@(c) A*c^(-beta-1)*sum([0:R].*(erf(([0:R]+1)/sqrt(c*t))-erf([0:R]/sqrt(c*t)))),cmin,cmax)^2 - 3*quad(@(c) A*c^(-beta-1)*sum([0:R].*(erf(([0:R]+1)/sqrt(c*t))-erf([0:R]/sqrt(c*t)))),cmin,cmax), first_guess);

# eigen value criterion
	
#	s(i) = fsolve(@(t) quad(@(c) A*c^(-beta-1)*sum([0:R].^2 .*(erf(([0:R]+1)/sqrt(c*t))-erf([0:R]/sqrt(c*t)))),cmin,cmax) - quad(@(c) A*c^(-beta-1)*sum([0:R].*(erf(([0:R]+1)/sqrt(c*t))-erf([0:R]/sqrt(c*t)))),cmin,cmax)^2 + quad(@(c) A*c^(-beta-1)*sum([0:R].*(erf(([0:R]+1)/sqrt(c*t))-erf([0:R]/sqrt(c*t)))),cmin,cmax) - 1, first_guess);
	
endfor

op=fopen('levyPerc5mr.dat','w');

for j=1:80
	fprintf(op,"%f %f \n", 0.025*j, s(j));
endfor

fclose(op);



 
