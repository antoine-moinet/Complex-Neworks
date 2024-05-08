clear

N = input("number of points = ");
max = input("maximum value = ");
min = input("minimum value = ");
filename = input("filename is ","s");

s = load('-ascii',[filename '.dat']);

integers = 0:N;
x = min*(max/min).^(integers/N);
x2 = s(:,1);
x3 = s(:,2);
[h,whichBin] = histc(x2,x);

for i = 1:N
    flagBinMembers = (whichBin == i);
    binMembers     = x3(flagBinMembers);
    binMean(i)     = mean(binMembers);
end

op=fopen([filename 'Log.dat'],'w');

for j=1:N
	fprintf(op,"%e %e \n", (x(j)+x(j+1))/2, binMean(j));
endfor
fclose(op);





 
