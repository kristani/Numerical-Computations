file = sprintf("zad3/u_%i.txt",i)
splot file u 1:2:3 w pm3d  title sprintf("t=%i",i)
i=i+1
if (i < n) reread