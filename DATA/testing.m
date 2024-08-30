%%
data = readmatrix("conf_MC.data", 'FileType', 'text');
histogram(data(1:500, 3));

%%
gdr = load("gdr.dat");
plot(gdr(:,1),gdr(:,2))

%%
pressure = load("pressure.dat");

