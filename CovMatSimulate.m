% will make a dummy MDS analysis before and after modulating the covariance
% matrix with typical modulations.
%%
figure(1);clf
ampc  = 1;
ampg  = 0;
freq  = 1;
phase = 0;
sigma = 0;

x0    = [ampc freq phase ampg sigma 0];
[y x] = meshgrid(linspace(0,2*pi-2*pi/8,8));
xdata = [y(:) x(:)];
%
model = reshape(CorrMatDecomposition(x0,xdata),8,8);
model = tril(model,0) + tril(model,-1)';
model = Scale(model) + eps;
model(1:9:end) = 1;
subplot(3,2,1);
imagesc(model)
[ori] = cmdscale(model);
subplot(3,2,2);
plot(zscore(ori(:,1)),zscore(ori(:,2)),'R.-')
hold on
text(zscore(ori(:,1))'+.1,zscore(ori(:,2))'+.1,num2str([1:8]'),'fontsize',15)
axis square
xx = xlim
yy = ylim

%% modulator
ampc  = 20;
ampg  = -50;
freq  = 1;
phase = 0;
sigma = 3;
x0    = [ ampc freq phase ampg sigma 0];
[y x] = meshgrid(linspace(0,2*pi-2*pi/8,8));
xdata = [y(:) x(:)];
modul = reshape(CorrMatDecomposition(x0,xdata),8,8);
modul = tril(modul,0) + tril(modul,-1)';
modul = Scale(modul) + eps;
subplot(3,2,3);
imagesc(modul)
%
new   = model + (modul-mean(modul(:)))*2;
new   = tril(new,0) + tril(new,-1)';
new   = Scale(new)  + eps;
new(1:9:end) = 1;
[D]   = cmdscale(new);
[d z] = procrustes(ori(:,1:2),D(:,1:2));
subplot(3,2,4);
plot(D(:,1),D(:,2),'o-')
hold on
text(D(:,1)',D(:,2)',num2str([1:8]'),'fontsize',15)
hold off
xlim(xx)
ylim(yy)
axis tight
%%
modul(1:3,1:3) = 0
modul(5:8,5:8) = 0
modul(5:8,1:3) = 0
modul(1:3,5:8) = 0
new   = model + -(modul-mean(modul(:)))*1;
subplot(3,2,5); 
imagesc(new);
subplot(3,2,6);
new   = tril(new,0) + tril(new,-1)';
new   = Scale(new)  + eps;
new(1:9:end) = 1;
[D]   = cmdscale(new)
%[d z] = procrustes(ori(:,1:2),D(:,1:2));
subplot(3,2,6);
plot(zscore(D(:,1)),zscore(D(:,2)),'o-')
xlim(xx)
ylim(yy)
%%
ampc  = 1;
ampg  = 1;
freq  = 1;
phase = pi;
sigma = 20;
x0    = [ ampc freq phase ampg sigma 0];
[y x] = meshgrid(linspace(0,2*pi-2*pi/8,8));
xdata = [y(:) x(:)];
modul = reshape(CorrMatDecomposition(x0,xdata),8,8);
modul = tril(modul,0) + tril(modul,-1)';
modul = Scale(modul) + eps;
new   = model.*(modul)*.15;
new   = tril(new,0) + tril(new,-1)';
new   = Scale(new)  + eps;
new(1:9:end) = 1;
subplot(3,2,5);
imagesc(modul)
[D] = cmdscale(new)
subplot(3,2,6);
plot(-zscore(D(:,1)),zscore(D(:,2)),'o-')
%