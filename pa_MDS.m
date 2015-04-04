function pa_MDS(betas,filename)
%given a VxCxS matrix M, with C conditions, V voxels and S subjects, this
%function will generate a 2D MDS-representation of different patterns stored
%in the columns of M. The error bars will be derived from resampling of
%different subjects.
B        = betas;
c        = GetFearGenColors;
dimen    = 2;
tsubject = size(B,3);
tbs      = 1000;
%%
D = zeros(8,8,tbs);
Y = zeros(8,2,tbs);
n = 1;
while n <= tbs
    subs     = randsample(1:tsubject,tsubject,1);
    Br       = mean(B(:,1:8,subs),3);
    D(:,:,n) = squareform( pdist(Br','euclidean'));
    Y(:,:,n) = mdscale(D(:,:,n),dimen,'start','cmdscale','criterion','metricstress');    
    if rem(n,100) == 0
        n
    end
    n        = n +1;
end

%%
for i =1:size(Y,3);[d, Z(:,:,i), tr] = procrustes(mean(Y,3),Y(:,:,i),'reflection',false);end
%
x = squeeze(Z(:,1,:));
y = squeeze(Z(:,2,:));
%
figure(2);clf;subplot(1,2,1);hold on;
for i =1:8;
    plot(x(i,:),y(i,:),'.','color',c(i,:));    
end
for i = 1:8
    plot(median(x(i,:)),median(y(i,:)),'s','color','k','markersize',20,'markerfacecolor',c(i,:));
end
hold off
axis tight;axis square
%%
subplot(1,2,2)
for i = 1:8;
    [mvar]=error_ellipse([x(i,:) ;y(i,:)]','color',c(i,:),'linewidth',2);
    plot(mean(x(i,:),2),mean(y(i,:),2),'o','color',c(i,:),'markersize',15,'markerfacecolor',c(i,:));
    hold on;
end
axis tight;axis square
SaveFigure(filename);