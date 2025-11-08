function fig = plot_intensity_grid(out, opts)
% opts: title, fontSize, cm ('turbo'默认), layout [rows,cols]
if nargin<2
    opts=struct(); 
end
if ~isfield(opts,'title')
    opts.title = ""; 
end
if ~isfield(opts,'fontSize')
    opts.fontSize = 12; 
end
if ~isfield(opts,'cm')
    opts.cm = 'turbo'; 
end
if ~isfield(opts,'layout') 
    opts.layout = [3,4];    %布局 三行四列
end

x=out.x; 
y=out.y; 
I=out.I;
M = size(I,3); 
R = out.meta.R; 
f = out.meta.f; 
z = out.meta.z; 
lambda = out.meta.lambda;


fig = figure('Color','w','Units','normalized','Position',[0.06 0.08 0.88 0.82]);
t = tiledlayout(opts.layout(1), opts.layout(2),'TileSpacing','compact','Padding','compact');

if strlength(opts.title)>0
    title(t, sprintf('%s  |  \\lambda=%.0f nm, f=%.1f mm, z=%.1f mm, R=%.2f mm', opts.title, lambda*1e9, f*1e3, z*1e3, R*1e3), 'FontSize',opts.fontSize+3,'FontWeight','bold');
end

for m=1:M
    nexttile;
    imagesc(x(1,:),y(:,1),I(:,:,m),[0 1]); 
    axis image off;
    title(sprintf('Set %d',m),'FontSize',opts.fontSize);
end
colormap(fig, opts.cm);
cb = colorbar; cb.Location='eastoutside'; cb.Label.String='强度(归一化)'; cb.Label.FontSize=opts.fontSize;
end
