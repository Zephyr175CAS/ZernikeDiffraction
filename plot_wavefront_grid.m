function fig = plot_wavefront_grid(out, opts)
if nargin<2
    opts=struct(); 
end
if ~isfield(opts,'title')
    opts.title = ""; 
end
if ~isfield(opts,'wrap')
    opts.wrap = true; 
end
if ~isfield(opts,'fontSize')
    opts.fontSize = 12; 
end
if ~isfield(opts,'cm')
    opts.cm = 'hsv';
end
if ~isfield(opts,'layout')
    opts.layout = [3,4];    %布局 三行四列 
end

x=out.x; 
y=out.y; 
pupil=out.pupil; 
wf=out.wf;
M = size(wf,3); 
R = out.meta.R; 
f = out.meta.f; 
z = out.meta.z; 
lambda = out.meta.lambda;

if opts.wrap
    W = angle(exp(1i*wf)); 
    clim = [-pi,pi]; 
else
    a=max(abs(wf),[],'all'); 
    W=wf; clim=[-a,a]; 
end

fig = figure('Color','w','Units','normalized','Position',[0.06 0.08 0.88 0.82]);
t = tiledlayout(opts.layout(1), opts.layout(2),'TileSpacing','compact','Padding','compact');

if strlength(opts.title)>0
    title(t, sprintf('%s  |  \\lambda=%.0f nm, f=%.1f mm, z=%.1f mm, R=%.2f mm', opts.title, lambda*1e9, f*1e3, z*1e3, R*1e3), 'FontSize',opts.fontSize+3,'FontWeight','bold');
end

for m=1:M
    nexttile;
    A = W(:,:,m); 
    A(~pupil)=NaN;
    imagesc(x(1,:),y(:,1),A,clim); 
    axis image off;
    title(sprintf('Set %d',m),'FontSize',opts.fontSize);
end
colormap(fig, opts.cm);
cb = colorbar; cb.Location='eastoutside'; cb.Label.String='相位'; cb.Label.FontSize=opts.fontSize;
end
