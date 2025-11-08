function out = simulate_zernike_imaging(varargin)
%校验传参
P = local_defaults();
if nargin==1 && isstruct(varargin{1})
    S = varargin{1};
else
    if mod(nargin,2)~=0
        error('错了');
    end
    S = struct(varargin{:});
end
fn = fieldnames(P);
for i=1:numel(fn)
    k = fn{i};
    if isfield(S,k) 
        P.(k) = S.(k); 
    end
end

req = {'f','z','R','lambda'};
for i=1:numel(req)
    assert(~isempty(P.(req{i})), '缺少必要参数：%s', req{i});
end
assert(P.N>0 && mod(P.N,2)==0, 'N需为正偶数(如512)。');
assert(P.R>0 && P.f>0 && P.lambda>0 && P.oversample>0, 'f/R/lambda/oversample需为正。');
%准备网格
N = P.N; 
R = P.R; 
lambda = P.lambda;
L = 2*R*P.oversample;              % 计算窗口边长
dx = L/N; 
x1 = (-N/2:N/2-1)*dx;
[y,x] = meshgrid(x1,x1);
r = hypot(x,y); 
theta = atan2(y,x);
rho = r/R; 
pupil = rho<=1;
k = 2*pi/lambda;
%传递函数
df = 1/L; 
f1 = (-N/2:N/2-1)*df;
[fy,fx] = meshgrid(f1,f1);
% %菲涅尔传递函数
% H = exp(1i*k*P.z) .* exp(-1i*pi*lambda*P.z*(fx.^2 + fy.^2));
% Hf = fftshift(H);
kx = 2*pi*fx;                
ky = 2*pi*fy;
kz = sqrt( (k)^2 - kx.^2 - ky.^2 );

H = exp(1i * P.z * kz);
Hf = fftshift(H);

%理想透镜透过率函数 相位部分
lens_phase = -k*(x.^2 + y.^2)/(2*P.f);

% Zernike系数
if isempty(P.coeffs)
    M = 12; 
    T = 8; %12组 8项
    waves_sigma = 0.25;                 % 系数高斯分布 随机
    C = (2*pi)*(waves_sigma*randn(M,T));
    C(:,1) = 0;                          % 去活塞
    C(1,:) = 0;                          % 第一组所有系数设为0 平面波
    
else
    C = P.coeffs;
    [M,T] = size(C);
    assert(T==8, 'coeffs需为 Mx8。');
    C(1,:) = 0;                          % 强制第一组为平面波
end

Z = cell(1,8);
for j=1:8
    Z{j} = zernike_noll(j, rho, theta); 
end

wf = zeros(N,N,M); %相位
I = zeros(N,N,M);  %强度
%Zernike项加起来成为了相位
for m=1:M
    phi = zeros(N,N);
    for j=1:8
        phi = phi + C(m,j).*Z{j}; 
    end
    U0 = pupil .* exp(1i*(phi + lens_phase)); %pupil外设置为0
    Uz = ifftshift(ifft2( fft2(fftshift(U0)).*Hf )); %fftshit 低频放中心
    It = abs(Uz).^2; 
    It = It./max(It(:)+eps); %模方归一化到 [0,1]
    wf(:,:,m) = phi.*pupil;
    I(:,:,m)  = It;
end

out = struct();
out.x=x; 
out.y=y; 
out.pupil=pupil; 
out.wf=wf; 
out.I=I;
out.meta = struct('N',N,'dx',dx,'L',L,'R',R,'lambda',lambda,'k',k, ...
    'fx',fx,'fy',fy,'H',H,'lens_phase',lens_phase,'coeffs',C, ...
    'z',P.z,'f',P.f,'oversample',P.oversample);

end
%子函数：获取默认参数设置
function P = local_defaults()
P = struct();
P.f = []; P.z = []; P.R = []; P.lambda = [];
P.N = 512;
P.oversample = 2;
P.coeffs = [];
end
