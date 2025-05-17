%% 清屏
clc;
clear;
close all; 


%% 导入数据

folderPath = 'E:\CODE\mat_app\TA\TR_TA_TRANS\MAPI';
new_deltaA = 'new-Hpvp-MAPI-DT-400ex-vispro-300pt-1spp-exp.dat';
samFile = 'Hpvp-MAPI-SAM.csv';
whiteFile = 'Hpvp-MAPI-WHITE.csv';
outFileName = 'deltaR-Hpvp-MAPI-DT-400ex-vispro-300pt-1spp-exp.dat';


fullFilePath = fullfile(folderPath, new_deltaA); % 拼接完整路径
if exist(fullFilePath, 'file')  % 判断文件是否存在
    rawData = importdata(fullFilePath);
else
    error('TA文件不存在！');
end

delaytime = rawData(1,2:end);  % 横轴时间（单位：ps）
wavelength = rawData(2:end,1);  % 纵轴波长（单位：nm）
deltaA_over_A = rawData(2:end,2:end);  % 矩阵形式：波长 × 时间
omega = 1240./wavelength;
omemin = 1.6; omemax = 2.5;
deomega = 0.001;
omega_seq = omemin:deomega:omemax;
deltaA_interp = zeros(length(omega_seq), length(delaytime));
for i = 1:length(delaytime)
    deltaA_interp(:,i) = interp1(omega, deltaA_over_A(:,i), omega_seq, 'spline');
end

%% 画图确认
figure
pcolor(omega_seq, delaytime, deltaA_interp');
shading interp;  % 平滑颜色过渡
colorbar;
xlabel('Delay Time (ps)');
ylabel('Photon Energy (eV)');
title('\DeltaA / A Interpolated');


%% 导入吸收谱
samPath = fullfile(folderPath, samFile);
whitePath = fullfile(folderPath, whiteFile);
if exist(samPath, 'file')
    samData = readmatrix(samPath);  % 导入样品透射光数据
else
    error('样品透射光文件不存在！');
end
if exist(whitePath, 'file')
    whiteData = readmatrix(whitePath);  % 导入入射光白光数据
else
    error('入射白光文件不存在！');
end

% 读取数据
samData = readmatrix(samPath);      % 1024x2，第一列波长，第二列强度
whiteData = readmatrix(whitePath);  % 1024x2，第一列波长，第二列强度

% 波长转能量 (eV)
omega_sam = 1240 ./ samData(:,1);
omega_white = 1240 ./ whiteData(:,1);

% 强度
intensity_sam = samData(:,2);
intensity_white = whiteData(:,2);

% 插值到统一能量序列 omega_seq
intensity_sam_interp = interp1(omega_sam, intensity_sam, omega_seq, 'pchip', 'extrap');
intensity_white_interp = interp1(omega_white, intensity_white, omega_seq, 'pchip', 'extrap');

% 计算吸收强度
absorption = (intensity_white_interp - intensity_sam_interp) ./ intensity_white_interp;

% 绘制吸收光谱（可选）
figure;
plot(omega_seq, absorption);
% plot(omega_sam, intensity_sam,omega_white, intensity_white, omega_seq, absorption);
xlabel('Photon Energy (eV)');
ylabel('Absorption Intensity');
title('Absorption Spectrum ( (white - sam) / white )');
grid on;

d = 400; % nm，样品厚度

% 波长序列，单位 nm（从之前的 omega_seq 转回波长）
lambda_seq = 1240 ./ omega_seq; % nm

% Δk = (lambda / (4*pi*d)) * ΔA
% deltaA_interp是矩阵，逐元素计算
delta_k = (lambda_seq ./ (4 * pi * d))' .* deltaA_interp;  % 注意lambda_seq转置成列向量
% k ≈ absorption * lambda / (4 * pi * d)
k = (absorption' .* lambda_seq') ./ (4 * pi * d);  % 全部转置为列向量

%% 绘制稳态 k（消光系数）
figure;
plot(omega_seq, k, 'b', 'LineWidth', 1.5);
xlabel('Photon Energy (eV)');
ylabel('k');
title('Steady-State Extinction Coefficient k');
grid on;

%% 绘制 Δk（二维色图）
figure;
pcolor(omega_seq, delaytime, delta_k');  % delta_k 尺寸为能量 × 时间，转置为时间 × 能量
shading interp;
colorbar;
xlabel('Delay Time (ps)');
ylabel('Photon Energy (eV)');
title('\Deltak (Interpolated from \DeltaA)');

n_static = kk_transform(k, omega_seq);
delta_n = kk_transform(delta_k, omega_seq);

%% Δn转换ΔR
n0 = 1; % 空气折射率

% 复折射率
n_complex = n_static + 1i * k;

R0 = abs((n_complex - n0) ./ (n_complex + n0)).^2;

% delta_n 和 delta_k 是矩阵，列为不同delay time，行对应能量点
delta_n_matrix = delta_n; 
delta_k_matrix = delta_k;

deltaR_over_R = zeros(size(delta_n_matrix)); % 初始化结果矩阵

for t = 1:length(delaytime)
    n_t = n_static + delta_n_matrix(:,t);
    k_t = k + delta_k_matrix(:,t);
    n_complex_t = n_t + 1i * k_t;
    
    R_t = abs((n_complex_t - n0) ./ (n_complex_t + n0)).^2;
    
    deltaR_over_R(:,t) = (R_t - R0) ./ R0; % 计算瞬态反射率变化
end

figure;
pcolor(omega_seq, delaytime, deltaR_over_R');
shading interp;
colorbar;
xlabel('Photon Energy (eV)');
ylabel('Delay Time (ps)');
title('\Delta R / R');

%% 导出
% 构造输出矩阵大小：(波长数+1) × (时间点数+1)
num_wavelengths = length(lambda_seq);
num_times = length(delaytime);

output_matrix = zeros(num_wavelengths + 1, num_times + 1);

% 第一行第一列设为0
output_matrix(1,1) = 0;

% 第一行第二列开始填delaytime
output_matrix(1, 2:end) = delaytime;

% 第一列第二行开始填lambda_seq
output_matrix(2:end, 1) = lambda_seq';

% 剩下填deltaR_over_R数据（波长 × 时间）
output_matrix(2:end, 2:end) = deltaR_over_R;

% 写文件
outFilePath = fullfile(folderPath, outFileName);

writematrix(output_matrix, outFilePath);

disp(['保存成功：', outFilePath]);

%% 函数 kk变换
function delta_n = kk_transform(delta_k, omega_seq)
% delta_k: 每列是一个 delay 时间的 Δk(ω)
% omega_seq: 行向量，单位 eV

    N = length(omega_seq);
    delta_n = zeros(N, size(delta_k,2)); % 初始化输出矩阵

    for t = 1:size(delta_k,2)  % 对每个时间点
        dk = delta_k(:,t);
        dn = zeros(N,1);
        for j = 1:N
            w0 = omega_seq(j);
            sum_k = 0;
            for i = 1:N
                if i ~= j
                    w = omega_seq(i);
                    sum_k = sum_k + dk(i) * w / (w^2 - w0^2);
                end
            end
            dn(j) = 2 / pi * sum_k * (omega_seq(2) - omega_seq(1));  % 近似积分
        end
        delta_n(:,t) = dn;
    end
end

