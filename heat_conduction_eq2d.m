% ������
clc; clear all;

% �O���[�o���ϐ��錾
global dt n ita L x y

% �p�����[�^�[
L = 10;% �v�Z�̈�
n = 30;% �v�f��
dL = L / n;% �v�f�T�C�Y
c = 1;% ��M
rho = 1;% ���x
k = 1;% �M�`����
alpha = k / (c * rho);% �M�g�U��

t_total = 10;% SIM�v�Z����
dt = 0.01;% �^�C���X�e�b�v
ita = t_total / dt;% ������

% ���W�̐���
x = 0 : L / n : L;
y = L : - L / n : 0;
[X, Y] = meshgrid(x, y);

% �z��m��
T2d = zeros(n, n, ita);
wall_2d = zeros(n, n);

% �M���E����
T1 = 100;% �k��
T2 = 0;% ����
T3 = 0;% �쑤
T4 = 0;% ����

% ���E������ݒ�ƍs��̍쐬
[T2d, wall_2d] = DefineWall(T2d, wall_2d, T1, T2, T3, T4);

% ���x�x�N�g���̍\�z
T1d = zeros(n * n, ita);
for i = 1 : ita
    T1d(:, i) = reshape(T2d(:, :, i).', [n * n, 1]);% ���V�F�C�v�̑O�ɓ]�u���Ă���
end
wall_1d = reshape(wall_2d.', [n * n, 1]);

% ���Ԕ��W�s��̍\�z
A = TransitionMatrix(wall_1d, alpha, dL);

% ���Ԕ��W
for i = 2 : ita
    
    % ���Ԕ��W
    T1d(:, i) = A * T1d(:, i - 1);
    T2d(:, :, i) = (reshape(T1d(:, i), [n, n])).';
    
    % �R�}���h�E�B���h�E�ւ̏o��
    txt = ['ita = ',num2str(i),' / ',num2str(ita)];
    disp(txt);
    
    % ����
    vis_contour('temp.gif', i, T2d(:, :, i), 0, 100, 1);
    
end

%% �ȉ��֐�

function[T2d, wall_2d] = DefineWall(T2d, wall_2d, T1, T2, T3, T4)

% �O���[�o���ϐ��Ăяo��
global n

% �Ǎs��
% 1|2|3
% 8|0|4
% 7|6|5

for j = 1 : n
    for i = 1 : n
        if i == 1 && j == 1% ����
            T2d(i, j, :) = T1;
            wall_2d(i, j) = 1;
        elseif i == 1 && j == n% �E��
            T2d(i, j, :) = T1;
            wall_2d(i, j) = 3;
        elseif i == n && j == n% �E��
            T2d(i, j, :) = T3;
            wall_2d(i, j) = 5;
        elseif i == n && j == 1% ����
            T2d(i, j, :) = T3;
            wall_2d(i, j) = 7;
        elseif i == 1
            T2d(i, j, :) = T1;
            wall_2d(i, j) = 2;
        elseif j == n
            T2d(i, j, :) = T3;
            wall_2d(i, j) = 4;
        elseif i == n
            T2d(i, j, :) = T2;
            wall_2d(i, j) = 6;
        elseif j == 1
            T2d(i, j, :) = T4;
            wall_2d(i, j) = 8;
        end
    end
end

end


function[A] = TransitionMatrix(wall_1d, alpha, dL)

% �O���[�o���ϐ��Ăяo��
global n dt

% �Ǎs��
% 1|2|3
% 8|0|4
% 7|6|5

A = zeros(n * n, n * n);
for i = 1 : n * n
    if wall_1d(i, 1) == 1% �p�͓�̋��E�����̕��ς��Ƃ�B
        A(i, i + 1) = 1 / 2;
        A(i, i + n) = 1 / 2;
    elseif wall_1d(i, 1) == 2% �k����
        A(i, i) = 1;
    elseif wall_1d(i, 1) == 3% �p�͓�̋��E�����̕��ς��Ƃ�B
        A(i, i - 1) = 1 / 2;
        A(i, i + n) = 1 / 2;
    elseif wall_1d(i, 1) == 4% ������
        A(i, i) = 1;
    elseif wall_1d(i, 1) == 5% �p�͓�̋��E�����̕��ς��Ƃ�B
        A(i, i - 1) = 1 / 2;
        A(i, i - n) = 1 / 2;
    elseif wall_1d(i, 1) == 6% �쑤��
        A(i, i) = 1;
    elseif wall_1d(i, 1) == 7% �p�͓�̋��E�����̕��ς��Ƃ�
        A(i, i + 1) = 1 / 2;
        A(i, i - n) = 1 / 2;
    elseif wall_1d(i, 1) == 8% ������
        A(i, i) = 1;
    elseif wall_1d(i, 1) == 0% �ǂłȂ��Ȃ��
        A(i, i) = 1 - 4 * dt * alpha / dL^2;
        A(i, i - 1) = dt * alpha / dL^2;
        A(i, i + 1) = dt * alpha / dL^2;
        A(i, i - n) = dt * alpha / dL^2;
        A(i, i + n) = dt * alpha / dL^2;
    end
end

end

function[] = vis_contour(filename, timestep, u, maxrange, minrange, fignum)

% �O���[�o���ϐ��Ăяo��
global dt L x y

figure(fignum);
h = imagesc(x, y, u);

h.AlphaData = isfinite(u); % NaN��Inf�𓧖��ɂ���
title(['time = ', num2str(timestep * dt, '%.3f')]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);
axis equal; axis tight; axis on;
xlabel('x')
ylabel('y')
view(0, 270); % ���_�̐ݒ�

xticks(0 : round(L / 5) : L);
yticks(0 : round(L / 5) : L);

colorbar;
caxis([maxrange minrange])
frame = getframe(fignum);
im = frame2im(frame);
[imind, cm] = rgb2ind(im, 256);
if timestep == 2
    imwrite(imind, cm, filename, 'gif', 'DelayTime', 0.001, 'Loopcount', inf);
elseif rem(timestep, 10) == 0
    imwrite(imind, cm, filename, 'gif', 'DelayTime', 0.001, 'WriteMode', 'append');
end

end