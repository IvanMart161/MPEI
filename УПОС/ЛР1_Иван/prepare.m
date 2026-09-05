clear; clc; close all;

% --- Исходные параметры ---
f0 = 1000;       % Несущая частота в МГц (1 ГГц)
F = 1;           % Частота модуляции в МГц
beta1 = 10;      % Первый индекс модуляции (девиация 10 МГц)
beta2 = 2;       % Второй индекс модуляции (девиация 2 МГц)

figure('Color', 'w', 'Position', [100, 100, 900, 700]);

% --- 1. Спектр гармонического сигнала (п. 2) ---
subplot(3, 1, 1);
stem(f0, 1, 'b', 'LineWidth', 2, 'MarkerFaceColor', 'b');
title('Пункт 2: Спектр гармонической несущей (f_0 = 1 ГГц)');
xlabel('Частота, МГц');
ylabel('Относительная амплитуда');
xlim([f0 - 20, f0 + 20]);
ylim([0, 1.2]);
grid on;

% --- 2. Спектр ЧМ при beta = 2 (п. 3) ---
n_max_2 = ceil(beta2 + 4); % По правилу Карсона берем значимые гармоники
n2 = -n_max_2:n_max_2;
f_ch2 = f0 + n2 * F;
amp2 = abs(besselj(n2, beta2));

subplot(3, 1, 2);
stem(f_ch2, amp2, 'r', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
title(['Пункт 3: Спектр ЧМ-сигнала (\beta = 2, \Delta f = 2 МГц, F = 1 МГц)']);
xlabel('Частота, МГц');
ylabel('|J_n(\beta)|');
xlim([f0 - 20, f0 + 20]);
ylim([0, 0.8]);
grid on;

% --- 3. Спектр ЧМ при beta = 10 (п. 3) ---
n_max_1 = ceil(beta1 + 5); % Значимые гармоники для beta = 10
n1 = -n_max_1:n_max_1;
f_ch1 = f0 + n1 * F;
amp1 = abs(besselj(n1, beta1));

subplot(3, 1, 3);
stem(f_ch1, amp1, 'm', 'LineWidth', 1.5, 'MarkerFaceColor', 'm');
title(['Пункт 3: Спектр ЧМ-сигнала (\beta = 10, \Delta f = 10 МГц, F = 1 МГц)']);
xlabel('Частота, МГц');
ylabel('|J_n(\beta)|');
xlim([f0 - 20, f0 + 20]);
ylim([0, 0.5]);
grid on;