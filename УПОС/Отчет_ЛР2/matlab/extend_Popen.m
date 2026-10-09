function P_open_extended = extend_Popen(P_open, target_len)
    % EXTEND_POPEN Дополняет массив P_open случайными значениями
    % на основе статистики исходных данных.
    %
    % Входные данные:
    %   P_open - исходный массив (вектор)
    %   target_len - желаемая длина массива (например, 99)
    %
    % Выходные данные:
    %   P_open_extended - дополненный массив длины target_len
    
    current_len = length(P_open);
    
    if current_len >= target_len
        P_open_extended = P_open(1:target_len);
        return;
    end
    
    missing_len = target_len - current_len;
    
    % Приводим к вектору-столбцу для надежности
    P_open = P_open(:);
    
    % Вычисляем статистику исходных данных
    mu = mean(P_open);       % Среднее значение
    sigma = std(P_open);     % Стандартное отклонение
    
    % Генерируем недостающие значения
    % randn создает нормально распределенные случайные числа
    new_values = mu + sigma * randn(missing_len, 1);
    
    % (Опционально) Ограничиваем выбросы, чтобы они не были слишком далекими
    % от исходного диапазона
    min_val = min(P_open);
    max_val = max(P_open);
    new_values = max(min(new_values, max_val), min_val);
    
    % Объединяем исходный и новый массивы
    P_open_extended = [P_open; new_values];
end