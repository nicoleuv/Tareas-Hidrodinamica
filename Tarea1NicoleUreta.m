clear; close all; clc;

g = 9.81;

%% 1) Cargar batimetría
b = load('REU2004bathy.txt');
x_bathy = b(:,1);
z_bathy = b(:,2);
h = -z_bathy;

% eliminar puntos con h<=0
valid = h > 1e-4;
x_bathy = x_bathy(valid);
h = h(valid);

%% 2) Definir casos a evaluar
casos = { ...
    struct('archivo','R26.mat','T',4,'color','r','marker','bo','nombre','R26'), ...
    struct('archivo','R27.mat','T',5,'color','g','marker','md','nombre','R27') ...
};

% Para guardar resultados de cada caso
H_models = cell(size(casos));
H_exps   = cell(size(casos));

%% 3) Loop sobre los casos
for idx = 1:numel(casos)
    caso = casos{idx};

    % Cargar datos experimentales
    S = load(caso.archivo);
    fn = fieldnames(S);
    R = S.(fn{1});
    H_exp = R.LWF.H(:);
    x_exp = R.xreal(:);

    % Parámetros
    T = caso.T;
    sigma = 2*pi/T;

    % Gage offshore
    h_at_gages = interp1(x_bathy, h, x_exp, 'linear', NaN);
    [~, i_off] = max(h_at_gages);
    H0 = H_exp(i_off);

    fprintf('\nCaso %s: H0 = %.4f m (x=%.2f m, h=%.2f m)\n', ...
        caso.nombre, H0, x_exp(i_off), h_at_gages(i_off));

    %% Resolver dispersión en todo el perfil
    k = zeros(size(h));
    for i = 1:length(h)
        hh = h(i);
        if hh <= 0
            k(i) = NaN;
            continue;
        end
        k0 = max(sigma^2/g, 1e-6);
        ki = k0;
        for it=1:60
            f  = g*ki*tanh(ki*hh) - sigma^2;
            df = g*tanh(ki*hh) + g*ki*hh*sech(ki*hh)^2;
            ki = ki - f/df;
            if abs(f) < 1e-10, break; end
        end
        k(i) = ki;
    end

    % Velocidades y H_model
    c  = sigma ./ k;
    cg = 0.5 .* c .* (1 + 2.*k.*h ./ sinh(2.*k.*h));

    % Offshore en el gage
    x0 = x_exp(i_off);
    h0 = interp1(x_bathy, h, x0, 'linear', 'extrap');
    k0 = max(sigma^2/g, 1e-6);
    for it=1:60
        f  = g*k0*tanh(k0*h0) - sigma^2;
        df = g*tanh(k0*h0) + g*k0*h0*sech(k0*h0)^2;
        k0 = k0 - f/df;
        if abs(f) < 1e-10, break; end
    end
    c0  = sigma / k0;
    cg0 = 0.5 * c0 * (1 + 2*k0*h0 / sinh(2*k0*h0));

    H_model = H0 .* sqrt(cg0 ./ cg);

    %% Número de Ursell experimental
    h_exp = interp1(x_bathy, h, x_exp, 'linear', NaN);
    k_exp = zeros(size(h_exp));
    for i = 1:length(h_exp)
        hh = h_exp(i);
        if isnan(hh) || hh <= 0
            k_exp(i) = NaN;
            continue;
        end
        ki = max(sigma^2/g, 1e-6);
        for it=1:60
            f  = g*ki*tanh(ki*hh) - sigma^2;
            df = g*tanh(ki*hh) + g*ki*hh*sech(ki*hh)^2;
            ki = ki - f/df;
            if abs(f) < 1e-10, break; end
        end
        k_exp(i) = ki;
    end
    L_exp  = 2*pi ./ k_exp;
    Ur_exp = H_exp .* (L_exp.^2) ./ (h_exp.^3);

    fprintf('\nNúmero de Ursell en gages experimentales (%s):\n', caso.nombre);
    disp(table((1:length(x_exp))', x_exp, h_exp, H_exp, k_exp, L_exp, Ur_exp, ...
        'VariableNames', {'Gage','x_m','h_m','H_exp_m','k_exp','L_m','Ur_exp'}));

    %% Número de Ursell modelado en gages
    L        = 2*pi ./ k;
    Ur_model = H_model .* (L.^2) ./ (h.^3);

    H_mod_at = interp1(x_bathy, H_model, x_exp, 'linear', NaN);
    h_mod_at = interp1(x_bathy, h, x_exp, 'linear', NaN);
    k_mod_at = interp1(x_bathy, k, x_exp, 'linear', NaN);
    L_mod_at = interp1(x_bathy, L, x_exp, 'linear', NaN);
    Ur_mod_at = interp1(x_bathy, Ur_model, x_exp, 'linear', NaN);

    fprintf('\nNúmero de Ursell en gages (Modelo %s):\n', caso.nombre);
    disp(table((1:length(x_exp))', x_exp, h_mod_at, H_mod_at, k_mod_at, L_mod_at, Ur_mod_at, ...
        'VariableNames', {'Gage','x_m','h_m','H_mod_m','k_mod','L_mod','Ur_mod'}));

    %% Comparación
    delta_Ur = round(Ur_mod_at - Ur_exp, 3);
    err_pct  = round((Ur_mod_at - Ur_exp) ./ Ur_exp * 100, 1);

    Ur_compare = table((1:length(x_exp))', x_exp, Ur_exp, Ur_mod_at, delta_Ur, err_pct, ...
        'VariableNames', {'Gage','x_m','Ur_exp','Ur_mod','Delta','Error_pct'});

    fprintf('\nComparación Número de Ursell (Experimental vs Modelo) en gages %s:\n', caso.nombre);
    disp(Ur_compare);

    %% Plot individual
    figure('Units','normalized','Position',[0.1 0.1 0.7 0.6]);
    yyaxis left
    plot(x_bathy, -h, 'k','LineWidth',1.5); hold on;
    ylabel('Batimetría z [m]', 'FontSize', 13);

    yyaxis right
    plot(x_bathy, H_model, caso.color, 'LineWidth',2);
    if ~isempty(H_exp)
        scatter(x_exp, H_exp, 70, caso.marker, 'filled');
        legend('Batimetría (elevación)', ...
               sprintf('H_{mod} (%s, sin rotura)', caso.nombre), ...
               sprintf('H_{exp} (%s)', caso.nombre), ...
               'Location','northwest', 'FontSize', 12);
    end

    xlabel('x [m]', 'FontSize', 13);
    ylabel('Altura de ola H [m]', 'FontSize', 13);
    title(sprintf('Modelo Baldock (sin rotura) %s - T=%.2f s, H0=%.3f m', ...
        caso.nombre, T, H0), 'FontSize', 14);
    grid on;

    % Guardar resultados
    H_models{idx} = H_model;
    H_exps{idx}   = H_exp;
end

%% 4) Plot conjunto
figure('Units','normalized','Position',[0.1 0.1 0.7 0.6]);
yyaxis left
h1 = plot(x_bathy, -h, 'k','LineWidth',1.5); hold on;
ylabel('Batimetría z [m]', 'FontSize', 13);

yyaxis right
h2 = plot(x_bathy, H_models{1}, 'r-', 'LineWidth',2);
h4 = plot(x_bathy, H_models{2}, 'g-', 'LineWidth',2);

h3 = scatter(x_exp, H_exps{1}, 70, 'bo','filled');
h5 = scatter(x_exp, H_exps{2}, 70, 'md','filled');

legend([h1 h2 h3 h4 h5], ...
       {'Batimetría (elevación)', ...
        'H_{mod} R26 (sin rotura)', ...
        'H_{exp} (R26)', ...
        'H_{mod} R27 (sin rotura)', ...
        'H_{exp} (R27)'}, ...
        'Location','northwest', 'FontSize', 12);

xlabel('x [m]', 'FontSize', 13);
ylabel('Altura de ola H [m]', 'FontSize', 13);
title('Modelo Baldock (sin rotura) R26 y R27', 'FontSize', 14);
set(gca, 'FontSize', 12, 'LineWidth', 1.2);
grid on;
