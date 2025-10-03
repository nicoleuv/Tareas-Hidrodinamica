%% Baldock 1998, solo rotura, sin roller
clear; close all; clc;

g = 9.81;

%% 1) Cargar batimetría
b = load('REU2004bathy.txt');
x_bathy = b(:,1);
z_bathy = b(:,2);
h = -z_bathy;

% eliminar puntos con h<=0 para trabajar solo con profundidades positivas
valid = h > 1e-4;
x_bathy = x_bathy(valid);
h = h(valid);

%% 2) Casos a evaluar RXX
casos = { ...
    struct('archivo','R26.mat','T',4,'color','r','marker','bo'), ...
    struct('archivo','R27.mat','T',5,'color','g','marker','md') ...
};

% añadir campo 'nombre''RXX' en casos 
for i = 1:numel(casos)
    [~,fname,~] = fileparts(casos{i}.archivo);
    casos{i}.nombre = fname;  
end

% Para guardar resultados
H_models = cell(size(casos));
H_exps   = cell(size(casos));

%% 3) Loop sobre los casos
for idx = 1:numel(casos)
    caso = casos{idx};

    % Cargar datos experimentales de RXX (gages)
    S = load(caso.archivo);
    fn = fieldnames(S);
    R = S.(fn{1});
    H_exp = R.LWF.H(:);
    x_exp = R.xreal(:);

    % Parámetros
    T = caso.T;
    omega = 2*pi/T;

    % Gage offshore (más profundo, mas mar adentro para condicion inicial)
    h_at_gages = interp1(x_bathy, h, x_exp, 'linear', NaN);
    [~, i_off] = max(h_at_gages);
    H0 = H_exp(i_off);
    x0 = x_exp(i_off);
    h0 = h_at_gages(i_off);

    fprintf('\nCaso %s: H0 = %.4f m (x=%.2f m, h=%.2f m)\n', ...
        caso.nombre, H0, x0, h0);

    %% Resolver dispersión en todo el perfil
    k = zeros(size(h));
    for i = 1:length(h)
        hh = h(i);
        if hh <= 0
            k(i) = NaN;
            continue;
        end

        % Newton-Raphson para resolver ecuación de dispersión: gk*tanh(kh)=ω²
        k0 = max(omega^2/g, 1e-6);
        ki = k0;
        for it=1:60
            f  = g*ki*tanh(ki*hh) - omega^2;
            df = g*tanh(ki*hh) + g*ki*hh*sech(ki*hh)^2;
            ki = ki - f/df;
            if abs(f) < 1e-10, break; end
        end
        k(i) = ki;
    end
    c  = omega ./ k;
    cg = 0.5 .* c .* (1 + 2.*k.*h ./ sinh(2.*k.*h));

    % Co-localizar condición inicial en el nodo más cercano al gage offshore
    [~, i0] = min(abs(x_bathy - x0));
    k0 = k(i0);
    c0 = omega/k0;
    cg0 = 0.5 * c0 * (1 + 2*k0*h0/sinh(2*k0*h0));

    %% 5) Baldock 1998 con rotura (marching onshore)
    rho = 1000; 
    fp  = 1/T; 
    Bbj = 1.0;
    A   = rho*g*Bbj*fp/4;
    gamma_b = 0.55; % condición de rompiente

    % inicializar
    H_model = NaN(size(h));
    H_model(i0) = H0;
    F = (1/8)*rho*g*H0^2 * cg0; % flujo de energía inicial

    % propagacion en direccion onshore
    for j = i0+1:length(x_bathy)
        dx   = x_bathy(j) - x_bathy(j-1);
        Hb   = gamma_b * h(j);
        Hrms = max(H_model(j-1)/sqrt(2), 1e-6);

        % disipación según Baldock <ε>
        epsj = A * trapz_linspace(Hb, max(6*Hrms,Hb+1e-6), ...
            @(H) H .* (2*H./Hrms.^2) .* exp(-(H./Hrms).^2));

        % actualiza flujo y energía
        F = max(F - epsj*dx, 0);
        E = F / max(cg(j), 1e-8);

        % altura nueva limitada con rompiente
        H_try      = sqrt(8*E/(rho*g));
        H_model(j) = min(H_try, Hb);
    end

    % propagacion en direccion offshore (sin disipacion relevante)
    for j = i0-1:-1:1
        F  = (1/8)*rho*g*H_model(j+1)^2 * cg(j+1);
        E  = F / max(cg(j),1e-8);
        H_model(j) = sqrt(8*E/(rho*g));
    end

    %% Guardar resultados
    H_models{idx} = H_model;
    H_exps{idx}   = H_exp;

    %% Plot individual
    %a la izquierda esta las unidades de la batimetria
    figure('Units','normalized','Position',[0.1 0.1 0.7 0.6]);
    yyaxis left
    h_bathy = plot(x_bathy, -h, 'k','LineWidth',1.5); hold on;
    ylabel('Batimetría z [m]', 'FontSize', 13);

    % a la derecha las unidades de altura 
    yyaxis right
    h_mod = plot(x_bathy, H_model, caso.color, 'LineWidth',2); hold on;
    h_break = plot(x_bathy, gamma_b*h, 'k--','LineWidth',1.5);
    if ~isempty(H_exp)
        h_exp = scatter(x_exp, H_exp, 70, caso.marker, 'filled');
        legend([h_bathy h_mod h_exp h_break], ...
               {'Batimetría (elevación)', ...
                sprintf('H_{mod} (%s, con rotura)', caso.nombre), ...
                sprintf('H_{exp} (%s)', caso.nombre), ...
                sprintf('Línea de rompiente (γ = %.2f)', gamma_b)}, ...
                'Location','northwest', 'FontSize', 12);
    else
        legend([h_bathy h_mod h_break], ...
               {'Batimetría (elevación)', ...
                sprintf('H_{mod} (%s, con rotura)', caso.nombre), ...
                sprintf('Línea de rompiente (γ = %.2f)', gamma_b)}, ...
                'Location','northwest', 'FontSize', 12);
    end

    xlabel('x [m]', 'FontSize', 13);
    ylabel('Altura de ola H [m]', 'FontSize', 13);
    title(sprintf('Modelo Baldock (con rotura) %s - T=%.2f s, H0=%.3f m', ...
        caso.nombre, T, H0), 'FontSize', 14);
    ylim([0 1.6]);
    grid on;
end

%% 4) Plot conjunto
%a la izquierda esta las unidades de la batimetria
figure('Units','normalized','Position',[0.1 0.1 0.7 0.6]);
yyaxis left
h1 = plot(x_bathy, -h, 'k','LineWidth',1.5); hold on;
ylabel('Batimetría z [m]', 'FontSize', 13);

% a la derecha las unidades de altura
yyaxis right
h2 = plot(x_bathy, H_models{1}, 'r-', 'LineWidth',2);
h4 = plot(x_bathy, H_models{2}, 'g-', 'LineWidth',2);

h3 = scatter(x_exp, H_exps{1}, 70, 'bo','filled');
h5 = scatter(x_exp, H_exps{2}, 70, 'md','filled');

% Línea de rompiente en plot conjunto
h_break = plot(x_bathy, gamma_b*h, 'k--','LineWidth',1.5);

legend([h1 h2 h3 h4 h5 h_break], ...
       {'Batimetría (elevación)', ...
        sprintf('H_{mod} %s (con rotura)', casos{1}.nombre), ...
        sprintf('H_{exp} (%s)', casos{1}.nombre), ...
        sprintf('H_{mod} %s (con rotura)', casos{2}.nombre), ...
        sprintf('H_{exp} (%s)', casos{2}.nombre), ...
        sprintf('Línea de rompiente (γ = %.2f)', gamma_b)}, ...
        'Location','northwest', 'FontSize', 12);

xlabel('x [m]', 'FontSize', 13);
ylabel('Altura de ola H [m]', 'FontSize', 13);
title(sprintf('Modelo Baldock (con rotura) %s y %s', casos{1}.nombre, casos{2}.nombre), ...
      'FontSize', 14);

ylim([0 1.6]);
set(gca, 'FontSize', 12, 'LineWidth', 1.2);
grid on;

%% Función auxiliar
function I = trapz_linspace(a,b,fun)
    n = 200;
    x = linspace(a,b,n);
    fx = fun(x);
    I = trapz(x,fx);
end


