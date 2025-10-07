%% Baldock 1998, SOLO ROTURA 

%para poder correr el codigo de roller, primero se debe correr este porque
%se guarda un archivo .mat (ver linea 167) del codigo de ROTURAFinaaaal

clear; close all; clc;

g   = 9.81;
rho = 1000;


gamma_b = 0.55;   
Bbj     = 1.0;    % coef Baldock 
eps_F   = 1e-8;   
eps_h   = 1e-4;   
Nquad   = 200;    
ylim_top = 1.6;   % tope eje Y derecho

%% 1) Cargar batimetría
b = load('REU2004bathy.txt');     
x_bathy = b(:,1);
z_bathy = b(:,2);
h = -z_bathy;     % profundidad positiva

% Solo h utiles y ordenar por x
valid   = h > eps_h & isfinite(x_bathy) & isfinite(h);
x_bathy = x_bathy(valid);
h       = h(valid);
[x_bathy, idx] = sort(x_bathy(:));
h = h(idx);

%% 2) Casos experimentales
casos = { ...
    struct('archivo','R26.mat','T',4, 'col',[1 0 0],   'mk','o', 'mksz',70, 'nombre','R26'), ...
    struct('archivo','R40.mat','T',4, 'col',[0 0.6 0], 'mk','d', 'mksz',70, 'nombre','R40') ...
};

H_models = cell(size(casos));
H_exps   = cell(size(casos));
X_exps   = cell(size(casos));

%% 3) Loop principal: simulacion por caso
for idx = 1:numel(casos)
    caso = casos{idx};

    %Datos experimentales, extrae datos de gage 
    S  = load(caso.archivo);
    fn = fieldnames(S);
    R  = S.(fn{1});
    H_exp = R.LWF.H(:);
    x_exp = R.xreal(:);
    X_exps{idx} = x_exp;  H_exps{idx} = H_exp;

    %Parametros de ola 
    T     = caso.T;
    omega = 2*pi/T;
    fp    = 1/T;

    %Determina condicion de borde offshore, con el gage mas offshore
    h_at_gages = interp1(x_bathy, h, x_exp, 'linear', NaN);
    [~, i_off] = max(h_at_gages);
    H0 = H_exp(i_off);
    x0 = x_exp(i_off);
    h0 = h_at_gages(i_off);
    [~, i0] = min(abs(x_bathy - x0));
    fprintf('\nCaso %s: H0=%.3f m (x=%.2f m, h=%.2f m), T=%.2f s\n', ...
        caso.nombre, H0, x0, h0, T);

    %Dispersion a lo largo del perfil para obtener k
    k  = solve_k_vector(omega, g, h);
    c  = omega ./ k;
    kh = k.*h;
    cg = 0.5 .* c .* (1 + (2.*kh)./sinh(2.*kh));

    A = rho*g*Bbj*fp/4;  %Constante de disipacion Baldock 

    %Integracon onshore con rotura
    H_model     = NaN(size(h));
    H_model(i0) = H0;                           %condicio inicial offshore
    F = (1/8)*rho*g*H0^2 * max(cg(i0), eps_F);  % flujo inicial

    
    % Integración ON-SHORE 
    for j = i0+1:length(x_bathy)
        dx = x_bathy(j) - x_bathy(j-1);
        if dx <= 0, H_model(j) = H_model(j-1); continue; end

        Hb   = gamma_b * h(j); %altura rompiente teorica
        Hrms = max(H_model(j-1)/sqrt(2), 1e-6);  % Altura RMS local
        Hmax = max(6*Hrms, Hb + 1e-6); %rango integracion para la distribución de alturas

        %<ε>
        epsj = A * trapz_linspace(Hb, Hmax, ...
               @(H) H .* (2*H./Hrms.^2) .* exp(-(H./Hrms).^2), Nquad);

        F = F - epsj*dx; %balance de energia:

        % si el flujo es 0 cero se termina el calculo
        if F < eps_F
            H_model(j:end) = 0; break;
        end

        % energ´a y la altura correspondiente
        Ej    = F / max(cg(j), eps_F);
        H_try = sqrt(8*Ej/(rho*g));
        H_model(j) = min(H_try, Hb); % Limite, no puede superar Hb
    end










    % Integración OFF-SHORE (sin disipacion) 
    for j = i0-1:-1:1
        F_up = (1/8)*rho*g*H_model(j+1)^2 * max(cg(j+1), eps_F);
        Ej   = F_up / max(cg(j), eps_F);
        H_model(j) = sqrt(8*Ej/(rho*g));
    end

    %Plot individual 
    figure('Units','normalized','Position',[0.1 0.1 0.72 0.60]);
    yyaxis left
    plot(x_bathy, -h, 'k','LineWidth',1.5); hold on; ylabel('Batimetría z [m]');

    yyaxis right
    plot(x_bathy, H_model, '-', 'Color', caso.col, 'LineWidth',2); hold on;
    scatter(x_exp, H_exp, caso.mksz, 'Marker', caso.mk, ...
            'MarkerFaceColor', caso.col, 'MarkerEdgeColor','k');
    ylim([0 ylim_top]);

    legend({'Batimetría (elevación)', ...
            sprintf('H_{mod} (%s, rotura)', caso.nombre), ...
            sprintf('H_{exp} (%s)', caso.nombre)}, 'Location','northwest');

    xlabel('x [m]'); ylabel('Altura de ola H [m]');
    title(sprintf('Baldock (1998) – Solo rotura | %s, T=%.2f s, H_0=%.3f m', ...
          caso.nombre, T, H0));
    grid on; box on;

    H_models{idx} = H_model;
end

%% 4) Plot conjunto
figure('Units','normalized','Position',[0.1 0.1 0.72 0.60]);
yyaxis left
plot(x_bathy, -h, 'k','LineWidth',1.5); hold on; ylabel('Batimetría z [m]');

yyaxis right
plot(x_bathy, H_models{1}, '-', 'Color', casos{1}.col, 'LineWidth',2); hold on;
plot(x_bathy, H_models{2}, '-', 'Color', casos{2}.col, 'LineWidth',2);
scatter(X_exps{1}, H_exps{1}, casos{1}.mksz, 'Marker', casos{1}.mk, ...
        'MarkerFaceColor', casos{1}.col, 'MarkerEdgeColor','k');
scatter(X_exps{2}, H_exps{2}, casos{2}.mksz, 'Marker', casos{2}.mk, ...
        'MarkerFaceColor', casos{2}.col, 'MarkerEdgeColor','k');
ylim([0 ylim_top]);

legend({'Batimetría (elevación)', ...
        sprintf('H_{mod} %s', casos{1}.nombre), sprintf('H_{exp} (%s)', casos{1}.nombre), ...
        sprintf('H_{mod} %s', casos{2}.nombre), sprintf('H_{exp} (%s)', casos{2}.nombre)}, ...
        'Location','northwest');

xlabel('x [m]'); ylabel('Altura de ola H [m]');
title(sprintf('Baldock (1998) – Solo rotura | %s y %s', casos{1}.nombre, casos{2}.nombre));
set(gca,'FontSize',12,'LineWidth',1.2); grid on; box on;

%Guardar resultados para comparar con el modelo roller
save('H_model_rotura.mat', 'x_bathy', 'casos', 'H_models');
fprintf('Resultados de rotura guardados en H_model_rotura.mat\n');



%% Funciones 

% para la relación de dispersion ω² = gk tanh(kh), con Newton-raphson
function k = solve_k_scalar(omega, g, h)
    if ~(isfinite(h) && h > 0), k = NaN; return; end
    k = max(omega^2/g, 1e-6); tol = 1e-12; maxit = 80;
    for it = 1:maxit
        kh = k*h;
        f  = g*k*tanh(kh) - omega^2;
        df = g*tanh(kh) + g*k*h*(1./cosh(kh)).^2;
        step = f/df; k = k - step;
        if abs(step) < tol, break; end
    end
end

% resolver k para todo el perfil
function k = solve_k_vector(omega, g, hvec)
    k = nan(size(hvec));
    for i = 1:numel(hvec)
        k(i) = solve_k_scalar(omega, g, hvec(i));
    end
end

%Integral de <ε>
function I = trapz_linspace(a,b,fun,n)
    if nargin<4, n=200; end
    x = linspace(a,b,n); I = trapz(x, fun(x));
end
