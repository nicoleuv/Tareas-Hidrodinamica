%% ROTURA (Baldock 1998) + ROLLER (Lippmann96) 

%para poder correr el codigo de roller, primero se debe correr este porque
%se guarda un archivo .mat (ver linea 167) del codigo de ROTURAFinaaaal

clear; close all; clc;

g   = 9.81;       
rho = 1000;       

gamma_b  = 0.55;               
Bbj      = 1.0;                
alpha    = 12 * pi/180;        % ángulo del roller 
eps_h    = 1e-4;               
eps_F    = 1e-8;               
ylim_top = 1.6;                % límite superior del eje derecho

%% 1) Cargar batimetría
b = load('REU2004bathy.txt');
x_bathy = b(:,1);
z_bathy = b(:,2);
h = -z_bathy;    % profundidad positiva


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
    X_exps{idx} = x_exp;  
    H_exps{idx} = H_exp;

    %Parámetros de ola
    T     = caso.T;
    omega = 2*pi/T;
    fp    = 1/T;

    %Dispersion a lo largo del perfil para obtener k
    k  = solve_k_vector(omega, g, h);
    c  = omega ./ k;
    kh = k.*h;
    cg = 0.5 .* c .* (1 + (2.*kh)./sinh(2.*kh));

    %Determina condicion de borde offshore, con el gage mas offshore
    h_at_gages = interp1(x_bathy, h, x_exp, 'linear', NaN);
    [~, i_off] = max(h_at_gages);
    H0 = H_exp(i_off);
    x0 = x_exp(i_off);
    h0 = h_at_gages(i_off);
    [~, i0] = min(abs(x_bathy - x0));
    fprintf('\nCaso %s: H0=%.3f m (x=%.2f m, h=%.2f m), T=%.2f s\n', ...
        caso.nombre, H0, x0, h0, T);

    %Constantes 
    A  = rho*g*Bbj*fp/4;   %Constante de disipacion Baldock 
    Lr_fun = @(hh) max(hh./tan(alpha), 1e-6);  % Longitud de roller 

    %--- Inicialización ---
    H_model     = NaN(size(h));
    H_model(i0) = H0;
    Ew0 = (1/8)*rho*g*H0^2;           % Energía inicial
    Fw  = Ew0 * max(cg(i0), eps_F);   % flujo de ola
    Fr  = 0;                          % flujo de roller inicial 

    %=== Vectores para diagnóstico ===%
    Fw_prof = nan(size(h));
    Fr_prof = nan(size(h));
    eps_prof = nan(size(h));
    Fw_prof(i0)=Fw; Fr_prof(i0)=Fr; eps_prof(i0)=0;

    
    % Integración ON-SHORE (rotura + roller)
    for j = i0+1:length(x_bathy)
        dx = x_bathy(j) - x_bathy(j-1);
        if dx <= 0, H_model(j)=H_model(j-1); continue; end

        hj   = h(j);
        Hb   = gamma_b * hj;   % Altura límite de rompiente
        Hrms = max(H_model(j-1)/sqrt(2), 1e-6);  % Altura RMS previa

        % Disipación por rotura (Baldock)
        Nquad = 200;
        Hmax = max(6*Hrms, Hb + 1e-6);
        epsj = A * trapz_linspace(Hb, Hmax, ...
               @(H) H .* (2*H./Hrms.^2) .* exp(-(H./Hrms).^2), Nquad);

        %Roller (Lippmann)
        Lr   = Lr_fun(hj);
        dFrdx = epsj - Fr / Lr;     % derivada de Fr
        Fr = max(Fr + dFrdx*dx, 0); % Integracion de Fr

        %Acoplamiento energético dFw/dx = -ε - dFr/dx
        dFw = -epsj*dx - dFrdx*dx;
        Fw = max(Fw + dFw, 0);

        %Energía de ola y altura 
        Ej    = Fw / max(cg(j), eps_F);
        H_try = sqrt(8*Ej/(rho*g));
        H_model(j) = min(H_try, Hb);

        % Guardar diagnóstico
        Fw_prof(j)=Fw; Fr_prof(j)=Fr; eps_prof(j)=epsj;
    end

    % Diagnóstico de clipping
    %clipped = sum(abs(H_model - gamma_b*h) < 1e-6 & ~isnan(H_model));
    %pct_clip = 100 * clipped / sum(~isnan(H_model));
    %fprintf('%% puntos con H = Hb (clipping): %.1f%%\n', pct_clip);


    % Integración OFF-SHORE (sin disipación)
    for j = i0-1:-1:1
        F_up = (1/8)*rho*g*H_model(j+1)^2 * max(cg(j+1), eps_F);
        Ej   = F_up / max(cg(j), eps_F);
        H_model(j) = sqrt(8*Ej/(rho*g));
    end


    % Gráfico individual (rotura vs rotura+roller)
    figure('Units','normalized','Position',[0.1 0.1 0.72 0.60]);
    yyaxis left
    plot(x_bathy, -h, 'k','LineWidth',1.5); hold on;
    ylabel('Batimetría z [m]');

    %yyaxis right
    %plot(x_bathy, H_model, '-', 'Color', caso.col, 'LineWidth',2); hold on;
    %scatter(x_exp, H_exp, caso.mksz, 'Marker', caso.mk, ...
            %'MarkerFaceColor', caso.col, 'MarkerEdgeColor','k');
    %ylim([0 ylim_top]);









    yyaxis right
    %Cargar resultados del modelo sin roller (para comparar) 
    if exist('H_model_rotura.mat','file')
        if ~exist('H_sin_roller','var')
            data_rotura = load('H_model_rotura.mat');
            H_sin_roller = data_rotura.H_models;
            x_rotura = data_rotura.x_bathy;
            casos_rotura = data_rotura.casos;
            disp('Modelos sin roller cargados correctamente.');
        end
    else
        warning('No se encontró H_model_rotura.mat. No se graficará la línea de rotura.');
        H_sin_roller = [];
    end

    %Graficar linea punteada: modelo sin roller
    if ~isempty(H_sin_roller)
        nombre_actual = caso.nombre;
        nombres_rotura = cellfun(@(c) c.nombre, casos_rotura, 'UniformOutput', false);
        idx_rot = find(strcmp(nombre_actual, nombres_rotura));
        if ~isempty(idx_rot)
            plot(x_rotura, H_sin_roller{idx_rot}, '--', 'Color', caso.col, ...
            'LineWidth',1.8); % línea punteada = rotura sola
        end
    end

    %Graficar modelo con roller (línea entera) 
    plot(x_bathy, H_model, '-', 'Color', caso.col, 'LineWidth',2); hold on;

    %Puntos experimentales
    scatter(x_exp, H_exp, caso.mksz, 'Marker', caso.mk, ...
        'MarkerFaceColor', caso.col, 'MarkerEdgeColor','k');
    ylim([0 ylim_top]);













    %legend({'Batimetría (elevación)', ...
            %sprintf('H_{mod} (%s, rotura+roller)', caso.nombre), ...
            %sprintf('H_{exp} (%s)', caso.nombre)}, 'Location','northwest');



    legend({'Batimetría (elevación)', ...
        sprintf('H_{mod} (%s, rotura)', caso.nombre), ...
        sprintf('H_{mod} (%s, rotura+roller)', caso.nombre), ...
        sprintf('H_{exp} (%s)', caso.nombre)}, ...
        'Location','northwest');


    xlabel('x [m]'); ylabel('Altura de ola H [m]');
    title(sprintf('Baldock + Lippmann – %s | T=%.2f s, H_0=%.3f m', ...
          caso.nombre, T, H0));
    grid on; box on;


    %Diagnostico energetico: permite ver se reparte la energía entre ola 
    % (Fw), roller (Fr) y disipación (ε) a lo largo del perfil.
    figure('Name',['Flujos ', caso.nombre]);
    subplot(3,1,1); plot(x_bathy, Fw_prof, 'b', x_bathy, Fr_prof, 'r--','LineWidth',1.3);
    legend('F_w (ola)','F_r (roller)'); ylabel('Flujo [W/m]');
    subplot(3,1,2); plot(x_bathy, Fw_prof+Fr_prof,'k','LineWidth',1.3); ylabel('F_{total}');
    subplot(3,1,3); plot(x_bathy, eps_prof,'m','LineWidth',1.3);
    ylabel('\epsilon'); xlabel('x [m]');

    H_models{idx} = H_model;
end

%% 4) Gráfico conjunto
%figure('Units','normalized','Position',[0.1 0.1 0.72 0.60]);
%yyaxis left
%plot(x_bathy, -h, 'k','LineWidth',1.5); hold on; ylabel('Batimetría z [m]');

%yyaxis right
%plot(x_bathy, H_models{1}, '-', 'Color', casos{1}.col, 'LineWidth',2); hold on;
%plot(x_bathy, H_models{2}, '-', 'Color', casos{2}.col, 'LineWidth',2);
%scatter(X_exps{1}, H_exps{1}, casos{1}.mksz, 'Marker', casos{1}.mk, ...
        %'MarkerFaceColor', casos{1}.col, 'MarkerEdgeColor','k');
%scatter(X_exps{2}, H_exps{2}, casos{2}.mksz, 'Marker', casos{2}.mk, ...
       % 'MarkerFaceColor', casos{2}.col, 'MarkerEdgeColor','k');
%ylim([0 ylim_top]);

%legend({'Batimetría (elevación)', ...
        %sprintf('H_{mod} %s (rotura+roller)', casos{1}.nombre), sprintf('H_{exp} (%s)', casos{1}.nombre), ...
        %sprintf('H_{mod} %s (rotura+roller)', casos{2}.nombre), sprintf('H_{exp} (%s)', casos{2}.nombre)}, ...
        %'Location','northwest');



%% 4) Gráfico conjunto (rotura vs rotura+roller)
% Muestra ambos casos (R26 y R40) con y sin roller en una misma figura.
figure('Units','normalized','Position',[0.1 0.1 0.72 0.60]);
yyaxis left
plot(x_bathy, -h, 'k','LineWidth',1.5); hold on;
ylabel('Batimetría z [m]');

yyaxis right

%resultados sin roller
if exist('H_model_rotura.mat','file')
    data_rotura = load('H_model_rotura.mat');
    H_sin_roller = data_rotura.H_models;
    x_rotura = data_rotura.x_bathy;
    casos_rotura = data_rotura.casos;
else
    warning('No se encontró H_model_rotura.mat. El conjunto mostrará solo las curvas con roller.');
    H_sin_roller = [];
end

%dibujar cada caso 
for idx = 1:numel(casos)
    caso = casos{idx};
    col = caso.col;
    nombre = caso.nombre;

    %Modelo sin roller (línea punteada)
    if ~isempty(H_sin_roller)
        nombres_rotura = cellfun(@(c) c.nombre, casos_rotura, 'UniformOutput', false);
        idx_rot = find(strcmp(nombre, nombres_rotura));
        if ~isempty(idx_rot)
            plot(x_rotura, H_sin_roller{idx_rot}, '--', 'Color', col, 'LineWidth',1.8);
        end
    end

    %Modelo con roller (línea entera)
    plot(x_bathy, H_models{idx}, '-', 'Color', col, 'LineWidth',2); hold on;

    %Datos experimentales
    scatter(X_exps{idx}, H_exps{idx}, caso.mksz, 'Marker', caso.mk, ...
        'MarkerFaceColor', col, 'MarkerEdgeColor','k');
end

ylim([0 ylim_top]);

legend({'Batimetría (elevación)', ...
        'H_{mod} R26 (rotura)', 'H_{mod} R26 (rotura+roller)', 'H_{exp} (R26)', ...
        'H_{mod} R40 (rotura)', 'H_{mod} R40 (rotura+roller)', 'H_{exp} (R40)'}, ...
        'Location','northwest');

xlabel('x [m]'); ylabel('Altura de ola H [m]');
title('Baldock (1998) + Lippmann et al. (1996): Comparación rotura vs rotura+roller');
set(gca,'FontSize',12,'LineWidth',1.2); grid on; box on;






xlabel('x [m]'); ylabel('Altura de ola H [m]');
title(sprintf('Baldock + Lippmann – %s y %s (rotura + roller)', ...
       casos{1}.nombre, casos{2}.nombre));
set(gca,'FontSize',12,'LineWidth',1.2); grid on; box on;


%% FUNCIONES 

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
    x = linspace(a,b,n);
    I = trapz(x, fun(x));
end

