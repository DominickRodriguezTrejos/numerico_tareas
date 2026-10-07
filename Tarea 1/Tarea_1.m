%% Ejercicio 1

function [M, cont] = division_multiples_secciones(f, a, b, n, tol)
    % DIVISION_MULTIPLES_SECCIONES Busca una raíz de f(x) = 0 dividiendo el
    % intervalo [a, b] en 'n' subintervalos por iteración.
    %
    % Entradas:
    %   f    - Function handle de la función objetivo, ej: @(x) exp(x) - 2*x - 1
    %   a    - Límite inferior del intervalo inicial a_0
    %   b    - Límite superior del intervalo inicial b_0
    %   n    - Número de subintervalos
    %   tol  - Tolerancia para detener el programa
    %
    % Salidas:
    %   M    - Matriz de resultados (Col 1: x_k, Col 2: e_k, Col 3: a_k, Col 4: b_k)
    %   cont - Número total de iteraciones
    
    % a. Verificación inicial de existencia de raíz (Teorema de Bolzano)
    if (f(a) * f(b) > 0)
        warning('El método puede no funcionar ya que f(a)*f(b) > 0.');
        M = [];
        cont = 0;
        return;
    end
    
    cont = 0;
    er = abs(b - a);
    
    % Bucle principal hasta alcanzar la tolerancia
    while er > tol
        h = (b - a) / n;
        a_nuevo = a;
        b_nuevo = b;
    
        % b. Recorrer los 'n' subintervalos para hallar en cuál está la raíz
        for i = 1:n
            s_i_prev = a + (i - 1) * h; % Puntos_s_{i-1}
            s_i      = a + i * h;       % Puntos_s_{i}
    
            % Si hay un cero en los extremos exactos del subintervalo
            if f(s_i) == 0
                a_nuevo = s_i;
                b_nuevo = s_i;
                break;
            end
    
            % Búsqueda del cambio de signo: f(s_{i-1}) * f(s_i) < 0
            if f(s_i_prev) * f(s_i) < 0
                a_nuevo = s_i_prev;
                b_nuevo = s_i;
                break;
            end
        end
    
        % Actualizar a y b con el nuevo intervalo hallado [a_{k+1}, b_{k+1}]
        a = a_nuevo;
        b = b_nuevo;
    
        % c, d, e. Asignación a las columnas de la matriz M
        x_k = (a + b) / 2;      % Columna 1: Aproximación x_k
        er  = abs(b - a);        % Columna 2: Error e_k = |a_{k+1} - b_{k+1}|
    
        cont = cont + 1;
        M(cont, 1) = x_k;
        M(cont, 2) = er;
        M(cont, 3) = a;         
        M(cont, 4) = b;         
    
        if a == b
            break;
        end
    end

end

f = @(x) exp(x) - 2*x - 1;
a = 1;
b = 2;
n = 4;        
tol = 1e-6;   

[M, cont] = division_multiples_secciones(f, a, b, n, tol);

disp('Matriz M [x_k, e_k, a_k, b_k]:');
disp(M);
fprintf('Total de iteraciones: %d\n', cont);

%% Ejercicio 2
disp("Ejercicio 2")
% Carga la matriz
load('bcsstk06.mat')
bcsstk = Problem.A;

% Defina las matrices de precondicionamento
D_probA = diag(diag(bcsstk)); L_probA = tril(bcsstk, -1); U_probA = triu(bcsstk, 1);
M1 = D_probA; M2 = (D_probA + L_probA)*(D_probA \ (D_probA + U_probA));

% Defina el vector de constantes b y x0
b = bcsstk*(1:length(bcsstk))';
x0 = zeros(length(bcsstk),1);

% Aproxima la solución del sistema Ax = b con la función grad_conj
[x_precond_M1, ers_precond_M1] = grad_conj(bcsstk, b, x0, M1, 10^-6);
[x_precond_M2, ers_precond_M2] = grad_conj(bcsstk, b, x0, M2, 10^-6);

% Grafica los errores 
figure
semilogy(ers_precond_M1, 'b-', 'LineWidth', 1.5); hold on
semilogy(ers_precond_M2, 'r-', 'LineWidth', 1.5);
xlabel('Iteración k'); ylabel('||b - A x_k||_2');
legend('M_1 = D', 'M_2 = (D+L)D^{-1}(D+U)'); grid on
title('Gradiente conjugado precondicionado, bcsstk06')

% Aproxima la solución del sistema Ax = b con gauss_seidel
[x_gauss, ers_gauss] = gauss_seidel(bcsstk, b, x0, 10^-6);
% Grafica los errores
figure
semilogy(ers_gauss, 'b-', 'LineWidth', 1.5); hold on;
xlabel('Iteración k'); ylabel('||b - A x_k||_2'); grid on
title('Gauss Seidel, bcsstk06')

%% Ejercicio 3

%% Ejercicio 4
disp("Ejercicio 4.d")
aproximarRaiz(1, 1e-6)
aproximarRaiz(1, 1e-4)
aproximarRaiz(1, 1e-10)

%% Ejercicio 5

%% Ejercicio 6

%% Ejercicio 7

% d)

function M = metodoNewtonModificadoMultiples(p_coef, df_coef, ddf_coef, x0, tol, iterMax)
% METODONEWTONMODIFICADORAIZMULTIPLE Implementa el algoritmo de Newton para encontrar raices multiples de un polinomio
%
% Entradas:
%   p_coef   - Vector de coeficientes del polinomio p(x)
%   df_coef  - Vector de coeficientes de la primera derivada p'(x)
%   ddf_coef - Vector de coeficientes de la segunda derivada p''(x)
%   x0       - Aproximacion inicial
%   tol      - Tolerancia 
%   iterMax  - Numero maximo de iteraciones 
%
% Salida:
%   M        - Matriz de resultados

k = 1;
xk_ant = x0;

px   = polyval(p_coef, xk_ant);
dpx  = polyval(df_coef, xk_ant);
ddpx = polyval(ddf_coef, xk_ant);

den = (dpx)^2 - px * ddpx;

if abs(den) < eps
    error('Error: El denominador es 0.');
end

M = zeros(iterMax + 1, 3);
M(1, :) = [0, xk_ant, inf];

error_abs = tol + 1;

while error_abs > tol && k <= iterMax
    xk_act = xk_ant - (px * dpx) / den;
    error_abs = abs(xk_act - xk_ant);

    M(k + 1, :) = [k, xk_act, error_abs];

    xk_ant = xk_act;
    px   = polyval(p_coef, xk_ant);
    dpx  = polyval(df_coef, xk_ant);
    ddpx = polyval(ddf_coef, xk_ant);

    den = (dpx)^2 - px * ddpx;

    if abs(den) < eps
        fprintf('Error: El denominador se anulo en la iteracion %d.\n', k);
        M = M(1:k+1, :);
        return;
    end

    k = k + 1;
end

M = M(1:k, :);
end

%% Ejercicio 8
disp("Ejercicio 8")
% (a)
f = @(x) 816 * x^3 - 3835 * x^2 + 6000 * x - 3125;
x0 = fzero(f, [1.4,1.7]);
x1 = fzero(f, [1.4,1.5]);
x2 = fzero(f, [1.5,1.6]);


disp("(a) Las raices de la funcion son:")
disp(x0)
disp(x1)
disp(x2)

% (b)
x_cero = linspace(1.4,1.7,10000);
y = [];
df = @(x) 2448 * x^2 - 7670 * x + 6000;
for x = x_cero
    resul = newtonR(f, df, x, 1e-5, 100);
    y = [y resul(length(resul), 1)];
end

scatter(x_cero,y, 8, "filled");
xlabel("x_0")
ylabel("Valor de convergencia")

%% Ejercicio 9

% c)

function suc  = iter(f, x0, tol)
    % iter: Calcula la sucesion de aproximaciones mediante el Metodo de Halley
    % Entradas:
    %   f   - Vector de coeficientes del polinomio [a_n, a_{n-1}, ..., a_1, a_0]
    %   x0  - Valor inicial de la sucesion (x_0)
    %   tol - Tolerancia 
    %
    % Salida:
    %   suc - Vector columna con la sucesion de valores [x0, x1, x2, ...]

    df = polyder(f);
    ddf = polyder(df);
    suc = x0;
    xk = x0;
    error = inf;

    while error > tol

        fx = polyval(f, xk);
        dfx = polyval(df, xk);
        ddfx = polyval(ddf,xk);

        xk_nuevo = xk - (fx/dfx) *(1/(1 - (fx/dfx)* (ddfx/(2*dfx))));
        error = abs(xk_nuevo - xk);

        suc = [suc; xk_nuevo];
        xk = xk_nuevo;
    end
end

%d)

function[M] = metodoNewton(f,df,x,tol)
% funcion para aproximar la solucion de f(x) = 0 usando el
% metodo de punto fijo
% Entradas:       f --- funcion 
%                df --- derivada de f
%                 x --- aproximacion inicial
%               tol --- tolerancia para detener el ciclo
%            iteMax --- iteraciones maxima
% Salidas:        x --- aproximacion final a la raiz de f
    cont = 1; 
    er = tol + 1; M(1,1) = x; 
    M(1,2) = inf;
    while er > tol 
        t = x;
        q = df(x);
        if (abs(q)<tol)
            return;
        end
        x = x - f(x)/q;
        er = abs(t-x);
    
        M(cont+1,1) = x;
        M(cont+1,2) = er;
        cont = cont + 1; 
    end
end

clear; clc; close all;

f_poly = [816, -3835, 6000, -3125]; 
x0  = 1.6;                           
tol = 1e-6;                         

f_handle  = @(x) polyval(f_poly, x);
df_handle = @(x) polyval(polyder(f_poly), x);

suc_halley = iter(f_poly, x0, tol); 
M_newton   = metodoNewton(f_handle, df_handle, x0, tol);
suc_newton = M_newton(:, 1);

c = 3125 / 2000; 

err_halley = abs(suc_halley(:) - c);
err_newton = abs(suc_newton(:) - c);


figure;
semilogy(0:length(err_halley)-1, err_halley, '-o', 'LineWidth', 1.8, 'MarkerSize', 6);
hold on;
semilogy(0:length(err_newton)-1, err_newton, '-s', 'LineWidth', 1.8, 'MarkerSize', 6);
grid on;

title('Comparación del Error: Halley vs. Newton');
xlabel('Iteración (k)');
ylabel('Error |x_k - c|');
legend('Método de Halley (Cúbico)', 'Método de Newton (Cuadrático)', 'Location', 'northeast');

% e)
clear; clc; close all;

f_poly = [816, -3835, 6000, -3125]; 
tol    = 1e-6;                         
xx     = linspace(1.4, 1.7, 500); 


converged_val = zeros(size(xx));

for i = 1:length(xx)
    suc_halley = iter(f_poly, xx(i), tol);
    converged_val(i) = suc_halley(end); 
end

figure;
plot(xx, converged_val, 'LineWidth', 1.8, 'Color', [0, 0.4470, 0.7410]);
hold on;

yline(1.5625, '--r', 'Raíz c = 1.5625', 'LineWidth', 1.5, ...
    'FontSize', 10, 'Interpreter', 'none', 'LabelHorizontalAlignment', 'left');

grid on;

title('Intervalo de Atracción: Método de Halley');
xlabel('Valor Inicial (xx)');
ylabel('Valor de Convergencia');

xlim([1.4, 1.7]);
ylim([1.3, 1.8]);

%% Ejercicio 10
disp("Ejercicio 10")
% Inicializa los valores constantes.
t_med  = [0 0.2 0.4 0.6 0.8]';
Ti_med = [37 36.72 36.41 36.12 35.90]';
Ta_med = 21; T0_med = 37;

% Defina la función dada en el enunciado y encuentre sus derivadas
f_med = @(k) Ta_med + (T0_med - Ta_med)*exp(-k*t_med) - Ti_med;
f_med_sum = @(k) sum(f_med(k)^.2);
df_med = @(k) sum(2*f_med(k).*(-(T0_med - Ta_med)*t_med.*exp(-k*t_med)));
df2_med = @(k) sum(2*((T0_med-Ta_med)*t_med.*exp(-k*t_med)).^2 + 2*f_med(k).*((T0_med-Ta_med)*t_med.^2.*exp(-k*t_med)));

% Encuentra el punto critico aproximando k tal que f'(k) = 0
pto_crit_med = fzero(df_med, 0.01);
fprintf("El punto crítico de la función es aproximadamente %.4f\n", pto_crit_med)

% Evalua el punto critico en la segunda derivada y clasifica min/max
df2_med_kc = df2_med(pto_crit_med);
if df2_med_kc > 0
    disp("En el punto crítico de la función f(k) se tiene un mínimo")
else
    disp("En el punto crítico de la función f(k) se tiene un máximo")
end

% Defina los valores para el caso de la persona fallecida en San Pedro
Ta_muerto = 31;
T0_muerto = 37;
T_muerto = 34;

t_muerte = -log((T_muerto - Ta_muerto)/(T0_muerto - Ta_muerto))/pto_crit_med;
%% Ejercicio 11

% c)
function [L, U, P, Q] = eliminacionGaussianaPivoteoTotal(A)
    % ELIMINACIONGAUSSIANAPIVOTEOTOTAL
    % Calcula la factorización PAQ = LU usando pivoteo total
    
    m = size(A, 1);
    U = A;
    L = eye(m);
    P = eye(m);
    Q = eye(m);
    
    for k = 1 : m-1
        
        subMatriz = abs(U(k:m, k:m));
        [maxVal, idxRow] = max(subMatriz, [], 1);
        [~, s_rel] = max(maxVal);
        r_rel = idxRow(s_rel);
    
        r = r_rel + k - 1; 
        s = s_rel + k - 1; 
    
        if r ~= k
            U([k, r], k:m) = U([r, k], k:m);
            P([k, r], :)   = P([r, k], :);
            if k > 1
                L([k, r], 1:k-1) = L([r, k], 1:k-1);
            end
        end
    
    
        if s ~= k
            U(:, [k, s]) = U(:, [s, k]);
            Q(:, [k, s]) = Q(:, [s, k]);
        end
    
    
        for i = k+1 : m
            if U(k, k) == 0
                error('El pivote es cero. La matriz puede ser singular.');
            end
            factor = U(i, k) / U(k, k);
            L(i, k) = factor;
            U(i, k:m) = U(i, k:m) - factor * U(k, k:m);
        end
    end
end

%% Ejercicio 12
disp("Ejercicio 12.f")
F = @(x) [3*x(1) - cos(x(2)*x(3)) - 0.5;
    x(1)^2-81*(x(2)+0.1)^2+sin(x(3))+1.06;
    exp(-x(1)*x(2))+20*x(3)-(3-10*pi)/(3)];
J = @(x) [3, x(3)*sin(x(2)*x(3)), x(2)*sin(x(2)*x(3));
    2*x(1), -162*(x(2) + 0.1), cos(x(3));
    -x(2)*exp(-x(1)*x(2)), -x(1)*exp(-x(1)*x(2)), 20];
x0 = [0.1 , 0.1 , -0.1];

NewtonSistemasNoLineales(F, J, x0, 1e-4)
%% Ejercicio 13

%b.1)

function A = generarMatrizA(m)
% generarMatrizA: Genera la matriz m x m pedida
% Entrada:
%   m - Dimensión de la matriz 
% Salida:
%   A - Matriz m x m con 1s en la diagonal y última columna y -1s debajo de la diagonal.


A = -tril(ones(m), -1);
A = A + eye(m);
A(:, end) = 1;

end

% b.2)
clear; clc; close all;

m = 50;
A = generarMatrizA(m);

[L, U, P] = lu(A);

rho_A = max(abs(U(:))) / max(abs(A(:)));

norma_error_LU = norm(A - L*U, 2);

fprintf('=== PARTE (b.2) para m = 50 ===\n');
fprintf('Factor de crecimiento rho(A): %.4e\n', rho_A);
fprintf('Norma ||A - L*U||_2:          %.4e\n\n', norma_error_LU);

% b.3)

m_vals = 10:60;
num_m = length(m_vals);

err_LU = zeros(num_m, 1);
err_QR = zeros(num_m, 1);

for idx = 1:num_m

    m_curr = m_vals(idx);
    A_curr = generarMatrizA(m_curr);
    xExact = rand(m_curr, 1);
    b = A_curr * xExact;


    [L, U, P] = lu(A_curr);
    y_lu = L \ (P * b);
    xApr_LU = U \ y_lu;


    [Q, R] = qr(A_curr);
    xApr_QR = R \ (Q' * b);

    err_LU(idx) = norm(xExact - xApr_LU, 2);
    err_QR(idx) = norm(xExact - xApr_QR, 2);
end


figure;
semilogy(m_vals, err_LU, '-o', 'LineWidth', 1.8, 'MarkerSize', 5);
hold on;
semilogy(m_vals, err_QR, '-s', 'LineWidth', 1.8, 'MarkerSize', 5);
grid on;

title('Error ||xExact - xApr||_2 en función de m');
xlabel('Dimensión m');
ylabel('Error ||xExact - xApr||_2 (Escala Logarítmica)');
legend('Factorización LU', 'Factorización QR', 'Location', 'northwest');



%% Ejercicio 14

%% Ejercicio 15

% c)

function [M, p_sig] = metodo_muller(p, p0, p1, p2, tol, N)
    % metodo_muller Implementa el algoritmo del metodo de Muller
    %
    % Entradas:
    %   p   - Vector con los coeficientes del polinomio [a_n, ..., a_1, a_0]
    %   p0, p1, p2 - Tres puntos iniciales dados
    %   tol - Tolerancia
    %   N   - Cantidad maxima de iteraciones
    %
    % Salidas:
    %   M     - Matriz de resultados
    %   p_sig - Ultima aproximacion obtenida
    
    h1 = p1 - p0;
    h2 = p2 - p1;
    
    f0 = polyval(p, p0);
    f1 = polyval(p, p1);
    f2 = polyval(p, p2);
    
    delta1 = (f1 - f0) / h1;
    delta2 = (f2 - f1) / h2;
    d = (delta2 - delta1) / (h2 + h1);
    
    M = [];
    i = 3;
    
    while i <= N
        b = delta2 + h2 * d;
        D = sqrt(b^2 - 4 * f2 * d);
    
        if abs(b - D) < abs(b + D)
            E = b + D;
        else
            E = b - D;
        end
    
        h = -2 * f2 / E;
        p_sig = p2 + h;
    
        error_abs = abs(p_sig - p2);
    
    
        M(end+1, :) = [p_sig, error_abs];
    
    
        if error_abs < tol
            break;
        end
    
    
        p0 = p1;
        p1 = p2;
        p2 = p_sig;
    
        f0 = polyval(p, p0);
        f1 = polyval(p, p1);
        f2 = polyval(p, p2);
    
        h1 = p1 - p0;
        h2 = p2 - p1;
        delta1 = (f1 - f0) / h1;
        delta2 = (f2 - f1) / h2;
        d = (delta2 - delta1) / (h2 + h1);
    
        i = i + 1;
    end
   end

% d)
clear; clc; close all;
p = [16, -40, 5, 20, 6];

p0 = 0.5;
p1 = -0.5;
p2 = 0;

tol = 1e-12; 
N = 6;     

[M, x8] = metodo_muller(p, p0, p1, p2, tol, N);
disp('==================================================');
disp('   Aproximación (x_k)       |   Error Absoluto (|x_k - x_{k-1}|)');
disp('==================================================');
disp(M);

fprintf('\nResultado Final:\n');
fprintf('La aproximación x_8 es: %.6f + %.6fi\n', real(x8), imag(x8));


%% FUNCIONES
function [M, cont] = division_multiples_secciones(f, a, b, n, tol)
% DIVISION_MULTIPLES_SECCIONES Busca una raíz f(x) = 0 dividiendo el intervalo [a, b] en 'n' subintervalos por iteración.
%
% Entradas:
%   f    - funcion
%   a    - Límite inferior del intervalo inicial a_0
%   b    - Límite superior del intervalo inicial b_0
%   n    - Número de subintervalos
%   tol  - Tolerancia 
%
% Salidas:
%   M    - Matriz de resultados
%   cont - Numero total de iteraciones

% Verificación inicial de existencia de raíz (Teorema de Bolzano)
if (f(a) * f(b) > 0)
    warning('El método puede no funcionar ya que f(a)*f(b) > 0.');
    M = [];
    cont = 0;
    return;
end

cont = 0;
er = abs(b - a);

while er > tol

    h = (b - a) / n;
    a_nuevo = a;
    b_nuevo = b;

    for i = 1:n
        s_i_prev = a + (i - 1) * h; 
        s_i      = a + i * h;       

        if f(s_i) == 0
            a_nuevo = s_i;
            b_nuevo = s_i;
            break;
        end

        if f(s_i_prev) * f(s_i) < 0
            a_nuevo = s_i_prev;
            b_nuevo = s_i;
            break;
        end
    end

    a = a_nuevo;
    b = b_nuevo;

    x_k = (a + b) / 2;     
    er  = abs(b - a);       

    cont = cont + 1;
    M(cont, 1) = x_k;
    M(cont, 2) = er;
    M(cont, 3) = a;         
    M(cont, 4) = b;         

    if a == b
        break;
    end
end

function[x, res] = grad_conj(A, b, x0, M, tol)
% Funcion del metodo iterativo de gradiente conjugado con precondicionador 
% Aplica el metodo de gradiente conjugado para aproximar la solucion x de
% la ecuacion Ax = b.
% Entradas: 
%         A -- matriz nxn, simetrica y definida positiva.
%         b -- vector de constantes.
%        x0 -- aproximacion inicial de la solucion.
%         M -- matriz que precondiciona el sistema.
%       tol -- tolerancia de convergencia
% Salidas: 
%          x -- aproximacion final de la solucion del sistema.
%        res -- vector columna con los errores absolutos de cada iteraccion.

% Verifica que A es simetrica.
if ~issymmetric(A)
    error("La matriz A no es simetrica")
end
% A es definida positiva si todos sus valores propios son positivos.
% Verifica que A es definida positiva.
if any(eig(A) <= 0)
    error("La matriz A no es definida positiva")
end

% Defina
r = b - A*x0; z = M\r;
p = z; alpha = 1; x = x0; cont = 0;
res = norm(r);

% Contruya el ciclo
while norm(alpha*p) >= tol && cont < 10*length(b)
    alpha = (r'*z)/(p'*A*p);
    x = x + alpha*p;
    r_new = r - alpha*A*p; 
    z_new = M\r_new;
    beta = (z_new'*r_new)/(z'*r);
    p = z_new + beta*p;
    r = r_new; z = z_new;
    res = [res; norm(r)];
    cont = cont + 1;
end

end

function [x, res] = gauss_seidel(A, b, x0, tol)
% Funcion del metodo de Gauss Seidel 
% Aplica el metodo de Gauss Seidel para aproximar la solucion x de la ecuacion Ax = b.
% Entradas: 
%         A -- matriz estrictamente diagonal dominante de dimensión n.
%         b -- vector de constantes.
%        x0 -- aproximacion inicial de la solucion.
%       tol -- tolerancia de convergencia.
% Salidas: 
%         x -- aproximacion final de la solucion del sistema. 
%       res -- vector columna con los errores absolutos de cada iteraccion.

% Defina la descomposición de A
DL = tril(A); U  = triu(A, 1);
x  = x0;

% Revisa convergencia (norma de (L+D)^-1*U < 1)
if norm(DL\U, inf) <= 1
    error("No se cumple el supuesto de convergencia")
end

res = norm(b - A*x);
cont = 0;
% Construya el ciclo
while norm(b - A*x) >= tol && cont < 10*length(b)
    x_new  = DL\b - DL\U*x;        
    x   = x_new;
    res = [res; norm(b - A*x)];
    cont = cont + 1;
end

end

function resultado = aproximarRaiz(x0, tolerancia)
% Funcion para aproximar la raiz cubica de veintiuno con una sucesion
% especifica
% Entradas:         
%           x0 -- estimacion inicial
%   tolerancia -- precision deseada por el usuario
%
% Salidas:   
%    resultado -- la estimacion del valor real

if (x0 > 0 && tolerancia > 0)
    contador = 1;
    er = 1;
    while (contador < 200 && er > tolerancia)
        resultado = (20*x0 + 21/x0^2)/21;
        er = abs(resultado - x0);
        x0 = resultado;
        contador = contador + 1;
    end
else
    error("el valor inicial y la tolerancia deben ser positivos")
end

end

function [M] = newtonR(f, df, x, tol, iterMax)
    % Encuentra raices de una funcion por medio del metodo de Newton-Raphson
    % Entradas:       
    %                 f -- una funcion continua
    %                df -- la derivada de f
    %                 x -- la estimacion inicial de la raiz
    %               tol -- tolerancia entre estimaciones para detener 
    %                       el algoritmo
    %           iterMax -- cantidad maxima de iteraciones del algoritmo
    %
    % Salidas:        
    %                 M -- matriz con: estimacion y error entre estimaciones
    
    er = tol + 1;
    contador = 1;
    M = [x inf];
    while er > tol && contador <= iterMax
        q = df(x);
        if abs(q) < tol
            return
        end
        x1 = x - f(x) / q;
        er = abs(x1 - x);
        contador = contador + 1;
        x = x1;
        M = [M; [x er]]; % concatenar por filas
    end

end

function x2 = NewtonSistemasNoLineales(F, J, x0, tol)
    % Devuelve una aproximacion para la solucion de un sistema no lineal por
    % medio del metodo de Newton
    %
    % Entradas:   
    %             F -- vector con las funciones del sistema
    %             J -- la matriz jacobiana
    %            x0 -- aproximacion inicial
    %           tol -- tolerancia para el criterio de parada
    %
    % Salidas:    
    %             x -- aproximacion a la solucion
    
    x1 = x0(:);
    error = inf;
    contador = 1;
    while error > tol && contador < 1000
        Fx = F(x1);
        Jx = J(x1);
        y = -Jx \ Fx;
        x2 = x1 + y;
        error = norm(x2 - x1, inf);
        x1 = x2;
        contador = contador + 1;
    end

end