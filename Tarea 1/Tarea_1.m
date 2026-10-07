%% Ejercicio 1

%% Ejercicio 2
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

%% Ejercicio 5

%% Ejercicio 6

%% Ejercicio 7

%% Ejercicio 8

%% Ejercicio 9

%% Ejercicio 10
% Inicializa los valores constantes.
t_med  = [0 0.2 0.4 0.6 0.8]';
Ti_med = [37 36.72 36.41 36.12 35.90]';
Ta_med = 21; T0_med = 37;

% Defina la función dada en el enunciado
f_med = @(k) Ta_med + (T0_med - Ta_med)*exp(-k*t_med) - Ti_med;
f_med_sum = @(k) sum(f_med(k)^.2);
df_med = @(k) sum(2*f_med(k).*(-(T0_med-Ta_med)*t_med.*exp(-k*t_med)));
df2_med = @(k) sum(2*((T0_med-Ta_med)*t_med.*exp(-k*t_med)).^2 + 2*f_med(k).*((T0_med-Ta_med)*t_med.^2.*exp(-k*t_med)));

% Defina la función g, que cumple el teorema de punto fijo
g = @(k) k - (1/500)*df(k);

% Aplique el método de punto fijo
pto_fijo_med = punto_fijo(g, 0.1, 10^-6);

%% Ejercicio 11

%% Ejercicio 12

%% Ejercicio 13

%% Ejercicio 14

%% Ejercicio 15

%% FUNCIONES
function[x, res] = grad_conj(A, b, x0, M, tol)
% Funcion del metodo iterativo de gradiente conjugado con precondicionador 
% Aplica el metodo de gradiente conjugado para aproximar la solucion x de
% la ecuacion Ax = b.
% Inputs: 
%         A -- matriz nxn, simetrica y definida positiva.
%         b -- vector de constantes.
%        x0 -- aproximacion inicial de la solucion.
%         M -- matriz que precondiciona el sistema.
%       tol -- tolerancia de convergencia
% Outputs: 
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
% Inputs: 
%         A -- matriz estrictamente diagonal dominante de dimensión n.
%         b -- vector de constantes.
%        x0 -- aproximacion inicial de la solucion.
%       tol -- tolerancia de convergencia.
% Outputs: 
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

function[x, ite, res] = punto_fijo(g, x0, tol)
% Funcion del metodo de punto fijo.
% Aproxima la solucion de g(x) = 0 por medio del metodo de punto fijo.
% Inputs: 
%         g  -- función que cumple el teorema de punto fijo.
%        x0  -- aproximación incial.
%       tol  -- tolerancia de convergencia.
% Outputs: 
%          x -- aproximacion final a la raiz de f.
%        ite -- cantidad de iteraciones realizadas
%        res -- matriz que contiene cada iteracion de x y su error.

ite = 1;
er = tol + 1;
res = [x0, er];

% Defina el ciclo
while er > tol && ite < 100
    t = x0;
    x = g(x0);
    er = abs(t - x);
    ite = ite + 1;
    res = [res; x, er];
    x0 = x;
end

end