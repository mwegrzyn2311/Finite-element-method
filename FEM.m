function FEM(n, c)
%Funkcja rozwi�zuje r�wnanie typu u''+u = c na przedziale [0,2]
%Z lewym warunkiem zerowym Dirichleta i prawym zerowym Cauchy'ego

%Dzielimy przedzia� (0,2) na n kawa�k�w o szeroko�ci h
h = 2/n;

%Tworzymy macierz potencja��w biliniowych B(u,v)...
matrix = zeros(n+1 ,n+1);
%...I wype�niamy j� warto�ciami
for i = 1:n+1
    for j = 1: n+1
        matrix(j,i) = B(i, j, h);
    end
end
%Tworzymy wektor dla potencja�u liniowego L(v)
L = zeros(n+1,1);
L(1)=(c/2)*(4/n^2);
for k = 2:n
    L(k) = c*(4/n^2);
end
L(n+1)=(c/2)*(4/n^2);
%Obliczamy uk�ad r�wna� otrzymuj�c wektor niewiadomych, z kt�rych
%skorzystamy w celu wyznaczenia przybli�enia u
result = matrix\L;

%I cz��, w kt�rej wyznaczamy pary (y,x) punkt�w na wykresie
m = 10*n;
h2 = 2/m;
x = 0: h2: 2;
y = 0: h2: 2;
for k = 1:m+1
    y(k) = u(x(k), n, h, result);
end

%Na koniec rysujemy wykres
plot(x,y);
end

%Funkcja u (Przybli�one rozwi�zanie zadania)
function u = u(x, n, h, result)
u = 0;
for k = 1:n+1
    u = u + result(k)*e(k,h,x);
end
end

%Potencja� biliniowy
function B = B(i, j, h)
%u = e(i)
%v = e(j)

B = (-e(i,h,2)*e(j,h,2) + myIntegral(i,j,h,0,2));
end

%Funkcja, kt�ra oblicza ca�k� z f, kt�ra jest funkcj� podca�kow� w
%potencjale biliniowym
function myIntegral = myIntegral(i,j,h,from,to)
res = 0;
%Z wi�kszym m zwi�ksza si� dok�adno�� wyniku, co pozwala uzyskiwa� dok�adny
%wykres dla wi�kszych n, ale zwi�ksza si� r�wnie� czas egzekucji programu
m = 12354;
g = (to - from)/m;
for k = 1:m
    res = res + f(i,j,h,g*k)*g;
end
myIntegral = res;
end

%Funkcja podca�kowa
function f = f(i, j, h, x)
f = de(i,h,x)*de(j,h,x) - e(i,h,x)*e(j,h,x);
end

%Funkcje bazowe do aproksymacji u
function e = e(k, h, x)
if (x < 0 || x < h*(k-2))
    e = 0;
elseif (x < h*(k-1))
    e = h*(2-k)+x;
elseif (x < h*k)
    e = h*k-x;
else
e = 0;
end
end

%Pochodne funkcji bazowych
function de = de(k, h, x)
if (x < 0 || x < h*(k-2))
    de = 0;
elseif (x < h*(k-1))
    de = 1;
elseif (x < h*k)
    de = -1;
else
    de = 0;
end
end