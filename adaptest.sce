// biblioteca de funções 


function [Rx1, Rx2, Rx3, P1, P2, P3] = comat3(x,d1, d2, d3, N1, N2, N3)
//N1
        px1=xcorr(x,N1, "biased");
            for i=1:(N1+1)
                Rx1(:,i)=px1((N1+2-i):(2*(N1+1)-i))
            end
//N2
        px2=xcorr(x,N2, "biased");
            for j=1:(N2+1)
                Rx2(:,j)=px2((N2+2-j):(2*(N2+1)-j))
            end
//N3
        px3=xcorr(x,N3, "biased");
            for k=1:(N3+1)
                Rx3(:,k)=px3((N3+2-k):(2*(N3+1)-k))
            end
//N1
        pdx1=xcorr(x,d1, N1, "biased");
        P1 = pdx1(N1+1:2*N1+1);
// N2
        pdx2=xcorr(x,d2, N2, "biased");
        P2 = pdx2(N2+1:2*N2+1);
// N3
        pdx3=xcorr(x,d3, N3, "biased");
        P3 = pdx3(N3+1:2*N3+1);
endfunction

//Calculo dos coeficientes de aprendizagem ótimo e tempo de convergencia em amostras.
function [mi, No] = copt(Rx, c) 
    
        aux = spec(Rx);
        
        lambdamax = max(aux);
        
        mi = c/((2^c)*lambdamax); 
        
        minimo = min(aux);
        
        No = 1/(4*mi*minimo) ;
        
endfunction

//Metodo de quase-newton para busca dos coeficientes ótimos do filtro

function [w, e, y]=quasenw(x, d, mi,R, N)
    xn=zeros(N+1,1);
    w = zeros(N+1,1);
    grad = zeros(N+1,N+1);
    Ri=inv(R); 
    for i=1:length(x)
        xn = [x(i);xn(1:N)];
        y(i) = w(:,i)'*xn;
        e(i) = d(i)-y(i);
        grad = [-2*e(i)*xn, grad(:,1:N)];
        sgrad = sum(grad,"c")/(N+1);
        w = [w w(:,i)-mi*Ri*sgrad];
    end
endfunction


//Metodo de quase-newton para busca dos coeficientes ótimos do filtro

function [wn, e, y] = metgradi(x, d, mi, Rx, P, N) 
    Sx = inv(Rx);
    wn = zeros(N+1, 1);
    xn = zeros(N+1,1);
    xn(1) = x(1);
    y = wn'*xn;
    e = d(1)-y;
    // calculo do grad estimado
    grade = -2*(P - Rx*wn);
    wn = [wn, wn-mi*Sx*grade];

for n=2:length(d)
    xn=[x(n);xn(1:N)];
    y(n) = wn(:,n)'*xn;
    e(n) = d(n)-y(n);   
    // calculo do grad estimado
    grade = -2*(P - Rx*wn(:,n));
    wn = [wn wn(:,n)-mi*Sx*grade];
    
end

endfunction


// Função de auxilio para plotar grafico para análise.

function [] = plotxd(wn1, wn2, wn3, en1, en2, en3 , b)
    if b == 1 then
        subplot(2,1,1);
        plot2d(wn1');
        title("Convergencia dos coeficientes do Filtro para  N = 1")
        xlabel("n (Amostras)");
        ylabel("Wn", "rotation", 0);
        
        subplot(2,1,2);
        plot2d(en1');
        title("Erro")
        xlabel("n (Amostras)");
        ylabel("e(n)", "rotation", 0);
    end
     if b == 2 then
        subplot(2,1,1);
        plot2d(wn2');
        title("Convergencia dos coeficientes do Filtro para  N = 2")
        xlabel("n (Amostras)");
        ylabel("Wn", "rotation", 0);
        
        subplot(2,1,2);
        plot2d(en2');
        title("Erro")
        xlabel("n (Amostras)");
        ylabel("e(n)", "rotation", 0);
    end
     if b == 3 then
        subplot(2,1,1);
        plot2d(wn3');
        title("Convergencia dos coeficientes do Filtro para  N = 3")
        xlabel("n (Amostras)");
        ylabel("Wn", "rotation", 0);
        
        subplot(2,1,2);
        plot2d(en3');
        title("Erro")
        xlabel("n (Amostras)");
        ylabel("e(n)", "rotation", 0);
    end
endfunction
