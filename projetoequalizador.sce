// FILTRAGEM ADAPTATIVA
// UFPE, RECIFE  xx/02/2019
// ALUNOS: GABRIEL VICTOR, RICARDO DANTAS, MATHEUS ROMA

// PRIMEIRO PROJETO DA DISCIPLINA: EQUALIZADOR


// (1)
//Definição dos filtros do canal e calculo da matriz Rx.
rand("normal");

exec ("adaptest.sce");


//preambulo sugerido
sn = [1 -1 1 -1 1 1 -1 -1]';

//parametros que podem ser mudados a escolha do usuario

c =5;   //constante utilizada para achar um coeficiente de aprendizagem adequado. C tem que ser inteiro positivo.

b1 =1;   //constante utilizado para plotagem de graficos, parametro da função plotxd. b pode assumir 3 valores, 1, 2 e 3.

n = 8;	// parametro utilizado para gerar um sinal bpsk de tamanho n.

k = 100;	// parametro utilizado para aumentar o tamanho da amostragem do sinal do item 3.

numero = 1000; // parametro utilizado para aumentar o tamanho da amostragem do sinal do item 4.

m = 0.01; // coeficiente de escala do ruido branco.


//Filtros sugeridos para a realização do projeto

H1 = [1 1/2];
H2 = [1 3/4 -1/4];


N1 = 1; N2 = 2; N3 = 3;

ni = rand(n, 1);

sn1 = filter(H1, 1, sn) + ni;
sn2 = filter (H2, 1, sn) + ni;

dn1 = filter ([0 1], 1, sn);
dn2 = filter ([0 0 1], 1, sn);
dn3 = filter([0 0 0 1], 1, sn); 

//Computação da Matriz Rx
// mi para H1: 0 < mii_1 < 1/lambdamaxi_1, o valor de mii_1 será escolhido da seguinte forma: mi = c/((2^(c+1))*lambdamaxi_1


//para H1

[Rx1_1, Rx2_1, Rx3_1, P1_1, P2_1, P3_1] = comat3(sn1, dn1, dn2, dn3, 1, 2, 3);


[mi1_1, No1_1] = copt(Rx1_1, c);
[mi2_1, No2_1] = copt(Rx2_1, c);
[mi3_1, No3_1] = copt(Rx3_1, c);

//mi1_1 = 0.5517;

//Para H2
// mi para H2 0 < mi2_i < 1/lambdamaxi_2,o valor de mii_2 será escolhido da seguinte forma: mi = c/((2^(c+1))*lambdamaxi_2


//para H1

[Rx1_2, Rx2_2, Rx3_2, P1_2, P2_2, P3_2] = comat3(sn2, dn1, dn2, dn3, 1, 2, 3);


[mi1_2, No1_2] = copt(Rx1_2, c);
[mi2_2, No2_2] = copt(Rx2_2, c);
[mi3_2, No3_2] = copt(Rx3_2, c);



//(2)
//Definir um preambulo sn de tamanho n;
//Dado um n para criar um vetor de tamanho n com probabilidade

s1= (-1).^(grand(2*ceil(No1_1)*k, 1, "bin", 1, 0.5));
s2= (-1).^(grand(2*ceil(No2_1)*k, 1, "bin", 1, 0.5));
s3= (-1).^(grand(2*ceil(No3_1)*k, 1, "bin", 1, 0.5));
s4= (-1).^(grand(2*ceil(No1_2)*k, 1, "bin", 1, 0.5));
s5= (-1).^(grand(2*ceil(No2_2)*k, 1, "bin", 1, 0.5));
s6 =(-1).^(grand(2*ceil(No3_2)*k, 1, "bin", 1, 0.5));


d1_1 = filter ([0 1], 1, s1);
d2_1 = filter ([0 0 1], 1, s2);
d3_1 = filter ([0 0 0 1], 1, s3);

d1_2 = filter ([0 1], 1, s4);
d2_2 = filter ([0 0 1], 1, s5);
d3_2 = filter ([0 0 0 1], 1, s6);

//para H1

s1_1 = filter(H1, 1, s1) + m*rand(k*2*ceil(No1_1), 1);
s2_1 = filter(H1, 1, s2) + m*rand(k*2*ceil(No2_1), 1);
s3_1 = filter(H1, 1, s3) + m*rand(k*2*ceil(No3_1), 1);

//para H2
s1_2 = filter(H2, 1, s4) + m*rand(k*2*ceil(No1_2), 1);
s2_2 = filter(H2, 1, s5) + m*rand(k*2*ceil(No2_2), 1);
s3_2 = filter(H2, 1, s6) + m*rand(k*2*ceil(No3_2), 1);


// (3) Metodo do gradiente estimado foi utilizado o lms para estimar o gradiente do erro.


[w1_1, e1_1, y1_1] = metgradi(s1_1, d1_1, mi1_1, Rx1_1, P1_1, 1);
[w2_1, e2_1, y2_1] = metgradi(s2_1, d2_1, mi2_1, Rx2_1, P2_1, 2);
[w3_1, e3_1, y3_1] = metgradi(s3_1, d3_1, mi3_1, Rx3_1, P3_1, 3);

plotxd(w1_1, w2_1, w3_1, e1_1, e2_1, e3_1, b1);
////para H2

//s2 = filter (H2, 1, s);

[w1_2, e1_2, y1_2] = metgradi(s1_2, d1_2, mi1_2, Rx1_2, P1_2, 1);
[w2_2, e2_2, y2_2] = metgradi(s2_2, d2_2, mi2_2, Rx2_2, P2_2, 2);
[w3_2, e3_2, y3_2] = metgradi(s3_2, d3_2, mi3_2, Rx3_2, P3_2, 3);

//plotxd(w1_2, w2_2, w3_2, e1_2, e2_2, e3_2, b1);


//(4)
// (Metodo Quasi Newton)
// Para Rx teorica calculada


Rxo1 = [2.25 0.5 ; 0.5 2.25];
Rxo2 = [2.25 0.5 0; 0.5 2.25 0.5 ;0 0.5 2.25];
Rxo3 = [2.25 0.5 0 0; 0.5 2.25 0.5 0; 0 0.5 2.25 0.5 ; 0 0 0.5 2.25];


// filtro ótimo teorico
// para H1

wo1_1 = inv(Rxo1)*P1_1; 
wo2_1 = inv(Rxo2)*P2_1;
wo3_1 = inv(Rxo3)*P3_1;

// para H2

wo1_2 = inv(Rxo1)*P1_2; 
wo2_2 = inv(Rxo2)*P2_2;
wo3_2 = inv(Rxo3)*P3_2;

// filtro ótimo utilizando o metodo quasi-newton

//para H1

[wqn1_1, eqn1_1, yqn1_1]= quasenw(s1_1, d1_1, mi1_1,Rxo1, 1);
[wqn2_1, eqn2_1, yqn2_1]= quasenw(s2_1, d2_1, mi2_1,Rxo2, 2);
[wqn3_1, eqn3_1, yqn3_1]= quasenw(s3_1, d3_1, mi3_1,Rxo3, 3);

//plotxd(wqn1_1, wqn2_1, wqn3_1, eqn1_1, eqn2_1, eqn3_1, 1);
//Para H2

[wqn1_2, eqn1_2, yqn1_2]=quasenw(s1_2, d1_2, mi1_2,Rxo1, 1);
[wqn3_2, eqn2_2, yqn2_2]=quasenw(s2_2, d2_2, mi2_2,Rxo2, 2);
[wqn3_2, eqn3_2, yqn3_2]=quasenw(s3_2, d3_2, mi3_2,Rxo3, 3);

//plotxd(wqn1_2, wqn2_2, wqn3_2, eqn1_2, eqn2_2, eqn3_2, 1);

// (5) COMPARAÇÃO PARA UM VETOR DE N BITS


preambulo = [1 -1 1 -1 1 1 -1 -1];
A = (-1).^(grand(numero , 1, "bin", 1, 0.5))
s_i5 = [preambulo A'];

d1_i5 = filter ([0 1], 1, s_i5);

d2_i5 = filter ([0 0 1], 1, s_i5);

d3_i5 = filter ([0 0 0 1], 1, s_i5);

mic = 1/(4*2.75) ;


// Para H1

s1_i5 = filter(H1, 1, s_i5) + m*rand(length(s_i5), 1)';

[wqnw1_1, e1_i5_1, grade1_1] = quasenw(s1_i5, d1_i5, mic ,Rxo1, 1);
[wqnw2_1, e2_i5_1, grade2_1] = quasenw(s1_i5, d2_i5, mic ,Rxo2, 2);
[wqnw3_1, e3_i5_1, grade3_1] = quasenw(s1_i5, d3_i5, mic ,Rxo3, 3);

//plotxd(wqnw1_1, wqnw2_1, wqnw3_1, e1_i5_1, e2_i5_1, e2_i5_1 , b1)

//para H2

s2_i5 = filter(H2, 1, s_i5) + rand(length(s_i5), 1)';
[wqnw1_2, e1_i5_2, grade1_2] = quasenw(s2_i5, d1_i5, mic ,Rxo1, 1);
[wqnw2_2, e2_i5_2, grade2_2] = quasenw(s2_i5, d2_i5, mic ,Rxo2, 2);
[wqnw3_2, e3_i5_2, grade3_2] = quasenw(s2_i5, d3_i5, mic ,Rxo3, 3);

//plotxd(wqnw1_2, wqnw2_2, wqnw3_2, e1_i5_2, e2_i5_2, e3_i5_2 , b1)
