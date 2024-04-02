import inspect, re
import numpy as np
import matplotlib.pyplot as plt
from prettytable import PrettyTable as pt
from sympy import symbols, diff, exp, lambdify

def eq_population(capac, popi, taxa, popt):
    """Função para tratamento da forma da equação utilizada.

    Parâmetros:
        capac: capacidade máxima que um determinado ambiente suporta
        popi: população inicial
        taxa: taxa de crescimento intrínseco da população
        popt: população total

    Retorno:
        fs: função em string formatada
        fn: função em formato para ser plotada
    """

    k, p0, r, t, pt = symbols('k, p0, r, t, pt')


    f = k / (1+((k-p0)/p0) * exp(-r*t)) - pt

    fs = f.subs({'k': capac, 'p0':popi, 'r': taxa, 'pt': popt})

    fn = lambdify(t, fs, "numpy")

    print(f'Equação particular: f(t) = {fs}')
    print('Solução usando o método da bisseção:')
    
    return (fs,fn)


capac, popi, taxa, popt = 1000000, 1000, 0.01, 10000

fs,fn = eq_population(capac, popi, taxa, popt)

t = np.linspace(1,500)

plt.plot(t, fn(t),t, 0*t);
plt.xlabel('Tempo (t)')
plt.ylabel('População (p)')
plt.savefig('crescimento_populacional_plot.png')
plt.show()

def bissecao(f,a,b,tol,N):
    """Algoritmo para determinação de raízes pelo método da bisseção.

    Parâmetros: 
        f: string dependendo de uma variável, i.e., a função-alvo
            (e.g., 'x**2 - 1', 'x**2*cos(x)', etc.) 
        a: estimativa inferior
        b: estimativa superior
        tol: erro desejado (tolerância)
        N: número máximo de iterações a repetir

    Retorno: 
        x: aproximação para a raiz da função   
    """
    
    # construtor de tabela
    table = pt()
    
    # substitui expressões da string pelas chamadas das funções do numpy
    f = re.sub('(sin|sinh|cos|cosh|tan|tanh|exp|log|sqrt|log10|arcsin|arccos|arctan|arcsinh|arccosh|arctanh)', r'np.\1', f)
    
    # identifica a variável independente
    var = re.search(r'([a-zA-Z][\w]*) ?([\+\-\/*]|$|\))', f).group(1)
    
    # cria função anônima
    f = eval('lambda ' + var + ' :' + f)
    
    # checa se a função é de uma variável, senão lança erro        
    if len(inspect.getfullargspec(f).args) - 1 > 0:    
        raise ValueError('O código é válido apenas para uma variável.')

    # calcula valor da função nos extremos
    fa = f(a) 
    fb = f(b)
    
    # verifica sinal da função para o intervalo passado     
    if fa*fb >= 0:
        raise ValueError('A função deve ter sinais opostos em a e b!')
    
    # flag usada para marcar caso f(xm) = 0
    done = False;
        
    # no. iterações mínimo
    niter = int(np.ceil(np.log((b-a)/tol)/np.log(2)))
    if N < niter:
        print(f'! São necessárias pelo menos {niter} iterações, mas N = {N}.\n')

    
    # cabeçalho de tabela
    table.field_names = ['i','xm','f(xm)','a','b','f(a)','f(b)','EA']

    # bisecta o intervalo
    xm = (a+b)/2
    
    # contador 
    i = 1 
    
    # loop 
    while abs(a-b) > tol and (not done and N != 0):    
        
        # avalia a função no ponto médio
        fxm = f(xm) 
                        
        # adiciona linha de tabela de resultado
        table.add_row([i,np.round(xm,8),np.round(f(xm),8),
                   np.round(a,4),np.round(b,4),
                   np.round(f(a),4),np.round(f(b),4),
                   f'{abs(a-b):e}'])
   
        if fa*fxm < 0:      # Raiz esta à esquerda de xm
            b = xm
            fb = fxm
            xm = (a+b)/2
        elif fxm*fb < 0:    # Raiz esta à direita de xm
            a = xm
            fa = fxm
            xm = (a+b)/2
        else:               # Achamos a raiz
            done = True            
    
        N -= 1              # Atualiza passo
        i += 1              # Atualiza contador
    
    # impressão de tabela
    table.add_row([i,np.round(xm,8),np.round(f(xm),8),
                   np.round(a,4),np.round(b,4),
                   np.round(f(a),4),np.round(f(b),4),
                   f'{abs(a-b):e}'])
    table.align = 'c'; print(table)
    
    return xm



bissecao(str(fs), 200, 300, 1e-5, 25)



def newton(x0,f,df,tol,N):
    """Algoritmo para determinação de raízes pelo método de Newton.

    Parâmetros: 
        x0: estimativa inicial
        f: string dependendo de uma variável, i.e., a função-alvo
            (e.g., 'x**2 - 1', 'x**2*cos(x)', etc.) 
        df: string dependendo de uma variável, i.e., a derivada da função-alvo
        tol: erro desejado (tolerância)
        N: número máximo de iterações a repetir

    Retorno: 
        x: aproximação para a raiz da função    
    """

    # construtor de tabela
    table = pt()
    
    # substitui expressões da string pelas chamadas das funções do numpy
    f = re.sub('(sin|sinh|cos|cosh|tan|tanh|exp|log|sqrt|log10|arcsin|arccos|arctan|arcsinh|arccosh|arctanh)', r'np.\1', f)
    df = re.sub('(sin|sinh|cos|cosh|tan|tanh|exp|log|sqrt|log10|arcsin|arccos|arctan|arcsinh|arccosh|arctanh)', r'np.\1', df)
    
    # identifica a variável independente em f
    var = re.search(r'([a-zA-Z][\w]*) ?([\+\-\/*]|$|\))', f).group(1)
    
    # cria função
    f = eval('lambda ' + var + ' :' + f)
    
    # checa se a função é de uma variável, senão lança erro        
    try: 
        len(inspect.getfullargspec(f).args) - 1 > 0
    except:
        raise ValueError('O código é válido apenas para uma variável.')
    finally:
        # cria função derivada
        df = eval('lambda ' + var + ' :' + df)
    
    it = 0 # contador de iteracoes
    
    # cabeçalho de tabela
    table.field_names = ['i','x','f(x)','f\'(x)','ER']

    # imprime estimativa inicial
    print(f'Estimativa inicial: x0 = {x0:.6f}\n')  

    # Loop 
    for i in range(0,N):
        
        x = x0 - f(x0)/df(x0) # funcao de iteracao 
        
        e = abs(x-x0)/abs(x) # erro
        
        # tabela
        # impressão de tabela
        table.add_row([i,np.round(x,8),np.round(f(x),8),np.round(df(x),4),f'{e:e}'])
        table.align = 'c';      
        
        if e < tol:
            break
        x0 = x                
        
    print('Solução obtida com o método de Newthon-Raphson:')
    print(table)
       
    if i == N:
        print(f'Solução não obtida em {N:d} iterações')
    else:
        print(f'Solução obtida: x = {x:.6f}')


    return x


t = symbols('t')
f = -10000 + 1000000/(1 + 999*exp(-0.01*t))
df = diff(f, t)


newton(200, str(fs), str(df), 1e-5, 30)
