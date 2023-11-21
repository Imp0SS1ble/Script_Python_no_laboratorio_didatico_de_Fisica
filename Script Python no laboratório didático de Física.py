#Created on Tue Sep 13 20:55:48 2022
#Autor: Alison
#Código da IC com o Prof. Dr. Vagner

# Bibliotecas:
import matplotlib.pyplot as plt           # Permite a plotagem de gráficos.
from scipy.optimize import curve_fit      # Permite o ajuste de curvas de funções.
from lmfit import Model                   # Permite o ajuste de curvas de funções mais complexas.
import numpy as np                        # Traz consigo funções matemáticas.
import sys                                # Permite a função sys.exit(0), que finaliza o processo se caso o usuário informar algum valor inválido.
import pandas as pd                       # Permite a criação de tabelas.
from sklearn.metrics import r2_score      # Permite o calculo de R².

#Constantes (cte) utilizadas nas equações dos experimentos
g = 9.81                                  # Aceleração gravitacional aproximada, utilizado nos experimentos: 3, 4 e 10.
cte_de_Stefan_Boltzmann = 5.67e-8         # Cte representada pela letra grega sigma, utilizada no experimento 6.
alpha = 4.5e-3                            # Valor apenas para o tungstênio, cte representada pela letra grega alpha, utilizada no experimento 6.
u0 = 4e-7*np.pi                           # Cte de permeabilidade magnética no vácuo, representado pela letra grega mi, utilizada no experimento 8.

# Seleção do tipo de experimento
print ('1 - Movimento retilíneo com velocidade constante; \n')
print ('2 - Movimento retilíneo com aceleração constante; \n')
print ('3 - Queda livre; \n')
print ('4 - Força elástica: Lei de Hooke; \n')
print ('5 - Movimento Harmônico Simples e a determinação da constante elástica de uma mola; \n')
print ('6 - Circuitos elétricos: Lei de OHM, Lei de Stefan-Boltzmann e Ebers-Moll; \n')
print ('7 - Circuito RC: Carga e descarga de um capacitor; \n')
print ('8 - Bobina de Helmholtz: Determinação da componente horizontal do campo magnético terrestre; \n')
print ('9 - Lei de Newton para o resfriamento; \n')
print ('10 - Pêndulo Simples; \n')
print ('11 - Refração: Lei de Snell-Descartes.\n')

# Usuário informa o experimento realizado para posterior análise.
tipo_do_experimento = int(input('Informe o número do experimento: '))

# Verificação do código para análise do experimento escolhido.
if tipo_do_experimento == 1:  
    # Criação das listas para posterior uso
    t = []
    x = []
    x_ajustado = []
    diferença_experimental_calculado = []
    # Usuário informa a quantidade de sensores usados, logo o número de dados sobre tempo
    n = int(input('Quantos sensores foram usados? '))
    # Usuário informa o tempo em cada sensor
    for i in range(0, n):
        t.append(float(input('Informe o tempo no ' +str(i+1)+'º sensor (em segundos)? ')))
    # A lista de t é invertida, pois estava dado erro sem a inversão da mesma
    t.reverse()
    # É subtraido o valor do 1° sensor, assim deixando os próximos sensores com tempos relativos ao primeiro sensor, e os dados são novamnetes invertidos
    for i in range(0, n):
        t[i] -= t[n-1]
    t.reverse()
    # A partir do valor informado pelo usuário sobre a quantidade de sensores usados, o 1° sensor é definido como ponto inicial (0) e as distâncias dos próximos é definida em relação ao 1° sensor
    x.append(0)
    n = n-1
    # Usuário informa a distância entre o 1° sensor até os demais sensores
    for i in range(0, n):
        x.append(float(input('Informe a distância do ' +str(1)+'º ao '+str(i+2)+'º sensor (em metros)? ')))
    # É criado uma função com a equação de posição em relação ao tempo, onde t são os valores de tempo, v a constante de velocidade e x0 a constante de posição inicial
    def y(t, v, x0):
       return t*v + x0
    # Define o tamanho da imagem do gráfico
    fig, ax = plt.subplots(figsize = (8,5))
    # Valores dos eixos x e y
    xData = np.array(t)
    yData = np.array(x)
    # Tamanho mínimo e máximo dos eixos x e y
    plt.axis(ymin=0, ymax=(x[n])*1.1, xmin=0, xmax=(t[n])*1.1)
    # Título do gráfico
    plt.title('Movimento retilíneo com velocidade constante')
    # Faz a plotagem no gráfico com os dados, informa o formato de ponto dos dados no gráfico e sua legenda
    plt.plot(xData, yData, 'bo', label='Dados')
    # Realiza o ajuste de curva do gráfico definindo os coeficientes v e x0
    popt, pcov = curve_fit(y, xData, yData)
      
    # Calcula a x ajustado com os valores de tempo informados e coeficientes calculados
    for i in range(0, n+1):
        x_ajustado.append(t[i]*popt[0] + popt[1])
        # Faz o calculo de diferença entre o valor obtido experimentalmente e o valor calculado com os coeficientes ajustados
        diferença_experimental_calculado.append(x[i] - x_ajustado[i])
        
    r2 = r2_score(x_ajustado, x) # A partir dos dados exeperimentais informados pelo usuário, compara-se estatisticamente com os dados gerados pelos parâmetros fixos estimados pelo script, obtendo-se o R², onde 1 representa uma compatibilidade perfeita entre dados experimentais e os estimados, já 0 representa uma incompatibilidade total entre eles.
    
    # Comprimento da linha do ajuste de curva
    comprimento_da_curva = t[n]*2
    # Define o valores de inicio, fim e de intervalo, para se calcular a função t*v + x0
    xFit = np.arange(0.0, comprimento_da_curva , 0.1)
    # Faz a plotagem no gráfico da linha, informa o formato de reta do ajuste no gráfico e sua legenda
    plt.plot(xFit, y(xFit, *popt), 'r', label=f'Parâmetros de ajuste: v=%5.5f m/s, x0=%5.5f m\nEquação: x = t*v + x0\nR² = {r2:.5f}' % tuple(popt))
    # Títulos dos eixos do gráfico
    plt.xlabel('t (s)')
    plt.ylabel('x (m)')
    # plotagem da legenda
    plt.legend()
    # Requisita ao código para mostar o gráfico com todas as informações
    plt.show()   
  
    # Cria um DataFrame com as listas
    df = pd.DataFrame({'Δt': t, 'x (m)': x, 'x (m) ajustado': x_ajustado, 'x - x ajustado (m)': diferença_experimental_calculado})
    # Exibe a tabela
    print(df)    
    
elif tipo_do_experimento == 2:
    # Criação das listas para posterior uso
    t = []
    x = []
    x_ajustado = []
    diferença_experimental_calculado = []
    # Usuário informa a quantidade de sensores usados, logo o número de dados sobre tempo
    n = int(input('Quantos sensores foram usados? '))
    # Usuário informa o tempo em cada sensor
    for i in range(0, n):
        t.append(float(input('Informe quantos segundos deu no ' +str(i+1)+'º sensor? ')))
    # A lista de t é invertida, pois estava dado erro sem a inversão da mesma
    t.reverse()
    # É subtraido o valor do 1° sensor, assim deixando os próximos sensores com tempos relativos ao primeiro sensor, e os dados são novamnetes invertidos
    for i in range(0, n):
        t[i] -= t[n-1]
    t.reverse()
    n = n-1
    # A partir do valor informado pelo usuário sobre a quantidade de sensores usados, o 1° sensor é definido como ponto inicial (0) e as distâncias dos próximos é definida em relação ao 1° sensor
    x.append(0)
    # Usuário informa a distância entre o 1° sensor até os demais sensores
    for i in range(0, n):
        x.append(float(input('Informe a distância do ' +str(1)+'º ao '+str(i+2)+'º sensor (em metros)? ')))
    # É criado uma função com a equação de posição em relação ao tempo, onde t são os valores de tempo, a é a aceleração constante, v0 é a velocidade inicial e x0 a constante de posição inicial
    def y(t, a, v0, x0):
       return  a*t**2 + v0*t + x0
    # Define o tamanho da imagem do gráfico
    fig, ax = plt.subplots(figsize = (10,6))
    # Valores dos eixos x e y
    xData = np.array(t)
    yData = np.array(x)
    # Tamanho mínimo e máximo dos eixos x e y
    plt.axis(ymin=0, ymax=(x[n])*1.1, xmin=0, xmax=(t[n])*1.1)
    # Título do gráfico
    plt.title('Movimento retilíneo com aceleração constante')
    # Faz a plotagem no gráfico com os dados, informa o formato de ponto dos dados no gráfico e sua legenda
    plt.plot(xData, yData, 'bo', label='Dados')
    # Realiza o ajuste de curva do gráfico definindo os coeficientes a, v0 e x0
    popt, pcov = curve_fit(y, xData, yData)
    # Comprimento da linha do ajuste de curva
        
    # Calcula a x ajustado com os valores de tempo informados e coeficientes calculados
    for i in range(0, n+1):
        x_ajustado.append(popt[0]*t[i]**2 + popt[1]*t[i] + popt[2])
        # Faz o calculo de diferença entre o valor obtido experimentalmente e o valor calculado com os coeficientes ajustados
        diferença_experimental_calculado.append(x[i] - x_ajustado[i])
        
    r2 = r2_score(x_ajustado, x) # A partir dos dados exeperimentais informados pelo usuário, compara-se estatisticamente com os dados gerados pelos parâmetros fixos estimados pelo script, obtendo-se o R², onde 1 representa uma compatibilidade perfeita entre dados experimentais e os estimados, já 0 representa uma incompatibilidade total entre eles.
    
    comprimento_da_curva = t[n]*2
    # Define o valores de inicio, fim e de intervalo, para se calcular a função a*t**2 + v0*t + x0
    xFit = np.arange(0.0, comprimento_da_curva , 0.1)
    # Faz a plotagem no gráfico da linha, informa o formato de reta do ajuste no gráfico e sua legenda
    plt.plot(xFit, y(xFit, *popt), 'r', label=f'Parâmetros de ajuste: a=%5.5f m/s², v0=%5.5f m/s, x0=%5.5f m\nEquação: x = a*t**2 + v0*t + x0\nR² = {r2:.5f}' % tuple(popt))
    # Títulos dos eixos do gráfico
    plt.xlabel('t (s)')
    plt.ylabel('x (m)')
    # plotagem da legenda
    plt.legend()
    # Requisita ao código para mostar o gráfico com todas as informações
    plt.show()

    # Cria um DataFrame com as listas
    df = pd.DataFrame({'Δt': t, 'x (m)': x, 'x (m) ajustado': x_ajustado, 'x - x ajustado (m)': diferença_experimental_calculado})
    # Exibe a tabela
    print(df)   
  
elif tipo_do_experimento == 3:
    # Criação das listas para posterior uso
    t = []
    y_eixo = []
    y_eixo_ajustado = []
    diferença_experimental_calculado = []
    # Usuário informa a quantidade de sensores usados, logo o número de dados sobre tempo e alturas de quedas
    n = int(input('De quantas alturas diferentes foi solto o objeto? '))
    # Usuário informa o tempo em cada sensor e a distância de queda
    for i in range(0, n):
        y_eixo.append(float(input('Informe a altura da ' +str(i+1)+'ª queda (em metros): \n')))
        t.append(float(input('Informe o tempo da ' +str(i+1)+'ª queda (em segundos): \n')))
    # É criado uma função com a equação de posição em relação ao tempo, onde t são os valores de tempo, a é a aceleração constante, v0 é a velocidade inicial e y0 é a posição inicial
    def y(t, a, v0, y0):
        return a*t**2 + v0*t + y0
    # Define o tamanho da imagem do gráfico
    fig, ax = plt.subplots(figsize = (8,5))
    # Valores dos eixos x e y
    xData = np.array(t)
    yData = np.array(y_eixo)
    # Tamanho mínimo e máximo dos eixos x e y
    plt.axis(ymin=0, ymax=(y_eixo[(n-1)])*1.1, xmin=0, xmax=(t[(n-1)])*1.1)
    # Título do gráfico
    plt.title('Gráfico de queda livre')
    # Faz a plotagem no gráfico com os dados, informa o formato de ponto dos dados no gráfico e sua legenda
    plt.plot(xData, yData, 'bo', label='Dados')
    # Define o valores de inicio, fim e de intervalo, para se calcular a função a*t**2 + v0*t + y0
    xFit = np.arange(0.0, ((t[(n-1)])*1.2), 0.000001)
    # Realiza o ajuste de curva do gráfico definindo os coeficientes a, v0 e x0
    popt, pcov = curve_fit(y, xData, yData)
      
    # Calcula a y ajustado com os valores de tempo informados e coeficientes calculados
    for i in range(0, n):
        y_eixo_ajustado.append(popt[0]*t[i]**2 + popt[1]*t[i] + popt[2])
        # Faz o calculo de diferença entre o valor obtido experimentalmente e o valor calculado com os coeficientes ajustados
        diferença_experimental_calculado.append(y_eixo[i] - y_eixo_ajustado[i])
    
    r2 = r2_score(y_eixo_ajustado, y_eixo) # A partir dos dados exeperimentais informados pelo usuário, compara-se estatisticamente com os dados gerados pelos parâmetros fixos estimados pelo script, obtendo-se o R², onde 1 representa uma compatibilidade perfeita entre dados experimentais e os estimados, já 0 representa uma incompatibilidade total entre eles.
    
    # Calcula o erro percentual da aceleração estimada em relação a aceleração gravitacional, uma vez que se trata de uma queda livre
    erro_percentual = (g - (popt[0]*2))/g*100
    # Garante que o valor do erro percentual é sempre positivo
    if erro_percentual <= 0:
        erro_percentual*-1
    # Faz a plotagem no gráfico da linha, informa o formato de reta do ajuste no gráfico e sua legenda
    plt.plot(xFit, y(xFit, *popt), 'r', label=f'Parâmetros de ajuste: a=%5.5f m/s²;\n v0=%5.5f m/s;\n y0=%5.5f m\nEquação: y = a*t**2 + v0*t + y0\nR² = {r2:.5f}' % tuple(popt))
    # Títulos dos eixos do gráfico
    plt.xlabel('t (s)')
    plt.ylabel('y (m)')
    # plotagem da legenda
    plt.legend()
    # Requisita ao código para mostar o gráfico com todas as informações
    plt.show()
    # Exibe o valor do erro percentual
    print('\nErro percentual comparado a aceleração gravitacional: e% = ' + str('%3.1f'%erro_percentual) + ' %\n')
  
    # Cria um DataFrame com as listas
    df = pd.DataFrame({'Δt': t, 'y (m)': y_eixo, 'y (m) ajustado': y_eixo_ajustado, 'y - y ajustado (m)': diferença_experimental_calculado})
    # Exibe a tabela
    print(df)
    
elif tipo_do_experimento == 4:
    # Criação das listas para posterior uso
    m = []
    F = []
    x = []
    F_ajustado = []
    diferença_experimental_calculado = []
    # Usuário informa a quantidade de molas usadas
    N = int(input('Quantas molas foram utilizadas? '))
    # Se o número de molas for superior a um há uma associação de molas, o usuário irá informar abaixo essa associação
    if N != 1:
        associação_de_molas = int(input('Qual foi a associação das molas (1 para série ou 2 para paralelo)? '))
        if associação_de_molas !=1 and associação_de_molas !=2:
            print('Valor inválido! Reinicie o código e digite um valor válido')
            sys.exit(0)
    # Usuário informa a quantidade de massas diferentes usadas, logo o número de dados sobre força e deformação
    n = int(input('Quantas massas diferentes foram utilizadas? '))
    for i in range(0, n):
        # Usuário informa a massa utilizada
        m.append(float(input('Informe a ' +str(i+1)+ 'ª massa utilizada (em quilogramas)? ')))
        # É calculado a força a partir da segunda lei de Newton
        F.append(m[i]*g)
        # Usuário informa a deformação da mola
        x.append(float(input('Qual foi a deformação da mola (em metros) com a '  +str(i+1)+ 'ª massa? ')))
    # É criado uma função com a equação de força em relação a elongação da mola, onde x são os valores de elongação, K a constante elática da mola e a é uma constante soma para ajuste
    def y(x,K,a):
        return x*K + a
    # Define o tamanho da imagem do gráfico
    fig, ax = plt.subplots(figsize = (8,5))
    # Valores dos eixos x e y
    xData = np.array(x)
    yData = np.array(F)
    # Tamanho mínimo e máximo dos eixos x e y
    plt.axis(ymin=0, ymax=(F[(n-1)])*1.1, xmin=0, xmax=(x[(n-1)])*1.1)
    # Define o título do gráfico de acordo com a configuração experimental
    if N == 1:
        plt.title('Força elástica: Lei de Hooke (Com uma mola)')
    elif associação_de_molas == 1:
        plt.title('Força elástica: Lei de Hooke (Associação em série de molas)')
    else:
        plt.title('Força elástica: Lei de Hooke (Associação em paralelo de molas)')
    # Faz a plotagem no gráfico com os dados, informa o formato de ponto dos dados no gráfico e sua legenda
    plt.plot(xData, yData, 'bo', label='Dados')
    # Realiza o ajuste de curva do gráfico definindo os coeficientes K e a
    popt, pcov = curve_fit(y, xData, yData)
    # Define o valores de inicio, fim e de intervalo, para se calcular a função x*K + a
    
    # Calcula a F ajustado com os valores de deformação informados e coeficientes calculados
    for i in range(0, n):
        # Calcula valores para a lista para posterior uso
        F_ajustado.append(x[i]*popt[0] + popt[1])
        # Faz o calculo de diferença entre o valor obtido experimentalmente e o valor calculado com os coeficientes ajustados
        diferença_experimental_calculado.append(F[i] - F_ajustado[i])
        
    r2 = r2_score(F_ajustado, F) # A partir dos dados exeperimentais informados pelo usuário, compara-se estatisticamente com os dados gerados pelos parâmetros fixos estimados pelo script, obtendo-se o R², onde 1 representa uma compatibilidade perfeita entre dados experimentais e os estimados, já 0 representa uma incompatibilidade total entre eles.
      
    xFit = np.arange(0.0, ((x[(n-1)])*1.2), 0.000001)
    # Faz de acordo com a configuração exeperimental a plotagem no gráfico da linha, informa o formato de reta do ajuste no gráfico e sua legenda
    if N == 1:
        plt.plot(xFit, y(xFit, *popt), 'r', label=f'Parâmetro de ajuste: K=%5.5f N/m, a=%5.5f N\nEquação: F = x*K + a\nR² = {r2:.5f}' % tuple(popt))
    else:
        plt.plot(xFit, y(xFit, *popt), 'r', label=f'Parâmetro de ajuste: K_eq.=%5.5f N/m, a=%5.5f N\nEquação: F = x*K_eq. + a\nR² = {r2:.5f}' % tuple(popt))  
    # Títulos dos eixos do gráfico
    plt.xlabel('x (m)')
    plt.ylabel('F (N)')
    # plotagem da legenda
    plt.legend()
    # Requisita ao código para mostar o gráfico com todas as informações
    plt.show()
      
    # Cria um DataFrame com as listas
    df = pd.DataFrame({'x (m)': x, 'm (kg)': m, 'F (N)': F, 'F (kg) ajustado': F_ajustado, 'F - F ajustado (N)': diferença_experimental_calculado})
    # Exibe o valor da aceleração gravitacional considerada
    print('\nValor da aceleração gravitacional utilizada : '+str(g)+ ' m/s²\n')
    # Exibe a tabela
    print(df)
    
elif tipo_do_experimento == 5:
    m = []
    t = []
    T = []
    T_ajustado = []
    diferença_experimental_calculado = []
    N = int(input('Quantas molas foram utilizadas? '))
    if N != 1:
        associação_de_molas = int(input('Qual foi a associação das molas (1 para série ou 2 para paralelo)? '))
        if associação_de_molas !=1 and associação_de_molas !=2:
            print('Valor inválido! Reinicie o código e digite um valor válido')
            sys.exit(0)
    n = int(input('Quantas massas diferentes foram utilizados? '))
    num_de_ciclos = int(input('Quantos ciclos foram realizados em um periodo de tempo? '))
    m.append(0)
    T.append(0)
    for i in range(0, n):
        m.append(float(input('Informe a ' +str(i+1)+ 'ª massa utilizada (em quilogramas)? ')))
        t.append(float(input('Qual foi o tempo de ' +str(num_de_ciclos)+ ' ciclos com a '  +str(i+1)+ 'ª massa? ')))
        T.append((t[i]/num_de_ciclos)**2)
    def y(m,K,a):
        return 4*(np.pi**2)*m/K + a
    fig, ax = plt.subplots(figsize = (8,5))
    xData = np.array(m)
    yData = np.array(T)
    plt.axis(ymin=0, ymax=(T[(n-1)])*1.1, xmin=0, xmax=(m[(n-1)])*1.1)
    if N == 1:
        plt.title('Movimento Harmônico Simples (Com uma mola)')
    elif associação_de_molas == 1:
        plt.title('Movimento Harmônico Simples (Associação em série)')
    else:
        plt.title('Movimento Harmônico Simples (Associação em paralelo)')
    plt.plot(xData, yData, 'bo', label='Dados')
    popt, pcov = curve_fit(y, xData, yData)
   
    for i in range(0, n):
        b = 4*(np.pi**2)*(1/popt[0])
        T_ajustado.append(4*(np.pi**2)*m[i]/popt[0] + popt[1])
        diferença_experimental_calculado.append(T[i] - T_ajustado[i])
        
    r2 = r2_score(T_ajustado, T) # A partir dos dados exeperimentais informados pelo usuário, compara-se estatisticamente com os dados gerados pelos parâmetros fixos estimados pelo script, obtendo-se o R², onde 1 representa uma compatibilidade perfeita entre dados experimentais e os estimados, já 0 representa uma incompatibilidade total entre eles.
       
    xFit = np.arange(0.0, ((m[(n-1)])*1.2), 0.000001)
    if N == 1:
        plt.plot(xFit, y(xFit, *popt), 'r', label=f'Parâmetro de ajuste: m/K={b:.5f} kg*m/N, a={popt[1]:.5f} s²\nEquação: T² = 4*π²*m/K + a\nR² = {r2:.5f} \nConstante da mola K={popt[0]:.5f} N/m')
    else:
        plt.plot(xFit, y(xFit, *popt), 'r', label=f'Parâmetro de ajuste: m/K={b:.5f} kg*m/N, a={popt[1]:.5f} s²\nEquação: T² = 4*π²*m/K_eq. + a\nR² = {r2:.5f} \nConstante equivalente da mola K={popt[0]:.5f} N/m')
    plt.xlabel('m (kg)')
    plt.ylabel('T² (s²)')
    plt.legend()
    plt.show()
  
    df = pd.DataFrame({'m (kg)': m, 'T² (s)': T, 'T² (s) ajustado': T_ajustado, 'T² - T² ajustado (s)': diferença_experimental_calculado})
    print(df)
    
elif tipo_do_experimento == 6:
    I = []
    V = [] 
    I_ajustado = []
    P_ajustado = []
    diferença_experimental_calculado = []
    
    grafico_da_lei_de_ohm = int(input('\nFoi feito o experimento com resistor(es) (1 para sim e 2 para não)? '))
    
    if grafico_da_lei_de_ohm != 1 and grafico_da_lei_de_ohm != 2:
        print('Valor inválido! Reinicie o código e digite um valor válido')
        sys.exit(0)
    
    grafico_da_lei_de_Stefan_Boltzmann = int(input('\nFoi feito o experimento com uma lâmpada (1 para sim e 2 para não)? '))
    
    if grafico_da_lei_de_Stefan_Boltzmann != 1 and grafico_da_lei_de_Stefan_Boltzmann != 2:
        print('Valor inválido! Reinicie o código e digite um valor válido')
        sys.exit(0)
    
    grafico_de_Ebers_Moll = int(input('\nFoi feito o experimento com diodo(s) (1 para sim e 2 para não)? '))
    
    if grafico_de_Ebers_Moll != 1 and grafico_de_Ebers_Moll != 2:
        print('Valor inválido! Reinicie o código e digite um valor válido')
        sys.exit(0)
    
    if grafico_da_lei_de_ohm == 1:
        Req = 0
        N = int(input('Quantos resistores foram utilizados? '))
        if N != 1:
            associação_de_resistores = int(input('Qual foi a associação dos resistores (1 para série ou 2 para paralelo)? '))
            if associação_de_resistores !=1 and associação_de_resistores !=2:
                print('Valor inválido! Reinicie o código e digite um valor válido')
                sys.exit(0)
            elif associação_de_resistores == 2:
                for i in range(0, N):
                    R = float(input('Qual o valor de resistência do '+str(i+1)+'º resitor em OHM (Ω)? '))
                    Req = Req + 1/R
                Req = 1/Req
            else: 
                for i in range(0, N):
                    R = float(input('Qual o valor de resistência do '+str(i+1)+'º resitor em OHM (Ω)? '))
                    Req = Req + R
        else:
            R = float(input('Qual o valor de resistência do resitor em OHM (Ω)? '))
            Req = R
            
        n = int(input('Quantos valores foram anotados de tensão e corrente sobre o(s) resistor(es)? '))
        for i in range(0, n):
            V.append(float(input('Informe qual a tensão no(s) resistor(es) na ' +str(i+1)+'ª anotação em volts (V)? ')))
            I.append(float(input('Informe qual a corrente no(s) resistor(es) na ' +str(i+1)+'ª anotação em ampère (A)? ')))
        
        def y(V,a):
            return V/a
        fig, ax = plt.subplots(figsize = (8,5))
        xData = np.array(V)
        yData = np.array(I)
        plt.axis(ymin=0, ymax=(I[(n-1)])*1.1, xmin=0, xmax=(V[(n-1)])*1.1)
        if N == 1:
            plt.title('Lei de OHM (Com um resitor)')
        elif associação_de_resistores == 1:
            plt.title('Lei de OHM (Associação em série de resistores)')
        else:
            plt.title('Lei de OHM (Associação em paralelo de resistores)')
    
        plt.plot(xData, yData, 'bo', label='Dados')
        popt, pcov = curve_fit(y, xData, yData)
 
        for i in range(0, n):
            I_ajustado.append(V[i]/popt[0])
            diferença_experimental_calculado.append(I[i] - I_ajustado[i])
        
        r2 = r2_score(I_ajustado, I) # A partir dos dados exeperimentais informados pelo usuário, compara-se estatisticamente com os dados gerados pelos parâmetros fixos estimados pelo script, obtendo-se o R², onde 1 representa uma compatibilidade perfeita entre dados experimentais e os estimados, já 0 representa uma incompatibilidade total entre eles.    
    
        xFit = np.arange(0.0, ((V[(n-1)])*1.2), 0.001)
        if N == 1:
            plt.plot(xFit, y(xFit, *popt), 'r', label=f'Parâmetro de ajuste: a=%5.5f Ω.\nEquação: I = V*a\nR² = {r2:.5f}' % tuple(popt))
        else:
            plt.plot(xFit, y(xFit, *popt), 'r', label=f'Parâmetro de ajuste: a=%5.5f Ω.\nEquação: I = V*a\nR² = {r2:.5f}' % tuple(popt))  
        plt.xlabel('V (V)')
        plt.ylabel('i (A)')
        plt.legend()
        plt.show()
       
        df = pd.DataFrame({'V (V)': V, 'I (A)': I, 'I (A) ajustado': I_ajustado, 'I - I ajustado (s)': diferença_experimental_calculado})
        print(df)
        print('\n')
        
        
    if grafico_da_lei_de_Stefan_Boltzmann == 1:
        P = []
        T = []
        R = []
        
        T_amb = float(input('Informe a temperatura ambiente do local do experimento em graus celcius (°C)? '))
        T_amb = T_amb + 273.15 # conversão de celcios para kelvin
        R0 = float(input('Informe a resistência da lâmpada quando não há passagem de corrente em ohm (Ω)? '))
        n = int(input('Quantos valores foram anotados de tensão, corrente e resistência sobre a lâmpada? '))
        for i in range(0, n):
            V.append(float(input('Informe qual a tensão sobre lâmpada na ' +str(i+1)+'ª anotação em volt (V)? ')))
            I.append(float(input('Informe qual a corrente sobre lâmpada na ' +str(i+1)+'ª anotação em ampère (A)? ')))
            P.append(V[i]*I[i])
            R.append(V[i]/I[i])
            T.append((R[i]/(alpha*R0) - (1/alpha) + T_amb)) # Equação de resistência elétrica do filamento de uma lâmpada
                     
        def y(T, a):
           return cte_de_Stefan_Boltzmann*a*T**4
        fig, ax = plt.subplots(figsize = (8,5))
        xData = np.array(T)
        yData = np.array(P)
        plt.axis(ymin=0, ymax=(P[n-1])*1.1, xmin=0, xmax=(T[n-1])*1.1)
        plt.title('Lei de Stefan Boltzmann')
        plt.plot(xData, yData, 'bo', label='Dados')
        popt, pcov = curve_fit(y, xData, yData)
        
        list.clear(diferença_experimental_calculado)
        for i in range(0, n):
            P_ajustado.append(cte_de_Stefan_Boltzmann*popt[0]*T[i]**4)
            diferença_experimental_calculado.append(P[i] - P_ajustado[i])
            
        r2 = r2_score(P_ajustado, P) # A partir dos dados exeperimentais informados pelo usuário, compara-se estatisticamente com os dados gerados pelos parâmetros fixos estimados pelo script, obtendo-se o R², onde 1 representa uma compatibilidade perfeita entre dados experimentais e os estimados, já 0 representa uma incompatibilidade total entre eles.
        coef = cte_de_Stefan_Boltzmann*popt[0]
        
        comprimento_da_curva = T[n-1]*1.5
        xFit = np.arange(0.0, comprimento_da_curva , 0.1)
        plt.plot(xFit, y(xFit, *popt), 'r', label=f'Parâmetro de ajuste: σ*a = {coef:.5e} \n a=%5.5f m².\nEquação: P(T) = σ*a*T**4\nR² = {r2:.5f}' % tuple(popt))
        plt.ylabel('P (W)')
        plt.legend()
        plt.show()   

        df = pd.DataFrame({'T (K)': T, 'P (W)': P, 'P (W) ajustado': P_ajustado, 'P - P ajustado (W)': diferença_experimental_calculado})
        print(df)
        print('\n')             
        
    if grafico_de_Ebers_Moll == 1:
        Ut = 25e-3
        tensao_fonte = float(input('Qual a tensão da fonte em volt (V)? '))
        resistor = float(input('\nQual a resistência usada junto com o diodo em ohm (Ω)? '))
        n = int(input('Quantos valores foram anotados de tensão e corrente sobre o diodo? '))
        for i in range(0, n):
            V.append(float(input('Informe qual a tensão sobre o diodo na ' +str(i+1)+'ª anotação em volt (V)? ')))
            I.append(float(input('Informe qual a corrente no diodo na ' +str(i+1)+'ª anotação em ampère (A)? ')))
        I_máx = ((tensao_fonte - V[n-1])/resistor)*1.25
        
        def y(V, a, b):
            return a*(np.exp(V/(b*Ut))-1)
        fig, ax = plt.subplots(figsize = (8,5))
        xData = np.array(V)
        yData = np.array(I)
        plt.axis(ymin=0, ymax=(I[n-1])*1.1, xmin=0, xmax=(V[n-1])*1.1)
        plt.title('Lei de Ebers Moll')
        plt.plot(xData, yData, 'bo', label='Dados')

        # Criar um modelo com a função de ajuste e limites para a e b
        model = Model(y)
        params = model.make_params(a=(I_máx/2), b=5)
        params['a'].set(min=0, max=I_máx)
        params['b'].set(min=1e-1, max=10)
   
        # Ajustar o modelo aos dados
        result = model.fit(I, params, V=V) 
       
        list.clear(diferença_experimental_calculado)
        for i in range(0, n):
            I_ajustado.append(result.params['a'].value*(np.exp(V[i]/(result.params['b'].value*Ut))-1))
            diferença_experimental_calculado.append(I[i] - I_ajustado[i])
            
        r2 = r2_score(I_ajustado, I) # A partir dos dados exeperimentais informados pelo usuário, compara-se estatisticamente com os dados gerados pelos parâmetros fixos estimados pelo script, obtendo-se o R², onde 1 representa uma compatibilidade perfeita entre dados experimentais e os estimados, já 0 representa uma incompatibilidade total entre eles.
        coef = 1/(result.params['b'].value*Ut)
        
        # Plotar os dados e a curva ajustada
        plt.scatter(V, I)
        a = result.params['a'].value
        eq_label = f'Parâmetros de ajuste: a= {a:.5e} A, 1/b*Ut = %5.5f\nb= %5.5f \nEquação: I(V) = a*(np.exp(V/(b*Ut))-1)\nR² = {r2:.5f}' % (coef, result.params['b'].value)
        plt.plot(V, result.best_fit,'r', label=eq_label)
        plt.xlabel("V (V)")
        plt.ylabel("I (A)")
        plt.legend()
        plt.show()
 
        df = pd.DataFrame({'V (V)': V, 'I (A)': I, 'I (A) ajustado': I_ajustado, 'I - I ajustado (s)': diferença_experimental_calculado})
        print(df)
    
elif tipo_do_experimento == 7:
    tc = []
    Vc = []
    td = []
    Vd = []
    Vc_ajustado = []
    Vd_ajustado = []
    diferença_experimental_calculado = []
    Ceq = 0
    
    Opcoes = int(input('Qual gráfico será gerado e ajustado (1 para carga de capacitor ou 2 para descarga de capacitor ou 3 para ambos)? '))
    if Opcoes != 1 and Opcoes != 2 and Opcoes != 3:       
        print('Valor inválido! Reinicie o código e digite um valor válido')
        sys.exit(0)
    N = int(input('Quantos capacitores foram utilizados? '))
    if N != 1:
        associação_de_capacitores = int(input('Qual foi a associação dos capacitores (1 para série ou 2 para paralelo)? '))
        if associação_de_capacitores != 1 and associação_de_capacitores != 2:
            print('Valor inválido! Reinicie o código e digite um valor válido')
            sys.exit(0)
        elif associação_de_capacitores == 1:
            for i in range(0, N):
                C = float(input('Qual o valor de capacitância do '+str(i+1)+'º capacitor em Microfarad (μF - 10^-6 F)? '))
                Ceq = Ceq + 1/C
            Ceq = 1/Ceq
        else: 
            for i in range(0, N):
                C = float(input('Qual o valor de capacitância do '+str(i+1)+'º capacitor em Microfarad (μF - 10^-6 F)? '))
                Ceq = Ceq + C
    else:
        C = float(input('Qual o valor de capacitância do capacitor em Microfarad (μF - 10^-6 F)? '))
        Ceq = C
 
    # Obter os dados de tensão da fonte, resistência do resistor usado e quantidade de valores de tempo e tensão anotados. 
    V = float(input('Qual o valor de tensão da fonte utilizada em volts (V)? '))
    R = float(input('Qual o valor de resistência do resistor utilizado no circuito em Quilo-ohm (kΩ)? '))
    R = R*1000
    if Opcoes == 1 or Opcoes == 3:
        
        # Obter a quantidade de dados e os dados de tensão e tempo de carga do capacitor do usuário 
        nc = int(input('Quantos valores foram anotados na carga do capacitor? '))
        for i in range(0, nc):
            Vc.append(float(input('Informe qual a tensão no(s) capacitor(es) na ' +str(i+1)+'ª anotação de carga em volt (V)? ')))
            tc.append(float(input('Informe o tempo na ' +str(i+1)+'ª tensão anotada de carga em segundos (s)? ')))
            
        # Definir a função de ajuste de carga de capacitor
        def y(tc, a, b):
            return a*(1-np.exp(-(b*tc)))
        fig, ax = plt.subplots(figsize = (8,5))
        xData = np.array(tc)
        yData = np.array(Vc)
        plt.axis(ymin=0, ymax=(Vc[nc-1])*1.1, xmin=0, xmax=(tc[nc-1])*1.1)
        if N == 1:
            plt.title('Circuito RC: Carga de um capacitor de ' +str(Ceq)+  ' μF')
        elif N != 1 and associação_de_capacitores == 1:
            plt.title(f'Circuito RC: Carga de capacitores em série com capacitância equivalente de {round(Ceq, 3)} μF')
        else:
            plt.title(f'Circuito RC: Carga de capacitores em paralelo com capacitância equivalente de {round(Ceq, 3)} μF')
        plt.plot(xData, yData, 'bo', label='Dados')

        # Criar um modelo com a função de ajuste e limites para a e b
        model = Model(y)
        params = model.make_params(a=1, b=(1/(R*Ceq*10**-6)))
        params['a'].set(min=Vc[nc-1], max=V)
        params['b'].set(min=(1/(R*Ceq*10**-6))*0.75, max=(1/(R*Ceq*10**-6))*1.25)
   
        # Ajustar o modelo aos dados
        result = model.fit(Vc, params, tc=tc) 
   
        for i in range(0, nc):
            Vc_ajustado.append(result.params['a'].value*(1-np.exp(-(result.params['b'].value*tc[i]))))
            diferença_experimental_calculado.append(Vc[i] - Vc_ajustado[i])
            
        r2 = r2_score(Vc_ajustado, Vc) # A partir dos dados exeperimentais informados pelo usuário, compara-se estatisticamente com os dados gerados pelos parâmetros fixos estimados pelo script, obtendo-se o R², onde 1 representa uma compatibilidade perfeita entre dados experimentais e os estimados, já 0 representa uma incompatibilidade total entre eles.
          
        # Plotar os dados e a curva ajustada
        plt.scatter(tc, Vc)
        eq_label = f'Parâmetros de ajuste: a= %5.5f V, b= %5.5f S*1/C\nEquação: a*(1-np.exp(-b*t))\nR² = {r2:.5f}' % (result.params['a'].value, result.params['b'].value)
        plt.plot(tc, result.best_fit ,'r', label=eq_label)
        plt.xlabel("t (s)")
        plt.ylabel("V (V)")
        plt.legend()
        plt.show()
       
        df = pd.DataFrame({'t (s)': tc, 'V (V)': Vc, 'V (V) ajustado': Vc_ajustado, 'V - V ajustado (V)': diferença_experimental_calculado})
        print(df)
        
        list.clear(Vc)
        list.clear(tc)
   
    if Opcoes == 2 or Opcoes == 3:
        
        # Obter  a quantidade de dados e os dados de tensão e tempo de descarga do capacitor do usuário 
        nd = int(input('Quantos valores foram anotados na descarga do capacitor? '))
        for i in range(0, nd):
            Vd.append(float(input('Informe qual a tensão no(s) capacitor(es) na ' +str(i+1)+'ª anotação de descarga em volt (V)? ')))
            td.append(float(input('Informe quantos segundos deu na ' +str(i+1)+'ª tensão anotada de descarga em segundos (s)? ')))
        
        # Definir a função de ajuste de descarga de capacitor
        def y(td, a, b):
            return a*(np.exp(-(b*td)))
        fig, ax = plt.subplots(figsize = (8,5))
        xData = np.array(td)
        yData = np.array(Vd)
        plt.axis(ymin=0, ymax=(Vd[0])*1.1, xmin=0, xmax=(td[nd-1])*1.1)
        if N == 1:
            plt.title('Circuito RC: Descarga de um capacitor de ' +str(Ceq)+  ' μF')
        elif N != 1 and associação_de_capacitores == 1:
            plt.title(f'Circuito RC: Descarga de capacitores em série com capacitância equivalente de {round(Ceq, 3)} μF')
        else:
            plt.title(f'Circuito RC: Descarga de capacitores em paralelo com capacitância equivalente de {round(Ceq, 3)} μF')
        plt.plot(xData, yData, 'bo', label='Dados')

        # Criar um modelo com a função de ajuste e limites para a e b
        model = Model(y)
        params = model.make_params(a=1, b=(1/(R*Ceq*10**-6)))
        params['a'].set(min=Vd[0], max=V)
        params['b'].set(min=(1/(R*Ceq*10**-6))*0.75, max=(1/(R*Ceq*10**-6))*1.25)
   
        # Ajustar o modelo aos dados
        result = model.fit(Vd, params, td=td) 
   
        list.clear(diferença_experimental_calculado)
        print('\n')
        
        for i in range(0, nd):
            Vd_ajustado.append(result.params['a'].value*(np.exp(-(result.params['b'].value*td[i]))))
            diferença_experimental_calculado.append(Vd[i] - Vd_ajustado[i])
            
        r2 = r2_score(Vd_ajustado, Vd) # A partir dos dados exeperimentais informados pelo usuário, compara-se estatisticamente com os dados gerados pelos parâmetros fixos estimados pelo script, obtendo-se o R², onde 1 representa uma compatibilidade perfeita entre dados experimentais e os estimados, já 0 representa uma incompatibilidade total entre eles.
         
        # Plotar os dados e a curva ajustada
        plt.scatter(tc, Vc)
        eq_label = f'Parâmetros de ajuste: a= %5.5f V, b= %5.5f S*1/C\nEquação: a*(np.exp(-b*t))\nR² = {r2:.5f}' % (result.params['a'].value, result.params['b'].value)
        plt.plot(td, result.best_fit,'r', label=eq_label)
        plt.xlabel("t (s)")
        plt.ylabel("V (V)")
        plt.legend()
        plt.show()
   
        df = pd.DataFrame({'t (s)': td, 'V (V)': Vd, 'V (V) ajustado': Vd_ajustado, 'V - V ajustado (V)': diferença_experimental_calculado})
        print(df)
    
elif tipo_do_experimento == 8:
    I = []
    ang_graus = []
    ang_rad = []
    B = []
    tan = []
    B_ajustado = []
    diferença_experimental_calculado = []
    
    N = int(input('Informe quantas voltas possui as bobinas? '))
    N = float(N)
    R = float(input('Informe o raio das bobinas (em metros)? '))
    n = int(input('Qual a quantidade de medições de correntes elétricas e ângulos? '))
    for i in range(0, n):
        I.append(float(input('Qual a ' +str(i+1)+ 'ª corrente elétrica medida (em ampère)? ')))
        ang_graus.append(float(input('Qual o ' +str(i+1)+ 'º ângulo medido (em graus)? ')))
        ang_rad.append(ang_graus[i]*(np.pi/180))
        B.append((8.0*u0*N*I[i])/(np.sqrt(5.0**3)*R))
        tan.append(np.tan(ang_rad[i]))
        
    def y(a,b,tan):
        return a*tan+b
    fig, ax = plt.subplots(figsize = (8,5))
    xData = np.array(tan)
    yData = np.array(B)
    plt.axis(ymin=0, ymax=(B[(n-1)])*1.1, xmin=0, xmax=(tan[(n-1)])*1.1)
    plt.title('Bobina de Helmholtz')
    plt.plot(xData, yData, 'bo', label='Dados')
    popt, pcov = curve_fit(y, xData, yData)
   
    for i in range(0, n):
        B_ajustado.append(popt[0]*tan[i]+popt[1])
        B[i] = B[i]*1e6
        B_ajustado[i] = B_ajustado[i]*1e6
        diferença_experimental_calculado.append(B[i] - B_ajustado[i])
        
    r2 = r2_score(B_ajustado, B) # A partir dos dados exeperimentais informados pelo usuário, compara-se estatisticamente com os dados gerados pelos parâmetros fixos estimados pelo script, obtendo-se o R², onde 1 representa uma compatibilidade perfeita entre dados experimentais e os estimados, já 0 representa uma incompatibilidade total entre eles.
     
    xFit = np.arange(0.0, ((tan[(n-1)])*1.2), 0.000001)
    plt.plot(xFit, y(xFit, *popt), 'r', label=f'Parâmetro de ajuste: a= {popt[0]:.3e}, b= {popt[1]:.3e}\nEquação: a*tan+b\nR² = {r2:.5f}')
    plt.xlabel('tan θ')
    plt.ylabel('B (T)')
    plt.legend()
    plt.show()
    campo_magnetico_da_terra = popt[0]*np.tan((np.pi/4))+popt[1]
    print(f'O campo magenético terrestre a partir da equação é de: B = {campo_magnetico_da_terra:.3e} T \n')
    
    df = pd.DataFrame({'I (A)': I, 'θ (°)': ang_graus, 'tan θ': I, 'B (μT)': B, 'B (μT) ajustado': B_ajustado, 'B - B ajustado (μT)': diferença_experimental_calculado})
    print(df)
    
elif tipo_do_experimento == 9:
    T = []
    t = []
    T_ajustado = []
    diferença_experimental_calculado = []
    
    T_amb = float(input('Informe a temperatura ambiente do local do experimento em graus celcius (°C)? '))
    n = int(input('Quantos valores foram anotados de temperatura do fluído e tempo? '))
    for i in range(0, n):
        T.append(float(input('Informe qual a temperatura do fluído na ' +str(i+1)+'ª anotação em graus celcius (°C)? ')))
        t.append(float(input('Informe qual o tempo da ' +str(i+1)+'ª anotação em minutos (min)? ')))
    T0 = T[0]
        
    def y(t, a):
        return (np.exp(a*t)*(T0 - T_amb)) + T_amb
    fig, ax = plt.subplots(figsize = (8,5))
    xData = np.array(t)
    yData = np.array(T)
    plt.axis(ymin=0, ymax=(T[0])*1.1, xmin=0, xmax=(t[n-1])*1.1)
    plt.title('Lei de Newton para o resfriamento')
    plt.plot(xData, yData, 'bo', label='Dados')
    popt, pcov = curve_fit(y, xData, yData)
    
    for i in range(0, n):
        T_ajustado.append((np.exp(popt[0]*t[i])*(T0 - T_amb)) + T_amb)
        diferença_experimental_calculado.append(T[i] - T_ajustado[i])
        
    r2 = r2_score(T_ajustado, T) # A partir dos dados exeperimentais informados pelo usuário, compara-se estatisticamente com os dados gerados pelos parâmetros fixos estimados pelo script, obtendo-se o R², onde 1 representa uma compatibilidade perfeita entre dados experimentais e os estimados, já 0 representa uma incompatibilidade total entre eles.
    
    comprimento_da_curva = t[n-1]*1.5
    xFit = np.arange(0.0, comprimento_da_curva , 0.001)
    plt.plot(xFit, y(xFit, *popt), 'r', label=f'Parâmetro de ajuste: a=%5.5f 1/s.\nEquação: T(t) = np.exp(a*t)*(T0 - T_amb) + T_amb\nR² = {r2:.5f}' % tuple(popt))
    plt.xlabel('t (s)')
    plt.ylabel('T (°C)')
    plt.legend()
    plt.show()  
    
    df = pd.DataFrame({'t (min)': t, 'T (°C)': T, 'T (°C) ajustado': T_ajustado, 'T - T ajustado (°C)': diferença_experimental_calculado})
    print(df)  
    
elif tipo_do_experimento == 10:
    T = []
    Tm = []
    x = []
    Tm_ajustado = []
    diferença_experimental_calculado = []
        
    num_ciclos = int(input('\nQuantos ciclos foram contados para a marcação do tempo? '))
    
    if num_ciclos <= 0:
        print('Valor inválido! Reinicie o código e digite um valor válido')
        sys.exit(0)
        
    n = int(input('\nQuantas anotações de tempo de ' +str(num_ciclos)+ ' ciclos foram feitas? '))
       
    print('\nInforme os valores em ordem crescente de comprimento')
    for i in range(0, n):
        T.append(float(input('Informe o ' +str(i+1)+ 'º tempo anotado em segundos (s)? ')))
        Tm.append(T[i]/num_ciclos)
        x.append(float(input('Qual foi o comprimento da linha do pêndulo no '  +str(i+1)+ 'º comprimento anotado em metros (m)? ')))
      
    def y(x,a,b):
        return (2*np.pi*(x**0.5/a**0.5)) + b
             
    fig, ax = plt.subplots(figsize = (8,5))
    xData = np.array(x)
    yData = np.array(Tm)
    plt.axis(ymin=0, ymax=(Tm[(n-1)])*1.1, xmin=0, xmax=(x[(n-1)])*1.1)
    plt.title('Pêndulo simples (Depêndencia da variável de comprimento da linha)')
    plt.plot(xData, yData, 'bo', label='Dados')
    popt, pcov = curve_fit(y, xData, yData)
     
    for i in range(0, n):
         Tm_ajustado.append((2*np.pi*(x[i]**0.5/popt[0]**0.5)) + popt[1])
         diferença_experimental_calculado.append(Tm[i] - Tm_ajustado[i])
         
    r2 = r2_score(Tm_ajustado, Tm) # A partir dos dados exeperimentais informados pelo usuário, compara-se estatisticamente com os dados gerados pelos parâmetros fixos estimados pelo script, obtendo-se o R², onde 1 representa uma compatibilidade perfeita entre dados experimentais e os estimados, já 0 representa uma incompatibilidade total entre eles.
    coef = 2*np.pi*np.sqrt(1/popt[0])
    
    xFit = np.arange(0.0, ((x[(n-1)])*1.2), 0.000001)
    plt.plot(xFit, y(xFit, *popt), 'r', label=f'Parâmetro de ajuste: 2*pi*sqrt(1/a) = {coef:.5f},  b=%5.5f \na=%5.5f m/s² \nEquação: T = 2*pi*sqrt(L/a) + b\nR² = {r2:.5f}' % (popt[1], popt[0]))
    
    plt.xlabel('L (m)')
    plt.ylabel('T (s)')
    plt.legend()
    plt.show()
    
elif tipo_do_experimento == 11:
    I = []
    r = []
    sen_I = []
    sen_r = []
    sen_I_ajustado = []
    diferença_experimental_calculado = []
    
    n = int(input('Quantos ângulos foram medidos? '))
    for i in range(0, n):
        I.append(float(input('Informe o ângulo de incidencia da luz na ' +str(i+1)+'ª anotação em graus (°)? ')))
        r.append(float(input('Informe o ângulo de refração da luz na ' +str(i+1)+'ª anotação em graus (°)? ')))
        sen_I.append((I[i]*np.pi)/180)
        sen_r.append((r[i]*np.pi)/180)
        
    def y(sen_r, a):
       return a*sen_r
    fig, ax = plt.subplots(figsize = (8,5))
    xData = np.array(sen_r)
    yData = np.array(sen_I)
    plt.axis(ymin=0, ymax=(sen_I[n-1])*1.1, xmin=0, xmax=(sen_r[n-1])*1.1)
    plt.title('Refração: Lei de Snell-Descartes')
    plt.plot(xData, yData, 'bo', label='Dados')
    popt, pcov = curve_fit(y, xData, yData)
    
    for i in range(0, n):
        sen_I_ajustado.append((2*np.pi*(popt[0]**0.5/g**0.5)) + popt[1])
        diferença_experimental_calculado.append(sen_I[i] - sen_I_ajustado[i])
        
    r2 = r2_score(sen_I_ajustado, sen_I) # A partir dos dados exeperimentais informados pelo usuário, compara-se estatisticamente com os dados gerados pelos parâmetros fixos estimados pelo script, obtendo-se o R², onde 1 representa uma compatibilidade perfeita entre dados experimentais e os estimados, já 0 representa uma incompatibilidade total entre eles.
    
    comprimento_da_curva = sen_I[n-1]*2
    
    xFit = np.arange(0.0, comprimento_da_curva , 0.000001)
    plt.plot(xFit, y(xFit, *popt), 'r', label=f'Parâmetros de ajuste: a=%5.5f \nEquação: sen i = a * sen r\nR² = {r2:.5f}' % tuple(popt))
    plt.xlabel('sen r')
    plt.ylabel('sen i')
    plt.legend()
    plt.show()
    
else:
    print ('Número de experimento inválido.')