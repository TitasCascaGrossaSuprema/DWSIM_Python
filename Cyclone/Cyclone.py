import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

##############################################################################################################
# Dados de Operação
##############################################################################################################
Temperatura_do_Ar= 35	                            #°C
Pressao_do_Ar= 1.00	                                #atm
Vazao_de_Ar= 85.00 	                                #m3/h
Massa_Especifica_do_Ar= 1.14E-03	                #g/cm3
Viscosidade_Cinematica_Ar= 0.166	                #cm2/s

Vazao_Massica_Particulas= 40.00	                    #g/s
Massa_Especifica_Particulas= 2.430	                #g/cm3

##############################################################################################################
# Dados do Cilcone
##############################################################################################################
Diametro_do_Ciclone= 0.19                           #m
Diametro_do_Ciclone_cm=Diametro_do_Ciclone*100      #cm
Numero_de_Ciclones_em_Paralelo= 1

##############################################################################################################
# Lê Dimensões Padrões dos Ciclones Lappe e Stairmand
##############################################################################################################
df = pd.read_excel('Cyclone_Dimentions.xlsx', sheet_name="Dimensões")
Cyclone_Configuration = 'Lappe'
b_D     = df.loc[df['Dimensão'] == 'b/D', Cyclone_Configuration].values[0]
De_D    = df.loc[df['Dimensão'] == 'De/D', Cyclone_Configuration].values[0]
a_D     = df.loc[df['Dimensão'] == 'a/D', Cyclone_Configuration].values[0]
h_D     = df.loc[df['Dimensão'] == 'h/D', Cyclone_Configuration].values[0]
H_D     = df.loc[df['Dimensão'] == 'H/D', Cyclone_Configuration].values[0]
S_D     = df.loc[df['Dimensão'] == 'S/D', Cyclone_Configuration].values[0]
B_D     = df.loc[df['Dimensão'] == 'B/D', Cyclone_Configuration].values[0]

df = pd.read_excel('Cyclone_Dimentions.xlsx', sheet_name="Faixas")
kf = df.loc[df['Configuração'] == Cyclone_Configuration, 'kf'].values[0]
K = df.loc[df['Configuração'] == Cyclone_Configuration, 'K'].values[0]
vmin = df.loc[df['Configuração'] == Cyclone_Configuration, 'Vmin'].values[0]
vmax = df.loc[df['Configuração'] == Cyclone_Configuration, 'Vmax'].values[0]

##############################################################################################################
# Lê arquivo com curva de distribuição de partículas
##############################################################################################################
sizes = pd.read_excel('SizeDistribution.xlsx')
x_dados = sizes['Di (mm)']
y_dados = sizes['Xi exp']

##############################################################################################################
# Fazer a interpolação cúbica para encontrar o valor do D632
##############################################################################################################
interpolacao = interp1d(y_dados, x_dados, kind='cubic')
D632 = interpolacao(0.632)

##############################################################################################################
# Determina o parâmetro n do Modelo de Rosin-Rammlet-Bennet (RRB)
##############################################################################################################
def Rosin_Rammlet_Bennet(x, n):
    return 1 - np.exp(-(x/D632)**n)

parametros_otimizados, _ = curve_fit(Rosin_Rammlet_Bennet, x_dados, y_dados)
n_otimizado = parametros_otimizados
n_otimizado = n_otimizado[0]

##############################################################################################################
# Cálculo do Ciclone
##############################################################################################################
a = Diametro_do_Ciclone_cm*a_D
print("a=", a)
b = Diametro_do_Ciclone_cm*b_D
print("b=", b)

Q_m3s = (Vazao_de_Ar/3600)/Numero_de_Ciclones_em_Paralelo
Q_cm3s = Q_m3s*1e6
v = Q_cm3s/(a*b)
print("v= ", v)

mf = Q_cm3s*Massa_Especifica_do_Ar
print("mf= ", mf)

Vazao_Volumetrica_Particulas=Vazao_Massica_Particulas/Massa_Especifica_Particulas
print("Vazão Volumétrica Partículas= ", Vazao_Volumetrica_Particulas)

Cv = Vazao_Volumetrica_Particulas/(Vazao_Volumetrica_Particulas+Q_cm3s)
print("Cv= ", Cv)

Cm = Vazao_Massica_Particulas/(Vazao_Massica_Particulas+mf)
print("Cm= ", Cm)

Beta_Cv = 1/np.power((4.8*np.power(1-Cv,2))-(3.8*(1-Cv)),1/2)
print("Beta_Cv= ", Beta_Cv)

D_estrela = Diametro_do_Ciclone_cm*K*np.power((Viscosidade_Cinematica_Ar*Diametro_do_Ciclone_cm)/
                                              (Q_cm3s*((Massa_Especifica_Particulas/Massa_Especifica_do_Ar)-1)),1/2)*Beta_Cv

print("D_estrela=", D_estrela )

uc_cms= 4.0*Q_cm3s/(np.pi*np.power(Diametro_do_Ciclone_cm,2.0))
print("uc_cms= ", uc_cms)

uc_ms= uc_cms/100.0
print("uc_ms= ", uc_ms)

Massa_Especifica_do_Ar_kgm3= Massa_Especifica_do_Ar*1000.0
DP =Massa_Especifica_do_Ar_kgm3*kf*np.power(uc_ms,2.0)/2.0
print("Dp= ", DP)

Eficiencia_do_Soprador = 0.50
DP_no_Ciclone=DP/(1+0.023*np.power(Cm,0.69))
print("DP no Ciclone=", DP_no_Ciclone)

Potencia_do_Soprador=Numero_de_Ciclones_em_Paralelo*Q_m3s*DP_no_Ciclone/Eficiencia_do_Soprador
print("Potencia do Soprador= ", Potencia_do_Soprador)

Velocidade_na_Entrada = v/100.0
Diametro_de_Corte = D_estrela*1e4
print("Diametro de Corte=", Diametro_de_Corte)

Eficiencia_de_Separacao=(((1.11*n_otimizado)/(0.118+n_otimizado))/(1.81-0.322*n_otimizado+(D632/Diametro_de_Corte)))*(D632/Diametro_de_Corte)
print("Eficiencia de Separacao=", Eficiencia_de_Separacao*100, '%')

sizes['Di/D*'] = sizes['Di (mm)']/Diametro_de_Corte
sizes['ni(Di.D*)'] =(np.power(sizes['Di/D*'],2)/(1+np.power(sizes['Di/D*'],2)))*100

##############################################################################################################
# Plota Gráficos
##############################################################################################################
y_ajustado = Rosin_Rammlet_Bennet(x_dados, n_otimizado)

fig, axs = plt.subplots(2, 1, figsize=(8, 10))  # 2 linhas, 1 coluna

# Primeiro subplot (Ajuste de Curva Exponencial)
axs[0].scatter(x_dados, y_dados, label='Dados', color='red')
axs[0].plot(x_dados, y_ajustado, label='Ajuste', color='blue')
axs[0].set_xlabel('X')
axs[0].set_ylabel('Y')
axs[0].legend()
axs[0].set_title('Ajuste de Curva Exponencial')

# Segundo subplot (Eficiência Individual de Coleta)
axs[1].plot(sizes['Di/D*'], sizes['ni(Di.D*)'], color='blue')
axs[1].set_xlabel('Di/D*')
axs[1].set_ylabel('%')
axs[1].set_title('Eficiência Individual de Coleta')

# Ajustar o layout para que os subplots não se sobreponham
plt.tight_layout()

# Mostrar a figura com os dois subplots
plt.show()