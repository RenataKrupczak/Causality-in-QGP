import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import make_interp_spline 

## As coisas que devem ser mudadas para cada evento (nome do arquivo), o tau (tempo) e o range das linhas (valor de i e valor a ser somando a i, grid e tamanho real do sistema)

## Leitura do arquivo de texto, cada coluna se torna um vetor
W_tautau, W_taux, W_tauy, W_taueta, W_xx, W_xy, W_xeta, W_yy, W_yeta, W_etaeta, pressao, som, entropy, energy, temperature, visc1, visc2 = np.loadtxt(
    'evolution_Wmunu.dat',
    skiprows = 0,
    unpack = True
)
Bulk=np.loadtxt(
    'evolution_bulk_pressure.dat',
    skiprows = 0,
    unpack = True
)
eixox_tempo = []
eixoy_energiacausal = []
eixoy_energianada = []
eixoy_energiaacausal = []
celulas_causal = []
celulas_nada = []
celulas_acausal = []

# Loop para os vários corte
i = 0
tempo = -0.13
while i<2038400:
    tempo = round((tempo + 0.5), 2)
    print('Time: ', tempo)
    
    # zerar todas as variáveis
    resultado = []
    densityEnergy = 0
    epsilonCausal = 0
    epsilonNada = 0
    epsilonAcausal = 0
    eqALPHAn = 0
    eqBETAn = 0
    eqGAMMAn = 0
    eqDELTAn = 0
    eqTHETAn = 0
    eqKAPPAn = 0
    eqRHOn = 0
    eqOMEGAn = 0
    eqA = 0
    eqB = 0
    eqC = 0
    eqD = 0
    eqE = 0
    eqF = 0
    eqALPHA = 0
    eqBETA = 0
    eqGAMMA = 0
    eqDELTA = 0
    eqTHETA = 0
    eqKAPPA = 0
    eqRHO = 0
    eqOMEGA = 0
    nada = 0
    acausal = 0
    causal = 0
    Nada = 0
    Acausal = 0
    Causal = 0
    
    # Loop para cada célula do fluído
    for linha in range(i, i+78400):
        ## Leitura do valor de cada variável para a célula nº linha
        a1 = W_tautau[linha]
        a2 = W_taux[linha]
        a3 = W_tauy[linha]
        a4 = W_taueta[linha]
        a5 = W_xx[linha]
        a6 = W_xy[linha]
        a7 = W_xeta[linha]
        a8 = W_yy[linha]
        a9 = W_yeta[linha]
        a10 = W_etaeta[linha]
        
        p = pressao[linha]
        cs2 = som[linha]
        etas = visc1[linha]  ## etas e zetas eh a viscosidade especifica: eta/s e zeta/s
        zetas = visc2[linha]
        s = entropy[linha]
        epsilon = energy[linha]
        T = temperature[linha]
        PI = Bulk[linha]
        
        eta = etas*s    ## Viscosidades não especificas
        zeta = zetas*s
        
        ## Com as componentes, formamos a matriz shear 4x4 para a célula
        matriz = np.array([[a1, a2, a3, a4], [a2, a5, a6, a7], [a3, a6, a8, a9], [a4, a7, a9, a10]])
        g = np.array([[-1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
        invariant = np.dot(g, matriz) 
        
        ## Tirar as células que nao tem fluido
        if T <= 0.151:
            condicao = 1
        else:
            densityEnergy = densityEnergy+epsilon
            ## Calculamos os autovalores dessa matriz e colocamos em um vetor chamado autovalor
            l = np.linalg.eigvals(invariant)
            autovalor = np.array(l)
            
            if np.trace(invariant)>=0.1 or np.trace(invariant)<=-0.1:
                print('erro')
                print(np.trace(invariant))
            naonulos = []
            ## Obtemos os autovalores dessa célula seguindo as restrições do artigo: Lambda0=0 e Lambda1<=Lambda2<=Lambda3, sendo Lambda1<=0<=Lambda3
            ## Primeiro tomamos o autovalor nulo como sendo Lambda0
            for r in range(0, 4):
                if autovalor[r]<=10**(-5) and autovalor[r]>=-10**(-5):
                    Lambda0 = 0 ## Se aparecer o erro que Lambda0 nao esta definido eh pq nao passou nessa linha
                else:
                    naonulos.append(autovalor[r])
            x = np.size(naonulos)
            ## passo necesario caso a gente tenha mais de um autovalor nulo
            while x < 3:
                naonulos.append(0)
                x = x + 1
            
            ## Segundo tomamos uma ordem crescente com os outros autovalores
            if naonulos[0] <= naonulos[1] <= naonulos[2]:
                Lambda1 = naonulos[0]
                Lambda2 = naonulos[1]
                Lambda3 = naonulos[2]
            elif naonulos[0] <= naonulos[2] <= naonulos[1]:
                Lambda1 = naonulos[0]
                Lambda2 = naonulos[2]
                Lambda3 = naonulos[1]
            elif naonulos[1] <= naonulos[0] <= naonulos[2]:
                Lambda1 = naonulos[1]
                Lambda2 = naonulos[0]
                Lambda3 = naonulos[2]
            elif naonulos[1] <= naonulos[2] <= naonulos[0]:
                Lambda1 = naonulos[1]
                Lambda2 = naonulos[2]
                Lambda3 = naonulos[0]
            elif naonulos[2] <= naonulos[1] <= naonulos[0]:
                Lambda1 = naonulos[2]
                Lambda2 = naonulos[1]
                Lambda3 = naonulos[0]
            elif naonulos[2] <= naonulos[0] <= naonulos[1]:
                Lambda1 = naonulos[2]
                Lambda2 = naonulos[0]
                Lambda3 = naonulos[1]
            else:
                print('erro')
            naonulos.clear()
            
            ## Calculamos os coeficientes de segunda ordem do manual MUSIC
            tau_pi = 5*etas
            delta_pipi = 4*tau_pi/3
            phi_7 = 9/(70*p)
            tau_pipi = 10*tau_pi/7
            lambda_piPI = 6/5
            tau_PI = zetas/(15*((1/3 - cs2)**2))
            delta_PIPI = 2*tau_PI/3
            lambda_PIpi = 8*tau_PI*(1/3 - cs2)/5
            
            ## Pre condicoes
            if tau_PI<=0 or tau_pi<=0:
                print('Essa célula tem tempo de relaxação negativo ou nulo, linha:', linha)
            if eta<0 or zeta<0 or tau_pipi<0 or delta_pipi<0 or lambda_PIpi<0 or delta_pipi<0 or lambda_piPI<0 or cs2<0:
                print('Essa célula tem coeficientes negativos, linha:', linha)
            q = epsilon+p+PI
            if epsilon<=0 or p<0 or q<=0 or q+Lambda1<=0 or q+Lambda2<=0 or q+Lambda3<=0:
                print('Essa célula tem problemas na pré condição, linha:', linha)
            if Lambda1>0 or Lambda3<0:
                print('Autovalores não estão na ordem correta')
            w = Lambda1+Lambda2+Lambda3
            if w>=0.1 or w<=-0.1:
                print('Autovalores não somam zero, soma: ', w)
            
            # Condições necessárias (equacao 4 do artigo)
            A  = 2*eta + lambda_piPI*PI - 0.5*tau_pipi*math.fabs(Lambda1)
            
            B  = epsilon + p + PI - (2*eta + lambda_piPI*PI)/(2*tau_pi) - tau_pipi*Lambda3/(4*tau_pi)
            
            C1 = (2*eta + lambda_piPI*PI)/(2*tau_pi) + tau_pipi*(Lambda2 + Lambda1)/(4*tau_pi)
            C2 = (2*eta + lambda_piPI*PI)/(2*tau_pi) + tau_pipi*(Lambda3 + Lambda1)/(4*tau_pi)
            C3 = (2*eta + lambda_piPI*PI)/(2*tau_pi) + tau_pipi*(Lambda2 + Lambda3)/(4*tau_pi)
            
            D1 = epsilon + p + PI + Lambda1 - (2*eta + lambda_piPI*PI)/(2*tau_pi) - tau_pipi*(Lambda2 + Lambda1)/(4*tau_pi)
            D2 = epsilon + p + PI + Lambda1 - (2*eta + lambda_piPI*PI)/(2*tau_pi) - tau_pipi*(Lambda3 + Lambda1)/(4*tau_pi)
            D3 = epsilon + p + PI + Lambda2 - (2*eta + lambda_piPI*PI)/(2*tau_pi) - tau_pipi*(Lambda1 + Lambda2)/(4*tau_pi)
            D4 = epsilon + p + PI + Lambda2 - (2*eta + lambda_piPI*PI)/(2*tau_pi) - tau_pipi*(Lambda3 + Lambda2)/(4*tau_pi)
            D5 = epsilon + p + PI + Lambda3 - (2*eta + lambda_piPI*PI)/(2*tau_pi) - tau_pipi*(Lambda1 + Lambda3)/(4*tau_pi)
            D6 = epsilon + p + PI + Lambda3 - (2*eta + lambda_piPI*PI)/(2*tau_pi) - tau_pipi*(Lambda2 + Lambda3)/(4*tau_pi)
            
            E1 = (2*eta + lambda_piPI*PI)/(2*tau_pi) + tau_pipi*Lambda1/(2*tau_pi) + (2*eta + lambda_piPI*PI + (6*delta_pipi - tau_pipi)*Lambda1)/(6*tau_pi) + (zeta + delta_PIPI*PI + lambda_PIpi*Lambda1)/tau_PI + (epsilon + p + PI + Lambda1)*cs2 
            E2 = (2*eta + lambda_piPI*PI)/(2*tau_pi) + tau_pipi*Lambda2/(2*tau_pi) + (2*eta + lambda_piPI*PI + (6*delta_pipi - tau_pipi)*Lambda2)/(6*tau_pi) + (zeta + delta_PIPI*PI + lambda_PIpi*Lambda2)/tau_PI + (epsilon + p + PI + Lambda2)*cs2
            E3 = (2*eta + lambda_piPI*PI)/(2*tau_pi) + tau_pipi*Lambda3/(2*tau_pi) + (2*eta + lambda_piPI*PI + (6*delta_pipi - tau_pipi)*Lambda3)/(6*tau_pi) + (zeta + delta_PIPI*PI + lambda_PIpi*Lambda3)/tau_PI + (epsilon + p + PI + Lambda3)*cs2
            
            F1 = epsilon + p + PI + Lambda1 - (2*eta + lambda_piPI*PI)/(2*tau_pi) - tau_pipi*Lambda1/(2*tau_pi) - (2*eta + lambda_piPI*PI + (6*delta_pipi - tau_pipi)*Lambda1)/(6*tau_pi) - (zeta + delta_PIPI*PI + lambda_PIpi*Lambda1)/tau_PI - (epsilon + p + PI + Lambda1)*cs2
            F2 = epsilon + p + PI + Lambda2 - (2*eta + lambda_piPI*PI)/(2*tau_pi) - tau_pipi*Lambda2/(2*tau_pi) - (2*eta + lambda_piPI*PI + (6*delta_pipi - tau_pipi)*Lambda2)/(6*tau_pi) - (zeta + delta_PIPI*PI + lambda_PIpi*Lambda2)/tau_PI - (epsilon + p + PI + Lambda2)*cs2
            F3 = epsilon + p + PI + Lambda3 - (2*eta + lambda_piPI*PI)/(2*tau_pi) - tau_pipi*Lambda3/(2*tau_pi) - (2*eta + lambda_piPI*PI + (6*delta_pipi - tau_pipi)*Lambda3)/(6*tau_pi) - (zeta + delta_PIPI*PI + lambda_PIpi*Lambda3)/tau_PI - (epsilon + p + PI + Lambda3)*cs2
            
            # Condições suficientes (equacao 5 do artigo)
            ALPHA = (epsilon + p + PI - math.fabs(Lambda1)) - (2*eta + lambda_piPI*PI)/(2*tau_pi) - tau_pipi*Lambda3/(2*tau_pi)
            
            BETA = 2*eta + lambda_piPI*PI - tau_pipi*math.fabs(Lambda1)
            
            GAMMA = tau_pipi - 6*delta_pipi
            
            DELTA = lambda_PIpi/tau_PI + cs2 - tau_pipi/(12*tau_pi)
            
            THETA = (4*eta + 2*lambda_piPI*PI + (3*delta_pipi+tau_pipi)*Lambda3)/(3*tau_pi) + (zeta + delta_PIPI*PI + lambda_PIpi*Lambda3)/(tau_PI) + math.fabs(Lambda1) + Lambda3*cs2 + (((12*delta_pipi - tau_pipi)/(12*tau_pi))*(lambda_PIpi/tau_PI + cs2 - tau_pipi/(12*tau_pi))*(Lambda3+math.fabs(Lambda1))**2)/(epsilon + p + PI - math.fabs(Lambda1) - (2*eta + lambda_piPI*PI)/(2*tau_pi) - tau_pipi*Lambda3/(2*tau_pi)) - (epsilon + p + PI)*(1 - cs2)  
            
            KAPPA = (2*eta + lambda_piPI*PI + (tau_pipi - 6*delta_pipi)*math.fabs(Lambda1))/(6*tau_pi) + (zeta + delta_PIPI*PI - lambda_PIpi*math.fabs(Lambda1))/tau_PI + (epsilon + p + PI - math.fabs(Lambda1))*cs2 
            
            RHO = 1-((12*delta_pipi - tau_pipi)/(12*tau_pi))*(lambda_PIpi/tau_PI + cs2 - tau_pipi/(12*tau_pi))*((Lambda3 + math.fabs(Lambda1))**2)/((2*eta + lambda_piPI*PI)/(2*tau_pi) - tau_pipi*math.fabs(Lambda1)/(2*tau_pi))**2
            
            OMEGA = (4*eta + 2*lambda_piPI*PI - (3*delta_pipi + tau_pipi)*math.fabs(Lambda1))/(3*tau_pi) + (zeta + delta_PIPI*PI - lambda_PIpi*math.fabs(Lambda1))/tau_PI + (epsilon + p + PI - math.fabs(Lambda1))*cs2 - ((epsilon + p + PI + Lambda2)*(epsilon + p + PI + Lambda3)/(3*(epsilon + p + PI - math.fabs(Lambda1)))*(1 + ((2*eta + lambda_piPI*PI + tau_pipi*Lambda3)/tau_pi)/(epsilon + p + PI - math.fabs(Lambda1))))             
            
            ## Condicoes necessárias satisfeitas:    
            ## if A>=0 and B>=0 and C1>=0 and C2>=0 and C3>=0 and D1>=0 and D2>=0 and D3>=0 and D4>=0 and D5>=0 and D6>=0 and E1>=0 and E2>=0 and E3>=0 and F1>=0 and F2>=0 and F3>=0:
            ## Condicao suficientes satisfeitas:
            ## if ALPHA>=0 and BETA>0 and GAMMA<=0 and DELTA>=0 and THETA<=0 and KAPPA>=0 and RHO>=0 and OMEGA>=0:
            
            ## Se satisfaz as equacoes necessaria e as equacoes suficientes: CAUSAL
            ## Se satisfaz as equacoes necessarias, mas viola alguma das equacoes suficientes: NADA A AFIRMAR
            ## Se nao satisfaz as equacoes necessarias e nao satisfaz alguma das equacoes suficientes: ACAUSAL
            ## Se nao satisfaz as equacoes necessarias e satisfaz as equacoes suficientes: erro pq isso nao deve existir
            if A>=0 and B>=0 and C1>=0 and C2>=0 and C3>=0 and D1>=0 and D2>=0 and D3>=0 and D4>=0 and D5>=0 and D6>=0 and E1>=0 and E2>=0 and E3>=0 and F1>=0 and F2>=0 and F3>=0:
                if ALPHA>=0 and BETA>0 and GAMMA<=0 and DELTA>=0 and THETA<=0 and KAPPA>=0 and RHO>=0 and OMEGA>=0:
                    condicao = 2 #CAUSAL, satisfaz tudo não precisa contar as eq violadas
                    epsilonCausal = epsilonCausal+epsilon
                    Causal = Causal + 1
                else:
                    condicao = 3 #Nada a afirmar, contar quais das condições suficientes são violadas
                    epsilonNada = epsilonNada+epsilon
                    Nada = Nada + 1
                    if ALPHA<0:
                        eqALPHAn = eqALPHAn+1
                    if BETA<=0:
                        eqBETAn = eqBETAn+1
                    if GAMMA>0:
                        eqGAMMAn = eqGAMMAn+1
                    if DELTA<0:
                        eqDELTAn = eqDELTAn+1
                    if THETA>0:
                        eqTHETAn = eqTHETAn+1
                    if KAPPA<0:
                        eqKAPPAn = eqKAPPAn+1
                    if RHO<0:
                        eqRHOn = eqRHOn+1
                    if OMEGA<0:
                        eqOMEGAn = eqOMEGAn+1
            else:
                if ALPHA<0 or BETA<=0 or GAMMA>0 or DELTA<0 or THETA>0 or KAPPA<0 or RHO<0 or OMEGA<0:
                    condicao = 4 #ACAUSAL, contar quais das condições suficientes e necessárias são violadas
                    epsilonAcausal = epsilonAcausal+epsilon
                    Acausal = Acausal +1
                    if A<0:
                        eqA = eqA+1
                    if B<0:
                        eqB = eqB+1
                    if C1<0 or C2<0 or C3<0:
                        eqC = eqC+1
                    if D1<0 or D2<0 or D3<0 or D4<0 or D5<0 or D6<0:
                        eqD = eqD+1
                    if E1<0 or E2<0 or E3<0:
                        eqE = eqE+1
                    if F1<0 or F2<0 or F3<0:
                        eqF = eqF+1
                    if ALPHA<0:
                        eqALPHA = eqALPHA+1
                    if BETA<=0:
                        eqBETA = eqBETA+1
                    if GAMMA>0:
                        eqGAMMA = eqGAMMA+1
                    if DELTA<0:
                        eqDELTA = eqDELTA+1
                    if THETA>0:
                        eqTHETA = eqTHETA+1
                    if KAPPA<0:
                        eqKAPPA = eqKAPPA+1
                    if RHO<0:
                        eqRHO = eqRHO+1
                    if OMEGA<0:
                        eqOMEGA = eqOMEGA+1
                else:
                    condicao = 5 #erro (nao deve passar aqui)
        
        ## Armazena no vetor resultado as cores que cada celula deve ser pintada
        resultado.append(condicao)
        #print(resultado)
        
    ## Converte o vetor resultado na MatrizResultado a qual posiciona cada célula no lugar certo do grid
    MatrizResultado = np.reshape(resultado, (280, 280)) ## linha, coluna
    
    ## Eq violadas para células com nada a afirmar
    eqsnada = eqALPHAn+eqBETAn+eqGAMMAn+eqDELTAn+eqTHETAn+eqKAPPAn+eqRHOn+eqOMEGAn
    if eqsnada != 0:
        porcentagemALPHAn = 100*eqALPHAn/eqsnada
        porcentagemBETAn = 100*eqBETAn/eqsnada
        porcentagemGAMMAn = 100*eqGAMMAn/eqsnada
        porcentagemDELTAn = 100*eqDELTAn/eqsnada
        porcentagemTHETAn = 100*eqTHETAn/eqsnada
        porcentagemKAPPAn = 100*eqKAPPAn/eqsnada
        porcentagemRHOn = 100*eqRHOn/eqsnada
        porcentagemOMEGAn = 100*eqOMEGAn/eqsnada
        nada = 1 # existem células roxas
        
    else:
        nada = 2 # não existem células roxas
    
    ## Eq violadas para células acausais
    eqsacausais4 = eqA+eqB+eqC+eqD+eqE+eqF
    eqsacausais5 = eqALPHA+eqBETA+eqGAMMA+eqDELTA+eqTHETA+eqKAPPA+eqRHO+eqOMEGA
    if eqsacausais4 !=0 and eqsacausais5 !=0:
        porcentagemA = 100*eqA/eqsacausais4
        porcentagemB = 100*eqB/eqsacausais4
        porcentagemC = 100*eqC/eqsacausais4
        porcentagemD = 100*eqD/eqsacausais4
        porcentagemE = 100*eqE/eqsacausais4
        porcentagemF = 100*eqF/eqsacausais4
        porcentagemALPHA = 100*eqALPHA/eqsacausais5
        porcentagemBETA = 100*eqBETA/eqsacausais5
        porcentagemGAMMA = 100*eqGAMMA/eqsacausais5
        porcentagemDELTA = 100*eqDELTA/eqsacausais5
        porcentagemTHETA = 100*eqTHETA/eqsacausais5
        porcentagemKAPPA = 100*eqKAPPA/eqsacausais5
        porcentagemRHO = 100*eqRHO/eqsacausais5
        porcentagemOMEGA = 100*eqOMEGA/eqsacausais5
        acausal = 1 # existem células vermelhas
        
    else:
        acausal = 2 # não existem células vermelhas
    ## Cálculo da porcentagem de energia
    #print(densityEnergy)
    #print(epsilonCausal)
    #print(epsilonNada)
    #print(epsilonAcausal)
    porcentagem_energyCausal = 100*epsilonCausal/densityEnergy
    porcentagem_energyNada = 100*epsilonNada/densityEnergy
    porcentagem_energyAcausal = 100*epsilonAcausal/densityEnergy
    ## Cálculo da porcentagem de células
    total = Causal+Acausal+Nada
    porcentagemCausal = 100*Causal/total
    porcentagemNada = 100*Nada/total
    porcentagemAcausal = 100*Acausal/total
    if porcentagem_energyCausal == 0:
        causal = 2 # não existem células causais
    else:
        causal = 1 # existem células causais

    if acausal == 1 and nada == 1 and causal == 1:
        ## Fazer a parte que pinta as células.
        mapa_cores = {
            1 : '#000000', # preto
            2 : '#069AF3',  # azul
            3 : '#7E1E9C',  # roxo
            4 : '#FF0000',  # vermelho
            #5 : '#FFFF00', # amarelo
        }
    
        #N = len(mapa_cores) 
        #valores = list(mapa_cores.keys())
        cores = list(mapa_cores.values())
    
        cmap = LinearSegmentedColormap.from_list('', cores)
         
        plt.imshow(MatrizResultado, cmap=cmap, extent=[-14,14,-14,14])
        #cbar = plt.colorbar()
    
        #largura_cor = (max(valores) - min(valores)) / N
        #posicoes = np.linspace(min(valores) + largura_cor/2, max(valores) - largura_cor/2, N)
        #cbar.set_ticks(posicoes)
        #cbar.set_ticklabels(valores)
        plt.xlabel('X (fm)', fontsize = 16) # Título dos eixos
        plt.ylabel('Y (fm)', fontsize = 16) 
        plt.title(fr'$\tau = $ {tempo} fm/c', fontsize = 16)
        plt.tight_layout()
        plt.subplots_adjust(hspace = 0, wspace = 0)
        plt.xlim(-14,14)
        plt.ylim(-14,14)
        plt.show()
    elif acausal == 1 and nada == 1 and causal == 2:
        for k in range(0,280):
            for j in range(0,280):
                if MatrizResultado[k][j] == 3:
                    MatrizResultado[k][j] = 2
                elif MatrizResultado[k][j] == 4:
                    MatrizResultado[k][j] = 3
        ## Fazer a parte que pinta as células.
        mapa_cores = {
            1 : '#000000', # preto
            2 : '#7E1E9C',  # roxo
            3 : '#FF0000',  # vermelho
            #5 : '#FFFF00', # amarelo
        }
    
        #N = len(mapa_cores) 
        #valores = list(mapa_cores.keys())
        cores = list(mapa_cores.values())
    
        cmap = LinearSegmentedColormap.from_list('', cores)
         
        plt.imshow(MatrizResultado, cmap=cmap, extent=[-14,14,-14,14])
        #cbar = plt.colorbar()
    
        #largura_cor = (max(valores) - min(valores)) / N
        #posicoes = np.linspace(min(valores) + largura_cor/2, max(valores) - largura_cor/2, N)
        #cbar.set_ticks(posicoes)
        #cbar.set_ticklabels(valores)
        plt.xlabel('X (fm)', fontsize = 16) # Título dos eixos
        plt.ylabel('Y (fm)', fontsize = 16) 
        plt.title(fr'$\tau = $ {tempo} fm/c', fontsize = 16)
        plt.tight_layout()
        plt.subplots_adjust(hspace = 0, wspace = 0)
        plt.xlim(-14,14)
        plt.ylim(-14,14)
        plt.show()
    elif acausal == 2 and nada == 1 and causal == 1:
        ## Fazer a parte que pinta as células.
        mapa_cores = {
            1 : '#000000', # preto
            2 : '#069AF3',  # azul
            3 : '#7E1E9C',  # roxo
            #4 : '#FF0000',  # vermelho
            #5 : '#FFFF00', # amarelo
        }
    
        #N = len(mapa_cores) 
        #valores = list(mapa_cores.keys())
        cores = list(mapa_cores.values())
        
        cmap = LinearSegmentedColormap.from_list('', cores)
        
        plt.imshow(MatrizResultado, cmap=cmap, extent=[-14,14,-14,14])
        #cbar = plt.colorbar()
    
        #largura_cor = (max(valores) - min(valores)) / N
        #posicoes = np.linspace(min(valores) + largura_cor/2, max(valores) - largura_cor/2, N)
        #cbar.set_ticks(posicoes)
        #cbar.set_ticklabels(valores)
        plt.xlabel('X (fm)', fontsize = 16) # Título dos eixos
        plt.ylabel('Y (fm)', fontsize = 16) 
        plt.title(fr'$\tau = $ {tempo} fm/c', fontsize = 16)
        plt.tight_layout()
        plt.subplots_adjust(hspace = 0, wspace = 0)
        plt.xlim(-14,14)
        plt.ylim(-14,14)
        plt.show()
    elif acausal == 2 and nada == 2 and causal == 1:
        ## Fazer a parte que pinta as células.
        mapa_cores = {
            1 : '#000000', # preto
            2 : '#069AF3',  # azul
            #3 : '#7E1E9C',  # roxo
            #4 : '#FF0000',  # vermelho
            #5 : '#FFFF00', # amarelo
        }
    
        #N = len(mapa_cores) 
        #valores = list(mapa_cores.keys())
        cores = list(mapa_cores.values())
        
        cmap = LinearSegmentedColormap.from_list('', cores)
        
        plt.imshow(MatrizResultado, cmap=cmap, extent=[-14,14,-14,14])
        #cbar = plt.colorbar()
    
        #largura_cor = (max(valores) - min(valores)) / N
        #posicoes = np.linspace(min(valores) + largura_cor/2, max(valores) - largura_cor/2, N)
        #cbar.set_ticks(posicoes)
        #cbar.set_ticklabels(valores)
        plt.xlabel('X (fm)', fontsize = 16) # Título dos eixos
        plt.ylabel('Y (fm)', fontsize = 16) 
        plt.title(fr'$\tau = $ {tempo} fm/c', fontsize = 16)
        plt.tight_layout()
        plt.subplots_adjust(hspace = 0, wspace = 0)
        plt.xlim(-14,14)
        plt.ylim(-14,14)
        plt.show()
    elif acausal == 2 and nada == 1 and causal == 2:
        for k in range(0,280):
            for j in range(0,280):
                if MatrizResultado[k][j] == 4:
                    MatrizResultado[k][j] = 2
        ## Fazer a parte que pinta as células.
        mapa_cores = {
            1 : '#000000', # preto
            2 : '#7E1E9C',  # roxo
            #5 : '#FFFF00', # amarelo
        }
    
        #N = len(mapa_cores) 
        #valores = list(mapa_cores.keys())
        cores = list(mapa_cores.values())
    
        cmap = LinearSegmentedColormap.from_list('', cores)
         
        plt.imshow(MatrizResultado, cmap=cmap, extent=[-14,14,-14,14])
        #cbar = plt.colorbar()
    
        #largura_cor = (max(valores) - min(valores)) / N
        #posicoes = np.linspace(min(valores) + largura_cor/2, max(valores) - largura_cor/2, N)
        #cbar.set_ticks(posicoes)
        #cbar.set_ticklabels(valores)
        plt.xlabel('X (fm)', fontsize = 16) # Título dos eixos
        plt.ylabel('Y (fm)', fontsize = 16) 
        plt.title(fr'$\tau = $ {tempo} fm/c', fontsize = 16)
        plt.tight_layout()
        plt.subplots_adjust(hspace = 0, wspace = 0)
        plt.xlim(-14,14)
        plt.ylim(-14,14)
        plt.show()
    else:
        print('Erro nas cores.')
    
    #Armazenar as informações para o gráfico das propriedades
    eixox_tempo.append(tempo)
    eixoy_energiacausal.append(porcentagem_energyCausal)
    eixoy_energianada.append(porcentagem_energyNada)
    eixoy_energiaacausal.append(porcentagem_energyAcausal)
    celulas_causal.append(porcentagemCausal)
    celulas_nada.append(porcentagemNada)
    celulas_acausal.append(porcentagemAcausal)

    resultado.clear()
    i = i+78400
