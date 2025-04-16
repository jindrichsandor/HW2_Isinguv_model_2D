##### Knihovny #####
import math
import numpy as np
import random as rand
import time
import os # zmena pracovniho adresare; vytvoreni slozky pro ukladani dat


##### Nastaveni simulace #####
n=16 # velikost mrizky: nxn
temperatures=[1.5+0.015*i for i in range(1,100)] # pole teplot, pro ktere budeme simulovat system (N.B.: 1.5 < T < 3.0)

N_thermalization=50 # pocet sweepu pouzity k termalizaci systemu
N=500 # pocet konfiguraci, ktere pouzijeme k vypoctu strednich hodnot merenych velicin

save_folder="data" # nazev slozky, do ktere se budou ukladat data ze simulace

### -------------------- ###
### Nastaveni pouzite k vygenerovani prezentovanych vysledku 
#n = 16
#temperatures=[1.5+0.01*i for i in range(1,51)]
#temperatures.extend([2.0+0.0005*i for i in range(1,601)]) 
#temperatures.extend([2.3+0.002*i for i in range(1,226)]) 
#temperatures.extend([2.6+0.01*i for i in range(1,40)]) 
#N=60000

#n = 32
#temperatures=[1.5+0.01*i for i in range(1,71)]
#temperatures.extend([2.1+0.0005*i for i in range(1,501)]) 
#temperatures.extend([2.35+0.002*i for i in range(1,151)]) 
#temperatures.extend([2.6+0.01*i for i in range(1,40)])  
#N=40000

#n = 64
#temperatures=[1.5+0.01*i for i in range(1,78)]
#temperatures.extend([2.15+0.0005*i for i in range(1,401)]) 
#temperatures.extend([2.35+0.002*i for i in range(1,76)]) 
#temperatures.extend([2.5+0.01*i for i in range(1,50)])  
#N=10000
### -------------------- ###

# Dalsi parametry simulace
J = 1 # vazbova konstanta
# B = 0 # vnejsi magneticke pole - ve vypoctech automaticky predpokladame, ze je nulove
mu = 1 # magneticky moment
kB = 1 # Boltzmannova konstanta

### Inicializace pocatecniho stavu systemu
config=[[1]*n for i in range(n)] # vygenerovani nxn mrizky se vsemi spiny nahoru
#config=[[-1]*n for i in range(n)] # vygenerovani nxn mrizky se vsemi spiny dolu
#config=[rand.choices([-1,1], k = n) for i in range(n)] # vygenerovani nxn mrizky s nahodnou konfiguraci spinu

##### Definice funkci #####
def neighbours(configuration, i,j): # funkce k urceni spinu sousedicich se spinem s_ij (funkce zajistuje splneni periodicke okrajove podminky)
    # spin na pozici (i-1,j) (up)
    s_U=configuration[i-1][j] # N.B.: pro i=0 mame [-1][j], coz Python automaticky chape jako element na pozici [n-1][j]
    # spin na pozici (i+1,j) (down)
    if i==n-1: s_D=configuration[0][j]  # N.B.: n je globalni promenna
    else: s_D=configuration[i+1][j] 
    # spin na pozici (i,j-1) (left)
    s_L=configuration[i][j-1] 
    # spin na pozici (i,j+1) (right)
    if j==n-1: s_R=configuration[i][0]
    else: s_R=configuration[i][j+1] 
    
    return [s_U, s_D, s_L, s_R]


##### Simulace #####
start_time = time.time() # zaznam o case zacatku simulace

N_s=n*n # pocet spinu na mrizce

m_data=[0]*N # pole pro ukladani namerenych hodnot energie H
H_data=[0]*N # pole pro ukladani namarenych hodnot magnetizace m
result=[[0,0,0,0,0] for i in range(len(temperatures))] # tabulka pro ulozeni vysledku simulace (tj. prumernych hodnot vnitrni energie U, merne tepelne kapacity c, magnetizace M a magneticke susceptibility chi pro danou teplotu systemu)

# pomocna pole pro ukladani mezivysledku pri prumerovani
x1=[0]*N
x2=[0]*N


for T_i in range (len(temperatures)): # cyklus pres vsechny teploty
    T=temperatures[T_i]
    beta=1/(kB*T)

    ##### Termalizace #####
    for sweep in range(N_thermalization): # cyklus pro sweepy
        for i in range(n): # cyklus pres "radky" mrizky
            for j in range(n): # cyklus pres "sloupce" mrizky
                s_ij=config[i][j]
                
                s_X=neighbours(config,i,j) # spiny sousedu
                
                dE=2*J*s_ij*(s_X[0] + s_X[1] + s_X[2] + s_X[3]) # zmena energie v dusledku otoceni spinu s_ij -> -s_ij
                
                if dE<0: 
                    config[i][j]=-s_ij # prijmuti nove konfigurace
                elif math.exp(-beta*dE) >= rand.random(): 
                    config[i][j]=-s_ij # prijmuti nove konfigurace
    

    ##### Mereni #####
    H=0 # energie systemu
    m=0 # magnetizace systemu
    # vypocet H a m
    for i in range (n): # cyklus pres "radky" mrizky
        for j in range(n): # cyklus pres "sloupce" mrizky
            # N.B.: nasobeni energie vazbovou konstantou J a nasobeni magnetizace momentem mu provadime az pozdeji
            H=H-config[i][j]*(config[i-1][j]+config[i][j-1]) # N.B.: prochazime celou mrizku, takze v sume pres sousedy staci vzdy zapocist pouze souseda vlevo od s_ij a nad s_ij (pokud bychom zapocetli pro kazdy spin vsechny sousedy, vysel by nam dvakrat vetsi energie systemu)
            m=m+config[i][j] 
    H=H*J
    
    for sweep in range(N): # cyklus pro sweepy
        for i in range(n): # cyklus pres "radky" mrizky
            for j in range(n): # cyklus pres "sloupce" mrizky
                s_ij=config[i][j]
                
                s_X=neighbours(config,i,j) # spiny sousedu
                # zmena energie v dusledku otoceni spinu s_ij -> -s_ij
                dE=2*J*s_ij*(s_X[0] + s_X[1] + s_X[2] + s_X[3]) 
                
                if dE<0: 
                    config[i][j]=-s_ij # prijmuti nove konfigurace
                    H=H+dE # nova hodnota energie
                    m=m-2*s_ij  # nova hodnota magnetizace
                elif math.exp(-beta*dE) >= rand.random(): 
                    config[i][j]=-s_ij # prijmuti nove konfigurace
                    H=H+dE # nova hodnota energie
                    m=m-2*s_ij  # nova hodnota magnetizace    
        
        H_data[sweep]=H
        m_data[sweep]=m*mu
        

          
    ##### Prumerovani #####         
    # vypocet vnitrni energie U, tepelne kapacity c, magnetizace M a magneticke susceptibility chi (vse na jednu castici)
    for i in range (N):
        x1[i]=m_data[i]**2
        x2[i]=H_data[i]**2
            
    y1=math.fsum(m_data)/N 
    y2=math.fsum(x1)/N
    
    M=y1/N_s # prumerna magnetizace (na castici)
    chi=beta*(y2-y1**2)/N_s # susceptibilita (na castici)
            
    y1=math.fsum(H_data)/N 
    y2=math.fsum(x2)/N
    
    U=y1/N_s # prumerna vnitrni energie castice
    c=(y2-y1**2)/(N_s*kB*T**2) # prumerna tepelna kapacita (na castici)
        
        
    result[T_i]=[T,M,chi,U,c] # ulozeni namerenych dat do souboru s vysledky 
    
end_time=time.time() # zaznam o case ukonceni simulace


##### Ulozeni dat ze simulace #####
file_directory=os.path.dirname(os.path.abspath(__file__)) # absolutni cesta k adresari s python skriptem
os.chdir(file_directory) # zmena pracovniho adresare

# vytvoreni slozky pro ukladani dat (pokud neexistuje)
if os.path.exists(save_folder)==False:
    os.mkdir(save_folder)
    
# Ulozeni vysledku simulace do datoveho souboru #####
file_name=save_folder+"/"+"grid_" + str(n) + "x" + str(n) + "_result.dat"
np.savetxt(file_name,result,fmt="%1.30s",delimiter="  ", header="T  M chi U c", comments="")    

# Ulozeni logu ze simulace do datoveho souboru #####
file_name=save_folder+"/"+"grid_" + str(n) + "x" + str(n) + "_log.dat"
log_file=open(file_name, "w")
log_file.write("N_thermalization = " + str(N_thermalization) + " ; pocet sweepu pouzity k termalizaci systemu\n")
log_file.write("N = " + str(N) + " ; pocet konfiguraci pouzitych k prumerovani\n")
log_file.write("t = " + str(end_time - start_time) + " s ; celkovy cas behu simulace")
log_file.close


