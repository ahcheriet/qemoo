##from numpy.random import binomial
from scipy.misc.common import *
from copy import *
from math import *
from random import *
MaxN = 0
Reuse = 0
fileprob = open("prob.txt","w")            
filen = open("probn.txt","w")            
fileN = open("probN.txt","w")            
m = 30

def observe(p) :
    r = random()	
    if r < p:
        return (1)
    else:
        return (0)
    
def Tofile( File , X ,ifile , AorW):
    if AorW == 1:
        fileData = open(File+str(ifile)+".txt","a")
    else:
        fileData = open(File+str(ifile)+".txt","w")            
    i = 0
    while i < len( X ):
        fonction1 = eval(fonctionObject[0])
        fonction2 = eval(fonctionObject[1])
        fileData.write(str(fonction1(X[i]))+'\t'+str(fonction2(X[i]))+'\n')
        i = i +1
    fileData.close()            


Contrainte03 = []
FonctionObject03 = []
Contrainte03.append('lambda x: True')
FonctionObject03.append('lambda x:x[0]')
FonctionObject03.append('lambda x:(1+9*sum(map(lambda a:a/(m-1),x[1:]))*(1-sqrt(x[0]/(1+9*sum(map(lambda a:a/(m-1),x[1:]))))-(x[0]/(1+9*sum(map(lambda a:a/(m-1),x[1:]))))*sin(10*pi*x[0])))')
fonction1 = eval(FonctionObject03[0])
fonction2 = eval(FonctionObject03[1])
qoctet = []
for i in range(30):
    qoctet.append([0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5])
    
bestp = deepcopy(qoctet)
octet1 = []
for i in range(30):
    octet1.append([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])

for Ind_gen in range(500):
    X = []
    x_array = []
    pop = []
    New_pop = []
    Dominated = []
    NonDominated = []
    Tind = 100
    Tbit = 20
    print "Je suis dans l'iteration ",Ind_gen
    for ind in range(Tind):
        X_x = []
        for i_x in range(30):
            for i in range(Tbit):
                octet1[i_x][i] =str( observe(qoctet[i_x][i])) # bien concevoir les contraintes
            float_ele = int("".join(octet1[i_x][1:]),2)* 1.0 / pow(10,6)
            int_ele = int("".join(octet1[i_x][:1]),2)    
            float_ele = float_ele + int_ele
            b = 0
            i = 1
            j = 1
            while i < Tbit:
                b = b + int(octet1[i_x][i])*pow(2,(j)*-1)     # conversion
                i = i + 1
                j = j + 1
            float_ele = b + int_ele;
            if ( octet1[i_x][0] == '1' ):
                float_ele = float_ele*1.0
            X_x.append(float_ele)
##        if ( float_ele <= 100 ) and (float_ele  >= -100 ): # contrainte
       ## if X_x not in X:
        tmpoctet = []
        tmpoctet = deepcopy(octet1)
        pop.append(tmpoctet) # il y a un problem au niveau de lecture
        X.append(X_x)
    for i in range(len(X)):
        if ( i not in Dominated) :
            for j in range(len(X)):
                if i != j and (j not in Dominated) :
                    if ( fonction1(X[i]) < fonction1(X[j]) ) and ( fonction2(X[i]) < fonction2(X[j]) ): # fonction objective
                        Dominated.append(j)
                    if ( fonction1(X[j]) < fonction1(X[i]) ) and ( fonction2(X[j]) < fonction2(X[i]) ): 
                        Dominated.append(i)
    for i in range(len(X)):
        if i not in Dominated:
            NonDominated.append(i)
    for i in range(len(NonDominated)):
        New_pop.append(deepcopy(pop[NonDominated[i]]))
    if ( MaxN < len(NonDominated) ): # c'est mieux de sauvgarder les meilleurs non la quantity
        Reuse = 0
        MaxN = len(NonDominated)
        bestp = deepcopy(qoctet)
    else:
        Reuse = Reuse + 1
    tt = zip(*New_pop) ##transposer des ind nondominated
    for i_x in range(m):
        ti = zip(*tt[i_x])
        if ( Ind_gen ==0 ):
            print ti
        for i in range(len(ti)):
            qoctet[i_x][i] = (int(ti[i].count('1'))*1.0/len(NonDominated))*1.0 # il faut recalculer qoctet
            if ( Ind_gen ==0 ):
                print qoctet[i_x][i], len(tt),len(ti)
    if ( Reuse >= 50 ):
            Reuse = 0
            qoctet = deepcopy(bestp)
    fileprob.write(str(Ind_gen)+'\t'+str(qoctet[10])+'\n')
    fileN.write(str(Ind_gen)+'\t'+str(len(NonDominated))+'\n')
    filen.write(str(Ind_gen)+'\t'+str(tt[10].count('1'))+'\n')
##    for i in range(len(NonDominated)):
##        print X[NonDominated[i]]   ## put one the functions test here with Tofile procedure
        
fileData = open("meilleur.txt","w")            
i = 0
while i < len( NonDominated ):
    fonction1 = pow(X[NonDominated[i]],2)
    fonction2 = pow(X[NonDominated[i]]-2,2)
    fileData.write(str(fonction1)+'\t'+str(fonction2)+'\n')
    i = i +1
fileData.close()            
fileData = open("donnee.txt","w")     ## you have to resolve the probalbilty variation       
i = 0
while i < len( X ):
    fonction1 = pow(X[i],2)
    fonction2 = pow(X[i]-2,2)
    fileData.write(str(fonction1)+'\t'+str(fonction2)+'\n')
    i = i +1

    
Dominated = []
NonDominated = []
Tind = 1000
Tbit = 20

print qoctet

    
for ind in range(Tind):
    for i in range(Tbit):
        octet1[i] =str( observe(qoctet[i])) # j'ai change sa qoctet par bestp
    float_ele = int("".join(octet1[3:]),2)* 1.0 / pow(10,6)
    int_ele = int("".join(octet1[1:3]),2)    
    float_ele = float_ele + int_ele
    b = 0
    i = 3
    j = 1
    while i < len(octet1):
        b = b + int(octet1[i])*pow(2,(j)*-1)     # conversion
        i = i + 1
        j = j + 1
    float_ele = b + int_ele;
    if ( octet1[0] == '1' ):
        float_ele = float_ele*-1.0
    #print float_ele
    if ( float_ele <= 3 ) and (float_ele  >= -1 ):
        if float_ele not in X:
            pop.append(deepcopy(octet1))
            X.append(float_ele)
for i in range(len(X)):
    for j in range(len(X)):
        if i != j and (j not in Dominated) :
            if ( pow(X[i],2) < pow(X[j],2) ) and ( pow(X[i]-2,2) < pow(X[j]-2,2) ):
                Dominated.append(j)
            if ( pow(X[j],2) < pow(X[i],2) ) and ( pow(X[j]-2,2) < pow(X[i]-2,2) ):
                Dominated.append(i)
for i in range(len(X)):
    if i not in Dominated:
        NonDominated.append(i)
fileData2 = open("meilleur2.txt","w")            
i = 0
while i < len( NonDominated ):
    fonction1 = pow(X[NonDominated[i]],2)
    fonction2 = pow(X[NonDominated[i]]-2,2)
    fileData2.write(str(fonction1)+'\t'+str(fonction2)+'\n')
    i = i +1
fileData2.close()            

fileData.close()            
fileprob.close()            