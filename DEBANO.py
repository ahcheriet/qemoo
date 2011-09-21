# -*- coding: cp1252 -*-
#import pdb
from math import *
from copy import *
from threading import *
from random import *
from Queue import *
from scipy.misc.common import *

filein = open('workfile', 'w')
fval = 1 
m = 30  #  nombre des sous-var x(x0,..xi)
class qbit:
    " Un qubit avec a est alpha et le b est beta"
    def __init__(self):
        self.a = 1/sqrt(2)
        self.b = 1/sqrt(2)
        self.Iterate = True
    def observe(self,Iterate) :
        r = random()	
        if r < (self.b * self.b ):
            return (1)
        else:
            return (0)
    def qgate( self,n,N ):
        angle = self.CalculAngle(n,N)
        #print angle
        e = 0.0001
        a_1 = cos(angle)* self.a - sin(angle) * self.b
        b_1 = sin(angle)* self.a + cos(angle)* self.b
        self.a = a_1
        self.b = b_1
	if (self.a)**2 <= e and (1- e <= (self.b)**2):
		self.b = sqrt(1-e)
	if (self.b)**2 <= e and (1- e <= (self.a)**2):
		self.b = sqrt(e)
		self.a = sqrt(1- e)

    def factoriel( self,n ):
        i = 0
        fact = 1
        while i < n:
            fact = fact *( i+1)
            i = i+1
        return fact
    def Combinaison( self,n , p):
        c = self.factoriel(n)/(self.factoriel(p) * self.factoriel(n-p))
        return c
    def binomial( self,x,N ,p,q):
	b1 = pow(p,x)
	b2 = pow(q,N-x)
	c = comb(N,x)
        b =  c * b1 * b2
        return b
    def CalculAngle(self,n,N):
        b = self.binomial(n,N,pow(self.a,2),pow(self.b,2))
        if ( N == 0 ):
            b = self.b
        else:
            print n, N
            b = n/N;
        Angle = atan(sqrt(b)/sqrt(1-b))-atan(self.a/self.b)
        #filein.write(str(Angle)+'\n')
        return Angle
    
            
class qoctet(qbit):
    " Un octet constitue de N qbits "
    def __init__( self, Nbits ):
        self.octet = []
        self.Nbits = Nbits
        i = 1
        while i <= self.Nbits :
            Qbit = qbit()
            self.octet.append(Qbit)
            i = i + 1
    def observe( self ):
		i = 0
		Iterate = False
		self.result = []
		while i < self.Nbits :
			if Iterate:
				Iterate = False
			else:
				Iterate = True
			self.result.append(self.octet[i].observe(Iterate))
			i = i + 1
		i = 0
		self.RString=''
		while i < self.Nbits:
			self.RString = self.RString + str(self.result[i])
			i = i+1
		return (self.result)

    def Update( self , T ,N):
        i = 0
        while i < self.Nbits :
            n = T[i]
            self.octet[i].qgate( n,N )
            i = i + 1
    
            
class QGA( qoctet ):
    def __init__( self , Nbits ):
        self.Population = []
        self.Gen = 0
        self.Population = qoctet(Nbits)

      
    def Update( self , Population ,Angle):
        Population.Update( Angle )

class QMOO(QGA):
    def __init__(self , Contrainte, FonctionObject,Nbits ):
        self.qga = QGA(  Nbits )
        self.ReadData(Contrainte, FonctionObject)
        self.Nbcontrainte = len(Contrainte)
        self.fileNonData = open("Nondata1.txt","w")
        self.fileData = open("data1.txt","w")
        self.CodedX = []
        self.ListDominated = []
        self.lock = Lock()
        self.lockList = Lock()
        self.lock2 = Lock()
	self.Nbits = Nbits
        
    def ReadData( self, Contrainte, FonctionObject ):
        self.Contrainte = []
        self.fonctionObject = []
        self.Population = []
        i = 0
        while i < len(Contrainte):
            self.Contrainte.append(eval(Contrainte[i]))
            i = i+1
        i = 0
        while i < len(FonctionObject):
            self.fonctionObject.append(FonctionObject[i])
            i = i +1   
        self.Population = []
        for i in range(m):
            self.Population.append(copy(self.qga.Population))

        
    def Decode( self, x ):
        i = 0
        String1=''
        String2=''
        if ( fval == 0 ):
            while i < len(x) - 10 :
                String1 = String1 + x[i]
                i = i+1
            while i < len(x):
                String2 = String2 + x[i]
                i = i+1
        else:
            while i < len(x) :
                String2 = String2 + x[i]
                i = i+1
            String1 = '0'
        Virgule = str(int(String2,2))
        #Virgule = '0.'+Virgule
        nmb = pow(2, self.Nbits)
        Result = int(String1,2)+ float(Virgule)
        #print nmb,  Result ,   float(Result)/nmb
        return float(Result)/nmb
    
    def IsDominate( self,Lfonction, x, y , MinOrMax):
        fonction = eval(Lfonction)
        if MinOrMax == 1:
            if min(fonction(x),fonction(y)) == fonction(x):
                return True
            else:
                return False
        else:
            if max(fonction(x),fonction(y)) == fonction(x):
                return True
            else:
                return False           
    def Dominate( self , queueD,i):
         j = 0
         while j < len(self.X):
            if ( i != j ):
		k = 0
	        CountI = 0
	        CountJ = 0		
		while k < len ( self.fonctionObject ):                    
			  if  self.IsDominate( self.fonctionObject[k], self.X[i], self.X[j],1): # minimisation problem
				    CountI = CountI + 1
			  if  self.IsDominate( self.fonctionObject[k], self.X[j], self.X[i],1):
				    CountJ = CountJ + 1                               
			  k = k +1
		if CountI == len( self.fonctionObject ):
			  self.ListDominated.append(copy(self.X[j]))
			  self.lockList.acquire()
			  queueD.append(copy(self.X[j]))
			  self.lockList.release()
		if CountJ == len( self.fonctionObject ):
			  self.lockList.acquire()
			  queueD.append(copy(self.X[i]))
			  self.lockList.release()
			  self.ListDominated.append(copy(self.X[i]))
            j = j + 1

    def RealisableSpace( self ,Population ):
        i = 0
        Dominated = 0
        Xcoded = []
        Ycoded = []
        Y = []
        Ystring = []
        Niteration = 100
       # X_ij = Population[0].observe()
        for i in range(m):
            Population[i].observe()
            Xdecoded = (self.Decode(copy(Population[i].RString)))
            Y.append(Xdecoded )
            Ystring.append( copy(Population[i].RString))
        if Y not in self.X and self.InContrainte(Y):
            self.X.append(Y)
            self.CodedX.append(Ystring)
    
    def Tofile( self , File , X ,ifile , AorW):
        if AorW == 1:
            fileData = open(File+str(ifile)+".txt","a")
        else:
            fileData = open(File+str(ifile)+".txt","w")            
        i = 0
        while i < len( X ):
            fonction1 = eval(self.fonctionObject[0])
            fonction2 = eval(self.fonctionObject[1])
            fileData.write(str(fonction1(X[i]))+'\t'+str(fonction2(X[i]))+'\n')
            i = i +1
        fileData.close()            
    def Evaluate( self ,queueD, queueCodeN , queueN, ifile ,Population ):
        i = 0
        self.X = []
        self.ListNonDominated = []
        self.lockList.acquire()
        self.ListCodedNonDominated =[]# copy(queueCodeN)
        self.ListDominated = copy(queueD)
        self.lockList.release()
        while i < 100:   #????
                self.RealisableSpace( Population )
                j = len(self.X)-1
                self.Dominate(queueD,j);
                i = i +1
        self.RealisableSpace( Population )
        j = len(self.X)-1
        self.Dominate(queueD,j);
        self.lock.acquire()
        self.Tofile('data',self.X,ifile , 1)
        self.lock.release()
        fonction1 = eval(self.fonctionObject[0])
        fonction2 = eval(self.fonctionObject[1])
        i = 0
        while i < len(self.X):
            if self.X[i] not in self.ListDominated:
                self.ListNonDominated.append(copy(self.X[i]))
                self.ListCodedNonDominated.append(copy(self.CodedX[i]))
                self.lockList.acquire()
                queueCodeN.append(copy(self.CodedX[i]))
                queueN.append(copy(self.X[i]))                
                self.lockList.release()
            i = i +1
    
    def InContrainte(self, x):
        i = 0
        while i < self.Nbcontrainte:
            if self.Contrainte[i](x) == False:
                return 0
            i = i +1
        return 1

    def Interferate( self, ListC , Pop ):
        N = len(ListC)
        JJ = 0
        for index in range(m):
            JJ = 0
            Tx = []
            while JJ < self.Nbits:
                Tx.append(0)
                JJ = JJ + 1
            i = 0 #ekteblo hnaya Tx et Ty kamel raaak taaraf hathi hiya sahiha
            j = 0 #egsem ko
            while i < N:
                j = 0
                while j < self.Nbits:
                    if ListC[i][index][j] == '1':
                        Tx[j] = Tx[j] + 1
                        #print 'j',j
                    j = j + 1    
                i = i + 1 #ziide eckteb
            if N != 0:
                Pop[index].Update(Tx,N)
        
class Qthread( QMOO):
    def __init__( self  , Nthreads , Contrainte, FonctionObject , Ngen , Nbits):
        self.Contrainte = Contrainte
        self.FonctionObject = FonctionObject
        self.Nbits = Nbits
        self.lock = Lock()
        self.lockFiltrate = Lock()
        self.dominated = []
        self.CodeNondominated = []
        self.Nondominated = []        
        self.Nthreads = Nthreads
        self.threads = []
        self.Qlist = []
        self.Qoom = Queue()
        self.Ngen = Ngen
        i = 0
        while i < self.Nthreads :
            qm = QMOO( Contrainte , FonctionObject ,  Nbits)
            self.Qlist.append(copy(qm))
            i = i +1
        i = 0
        while i < self.Nthreads :
            thred = Thread(target = self.RunOnethread,args =() )
            self.threads.append(copy(thred))
            i = i+1
        i = 0
        while i < len( self.Qlist ):
            self.Qlist[i].Evaluate( self.dominated, self.CodeNondominated,self.Nondominated,1 ,self.Qlist[i].Population )
            self.Qoom.put(copy(self.Qlist[i]))
            i = i+1
            
        
    def RunOnethread( self  ):
        i = 0
        Qlist = self.Qoom.get()
        Ngen = self.Ngen
        while i < Ngen :
            self.Filterate()
            Qlist.Interferate( Qlist.ListCodedNonDominated , Qlist.Population  )
            Qlist.Evaluate( self.dominated, self.CodeNondominated,self.Nondominated, 1,Qlist.Population )
            i = i + 1
    def Filterate( self ):
        i = 0
        self.lockFiltrate.acquire()
        while i < len(self.Nondominated)-1:
                j = i+1
                while j < len(self.Nondominated):
                        CountI = 0
                        CountJ = 0
                        k = 0
                        while k < len ( self.Qlist[0].fonctionObject ):                    
                            if  self.Qlist[0].IsDominate( self.Qlist[0].fonctionObject[k], self.Nondominated[i], self.Nondominated[j],1):
                                CountI = CountI + 1
                            if  self.Qlist[0].IsDominate( self.Qlist[0].fonctionObject[k], self.Nondominated[j], self.Nondominated[i],1):
                                CountJ = CountJ + 1                               
                            k = k +1
                        if CountI == len( self.Qlist[0].fonctionObject ):
                                x = self.Nondominated[j]
                                j = j - 1
                                self.Nondominated.remove(x)
                        if CountJ == len( self.Qlist[0].fonctionObject ):
                                x = self.Nondominated[i]
                                j = j - 1
                                i = i - 1
                                self.Nondominated.remove(x)
                                break
                        j = j + 1
                i = i +1
        self.lockFiltrate.release()

        
    def RunThreads( self ):
        i = 0
        while i < self.Nthreads :
            self.threads[i].start()
            i = i+1
        i = 0
        while i < self.Nthreads :
            self.threads[i].join()
            i = i+1
        self.Filterate()
        self.Qlist[0].Tofile('Nondata',self.Nondominated,1 , 2)
        
        
       
##Contrainte = []
##FonctionObject = []
##Contrainte.append('lambda x: x[0]>=-5 and x[0] <=10 and x[1]>=-5 and x[1] <=10')
##FonctionObject.append('lambda x:x[0]*x[0]+x[1]*x[1]')
##FonctionObject.append('lambda x:(x[0]-5)*(x[0]-5)+(x[1]-5)*(x[1]-5)')
##
##Contrainte1 = []
##FonctionObject1 = []
##Contrainte1.append('lambda x: x[0]>= x[1] + 0.04 and x[0] <=0.1')
##FonctionObject1.append('lambda x:x[0]*x[0]-x[1]*x[1]')
##FonctionObject1.append('lambda x:1000 + ( 0.0001 / ( 192+200000+ (x[0]**4 - x[1]**4) /12 ) )')
##
##Contrainte2 = []
##FonctionObject2 = []
##Contrainte2.append('lambda x: x[0]>= 0 and x[1] >= 0 and x[0] <=1 and x[1] <=1')
##FonctionObject2.append('lambda x:x[0]')
##FonctionObject2.append('lambda x:(2- exp(-((x[1]-0.2)/0.04)**2)- 0.8*exp(-((x[1]-0.6)/0.4)**2))* (1/x[0])')
##
##Contrainte4 = []
##FonctionObject4 = []
##Contrainte4.append('lambda x: x[0]>= -100  and x[0] <=100')
##FonctionObject4.append('lambda x:x[0]*x[0]')
##FonctionObject4.append('lambda x:(x[0]-2)*(x[0]-2)')
##
##
##
##
##
##Contrainte3 = []
##FonctionObject3 = []
##Contrainte3.append('lambda x: sum(map(lambda a,y: a*y, x[0],W1)) <= C1')
##Contrainte3.append('lambda x: sum(map(lambda a,y: a*y, x[0],W2)) <= C2')
##FonctionObject3.append('lambda x: sum(map(lambda a,y: a*y, x[0],P1))')
##FonctionObject3.append('lambda x: sum(map(lambda a,y: a*y, x[0],P2))')
##
Contrainte01 = []
FonctionObject01 = []
Contrainte01.append('lambda x: True')
FonctionObject01.append('lambda x:x[0]')
FonctionObject01.append('lambda x:(1+9*sum(map(lambda a:a/(m-1),x[1:]))*(1-sqrt(x[0]/(1+9*sum(map(lambda a:a/(m-1),x[1:]))))))')

Contrainte02 = []
FonctionObject02 = []
Contrainte02.append('lambda x: True')
FonctionObject02.append('lambda x:x[0]')
FonctionObject02.append('lambda x:(1+9*sum(map(lambda a:a/(m-1),x[1:]))*(1-pow((x[0]/(1+9*sum(map(lambda a:a/(m-1),x[1:])))),2)))')


Contrainte03 = []
FonctionObject03 = []
Contrainte03.append('lambda x: True')
FonctionObject03.append('lambda x:x[0]')
FonctionObject03.append('lambda x:(1+9*sum(map(lambda a:a/(m-1),x[1:]))*(1-sqrt(x[0]/(1+9*sum(map(lambda a:a/(m-1),x[1:]))))-(x[0]/(1+9*sum(map(lambda a:a/(m-1),x[1:]))))*sin(10*pi*x[0])))')

Contrainte04 = []
FonctionObject04 = []
Contrainte04.append('lambda x: True')
FonctionObject04.append('lambda x:x[0]')
FonctionObject04.append('lambda x:(1+10*(m-1)+sum(map(lambda a:a*a-10*cos(4*pi*a),x[1:])))*(1-sqrt(x[0]/(1+10(m-1)+sum(map(lambda a:a*a-10*cos(4*pi*a),x[1:])))))')

Contrainte05 = []
FonctionObject05 = []
Contrainte05.append('lambda x: True')
FonctionObject05.append('lambda x:x[0]')
FonctionObject05.append('lambda x:(1+10*(m-1)+sum(map(lambda a:a*a-10*cos(4*pi*a),x[1:])))*(1-sqrt(x[0]/(1+10(m-1)+sum(map(lambda a:a*a-10*cos(4*pi*a),x[1:])))))')

Contrainte06 = []
FonctionObject06 = []
Contrainte06.append('lambda x: True')
FonctionObject06.append('lambda x:1-exp(-4*x[0])*pow(sin(6*pi*x[0]),6)')
FonctionObject06.append('lambda x:(1+9*pow(sum(x[1:])/(m-1),0.25))*(1-pow( ( 1-exp(-4*x[0])*pow(sin(6*pi*x[0]),6) )/( 1+9*pow(sum(x[1:])/(m-1),0.25 )) ,2))')

#    Qthread( Nthreads , Contrainte, FonctionObject , Ngen , Nbits)
#m = 10  sur fnct 6
qm = Qthread( 1      ,Contrainte03 , FonctionObject03 ,10,  60)   # changer la valeur de m
qm.RunThreads()

