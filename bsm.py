import numpy as np
import math
from math import *
from scipy.stats import norm

class BSM:
#class constructor; input: S,vol,r,q
    def __init__(self,S,vol,r,q):
        self.S, self.vol, self.r, self.q = S, vol, r, q

        
        """class constructor; input: S,vol,r,q"""
# alternative constructor; take data from dictionary
    @classmethod
    def from_dict( cls_, d ):
        S=d['S'] if 'S' in d else 100
        vol=d['vol'] if 'vol' in d else 0.3
        r=d['r'] if 'r' in d else 0
        q=d['q'] if 'q' in d else 0
        return cls_(S,vol,r,q)
    
    def d(self,K,T):
        lg =  np.log(self.S / K)
        volT = (self.vol * np.sqrt(T));
        diff = (self.r-self.q)*T
        d2 = (lg + diff)/volT - volT/2; 
        d1 = (lg + diff)/volT + volT/2;
        d3 = (lg - diff)/volT + volT/2;
        return [d1, d2,d3,volT]
            



    def Call(self,K,T):
        """Calculate the price of Call Option. Inputs: Strike, Time to maturity."""
        d = self.d(K,T)
        return np.exp(-self.q * T)*self.S*norm.cdf(d[0]) - K*np.exp(-self.r*T)*norm.cdf(d[1])
    
    
    def Call_Vol(self,vol,K,T):
        """Calculate the price of Call Option. Inputs: Strike, Time to maturity."""
        self.vol = vol
        return self.Call(K,T)
    

    def Put(self,K,T):
        """Calculate the price of Put Option. Inputs: Strike, Time to maturity."""
        d = self.d(K,T)
        return -np.exp(-self.q * T) *self.S * norm.cdf(-d[0]) + K * np.exp(-self.r * T) * norm.cdf(-d[1])

    
    def Put_Vol(self,vol,K,T):
        """Calculate the price of Call Option. Inputs: Strike, Time to maturity."""
        self.vol = vol
        return self.Put(K,T)    

    def DeltaCall(self,K,T):
        d = self.d(K,T)
        return np.exp(-self.q * T) * norm.cdf(d[0])

    def DeltaPut(self,K,T):
        d = self.d(K,T)
        return np.exp(-self.q * T) * (norm.cdf(d[0])-1)
    

    def Gamma(self,K,T):
        d = self.d(K,T)
        return norm.pdf(d[0])/(self.S*d[2])
    
    
    
    def ThetaCall(self,K,T):
        d = self.d(K,T)
        return -self.S*norm.pdf(d[0])*self.vol/(2*np.sqrt(T)) - self.r*K*np.exp(-self.r * T)*norm.cdf(d[1])
    
    
        
    def ThetaPut(self,K,T):
        d = self.d(K,T)
        return -self.S*norm.pdf(d[0])*self.vol/(2*np.sqrt(T)) + self.r*K*np.exp(-self.r * T)*norm.cdf(-d[1])


    def Pr(self,K,T):
        """Calculate the probability of S(T)<K. """
        d = self.d(K,T)
        return 1-norm.cdf(d[1])
    
    def BarrierHit(self,K,T):
        """Calculate the probability of S(T)<K before T. """
        Up = True if K>self.S else False
        v = 2*self.r/self.vol**2-1
        d = self.d(K,T)
        term1 = pow(K/self.S,v)
        term2 = norm.cdf(d[1])
        term3 = norm.cdf(d[2])
#        res = norm.cdf(d[1])*(pow(K/self.S,v)+1)
        res  = (term2 + term3*term1) if Up else ((1-term2) + (1-term3)*term1)
        return res


    def Vega(self,K,T):
        """Calculate the vega of Call option. """
        d1=(np.log(self.S/K) + (self.r -self.q+ self.vol**2 / 2) * T)/(self.vol * np.sqrt(T))
        vega = np.exp(-self.q * T)*self.S*np.sqrt(T)*norm.pdf(d1);
        return [vega, d1] 
    
    def ImpliedVol(self, option_price, C_P, K, T):
        """Calculate Implied Vol from option price """
        price = {"Call": self.Call_Vol, "Put": self.Put_Vol}
        l_Vol = 0
        r_Vol = .01
        while price[C_P](r_Vol,K,T) < option_price:
            l_Vol = r_Vol
            r_Vol *= 2
#        print(l_Vol,r_Vol, price[C_P](r_Vol,K,T), option_price)
    
        m_Vol = (l_Vol + r_Vol) / 2
        while (r_Vol - l_Vol) > 0.00001:
            if price[C_P](m_Vol,K,T) < option_price:
                l_Vol = m_Vol
            else:
                r_Vol = m_Vol
            m_Vol = (l_Vol + r_Vol) / 2
#            print(m_Vol)
        return m_Vol

    
    
    
    

def nCr(n,r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)


class Binomial:
    
# Class constructor    
    def __init__(self,S,T,r,q,vol,n):
#variables
        self.S = S
        self.n = n
        self.h = T/n
        self.u = math.exp(sqrt(self.h)*vol)
        self.d = 1/self.u
        self.df = exp(-r*self.h)
        self.dd = exp(-q*self.h)
        self.p = (exp((r-q)*self.h)-self.d)/(self.u-self.d)
#matrices = trees        
        self.TreeS    = [[S*self.u**(j)*self.d**(i-j) for j in range(n+1)] for i in range(n+1)]
        self.TreeProb = [[nCr(i,j)*self.p**(j)*(1-self.p)**(i-j) for j in range(n+1) if j<=i] for i in range(n+1)]
        self.TreeValue     = [[0 for j in range(n+1) if j<=i] for i in range(n+1)]
        self.TreePositionValue     = [[0 for j in range(n+1) if j<=i] for i in range(n+1)]
# paths         
        self.paths     = [[]]
        self.pathsProb = []
        # Build Tree for S    

    def PositionAdd(self,plusMinus):
        for i in range(self.n+1):
            for j in range(self.n+1):
                if j<=i:
                    if plusMinus == 'plus':
                        self.TreePositionValue[i][j]+=  self.TreeValue[i][j]     
                    else: 
                        self.TreePositionValue[i][j]-=  self.TreeValue[i][j] 
                        
        
    def ShowTree(self,name):
        trees = {'S': self.TreeS, 'Prob': self.TreeProb, 'value': self.TreeValue, 'PositionValue': self.TreePositionValue}
        print('Tree for ', name)
        print('\n'.join([''.join(['{0:.2f}   '.format(item) for j,item in enumerate(row) if i>=j ]) 
            for i,row in enumerate(trees[name])]))

    def ShowTreeValue(self):
        print('Tree for S')
        print('\n'.join([''.join(['{0:.2f}   '.format(item) for j,item in enumerate(row) if i>=j ]) 
            for i,row in enumerate(self.TreeValue)]))        
        

    def ShowTreeS(self):
        print('Tree for S')
        print('\n'.join([''.join(['{0:.2f}   '.format(item) for j,item in enumerate(row) if i>=j ]) 
            for i,row in enumerate(self.TreeS)]))
        
    def ShowTreeProb(self):
        print('Tree for Probs')
        print('\n'.join([''.join(['{0:.2f}   '.format(item) for j,item in enumerate(row) if i>=j ]) 
            for i,row in enumerate(self.TreeProb)]))    
        
       
         
# set payoff of derivative, i: maturity in periods    
    def payoffCall(self,K,i):
        self.TreeValue[i]=[max(ST-K,0) for ST in self.TreeS[i]]
        return self.TreeValue[i]


    def payoffPut(self,K,i):
        self.TreeValue[i]=[max(K-ST,0) for ST in self.TreeS[i]]
        return self.TreeValue[i]    

    
    
# list of all paths starting with note (0,0)   
    def allPaths(self):
        l = [[0],[1]]
        for i in range(2,self.n+1):
            g=list(l)
            l = l+g
#    print(l)
            for k in range(int(len(l)/2)):
                g=list(l[k])
                g.append(0)
#        print(g)
                l[k]=g
            for k in range(int(len(l)/2),int(len(l))):
                g=list(l[k])
                g.append(1)
#        print(g)
                l[k]=g  
        
        def f(x):
            res = [self.S]
            for i in range(len(x)):
                res.append(res[i]*self.u**x[i]*self.d**(1-x[i]))
            return res    
        
        res1 = [f(e) for e in l]
        jumps = [sum(e) for e in l]
        probs = [self.p**(e)*(1-self.p)**(self.n-e) for e in jumps]
        self.paths = res1
        self.jumps = l
        self.pathsProbs = probs
        

        

# list of paths starting with note (i,j)   
    def Paths(self,i,j):
        res = [k for k in range(len(self.paths)) if sum(self.jumps[k][:i])==j]
        return res
    
    def asian(self,K,i,j):
        self.TreeValue[i]=[max(K-ST,0) for ST in self.TreeS[i]]
        return self.TreeValue[i]    
    
    
    def price(self,K,i):
        payoff=self.payoffPut(K,i)
        probs=self.TreeProb[i]
        price=self.df**i*sum([x1*x2 for (x1,x2) in zip(payoff,probs)])
        return price

    
# American option by backward induction    
    def AmericanPut(self,K,i):
# check if final note, calculate final payoffs
        if i==self.n: 
            self.payoffPut(K,i) 
#otherwise
        else:
#calculate intrinsic values IV       
            IV=self.payoffPut(K,i)
# and value is we wait    
            W=[self.df*(self.p*self.TreeValue[i+1][j+1]+(1-self.p)*self.TreeValue[i+1][j]) for j in range(i+1)]
# choose max metween them
            self.TreeValue[i]=[max(x1,x2) for (x1,x2) in zip(IV,W)]
# if i=0, stop
        if i==0:
            return self.TreeValue[i]
# next step otherwise
        else:
            return self.AmericanPut(K,i-1) 
        
        
# Europian option by backward induction    
    def TreeCall(self,K,i):
# check if final note, calculate final payoffs

        if i==self.n: 
            self.payoffCall(K,i) 
        else:   
            W=[self.df*(self.p*self.TreeValue[i+1][j+1]+(1-self.p)*self.TreeValue[i+1][j]) for j in range(i+1)]
            self.TreeValue[i] = W
# if i=0, stop
        if i==0:
            return self.TreeValue[i]
# next step otherwise
        else:
            return self.TreeCall(K,i-1)        

# Europian option by backward induction    
    def TreeCallNew(self,K,i):
# check if final note, calculate final payoffs
        n = self.n
        self.TreeValue = [[0 for j in range(n+1) if j<=i] for i in range(n+1)]
        if i==self.n: 
            self.payoffCall(K,i) 
        else:   
            W=[self.df*(self.p*self.TreeValue[i+1][j+1]+(1-self.p)*self.TreeValue[i+1][j]) for j in range(i+1)]
            self.TreeValue[i] = W
# if i=0, stop
        if i==0:
            return self.TreeValue[i]
# next step otherwise
        else:
            return self.TreeCall(K,i-1)  

        
# Europian option by backward induction    
    def TreePut(self,K,i):
# check if final note, calculate final payoffs
        if i==self.n: 
            self.payoffPut(K,i) 
        else:   
            W=[self.df*(self.p*self.TreeValue[i+1][j+1]+(1-self.p)*self.TreeValue[i+1][j]) for j in range(i+1)]
            self.TreeValue[i] = W
# if i=0, stop
        if i==0:
            return self.TreeValue[i]
# next step otherwise
        else:
            return self.TreePut(K,i-1)        

        
def CallPrice(S,r,q,vol,K,T):
    """Calculate the price of Call Option. Inputs: Strike, Time to maturity."""
    if T.any()>0:
        d2=(np.log(S / K) + (r-q - vol**2 / 2) * T) / (vol * np.sqrt(T));
        d1=(np.log(S/K) + (r -q+ vol**2 / 2) * T)/(vol * np.sqrt(T))
        price = np.exp(-q * T) *S * norm.cdf(d1) - K * np.exp(-r * T) * norm.cdf(d2)
    else:     
        price = np.max(S-K,0)
        delta = 'n/a'
    return price

def CallDelta(S,r,q,vol,K,T):
    """Calculate the price of Call Option. Inputs: Strike, Time to maturity."""
    if T.any()>0:
        d2=(np.log(S / K) + (r-q - vol**2 / 2) * T) / (vol * np.sqrt(T));
        d1=(np.log(S/K) + (r -q+ vol**2 / 2) * T)/(vol * np.sqrt(T))
        price = np.exp(-q * T) *S * norm.cdf(d1) - K * np.exp(-r * T) * norm.cdf(d2)
        delta = np.exp(-q * T) * norm.cdf(d1)
    elif S.any()>K:     
        delta = 1
    else: 0    
    return delta


def Gamma(S,r,q,vol,K,T):
    d1=(np.log(S/K) + (r -q+ vol**2 / 2) * T)/(vol * np.sqrt(T))
    return norm.pdf(d1)/(S*vol*np.sqrt(T))



def ST(t,B,param):
    yt = (param['r']-param['sigma']**2/2)*t + B*param['sigma']
    res = param['S0']*np.exp(yt)
    return res
