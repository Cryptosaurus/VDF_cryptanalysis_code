import math,random
from math import log


#compute log2(2^x +2^y)
def log_add(x,y):
    if (y>x):
        (x,y)=(y,x)
    if(y+100<x):
        return x
    res = x + log(1+math.exp((y-x)*log(2.0)))/log(2.0)
    return res

def lgRho(x):
    if x<=1:
        return 0
    else:
        return -x*(log(x))/log(2.0)

def lgPiB(lgB):
    return  lgB-log(lgB*log(2.0))/log(2.0)

def lgProbAlmSmooth(lgB,lgBPrime,lgPrime):  #ProbAlmSmooth(B,B',x) = sum_(log2(B)<=i<=log2(B')) ProbSmooth(B,x/2^i)/ln(2)*((i-1)/(i*(i+1))))
    lgCumul=-10000 
    for i in range(lgB,lgBPrime):
        logval = lgRho((lgPrime-i)/lgB)+(log(i-1)-log(i)-log(i+1)-log(log(2.0)))/log(2.0)
        lgCumul = log_add(logval,lgCumul)
    return lgCumul


#logarithm of probability that a random number below 2^lgPrime has a 2^lgB0 smooth factor bigger than 2^lgT
def lgProbNonSmoothFactor(lgB0,lgT,lgPrime): #ProbPartSmooth(2^x,2^y,2^z)  = 1/y*\sum_{0\leq i<{x-z} }P_{sm,2^{x-i}}(2^y)
    lgCumul=-10000
    for i in range(0,lgPrime-lgT):
        if(i>=lgPrime-lgB0):
            logval = 0
        else:
            logval = lgRho((lgPrime-i)/lgB0)
        lgCumul = log_add(logval,lgCumul)
    lgCumul = lgCumul -log(lgT)/log(2.0)+0.23  #divide prob by y
    return lgCumul

def prob241(lgB,lgBPrime,lgR,lgPrime):
    lgAlmSmooth=lgProbAlmSmooth(lgB,lgBPrime,lgPrime)
    lgSuccess = lgAlmSmooth+lgR
    return lgSuccess

def findOpt241(lgPrime):
    mem = -10000
    mincost = 10000
    optLgB = 1
    optLgBPrime=2
    optLgR = 1
    for lgB in range(2,min(100,lgPrime-1)):
        for lgBPrime in range(lgB+1,lgB+20):
            for lgR in range(1,min(100,lgPrime-1)):
                val = prob241(lgB,lgBPrime,lgR,lgPrime)
                piB = lgPiB(lgB)
                cost = piB+lgR # group size * num of groups
                if (val>0 and cost<mincost):
                    optLgB=lgB
                    optLgBPrime=lgBPrime
                    optLgR = lgR
                    mem= lgPiB(lgBPrime)
                    mincost=cost
    return (mem,mincost,optLgB,optLgBPrime,optLgR)

def prob244(lgB,lgBPrime,lgB0,lgT,lgR,lgPrime):
    lgFirstFilter = lgPiB(lgB)-lgPiB(lgB0)+2*lgProbNonSmoothFactor(lgB0,lgT,lgPrime//2)
    lgSecondFilter= 2*lgProbAlmSmooth(lgB,lgBPrime,lgPrime/2-lgT)+lgR
    return min(lgFirstFilter,lgSecondFilter)

def findOpt244(lgPrime):
    mem = -10000
    mincost = 10000
    optLgB = 1
    optLgBPrime=2
    optLgR = 1
    optLgT = 1
    optLgB0 = 2
    for lgB in range(2,min(100,lgPrime-1)):
        for lgBPrime in range(lgB+1,lgB+20):
            for lgR in range(1,min(100,lgPrime-1)):
                for lgB0 in range(2,lgB,4):
                    for lgT in range(1,100,50):
                        val = prob244(lgB,lgBPrime,lgB0,lgT,lgR,lgPrime)
                        piB = lgPiB(lgB)
                        cost = piB+lgR # group size * num of groups
                        if (val>0 and cost<mincost):
                            optLgB=lgB
                            optLgBPrime=lgBPrime
                            optLgR = lgR
                            optLgT = lgT
                            optLgB0 = lgB0
                            mem= lgPiB(lgBPrime)
                            mincost=cost
    return (mem,mincost,optLgB,optLgBPrime,optLgB0,optLgT,optLgR)

def prob245(lgB,lgBPrime,lgR,lgPrime):
    lgAlmSmooth=lgProbAlmSmooth(lgB,lgBPrime,lgPrime+lgR)
    lgSuccess = lgAlmSmooth+lgR
    return lgSuccess


def findOpt245(lgPrime):
    minval = -10000
    mem = -10000
    mincost = 10000
    optLgB = 1
    optLgBPrime=2
    optLgR = 1
    for lgB in range(2,100):
        for lgBPrime in range(lgB+1,lgB+20):
            for lgR in range(1,100):
                val = prob245(lgB,lgBPrime,lgR,lgPrime)
                piB  = lgPiB(lgB)
                cost = max(piB,log(log(lgB*log(2.0))+1.2)/log(2.0)+lgR) #sum_{q<B} R/q
                if (val>0 and cost<mincost):
                    optLgB=lgB
                    optLgBPrime=lgBPrime
                    optLgR = lgR
                    minval=val                    
                    mem= lgPiB(lgBPrime)
                    mincost=cost
    return (mem,mincost,optLgB,optLgBPrime,optLgR)
 
def prob246(lgB,lgBPrime,lgR,lgPrime):
    lgAlmSmooth=lgProbAlmSmooth(lgB,lgBPrime,lgPrime/2+lgR)
    lgSuccess = 2*lgAlmSmooth+lgR
    return lgSuccess


def findOpt246(lgPrime):
    mem = -10000
    mincost = 10000
    optLgB = 1
    optLgBPrime=2
    optLgR = 1
    for lgB in range(2,100):
        for lgBPrime in range(lgB+1,lgB+20):
            for lgR in range(1,100):
                val = prob246(lgB,lgBPrime,lgR,lgPrime)
                piB  = lgPiB(lgB)
                cost = 1+max(piB,log(log(lgB*log(2.0))+1.2)/log(2.0)+lgR) #sum_{q<B} R/q
                if (val>0 and cost<mincost):
                    optLgB=lgB
                    optLgBPrime=lgBPrime
                    optLgR = lgR            
                    mem= lgPiB(lgBPrime)
                    mincost=cost
    return (mem,mincost,optLgB,optLgBPrime,optLgR)
 

def prob233(lgB,lgR,lgPrime):
    lgProbSmooth=lgRho(lgPrime/lgB)
    lgSuccess = lgProbSmooth+lgR
    return lgSuccess

def findOpt233(lgPrime):
    minval = -10000
    mincost = 10000
    optLgB = 1
    optLgR = 1
    for lgB in range(1,100):
        for lgR in range(1,100):
            val = prob233(lgB,lgR,lgPrime)
            piB = lgB-log(lgB*log(2.0))/log(2.0)
            cost = piB+lgR
            if (val>0 and cost<mincost):
                optLgB=lgB
                optLgR = lgR
                minval=val
                mincost=cost
    return (minval,mincost,optLgB,optLgR)

 

for prime in (96,112,128,192,256,384,512,768,1024):
    (minval,mincost,optLgB,optLgR) = findOpt233(prime);
    print("Alg 2.3.3 for prime {}: ",prime)
    print("logProb: ",minval);
    print("logCPU: ",mincost);
    print("logB: ",optLgB);
    print("logR: ",optLgR);
    (mem,mincost,optLgB,optLgBPrime,optLgR) = findOpt241(prime)
    print("Alg 2.4.1 for prime {}: ",prime)
    print("memory: ",mem);
    print("logCPU: ",mincost);
    print("logB: ",optLgB);
    print("logB': ",optLgBPrime);
    print("logR: ",optLgR);
    (mem,mincost,optLgB,optLgBPrime,optLgR) = findOpt245(prime);
    print("Alg 2.4.5 for prime: ", prime)
    print("memory: ",mem);
    print("logCPU: ",mincost);
    print("logB: ",optLgB);
    print("logB': ",optLgBPrime);
    print("logR: ",optLgR); 
    (mem,mincost,optLgB,optLgBPrime,optLgR) = findOpt246(prime);
    print("Alg 2.4.6 for prime: ", prime)
    print("memory: ",mem);
    print("logCPU: ",mincost);
    print("logB: ",optLgB);
    print("logB': ",optLgBPrime);
    print("logR: ",optLgR); 
 