# methods required to perform the transverse coordinate transformations

from copy import copy as cp
from numpy import sqrt as Sqrt
import numpy as np

#represents iterative x transformation
##starts at n-E
##iterates (c+d) times, this is x order
##returns a queue with X_{+/-}^1(n) for each n
#n start is n-e

class Item:
    def __init__(self, n,N,coeff):
        self.n = n
        self.N = N
        self.coeff = coeff

# represented by p_order
def transform_x (E,x_order, a,w,gouy,n, w0,z,zr):
    q1 = []
    q2 = []
    #queue of items which are passed each iteration
    coeff = 1 #coeff can be worked in this fxn, then multiplied to overall coupling after
    
    N = n-E #start at n-E
    
    start_item = Item(n,N,coeff)
    q1.append(start_item) #just for printing, need 3-array. to track inc/dec
    
    x_count = x_order  #counter to track iterations
    #print("n,E,N",n,E,N)
    
    #empty q1 iteratively and raise and lower
    while(x_count!=0):
        while(q1):
            #get from q1
            item = cp(q1.pop())
            
            #raise n, factor X_+^1(n)  **raising doesn't involve factor n**
            item_raise = cp(item)
            
            item_raise.coeff = item_raise.coeff*x_plus_1(w0,z,zr,item_raise.N)
            item_raise.N += 1 #overall n mode up
            
            
            q2.append(item_raise) # put into temp q
            
            
            #lower n, factor X_-^1(n)
            item_lower = cp(item)
            
            if(item_lower.N>0):
                item_lower.coeff = item_lower.coeff*x_minus_1(w0,z,zr,item_lower.N)
                item_lower.N -= 1 #overall n mode down

                q2.append(item_lower) #put into temp q
            
        #empty q2 back into q1 to re-iterate
        while(q2):
            q1.append(q2.pop())
            
        x_count-=1 #iteration done, decrement count
        
    
    #all x factors transformed, return full queue
    return(q1)

#x dep
def x_plus_1(w0,z,zr,n):    
    factor = (w0/2)*( ( 1-(1j)*(z/zr) )*np.sqrt(n+1))
    return(factor)

#x dep
def x_minus_1(w0,z,zr,n):
    factor = (w0/2)*( np.sqrt(n)*(1+(1j)*(z/zr)) )
    return(factor)

class Symbolic_Item:
    def __init__(self, n,N, p_count,p_coeff,wover2_count,wover2_coeff,n_coeff,overall):
        self.n = n #start n
        self.N = N #final n
        self.p_count=p_count # exponenent p is raised to. X+ -> p^-1, X- -> p^+1
        self.p_coeff=p_coeff # p term string ("p**(p_count)")
        self.wover2_count = wover2_count # exp. w/2 raised to. X+ or X- -> (w/2)**(+1)
        self.wover2_coeff = wover2_coeff # w/2 term. (w/2)**(wover2_count)
        self.n_coeff = n_coeff # n-dep. term. X+ -> Sqrt(n+1)* X- -> Sqrt(n)*   ... builds up
        self.overall = overall # overall term which builds up on n_coeff. (w/2)**(wover2_count)*p**(p_count)*n_coeff
        
        

def symbolic_transform_x (p_order,x_order,n):
    #queues of items which are passed each iteration
    q1 = []
    q2 = []
    
    
    p_count=0
    p_coeff = ""
    wover2_count=0
    wover2_coeff = ""
    n_coeff = "" #coeff can be worked in this fxn, then multiplied to overall coupling after
    overall=""
    
    N = n-p_order #start at n-p_order
    
    start_item = Symbolic_Item(n,N, p_count,p_coeff,wover2_count,wover2_coeff,n_coeff,overall) #WHAT IS COEFF
    q1.append(start_item) #just for printing, need 3-array. to track inc/dec
    
    x_count = x_order  #counter to track iterations
    
    #empty q1 iteratively and raise and lower
    while(x_count!=0):
        while(q1):
            #get from q1
            item = cp(q1.pop())
            
            #Perform X+ on this item and stack      
            item_raise = cp(item) 
            item_raise = symbolic_raise(item_raise)
            q2.append(item_raise) # put into temp q
            
            #lower n, factor X_-^1(n)
            item_lower = cp(item) 
            if(item_lower.N>0):
                item_lower = symbolic_lower(item_lower)
                q2.append(item_lower) #put into temp q
            
        #empty q2 back into q1 to re-iterate
        while(q2):
            q1.append(q2.pop())
            
        x_count-=1 #iteration done, decrement count

    #all x factors transformed, return full queue
    return(q1)

def symbolic_raise(item):
    item.p_count -= 1
    item.p_coeff = "p**"+str(item.p_count)
    item.wover2_count += 1
    item.wover2_coeff = "(w/2)**"+str(item.wover2_count)
    item.n_coeff += "*"+str(Sqrt(item.N+1))
    item.overall = item.p_coeff+"*"+item.wover2_coeff+item.n_coeff
    item.N += 1 #overall n mode up
    return(item)
   
def symbolic_lower(item):
    item.p_count += 1
    item.p_coeff = "p**"+str(item.p_count)
    item.wover2_count += 1
    item.wover2_coeff = "(w/2)**"+str(item.wover2_count)
    item.n_coeff += "*"+str(Sqrt(item.N))
    item.overall = item.p_coeff+"*"+item.wover2_coeff+item.n_coeff
    item.N -= 1 #overall n mode up
    return(item)

def build_symbolic_dictionary(n_max,p_max,x_max):
    d = {}
    for n_ind in range(0,n_max+1):
        for p_ind in range(0,p_max+1):
            for x_ind in range(0,x_max+1):
                q = symbolic_transform_x(p_ind,x_ind,n_ind)
                key = 'n'+str(n_ind)+'p'+str(p_ind)+'x'+str(x_ind)

                if (p_ind<=n_ind):
                    while(q):              
                        item = cp(q.pop())
                        if (x_ind==0):
                            item.overall=1
                        d.setdefault(key, []).append(item)
                    
    return(d)
