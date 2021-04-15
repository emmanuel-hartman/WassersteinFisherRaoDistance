import PLShapes
import numpy as np
import numpy.matlib
import torch
import open3d as o3d
import numpy as np


device = torch.device("cuda:0")

class Measure:
    """
    Init        
    ----------
    Measure Object
        returns a Measure Object itialized from a np array representing the supports and the corresponding masses
    """
    def __init__(self, supports, masses):
        self.supports = supports
        self.masses = masses 


###########################################################################
###################--Distance Computation--################################
###########################################################################

    
def measureDistance(measure1,measure2,NoIterations,eps):
    """
    Parameters
    ----------
    measure1 : Measure Object
        source measure
    ----------
    measure2 : Measure Object
        target measure
    ----------
    NoIterations : int
        maximum number of iterations for the WFR algorithm
    ----------
    eps : float
        minimum increase in the cost function F before the algorithm terminates

    Returns        
    ----------
    dist: Float
        distance of between self and shape2      
    ----------
    cost: np.array
        evolution of the cost function F        
    ----------
    ind: int
        number of iterations until convergence        
    ----------
    P,Q: two tensors
        pair of tensors representing the  square root of optimal semi-coupling
    """
    a = torch.from_numpy(measure1.masses).to(device)
    b = torch.from_numpy(measure2.masses).to(device)    
    m = a.size()[0]
    n = b.size()[0]
    u = torch.from_numpy(measure1.supports).to(device)
    v = torch.from_numpy(measure2.supports).to(device)
    Omega = calcOmega(u,v).to(device)
    
    P = Omega
    P[P<0]=0   
    Q = P
    
    cost=np.zeros((NoIterations+1,1))
    for k in range(0,NoIterations):
        P,Q=contractionRowCol(P,Q,Omega,a,b,m,n)
        cost[k+1,:]=calcF(P,Q,Omega).cpu()
        ind=k+1
        if (cost[k+1]-cost[k,:])/cost[k+1]<eps:
            break   
    dist=torch.sqrt(sum(a.cpu())+sum(b.cpu())-2*calcF(P,Q,Omega).cpu())
    return dist,cost,ind,P,Q 
    
def measureDistanceModSO3(measure1,measure2,RotIterations,RotDepth,NoIterations,eps):
    """
    Parameters
    ----------
    measure1 : Measure Object
        source measure
    ----------
    measure2 : Measure Object
        target measure
    ----------
    RotIterations : int
        maximum number of iterations for determining the optimal rotation
    ----------
    RotDepth : int
        maximum number of iterations for each optimal rotation step
    ----------
    NoIterations : int
        maximum number of iterations for the WFR algorithm
    ----------
    eps : float
        minimum increase in the cost function F before the algorithm terminates

    Returns        
    ----------
    meas: Measure Object
        measure at the midpoint between measure1 and measure2 with respect to the WFR metric
    """
    a = torch.from_numpy(measure1.masses).to(device)
    b = torch.from_numpy(measure2.masses).to(device)    
    m = a.size()[0]
    n = b.size()[0]
    u = torch.from_numpy(measure1.supports).to(device)
    v = torch.from_numpy(measure2.supports).to(device)
    R = torch.eye(3).type('torch.DoubleTensor').to(device)
    for i in range(0,RotIterations):    
        Omega=calcOmegaWithRotation(u,v,R).to(device)
        P,Q = getAssignment(a,b,m,n,Omega,RotDepth,eps)
        Rn = getRotation(P,Q,u,v,m,n).to(device)
        Rd = (R-Rn)*(R-Rn)
        s = torch.sum(Rd)
        if s<eps:
            R=Rn
            break
        else:
            R=Rn  
    Omega = calcOmegaWithRotation(u,v,R).to(device)
    P,Q = getAssignment(a,b,m,n,Omega,NoIterations,eps)
    dist = torch.sqrt(sum(a.cpu())+sum(b.cpu())-2*calcF(P,Q,Omega).cpu())
    return dist,P,Q,R


def midpointMeasure(measure1,measure2,NoIterations,eps):
    """
    Parameters
    ----------
    measure1 : Measure Object
        source measure
    ----------
    measure2 : Measure Object
        target measure
    ----------
    NoIterations : int
        maximum number of iterations for the WFR algorithm
    ----------
    eps : float
        minimum increase in the cost function F before the algorithm terminates

    Returns        
    ----------
    dist: Float
        distance of between measure1 and measure2 mod SO3      
    ----------
    P,Q: two tensors
        pair of tensors representing the square root of optimal semi-coupling
    ----------
    R : tensor
        tensor representing the optimal rotation
    """
    dist,cost,ind,P,Q = measureDistance(measure1,measure2,NoIterations,eps)
    t=.5
    
    P = P.cpu().numpy()
    Q = Q.cpu().numpy()
    
    m = measure1.masses.shape[0]
    n = measure2.masses.shape[0]
    
    u = measure1.supports
    v = measure2.supports 
    
    norms=[]
    areas=[]
    for i in range(0,m):
        for j in range(0,n):
            inv = t*P[i,j]*u[i]+(1-t)*Q[i,j]*v[j]
            np.expand_dims(inv,1)
            ninv =  np.sqrt(np.sum(np.power(inv,2)))
            vec=ninv*inv
            nvec=np.sqrt(np.sum(np.power(vec,2)))
            if nvec>0:
                norms.append(vec/nvec)
                areas.append(nvec)
    return Measures.Measure(np.array(norms), np.expand_dims(np.array(areas),1))
    

###########################################################################
###################--Helper Functions--####################################
###########################################################################

def calcF(P,Q,Omega):
    """
    Parameters
    ----------
    P,Q: two tensors
    pair of tensors representing the square root of current semi-coupling
    ----------
    Omega: tensor 
    tensor representing the cost matrix
    ----------        

    Returns
    ----------
    cost: float
    value of the cost function F for the semi-coupling (P,Q)
    """
    cost=torch.sum(P*Q*Omega)
    return cost

def contractionRowCol(P,Q,Omega,a,b,m,n):
    """
    Parameters
    ----------
    P,Q: two tensors
    pair of tensors representing the square root of current semi-coupling
    ----------
    Omega: tensor 
    tensor representing the cost matrix
    ----------   
    a: tensor 
    tensor representing the masses of measure1
    ----------   
    b: tensor 
    tensor representing the masses of measure2
    ----------      
    m: int
    number of supports of measure1
    ----------   
    n: int 
    number of supports of measure2
    ----------        

    Returns
    ----------
    P,Q: two tensors
    pair of tensors representing the updated square root of current semi-coupling
    """    
    P = rowNormalize(Q*Omega,a,n)
    Q = colNormalize(P*Omega,b,m)
    return P,Q  

def rowNormalize(Pnew,a,n):
    """
    Parameters
    ----------
    Pnew: two tensors
    pair of tensors representing the square root of current semi-coupling
    ----------
    Omega: tensor 
    tensor representing the cost matrix
    ----------   
    a: tensor 
    tensor representing the masses of measure1
    ----------  
    n: int 
    number of supports of measure2
    ----------        

    Returns
    ----------
    PnewNormalized: tensor
    tensor representing the updated square root of current semi-coupling 
    """  
    RowNormPnew = torch.sqrt(torch.sum(Pnew*Pnew,dim=1)/a.transpose(0,1))
    RowNormMatrix = RowNormPnew.repeat([n,1]).transpose(0,1)
    PnewNormalized = Pnew/RowNormMatrix
    return PnewNormalized
    
def colNormalize(Qnew,b,m):
    """
    Parameters
    ----------
    Qnew: two tensors
    tensor representing the square root of current semi-coupling
    ----------   
    b: tensor 
    tensor representing the masses of measure2
    ---------- 
    m: int 
    number of supports of measure1
    ----------        

    Returns
    ----------
    QnewNormalized: tensor
    tensor representing the updated square root of current semi-coupling 
    """  
    ColumnNormQnew = torch.sqrt(torch.sum(Qnew*Qnew,dim=0)/b.transpose(0,1))
    ColumnNormMatrix = ColumnNormQnew.repeat([m,1])
    QnewNormalized = Qnew/ColumnNormMatrix
    return QnewNormalized
    

def getAssignment(a,b,m,n,Omega,NoIterations,eps):
    """
    Parameters
    ----------
    a: tensor 
        tensor representing the masses of measure1
    ----------   
    b: tensor 
        tensor representing the masses of measure2
    ----------      
    m: int
        number of supports of measure1
    ----------   
    n: int 
        number of supports of measure2
    ----------
    Omega: tensor 
        tensor representing the cost matrix
    ----------
    NoIterations : int
        maximum number of iterations for the WFR algorithm
    ----------
    eps : float
        minimum increase in the cost function F before the algorithm terminates           

    Returns
    ----------
    P,Q: two tensors
        pair of tensors representing the updated square root of optimal semi-coupling
     """    
    P = Omega
    P[P<0]=0
    Q = P
    for j in range(0,NoIterations):
        P,Q=contractionRowCol(P,Q,Omega,a,b,m,n)
    return P,Q   

def getRotation(P,Q,u,v,m,n):
    """
    Parameters
    ----------
    P,Q: two tensors
        pair of tensors representing the updated square root of current semi-coupling
    ----------
    u: tensor 
        tensor representing the supports of measure1
    ----------   
    v: tensor 
        tensor representing the supports of measure2
    ----------      
    m: int
        number of supports of measure1
    ----------   
    n: int 
        number of supports of measure2     

    Returns
    ----------
    R : tensor
        tensor representing the optimal rotation
     """ 
    Y=torch.einsum('ij,ij,ik,jl->kl',P,Q,u,v)         
    (U,S,V)=torch.svd(Y)
    if S[2]>=0:
        R=torch.matmul(U,torch.transpose(V,0,1))
    else:
        I=torch.eye(3).type('torch.DoubleTensor').to(device)
        I[2,2]=-1;
        R=torch.matmul(torch.matmul(U,I),torch.transpose(V,0,1))
        
    return R

def calcOmega(u,v):
    """
    Parameters
    ----------
    u: tensor 
        tensor representing the supports of measure1
    ----------   
    v: tensor 
        tensor representing the supports of measure2
    ----------

    Returns
    ----------
    Omega : tensor
        tensor representing the cost matrix
     """ 
    R = torch.eye(3).type('torch.DoubleTensor').to(device)
    return calcOmegaWithRotation(u.to(device),v.to(device),R) 

def calcOmegaWithRotation(u,v,R):
    """
    Parameters
    ----------
    u: tensor 
        tensor representing the supports of measure1
    ----------   
    v: tensor 
        tensor representing the supports of measure2
    ----------
    R : tensor
        tensor representing a rotation of measure1

    Returns
    ----------
    Omega : tensor
        tensor representing the cost matrix
     """ 
    return torch.einsum('kl,ik,jl->ij',R,u,v)    
    
