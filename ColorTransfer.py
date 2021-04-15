import open3d as o3d
import numpy as np
import matplotlib.image as mpimg 
import matplotlib.pyplot as pyplot
import math
from skimage.segmentation import slic
import WFR


class Image:
    def __init__(self, filename=None, img=None,n_segments=1024):
        """
        Init        
        ----------
        PLShape Object
            returns a Shape Object itialized from a PNG file or a numpy array
        ----------
        Parameters
        ----------
        n_segments : int
            maximum number of supports for the associated measure
            
        """
        if filename==None:
            self.img=img
        else:
            rbga=mpimg.imread(filename)
            self.img=np.zeros((rbga.shape[0],rbga.shape[1],3))
            for i in range(rbga.shape[0]):
                for j in range(rbga.shape[1]):
                    if rbga.shape[2]>3:
                        self.img[i,j] = rbga[i,j,0:3]*rbga[i,j,3]
                    else:
                        self.img[i,j] = rbga[i,j,0:3]
        self.segments=slic(self.img, n_segments=n_segments, sigma=5)    
        self.meas=self.getMeasureFromImage()
        
    def getMeasureFromImage(self):
        """
        Returns
        ----------
        measure : Measure Object
            measure associated with the Image
            
        """
        segment_vectors = np.zeros((np.unique(self.segments).size,3))
        for i in range(self.img.shape[0]):
            for j in range(self.img.shape[1]):
                segmentid=self.segments[i,j]
                segment_vectors[segmentid,:]=segment_vectors[segmentid,:]+(self.img[i,j,:]/np.sqrt(np.linalg.norm(self.img[i,j,:])))
                
        supports = []
        norms = []        
        for i in range(len(segment_vectors)):
            norm=np.linalg.norm(segment_vectors[i])
            if norm>0:
                norms.append(norm)
                supports.append(segment_vectors[i]/norm)
            else:
                norms.append(.0000000001)
                supports.append(np.repeat(.0001,3))
            
        supports=np.array(supports).astype(np.float64)
        norms=np.transpose(np.array(norms,ndmin=2).astype(np.float64))
        
        return WFR.Measure(supports,norms)    
    
    def colorizeFromMap(self,P,img2):
        """
        Parameters
        ----------
        P : tensor
            tesor representing the square root of the optimal semi-coupling 
        ----------
            
        Returns
        ----------
        image : Image Object
          Image Object representing self recolored by the map from self to img2
            
        """
        norms2d=np.linalg.norm(self.img, axis=2)   
        nimg=np.zeros((self.img.shape[0],self.img.shape[1],3))
        
        newSupports=np.zeros((P.shape[0],3))
        for i in range(0,P.shape[0]):
            support=img2.meas.supports
            masses=P[i,:].cpu().numpy()
            masses*=masses
            newSupports[i,:]=np.sum(support*masses[:,None],axis=0)
            if np.linalg.norm(newSupports[i,:])>0:
                newSupports[i,:]/=np.linalg.norm(newSupports[i,:])        
        
        for i in range(self.img.shape[0]):
            for j in range(self.img.shape[1]):
                segmentid=self.segments[i,j]
                nimg[i,j,:]=norms2d[i,j]*norms2d[i,j]*newSupports[segmentid]
                    
        return Image(img=nimg)
                
    def exportToFile(self,filename):
        """
        Parameters
        ----------
        filename : string
            path to save the Image Object as a PNG file
        ----------
        """
        pyplot.imshow(self.img)
        pyplot.savefig(filename)
        
    def visualizeImage(self):
        """
        Visualizes the Image Object
        """     
        pyplot.imshow(self.img)
        pyplot.show()
        
    def colorDistance(self,image2,NoIterations,eps):
        """
        Parameters
        ----------
        image2 : Image Object
            target image
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
        dist,cost,ind,P,Q=WFR.measureDistance(self.meas,image2.meas,NoIterations,eps)    
        return dist,cost,ind,P,Q
    
    def colorDistanceModSO3(self,image2,RotIterations,RotDepth,NoIterations,eps):
        """
        Parameters
        ----------
        shape2 : Image Object
            target image
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
        ----------
            
        Returns        
        ----------
        dist: Float
            distance of between  self and shape2 mod SO3     
        ----------
        P,Q: two tensors
            pair of tensors representing the  square root of optimal semi-coupling
        ----------
        R : tensor
            tensor representing the optimal rotation
        """
        dist,P,Q,R=WFR.measureDistanceModSO3(self.meas,image2.meas,RotIterations,RotDepth,NoIterations,eps)
        return dist,P,Q,R
    