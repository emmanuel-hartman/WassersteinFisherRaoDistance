import open3d as o3d
import WFR
import numpy as np

class PLShape:
    def __init__(self, filename=None, mesh=None, meas=None): 
        """
        Init        
        ----------
        PLShape Object
            returns a Shape Object itialized from a PLY file, a Open3D mesh, or a Measure Object
        """
        if mesh!=None:
            self.mesh=mesh
            self.meas=self.getMeasureFromShape()
        elif filename!=None:
            self.mesh=o3d.io.read_triangle_mesh(filename)
            self.meas=self.getMeasureFromShape()
        elif meas!=None:            
            self.meas=meas
            self.mesh=self.getShapeFromMeasure()
            
    def getMeasureFromShape(self):
        """
        Returns        
        ----------
        measure: Measure Object
            returns a Measure Object that corresponds to the PL shape opject

        """
        measure = WFR.Measure(self.getNormalsAsArray(),self.getAreasAsArray())
        return measure
    
    def getShapeFromMeasure(self):
        """
        Returns        
        ----------
        shape: PLShape Object
            returns a Shape Object with the mesh of the unique convex PL surface that corresponds to the prescribed measure

        """
        u = [self.meas.supports[i,:] for i in range(0, self.supports.shape[0])] 
        a = np.squeeze(self.meas.masses).tolist()        
        P = polyhedrec.reconstruct(u,a)
        numV = len(P.vertices)
        numF = len(P.faces)

        file = open("temp.ply", "w")
        lines = ("ply","\n","format ascii 1.0","\n", "element vertex {}".format(numV),"\n","property float x","\n","property float y","\n","property float z","\n","element face {}".format(numF),"\n","property list uchar int vertex_index","\n","end_header","\n")

        file.writelines(lines)
        lines=[]
        for i in range(0,numV):
            for j in range(0,3):
                lines.append(str(P.vertices[i][j]))
                lines.append(" ")
            lines.append("\n")
        for i in range(0,numF):
            l=len(P.faces[i].vertices)
            lines.append(str(l))
            lines.append(" ")

            for j in range(0,l):
                lines.append(str(P.faces[i].vertices[j]))
                lines.append(" ")
            lines.append("\n")

        file.writelines(lines)
        file.close()
        cat1=o3d.io.read_triangle_mesh("temp.ply")
        os.remove("temp.ply")
        return cat1
    
    def exportToFile(self,filename):
        """
        Parameters
        ----------
        filename : String
             path to save location for self.mesh as a PLY file
        
        Save       
        ----------
        filename.ply : PLY File
            save self.mesh as a PLY file

        """
        
        o3d.io.write_triangle_mesh(filename,self.mesh)
        
    def visualizeShape(self):
        """
        Visualizes the PL shape
        """
        o3d.visualization.draw_geometries([self.mesh])
    
    def normalColorize(self, visualize):
        """
        Colors the PL shape so that each face is colored according to its unit normal vector.
        
        Parameters
        ----------
        visualize : boolean
            whether or not to visualize the shape ofter recoloring the shape
            
        """
        x = np.asarray(self.mesh.vertices)
        G = np.asarray(self.mesh.triangles)
        N = self.getNormalsAsArray()
        m=len(N)

        newmesh=self.mesh
        newmesh.clear()

        lTri=[]
        for i in range(0, m):

            Gi=np.array([[0,1,2]])
            xi=np.zeros((3,3))


            xi[0,:]=x[G[i,0]]
            xi[1,:]=x[G[i,1]]
            xi[2,:]=x[G[i,2]]



            xi= o3d.open3d_pybind.utility.Vector3dVector(xi)
            Gi= o3d.open3d_pybind.utility.Vector3iVector(Gi)
            iTri= o3d.geometry.TriangleMesh(xi,Gi)
            iTri.paint_uniform_color([(N[i,0]+1)/2,(N[i,1]+1)/2,(N[i,2]+1)/2])
            newmesh=newmesh+iTri

        self.mesh=newmesh
        if visualize:
            self.visualizeShape()
        
    def shadeColorize(self,Base,visualize):
        """
        Colors the PL shape so that each face is shaded according to its unit normal vector
        
        Parameters
        ----------
        Base : np.array 
            The base color of the mesh
        ----------
        visualize : boolean
            whether or not to visualize the shape ofter recoloring the shape
        """
        x = np.asarray(self.mesh.vertices)
        G = np.asarray(self.mesh.triangles)
        N = self.getNormalsAsArray()
        m=len(N)
        S=np.array([np.sqrt(3),np.sqrt(3),np.sqrt(3)])

        newmesh=self.mesh
        newmesh.clear()

        lTri=[]
        for i in range(0, m):
            Gi=np.array([[0,1,2]])
            xi=np.zeros((3,3))


            xi[0,:]=x[G[i,0]]
            xi[1,:]=x[G[i,1]]
            xi[2,:]=x[G[i,2]]
            
            Int=N[i]@S/3

            xi= o3d.open3d_pybind.utility.Vector3dVector(xi)
            Gi= o3d.open3d_pybind.utility.Vector3iVector(Gi)
            iTri= o3d.geometry.TriangleMesh(xi,Gi)
            iTri.paint_uniform_color(Base*(.5+.5*Int))
            newmesh=newmesh+iTri

        self.mesh=newmesh
        if visualize:
            self.visualizeShape()
    
    def colorizeFromMap(self,P,N,visualize):
        """
        Colors the PL shape so that each face is shaded according semi-coupling P,
        
        Parameters
        ----------
        P : tensor 
            The semi-coupling that defines the map from the PL shape to some other PL shape
        ----------
        N : np.array 
            The unit normal vectors of the other PL shape
        ----------
        visualize : boolean
            whether or not to visualize the shape ofter recoloring the shape
        """

        x1 = np.asarray(self.mesh.vertices)
        G1 = np.asarray(self.mesh.triangles)
        
        cats=self.mesh
        cats.clear()
        
        n = P.shape[1]

        for j in range(0,n):

            Gj=np.array([[0,1,2]])
            xj=np.zeros((3,3))


            xj[0,:]=x1[G1[j,0]]
            xj[1,:]=x1[G1[j,1]]
            xj[2,:]=x1[G1[j,2]]

            xj= o3d.open3d_pybind.utility.Vector3dVector(xj)
            Gj= o3d.open3d_pybind.utility.Vector3iVector(Gj)
            jTri= o3d.geometry.TriangleMesh(xj,Gj)

            i=np.argmax(P[:,j].cpu().numpy())


            jTri.paint_uniform_color([(N[i,0]+1)/2,(N[i,1]+1)/2,(N[i,2]+1)/2])
            cats=cats+jTri
        self.mesh=cats
        if visualize:
            self.visualizeShape()
        

    def getAreasAsArray(self):
        """
        Returns        
        ----------
        areas: np.array
            np array of the areas of the faces of the PL surface

        """
        x = np.asarray(self.mesh.vertices)
        G = np.asarray(self.mesh.triangles)
        n= G.shape[0]
        areas= np.zeros((n,1))

        for i in range(0,n):
            P=x[G[i,0],:]
            Q=x[G[i,1],:]
            R=x[G[i,2],:]

            PQ=Q-P
            PR=R-P
            areas[i,:]=.5*np.linalg.norm(np.cross(PQ,PR))


        return areas

    def getNormalsAsArray(self):
        """
        Returns        
        ----------
        normals: np.array
            np array of the unit normals of the faces of the PL surface

        """
        self.mesh.compute_triangle_normals(normalized=True)
        normals = np.asarray(self.mesh.triangle_normals)
        return normals
    
    
    def getRotatedShape(self,R):
        """
        Parameters
        ----------
        R : tensor 
            tensor representing a rotation
        Returns        
        ----------
        shape: PLShape
            PLShape of the rotated object

        """
        trimesh=self.mesh.rotate(R.cpu().numpy(),  center=self.mesh.get_center())
        return PLShape(mesh=trimesh)    
    
    def downsampleShape(self,n):
        """
        Parameters
        ----------
        n : int 
            The number of desired faces for the PL shape
            
            
        Returns        
        ----------
        shape: PLShape
            PLShape of the downsampled mesh

        """
        trimesh=self.mesh.simplify_quadric_decimation(target_number_of_triangles=n)
        return PLShape(mesh=trimesh)
    
    def shapeDistance(self,shape2,NoIterations,eps):
        """
        Parameters
        ----------
        shape2 : PLShape Object
            target shape
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
        dist,cost,ind,P,Q=WFR.measureDistance(self.meas,shape2.meas,NoIterations,eps)    
        return dist,cost,ind,P,Q

    def shapeDistanceModSO3(self,shape2,RotIterations,RotDepth,NoIterations,eps):
        """
        Parameters
        ----------
        shape2 : PLShape Object
            target shape
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
        dist: Float
            distance of between  self and shape2 mod SO3     
        ----------
        P,Q: two tensors
            pair of tensors representing the  square root of optimal semi-coupling
        ----------
        R : tensor
            tensor representing the optimal rotation
        """
        dist,P,Q,R=WFR.measureDistanceModSO3(self.meas,shape2.meas,RotIterations,RotDepth,NoIterations,eps)
        return dist,P,Q,R


        