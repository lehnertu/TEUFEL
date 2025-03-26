import numpy as np

class LorentzTransform:
    """
    Transform 4-vectors from one inertial system to another.
    Internally represented by a 4x4 numpy array.
    """

    def __init__(self,gamma,beta):
        """
        Compute the transform matrix for transforming 4-vectors
        from the reference (lab) system into a system moving with velocity
        beta=v/c and relativistic factor gamma in the direction given by beta.
        No checks are performed whether gamma and beta match
        """
        g_fac = gamma**2/(1+gamma)
        self.matrix = np.eye(4)
        self.matrix[0,0] = gamma
        self.matrix[0,1:4] = -gamma*beta
        self.matrix[1:4,0] = -gamma*beta
        self.matrix[1,1:4] += g_fac*beta[0]*beta
        self.matrix[2,1:4] += g_fac*beta[1]*beta
        self.matrix[3,1:4] += g_fac*beta[2]*beta
        
    @classmethod
    def from_betagamma(cls,BG):
        """
        Construct the transfor from a 3-vector beta*gamma.
        This is the preferred method as it is least prone to numeric inaccuracies.
        """
        gamma = np.sqrt(np.dot(BG,BG)+1.0)
        beta = BG/gamma
        return cls(gamma,beta)
        
    def __call__(self,X):
        """
        Apply the transformation to a 4-vector.
        Internally executed as a mmatrix-vector multiplication.
        """
        return np.dot(self.matrix,X)