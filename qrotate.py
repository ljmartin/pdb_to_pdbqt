import numpy as np


class Align(object):
    def __init__(self):
        pass

    def fast_normalise(self, q):
        """
        https://github.com/KieranWynn/pyquaternion/blob/99025c17bab1c55265d61add13375433b35251af/pyquaternion/quaternion.py#L513
        Normalise the object to a unit quaternion using a fast approximation method if appropriate.
        Object is guaranteed to be a quaternion of approximately unit length
        after calling this operation UNLESS the object is equivalent to Quaternion(0)
        """
    
        mag_squared = np.dot(q, q)
        if (mag_squared == 0):
            return
        if (abs(1.0 - mag_squared) < 2.107342e-08):
            mag =  ((1.0 + mag_squared) / 2.0) # More efficient. Pade approximation valid if error is small
        else:
            mag =  np.sqrt(mag_squared) # Error is too big, take the performance hit to calculate the square root properly
        
        return q / mag


    def quaternion_rotation_matrix(self, Q):
        """
        https://automaticaddison.com/how-to-convert-a-quaternion-to-a-rotation-matrix/
        Covert a quaternion into a full three-dimensional rotation matrix.
        
        Input
        :param Q: A 4 element array representing the quaternion (q0,q1,q2,q3) 
        
        Output
        :return: A 3x3 element matrix representing the full 3D rotation matrix. 
             This rotation matrix converts a point in the local reference 
             frame to a point in the global reference frame.
        """
        # Extract the values from Q
        q0 = Q[0]
        q1 = Q[1]
        q2 = Q[2]
        q3 = Q[3]
        
        # First row of the rotation matrix
        r00 = 2 * (q0 * q0 + q1 * q1) - 1
        r01 = 2 * (q1 * q2 - q0 * q3)
        r02 = 2 * (q1 * q3 + q0 * q2)
        
        # Second row of the rotation matrix
        r10 = 2 * (q1 * q2 + q0 * q3)
        r11 = 2 * (q0 * q0 + q2 * q2) - 1
        r12 = 2 * (q2 * q3 - q0 * q1)
        
        # Third row of the rotation matrix
        r20 = 2 * (q1 * q3 - q0 * q2)
        r21 = 2 * (q2 * q3 + q0 * q1)
        r22 = 2 * (q0 * q0 + q3 * q3) - 1
        
        # 3x3 rotation matrix
        rot_matrix = np.array([[r00, r01, r02],
                               [r10, r11, r12],
                               [r20, r21, r22]]) 
        return rot_matrix


    def from_axis_angle(self, axis, angle):
        """Initialise from axis and angle representation
        Create a Quaternion by specifying the 3-vector rotation axis and rotation
        angle (in radians) from which the quaternion's rotation should be created.
        Params:
        axis: a valid numpy 3-vector
        angle: a real valued angle in radians
        """
        mag_sq = np.dot(axis, axis)
        if mag_sq == 0.0:
            raise ZeroDivisionError("Provided rotation axis has no length")
        # Ensure axis is in unit vector form
        if (abs(1.0 - mag_sq) > 1e-12):
            axis = axis / np.sqrt(mag_sq)
        theta = angle / 2.0
        r = np.cos(theta)
        i = axis * np.sin(theta)
        
        return np.array([r, i[0], i[1], i[2]])
    
    def compute_com(self, pts, masses=None):
        """Computes the centre of mass of a point cloud, assuming mass=1"""
        
        if masses is None:
            masses = np.ones(pts.shape[0])
            masses /= masses.sum()
            
        return pts.T.dot(masses).astype('float64')

    def calc_m_i(self, pcl):
        """
        Computes the moment of inertia tensor.
        a more convoluted but easier to understand alternative is in here:
        https://github.com/jwallen/ChemPy/blob/master/chempy/geometry.py
        """
        A = np.sum((pcl**2) * np.ones(pcl.shape[0])[:,None],0).sum()
        B = (np.ones(pcl.shape[0]) * pcl.T).dot(pcl)
        eye = np.eye(3)
        return A * eye - B


    def get_pmi(self, coords):
        """
        Calculate principal moment of inertia
        This code is a re-write of the function from MDAnalysis:
        https://github.com/MDAnalysis/mdanalysis/blob/34b327633bb2aa7ce07bbd507a336d1708871b6f/package/MDAnalysis/core/topologyattrs.py#L1014
        """
        momint = self.calc_m_i(coords)
        eigvals, eigvecs = np.linalg.eig(momint)
        
        # Sort
        indices = np.argsort(-eigvals) #sorting it returns the 'long' axis in index 0.
        # Return transposed which is more intuitive format
        return eigvecs[:, indices].T

    def unit_vector(self, vector):
        """ Returns the unit vector of the vector.  """
        return vector / np.linalg.norm(vector)
    
    def angle_between(self, v1, v2):
        """ 
        Returns the angle in radians between vectors 'v1' and 'v2'::
        
            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
        """
        v1_u = self.unit_vector(v1)
        v2_u = self.unit_vector(v2)
        return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


    def rotaxis(self, a, b):
        """Calculate the vector required in order to rotate `a` to align with `b`
        This is just the normal vector between the two vectors, normalized to unit length"""
        if np.allclose(a, b):
            return np.array([1, 0, 0])
        c = np.cross(a, b)
        return c / np.linalg.norm(c)

    def align_pcl(self, pcl, get_angles=False):
        pcl_aligned = pcl.copy()

        angles = list()
        
        for axidx in range(2):
            #choose an axis to align:
            axis_vector = np.zeros(3)
            axis_vector[axidx] = 1

            ##get pmi:
            pmi = self.get_pmi(pcl_aligned)
            
            #get angle to that vector
            angle = self.angle_between(pmi[axidx], axis_vector)
            
            #get axis around which to rotate
            ax = self.rotaxis(pmi[axidx], axis_vector)
            
            q = self.from_axis_angle(ax, -angle) 
            nq = self.fast_normalise(q)
            rotmat = self.quaternion_rotation_matrix(nq)
            
            #rotmat = rotation_matrix(angle, ax)[:3,:3].T
            pcl_aligned = pcl_aligned.dot(rotmat)
            if get_angles:
                angles.append(rotmat)                

        if get_angles:
            return angles
        else:
            return pcl_aligned
