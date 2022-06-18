import numpy as np

class Quaternion():

    def __init__(self, real, ipart, jpart, kpart):
        self.q = np.array([real, ipart, jpart, kpart])
    
    def __str__(self):
        return f"[{self.q[0]}, {self.q[1]} i + {self.q[2]} j + {self.q[3]} k]"

    def __add__(self, other):
        return Quaternion(*(self.q + other.q))
    
    def __sub__(self, other):
        return Quaternion(*(self.q - other.q))

    def __mul__(self, other):
        
        if isinstance(other, Quaternion):
            real = self.q[0]*other.q[0] - np.dot(self.q[1:], other.q[1:])
            cross = np.cross(self.q[1:], other.q[1:])
            vector = self.q[0]*other.q[1:] + other.q[0]*self.q[1:] + cross
            return Quaternion(real, vector[0], vector[1], vector[2])

        if isinstance(other, int) or isinstance(other, float):
            return Quaternion(*self.q * other)
    
    def conjugate(self):
        return Quaternion(self.q[0], -self.q[1], -self.q[2], -self.q[3])

    def inverse(self):
        conj = self.conjugate()
        return Quaternion(conj.q[0]/self.norm()**2, conj.q[1]/self.norm()**2,
                             conj.q[2]/self.norm()**2, conj.q[3]/self.norm()**2)
    
    def normalize(self):
        return Quaternion(self.q[0]/self.norm(), self.q[1]/self.norm(), 
                            self.q[2]/self.norm(), self.q[3]/self.norm())

    def norm(self):
        return np.sqrt(np.sum(self.q**2))

    def rotate(self, other):
        q = other.normalize()
        inv = other.inverse()
        return q * self * inv

    def tovec(self):
        return(np.round(self.q[1:], 8))

    def liexy(self):
        return np.arctan2(-self.q[3], self.q[2])

    def liex(self):
        ang = self.liexy()
        tmp = -((np.cos(ang) * self.q[2]) - (np.sin(ang) * self.q[3]))
        return ang, np.arctan2(tmp, self.q[1])

    @staticmethod
    def angaxis(angle, axis):
        '''
        Create the quaternion needed to rotate by a given angle for a
        given arbitrary axis.
        '''
        axis = np.array(axis)
        axis = axis / np.linalg.norm(axis)
        axis = np.sin(angle / 2) * axis
        return Quaternion(np.cos(angle / 2), *axis)

    @staticmethod
    def vec2quat(vec):
        '''
        Convert vector to quaternion.
        '''
        return Quaternion(*np.array([0,*vec]))