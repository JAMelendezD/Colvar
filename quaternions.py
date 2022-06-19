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
        return other * self * other.conjugate()

    def tovec(self):
        return np.round(self.q[1:], 8)

    def liexy(self):
        return np.arctan2(-self.q[3], self.q[2])

    def liex(self):
        ang = self.liexy()
        tmp = -((np.cos(ang) * self.q[2]) - (np.sin(ang) * self.q[3]))
        return ang, np.arctan2(tmp, self.q[1])

    def half_angle(self, other):
        cos_half_ang = np.dot(self.normalize().q, other.normalize().q)
        return np.arccos(cos_half_ang), cos_half_ang 

    def interpolate(self, other, t, half_ang, cos_half_ang):
        t = t/2.0 # Not sure why this fixes it
        sin_half_ang = np.sqrt(1-cos_half_ang*cos_half_ang)
        numerator = self*np.sin((1.0-t)*half_ang) + other*np.sin(t*half_ang)
        denominator = 1.0/sin_half_ang
        return numerator*denominator

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