class Vertex(object):
    x = 0.0
    y = 0.0
    z = 0.0
    
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
    
    def vertex_to_vector(a):
        return Vector(a.x, a.y, a.z)
    
class Vector(object):
    x = 0.0
    y = 0.0
    z = 0.0
    
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
    
    def dot_product(a, b):
        return (a.x*b.x + a.y*b.y + a.z*b.z)
    
    def normalize_vector(a):
        a_mag = sqrt(a.x**2 + a.y**2 + a.z**2)
        return Vector(a.x/a_mag, a.y/a_mag, a.z/a_mag)
    
    def subtract_vectors(a, b):
        return Vector(a.x-b.x, a.y-b.y, a.z-b.z)
    
    def add_vectors(a, b):
        return Vector(a.x+b.x, a.y+b.y, a.z+b.z)
    
    def cross_product(a, b):
        return Vector(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x)
    
    

class Ray(object):
    origin = Vertex(0, 0, 0) #vertex
    vert = Vertex(0, 0, 0) #vertex
    direction = 0
    def __init__(self, origin, vert):
        self.origin = origin
        self.vert = vert
        self.direction = Vector(vert.x - origin.x, vert.y - origin.y, vert.z - origin.z) 

class Sphere(object):
    type = "sphere"
    vert = Vertex(0, 0, 0)
    radius = 0
    surface = 0
    
    def __init__(self, vert, radius, surface):
        self.vert = vert
        self.radius = radius
        self.surface = surface
        
class Triangle(object):
    type = "triangle"
    v1 = 0
    v2 = 0
    v3 = 0
    vert = 0
    surface = 0
    
    def __init__(self, v1, v2, v3, surface):
        self.v1 = Vertex.vertex_to_vector(v1)
        self.v2 = Vertex.vertex_to_vector(v2)
        self.v3 = Vertex.vertex_to_vector(v3)
        self.surface = surface
        self.vert = v1
        
class Color(object):
    r = 0
    g = 0
    b = 0
    
    def __init__(self, r, g, b):
        self.r = r
        self.g = g
        self.b = b
        
    def color_to_vector(a):
        return Vector(a.r, a.g, a.b)
        

class Light(object):
    c = Color(0, 0, 0)
    vert = Vertex(0, 0, 0)
    
    def __init__(self, c, vert):
        self.vert = vert
        self.c = c

class SurfaceMaterial(object):
    Cd = 0
    Ca = 0
    Cs = 0
    P = 0
    Krefl = 0
    def __init__(self, Cd, Ca, Cs, P, Krefl):
        self.Cd = Cd #diffuse
        self.Ca = Ca #ambient
        self.Cs = Cs #specular
        self.P = P #Specular Power - Phong
        self.Krefl = Krefl #reflection coefficient

    