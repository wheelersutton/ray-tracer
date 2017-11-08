# David Wheeler Sutton
# The most important part of this code is the interpreter, which will
# help you parse the scene description (.cli) files.

from components import *
from components import Vertex
from components import Vector

#scene variables
fov = 0 #field of view
bg = Color(0, 0, 0) #background color
shapeidx = 0
lightsources = []
currentSurface = 0
shapes = []
vertices = []
ray_shape_intersections = [] #[t value, Vertex of intersection, reference to the shape]
loops = 0

def setup():
    size(500, 500) 
    noStroke()
    colorMode(RGB, 1.0)  # Processing color values will be in [0, 1]  (not 255)
    background(0, 0, 0)

# read and interpret the appropriate scene description .cli file based on key press
def keyPressed():
    if key == '1':
        interpreter("i1.cli")
    elif key == '2':
        interpreter("i2.cli")
    elif key == '3':
        interpreter("i3.cli")
    elif key == '4':
        interpreter("i4.cli")
    elif key == '5':
        interpreter("i5.cli")
    elif key == '6':
        interpreter("i6.cli")
    elif key == '7':
        interpreter("i7.cli")
    elif key == '8':
        interpreter("i8.cli")
    elif key == '9':
        interpreter("i9.cli")
    elif key == '0':
        interpreter("i10.cli")

def interpreter(fname):
    resetScene()
    fname = "data/" + fname
    # read in the lines of a file
    with open(fname) as f:
        lines = f.readlines()

    # parse each line in the file in turn
    for line in lines:
        words = line.split()  # split the line into individual tokens
        if len(words) == 0:   # skip empty lines
            continue
        if words[0] == 'sphere':
            radius = float(words[1])
            x = float(words[2])
            y = float(words[3])
            z = float(words[4])
            # call your sphere creation routine here
            # for example: create_sphere(radius,x,y,z)
            s = Sphere(Vertex(x, y, z), radius, currentSurface)
            shapes.append(s)
        elif words[0] == 'fov':
            global fov
            fov = float(words[1])
        elif words[0] == 'background':
            global bg
            r = float(words[1])
            g = float(words[2])
            b = float(words[3])
            bg = Color(r, g, b)
        elif words[0] == 'light':
            x = float(words[1])
            y = float(words[2])
            z = float(words[3])
            r = float(words[4])
            g = float(words[5])
            b = float(words[6])
            lightsources.append(Light(Color(r, g, b), Vertex(x, y, z)))
        elif words[0] == 'surface':
            #may need to make the values floats, not strings
            Cd = Color(float(words[1]), float(words[2]), float(words[3]))
            Ca = Color(float(words[4]), float(words[5]), float(words[6])) #maybe just zero this out for now
            Cs = Color(float(words[7]), float(words[8]), float(words[9])) #maybe just zero this out for now
            P = float(words[10])
            Krefl = float(words[11])
            currentSurface = SurfaceMaterial(Cd, Ca, Cs, P, Krefl)
        elif words[0] == 'begin':
            del vertices[:] #clear out vertices for new shape
        elif words[0] == 'vertex':
            tVertex = Vertex(float(words[1]), float(words[2]), float(words[3]))
            vertices.append(tVertex)
        elif words[0] == 'end':
            tri = Triangle(vertices[0], vertices[1], vertices[2], currentSurface)
            shapes.append(tri)
        elif words[0] == 'write':
            render_scene()    # render the scene
            save(words[1])  # write the image to a file
            pass

# render the ray tracing scene
def render_scene():
    k = float(tan(radians(fov)/2))
    ray = 0
    count = 0
    for j in range(height):
        for i in range(width):
            # create an eye ray for pixel (i,j) and cast it into the scene
            z = -1
            x_prime = float( i / abs(z) )
            y_prime = float( (height - j) / abs(z) )
            x = (x_prime - float(width/2)) * float((2*k)/width)
            y = (y_prime - float(height/2)) * float((2*k)/height)
            ray = Ray(Vertex(0, 0, 0), Vertex(x, y, -1)) #change first vert to eye position
            
            for s in shapes:
                if (s.type == "sphere"):
                    intersect_ray_sphere(ray, s)
                elif (s.type == "triangle"):
                    intersect_ray_triangle(ray, s)
            
            smallest_t = 100000 #arbitrary large value that t would never be
            smallest_intersection = 0
            smallest_shape = 0
            for intersection in ray_shape_intersections:
                if (intersection[0] < smallest_t and intersection[0] >= 0):
                    smallest_t = intersection[0]
                    smallest_intersection = intersection[1]
                    smallest_shape = intersection[2]

            #calculate color 
            pix_color = 0
            if (smallest_t != 100000): #again, arbitrary large value
                pix_color = calculate_color(smallest_intersection, smallest_shape, 0)
            else:
                pix_color = bg
            
            pix_color = color(pix_color.r, pix_color.g, pix_color.b) #convert Color to color (why do i do this to myself?)
            set (i, j, pix_color)         # fill the pixel with the calculated color
            
            

def intersect_ray_sphere(r, s):
    dx = r.vert.x - r.origin.x #slope of x
    dy = r.vert.y - r.origin.y #slope of y
    dz = r.vert.z - r.origin.z #slope of z
    x0 = r.origin.x
    y0 = r.origin.y
    z0 = r.origin.z
    cx = s.vert.x
    cy = s.vert.y
    cz = s.vert.z
    r = s.radius
    a = dx**2 + dy**2 + dz**2
    b = 2*((x0*dx - cx*dx) + (y0*dy - cy*dy) + (z0*dz - cz*dz))
    c = (x0 - cx)**2 + (y0 - cy)**2 + (z0 - cz)**2 - r**2
    discr = b**2 - 4*a*c
    
    #figure out # of intersections based on discriminant value
    if (discr > 0):
        t1 = (-b + sqrt(discr)) / (2*a)
        t2 = (-b - sqrt(discr)) / (2*a)
        ray_shape_intersections.append([t1, Vertex(x0+t1*dx, y0+t1*dy, z0+t1*dz), s])
        ray_shape_intersections.append([t2, Vertex(x0+t2*dx, y0+t2*dy, z0+t2*dz), s])
    elif (discr == 0):
        t1 = (-b + sqrt(discr)) / (2*a)
        ray_shape_intersections.append([t1, Vertex(x0+t1*dx, y0+t1*dy, z0+t1*dz), s])
        
def intersect_ray_triangle(r, s):
    #A = v1, B = v2, C = v3
    AB = Vector.subtract_vectors(s.v2, s.v1)
    BC = Vector.subtract_vectors(s.v3, s.v2)
    AC = Vector.subtract_vectors(s.v3, s.v1)
    vec = Vector.normalize_vector(Vector.cross_product(AB, BC)) #N
    
    #prevents a division by 0 error when ray and plane are parallel
    raydir = Vector.dot_product(vec, r.direction)
    if (abs(raydir) < 0.0000001):
        return 0
    
    d = -Vector.dot_product(vec, s.v1)
    
    t = -(Vector.dot_product(vec, Vertex.vertex_to_vector(r.origin)) + d) / Vector.dot_product(vec, r.direction)
    
    intersect = Vector.add_vectors(Vertex.vertex_to_vector(r.origin), Vector(t*r.direction.x, t*r.direction.y, t*r.direction.z))
    
    #half-plane test
    
    edge1 = Vector.subtract_vectors(s.v2, s.v1)
    p1 = Vector.subtract_vectors(intersect, s.v1)
    C1 = Vector.cross_product(edge1, p1)
    
    edge2 = Vector.subtract_vectors(s.v3, s.v2)
    p2 = Vector.subtract_vectors(intersect, s.v2)
    C2 = Vector.cross_product(edge2, p2)
    
    edge3 = Vector.subtract_vectors(s.v1, s.v3)
    p3 = Vector.subtract_vectors(intersect, s.v3)
    C3 = Vector.cross_product(edge3, p3)
    
    if (Vector.dot_product(vec, C1) > 0 and Vector.dot_product(vec, C2) > 0 and Vector.dot_product(vec, C3) > 0):
        ray_shape_intersections.append([t, Vertex(intersect.x, intersect.y, intersect.z), s])
    
        
def calculate_color(intersect, obj, depth):
    s = obj.surface
    sum_r = s.Ca.r
    sum_g = s.Ca.g
    sum_b = s.Ca.b
    
    if (obj.type == "sphere"):
        surface_normal = Vector.normalize_vector(Vector.subtract_vectors(Vertex.vertex_to_vector(intersect), Vertex.vertex_to_vector(obj.vert)))
    elif (obj.type == "triangle"):
        AB = Vector.subtract_vectors(obj.v1, obj.v2)
        AC = Vector.subtract_vectors(obj.v1, obj.v3)
        normalvec = Vector.cross_product(AB, AC)
        surface_normal = Vector.normalize_vector(normalvec)
        if(surface_normal.z <= 0):
            surface_normal = Vector.subtract_vectors(Vector(0, 0, 0), surface_normal)
    
    for light in lightsources:
        l = Vertex.vertex_to_vector(light.vert)
        light_vector = Vector.normalize_vector(Vector.subtract_vectors(l, intersect))
        light_color = light.c
        
        offset = Vector(light_vector.x * 0.001, light_vector.y * 0.001, light_vector.z * 0.001)
        
        intersect_to_light = Ray(Vertex(intersect.x + offset.x, intersect.y + offset.y, intersect.z + offset.z), light.vert)
        
        del ray_shape_intersections[:]
        for shp in shapes:
            if (shp.type == "sphere"):
                intersect_ray_sphere(intersect_to_light, shp)
            elif (shp.type == "triangle"):
                intersect_ray_triangle(intersect_to_light, shp)
        
        smallest_t = 100000 #arbitrary large value that t would never be
        smallest_intersection = 0
        smallest_shape = 0
        for intersection in ray_shape_intersections:
            if (intersection[0] < smallest_t and intersection[0] >= 0):
                smallest_t = intersection[0]
                smallest_intersection = intersection[1]
                smallest_shape = intersection[2]
                
        shade = 1
        if smallest_intersection != 0:
            shade = 0
        
        #diffuse stuff
        Cd_max = max(0, Vector.dot_product(surface_normal, light_vector))
        
        #phong stuff
        v = Vector.subtract_vectors(Vector(0, 0, 0), intersect)
        h = Vector.normalize_vector(Vector.add_vectors(v, light_vector))
        phong = float(max(0, pow(Vector.dot_product(h, surface_normal), obj.surface.P)))
        
        sum_r += (s.Cd.r * light.c.r * Cd_max * shade) + (s.Cs.r * light.c.r * phong * shade)
        sum_g += (s.Cd.g * light.c.g * Cd_max * shade) + (s.Cs.g * light.c.g * phong * shade)
        sum_b += (s.Cd.b * light.c.b * Cd_max * shade) + (s.Cs.b * light.c.b * phong * shade)
    
    #if (s.Krefl > 0 and depth < 5):
        #calculate_color(intersect, obj, depth + 1)

        del ray_shape_intersections[:]
    return Color(min(1, sum_r), min(1, sum_g), min(1, sum_b))

#before rendering a new image
def resetScene():
    global shapeidx, bg
    del lightsources[:]
    del shapes[:]
    currentSurface = 0
    shapeidx = 0
    bg = Color(0, 0, 0)
    
# should remain empty for this assignment
def draw():
    pass