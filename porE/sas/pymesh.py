
import math
import numpy as np 
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay,ConvexHull,SphericalVoronoi
from ase.visualize import view 
from ase.atoms import Atoms 
from ase.data import atomic_numbers, covalent_radii, vdw_radii
from ase.io import read 

# Shoelace method 
#- https://ysjournal.com/tetrahedral-shoelace-method-calculating-volume-of-irregular-solids/

# Volume 
# - https://mathematica.stackexchange.com/questions/25811/how-to-calculate-the-volume-of-a-convex-hull

# Surface area, volume of convex hull 
# - https://mathematica.stackexchange.com/questions/6908/how-to-calculate-volume-of-convex-hull-and-volume-of-a-3d-object
# - https://stackoverflow.com/questions/24733185/volume-of-convex-hull-with-qhull-from-scipy

# Plotting 
# - https://stackoverflow.com/questions/27270477/3d-convex-hull-from-point-cloud

def get_radius(atom):
    """ 
        Get the radii for a given ase.atom 
    """
    #tmp = vdw_radii.tolist() 
    #tmp.append([40,2.36])
    #vdw_radii = numpy.array(tmp)
    tmp = vdw_radii.copy()
    tmp[40] = 2.36
    #print(tmp[40])
    return tmp[atomic_numbers[atom.symbol]]

def get_points_on_sphere(n, radius=1, center=[0,0,0], verbose=3):
    """
        Get points on a sphere 
        using the Golden-Section Spiral algorithm.
        Input: 
            n: number of points 
        Output: 
            points: coordinates for points on a sphere 
    """
    points = []
    inc = math.pi * (3 - math.sqrt(5))
    offset = 2 / float(n)
    for k in range(int(n)):
        y = k * offset - 1 + (offset / 2)
        #print('y : {}'.format(y))
        r = math.sqrt(1 - y*y)
        phi = k * inc
        points.append([math.cos(phi)*r, y, math.sin(phi)*r])
    com = np.array(points).sum(axis=0)/len(points)
    points = np.array(points)*radius+np.array(center)
    com = np.array(points).sum(axis=0)/len(points)
    r0 = np.linalg.norm(points[0]-center)
    if verbose >= 4: 
        print('COM = {}'.format(com))
        print('r = {}'.format(r0))
    return points, r0, com

def PolygonArea(corners):
    # 2D 
    # https://stackoverflow.com/questions/24467972/calculate-area-of-polygon-given-x-y-coordinates
    # https://en.wikipedia.org/wiki/Shoelace_formula
    n = len(corners) # of corners
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += corners[i][0] * corners[j][1]
        area -= corners[j][0] * corners[i][1]
    area = abs(area) / 2.0
    return area

def TriSurfArea(tri_surf_points):
    A = 0
    for p_tri in tri_surf_points:
        A += PolygonArea(p_tri)
    print('SA = {} '.format(A))

def tetrahedron_volume(a, b, c, d):
    return np.sum(np.abs(np.einsum('ij,ij->i', a-d, np.cross(b-d, c-d))) / 6)

def TetraSurfVol(pts):
     V = 0 
     for pt in pts:
         V += tetrahedron_volume(np.array([pt[0]]),np.array([pt[1]]),np.array([pt[2]]),np.array([pt[3]])) 
     print('V(TetraSurfVol) = {} '.format(V))

def tetrahedron_area(a,b,c):
    return 1/2.*np.linalg.norm(np.cross(b-a,c-a))

def TetraSurfA(pts):
     A = 0 
     for pt in pts: 
         A += tetrahedron_area(pt[0],pt[1],pt[2])
         #A += tetrahedron_area(pt[0],pt[1],pt[3])
         #A += tetrahedron_area(pt[1],pt[2],pt[3])
     print('A(TetraSurfA) = {} '.format(A))

def TetraVol(a):
    """ Analytical tetrahedron volume """
    # Ref.: https://en.wikipedia.org/wiki/Tetrahedron
    V = a**3/(6*np.sqrt(2))
    SA = np.sqrt(3)*a**2 
    print('V(tetrahedron,anayltically) = {} '.format(V))
    print('SA(cube,analytical) = {}'.format(SA))
    return V

def CubeVol(a):
    """ Analytical tetrahedron volume """
    # Ref.: https://en.wikipedia.org/wiki/Tetrahedron
    V = a**3
    SA = 6*a**2 
    print('V(cube,analytical) = {} '.format(V))
    print('SA(cube,analytical) = {}'.format(SA))
    return V

def SphereVol(r): 
    V= 4/3.*np.pi*r**3 
    SA = 4*np.pi*r**2 
    print('V(sphere,analytical) = {} '.format(V))
    print('SA(sphere,analytical) = {} '.format(SA))

def test_area():
    points = np.array([[0, 0], [0, 2], [2, 0], [2, 2]])
    tri = Delaunay(points)
    # strange 
    hull = ConvexHull(points)
    print(hull.area)
    print(hull.volume)
    
    idx = tri.simplices
    tri_surf_points = points[idx]
    TriSurfArea(tri_surf_points)
    plt.triplot(points[:,0], points[:,1], tri.simplices)
    plt.plot(points[:,0], points[:,1], 'o')
    plt.show()

def get_delaunay(pts): 
    # - each set is a tetrahedra 
    # - sum of tetrahedra volumes give volume 
    tri = Delaunay(pts)
    idx = tri.simplices
    tetra_pts = tri.points[tri.simplices]
    TetraSurfVol(tetra_pts)


def get_convexhull(pts): 
    # - each set is a triangle 
    # - some of all triangles of a convex hull gives surface area 
    hull = ConvexHull(pts)
    print('SA(convexhull) : {}'.format(hull.area))
    print('V(convexhull) : {}'.format(hull.volume))    
    pts_hull = pts[hull.simplices]
    TetraSurfA(pts_hull)
    simplices = np.column_stack((np.repeat(hull.vertices[0], hull.nsimplex),hull.simplices))
    tetra_pts = hull.points[simplices]
    TetraSurfVol(tetra_pts)

def get_sphericalvoronoi(pts,radius,com): 
    sv = SphericalVoronoi(points=pts) #,radius=radius,center=com)
    SA = sv.calculate_areas().sum()
    print('SA(sphericalvoronoi) : {}'.format(SA))


def test_tetrahedron():
    print('test: tetrahedron')
    # Tetrahedron 
    tetra = np.array([[1,1,1], [1,-1,-1],[-1,1,-1], [-1,-1,1]])
    pts = tetra
    a_tetra = 2*np.sqrt(2)
    TetraVol(a_tetra)
    get_delaunay(pts)
    get_convexhull(pts)

def test_cube(): 
    print('test: cube')
    # Cube 
    a = [-1, -1, -1]
    b = [-1, -1, +1]
    c = [-1, +1, -1]
    d = [-1, +1, +1]
    e = [+1, -1, -1]
    f = [+1, -1, +1]
    g = [+1, +1, -1]
    h = [+1, +1, +1]
    cube = np.array([a, b, c, d, e, f, g, h])
    # struct = Atoms(8*'X',cube) 
    # view(struct)
    pts = cube
    a_cube = 2
    CubeVol(a_cube)
    get_delaunay(pts)
    get_convexhull(pts)

def test_sphere(): 
    print('Test: one sphere')
    pts, r, com  = get_points_on_sphere(96)
    print('r : {}'.format(r))
    SphereVol(r)
    get_delaunay(pts)
    #get_convexhull(pts)
    hull = ConvexHull(pts)
    print('SA(convexhull) : {}'.format(hull.area))
    print('V(convexhull) : {}'.format(hull.volume))
    # struct = Atoms(len(pts)*'X'+'X',pts.append([0,0,0])) 
    # view(struct)
    #get_sphericalvoronoi(pts,r,com)


def test_two_spheres():
    print('Test: two spheres') 
    pts = []
    A_pts, A_r, A_com  = get_points_on_sphere(960,radius=1.2,center=[0,0,0])
    B_pts, B_r, B_com  = get_points_on_sphere(960,radius=2.02,center=[0,0,3])
    pts.extend(A_pts) 
    pts.extend(B_pts) 
    
    #SphereVol(r)
    get_delaunay(pts)
    ##get_convexhull(pts)
    hull = ConvexHull(pts)
    print('SA(convexhull) : {}'.format(hull.area))
    print('V(convexhull) : {}'.format(hull.volume))
    ## struct = Atoms(len(pts)*'X'+'X',pts.append([0,0,0])) 
    ## view(struct)
    #get_sphericalvoronoi(pts)

def vdw(atoms,n_points,verbose=3):
    """ 
        Calculate vdW surface using convex hull 
        Input: 
            atoms: ase.atoms 
            n_points: number of points on a sphere 
            verbose: verbosity 
        Output: 
            SA: vdW surface area [AA**2]
            V : vdW volume [AA**3] 
    """
    # NOTE: it seems to work for none-pbc structues  
    pts = []  
    for atom in atoms:
        tmp_pts, tmp_r, tmp_com  = get_points_on_sphere(n_points,radius=get_radius(atom),center=atom.position)
        pts.extend(tmp_pts)
    hull = ConvexHull(pts)
    SA = hull.area
    V = hull.volume  
    if verbose >= 3: 
        print('Surface area SA(convexhull,n:{}) : {} AA**2'.format(n_points,SA))
        print('Volume V(convexhull,n:{}) : {} AA**3'.format(n_points,V))
    return SA,V

def test_volume_surface(): 
    test_tetrahedron() 
    print('\n') 
    test_cube() 
    print('\n') 
    test_sphere()
    print('\n')
    test_two_spheres()

def test_C6H6():
    f_file ='236.mol'
    atoms=read(f_file)
    vdw(atoms,n_points=96)


if __name__ == '__main__':
    test_volume_surface()
    test_C6H6()
