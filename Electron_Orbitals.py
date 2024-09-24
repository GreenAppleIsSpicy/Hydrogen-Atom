import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib as mpl
import sympy as sym
import open3d as o3d

a=1 #Bohr Radius, the radius of the heighest electron probability density of the 1s Hydrogen atom, defines the unit of distance in this program

def N(n, l):
    """
    The Normalization for the radial wavefunction.
    
    Takes input:
    n: the principal quantum number
    l: the angular momentum quantum number

    Returns output:
    normalization for the radial wavefunction.
    """
    ro = 2/(n*a)
    return np.sqrt(ro**3 * math.factorial(n-l-1)/(2*n*math.factorial(n+l)**3))

def L(r, s, x):
    """
    The associated Laguerre polynomials.

    Takes input:
    r: distance from center (in Bohr radii) 
    s: Laguerre index 1
    x: Polynomial variable in the Real numbers as a numpy array

    Returns output:
    The associated Laguerre polynomial evaluated at x as a numpy array
    """
    k = 0
    for q in range(0, s+1):
        k += (-1)**q * (math.factorial(s+r)) / (math.factorial(s-q)*math.factorial(r+q)) * x**q/math.factorial(q)
    return k

def F(n, l, r):
    """
    Unnormalized radial wavefunction.

    Takes input:
    n: principal quantum number
    l: angular momentum quantum number
    r: distance from nucleus (in Bohr radii) as numpy array must be >= 0
    
    Return output:
    unnormalized radial wavefunction as numpy array
    """
    ro = 2/(n * a)
    return (ro * r)**l * L(2*l+1, n-l-1, ro*r) * np.exp(-ro * r/2)

def R(n, l, r):
    """
    The radial wavefunction.
    
    Takes input:
    n: principal quantum number
    l: angular momentum quantum number
    r: distance from nucleus (in Bohr radii) as numpy array must be >= 0

    Returns output:
    radial wavefunction as numpy array
    """
    return N(n, l) * F(n, l, r)

def P(l, m, x):
    """
    The associated Legendre polynomials.
    
    Takes input:
    l: angular momentum quantum number
    m: magnetic quantum number
    x: continuous variable on the bounds [0, 1] as numpy array

    Returns output:
    The value of the associated Legendre polynomials at x as a numpy array
    """
    y = sym.Symbol('y')
    d = sym.diff((y**2-1)**l, y, l+np.abs(m))
    D = sym.lambdify(y, d, "numpy")
    if m<0:
        q = (-1)**m * math.factorial(l-m)/math.factorial(l+m)
    else:
        q=1
    return q * (-1)**m / (2**l * math.factorial(l)) * (1-x**2)**(np.abs(m)/2) * D(x)


def Y(l, m, theta, phi):
    """
    The spherical harmonics/spherical wavefunction.

    Takes input:
    l: angular momentum quantum number
    m: magnetic quantum number
    theta: zenith angle
    phi: azimuthal angle
    """
    k = np.sqrt((2 * l + 1) / (4 * np.pi) * math.factorial(l-m)/math.factorial(l+m))
    Plm = P(l, m, np.cos(theta))
    e = np.exp(1j * m * phi)
    return k * Plm * e

def u(n, l, r):
    """
    The radial probability density.

    Takes input:
    n: principal quantum number
    l: angular momentum quantum number
    r: distance from nucleus (in Bohr radii) as numpy array must be >= 0

    Returns output:
    The value of the radial probability density at r as a numpy array
    """
    return np.abs(r * R(n, l, r))**2

def uy(l, m, theta, phi):
    """
    The spherical probability density.

    Takes input:
    l: angular momentum quantum number
    m: magnetic quantum number
    theta: zenith angle as a numpy array on [0, π]
    phi: azimuthal angle as a numpy array on [0, 2π)

    Returns output:
    The spherical probability density
    """
    return np.real(Y(l, m, theta, phi))**2

def up(m, phi):
    """
    The azimuthal probability density.

    Takes input:
    m: magnetic quantum number
    phi: azimuthal angle as a numpy array on [0, 2π)

    Returns output:
    The azimuthal probability density
    """
    return np.real(np.exp(1j * m * phi))**2

def Norm_integrate(func, args):
    """
    Integrates over func with variable kwargs and normalizes.

    Takes input:
    func: the function being integrated over
    args: tuple of function inputs

    Returns output:
    Normalized integral of the specified function over the given range as a numpy array
    
    Example: If the function was f(arg1, arg2, arg3), then func = f and arg = (arg1, arg2, arg3,)
    Note: If func has a single argument then arg should be input as arg=(arg1,) i.e. with the comma, otherwise will return an error.
    """
    I = np.zeros_like(func(*args))
    I[0] = func(*args)[0]
    for i in range(1, len(func(*args))):
        I[i] = I[i-1] + func(*args)[i]
    return I/np.max(I)

def MonteCarlo(func, n, l, m, count=250, res=1000, rmax=10):
    """
    Shows plot where points are chosen based on chosen probability distribution.

    Takes input:
    func:
    n: principal quantum number
    l: angular momentum quantum number
    m: magnetic quantum number

    (Optional)
    count: number of scatter points [by default this value is set to 250]
    res: the number of divistions of the specified variable [by default this value is set to 1000]
    rmax: the maximum value for r (in Bohr radii) [by default this value is set to 10]

    Returns output:
    Plot of the probability distribution and the chosen points
    """
    plt.close('all')
    checking = plt.figure()
    axes = checking.add_subplot()

    if func==up:
        phi = np.linspace(0, 2 * np.pi, res)
        p = np.random.rand(count)
        Ps = np.interp(p, Norm_integrate(up, (m, phi,)), phi)
    
        axes.plot(phi, up(m, phi))
        # axes.plot(phi, Norm_integrate(up, (m, phi,)))
        plt.scatter(Ps, np.zeros_like(Ps))

    elif func==uy:
        theta = np.arccos(2 * np.linspace(0, 1, res) -1)
        t = np.random.rand(count)
        Ts = np.interp(t, Norm_integrate(uy, (l, m, theta, 0,)), theta)
    
        axes.plot(theta, uy(l, m, theta, 0))
        # axes.plot(theta, Norm_integrate(uy, (l, m, theta, 0,)))
        plt.scatter(Ts, np.zeros_like(Ts))

    elif func==u:
        r = np.linspace(0, rmax, res)
        k = np.random.rand(count)
        Rs = np.interp(k, Norm_integrate(u, (n, l, r,)), r)

        axes.plot(r, u(n, l, r))
        # axes.plot(r, Norm_integrate(u, (n, l, r,)))
        plt.scatter(Rs, np.zeros_like(Rs), s=0.1, alpha=0.1)

    else:
        return 'This function is not yet impletmented please choose one of the defined functions from this list, [u, uy, up]'

def sph2crt(r, a, ze):
    """
    Converts from spherical coordinates to Cartesian coordinates.

    Takes input:
    r: radial coordinate
    a: azimuthal coordinate (radians)
    ze: zenith coordinate (radians)

    Returns output:
    x: the x-axis is the axis from which azimuthal angle is measured.
    y: the axis \pi radians along the azimuth from the x-axis.
    z: the axis perpendicular to the azimuthal plane, in the direction of zenith = 0 radians.
    """
    x = r * np.cos(a) * np.sin(ze)
    y = r * np.sin(a) * np.sin(ze)
    z = r * np.cos(ze)
    return x, y, z

def crt2sph(x, y, z):
    """
    Converts from Cartesian coordinates to spherical coordinates.

    Takes input:
    x: x coordinate
    y: y coordinate
    z: z coordinate

    Returns output:
    r: the radial coordinate, defines the distance from the origin (x, y, z) = (0, 0, 0)
    a: the azimuthal coordinate, defined by the angle from the x axis in the xy-plane. (radians)
    ze: the zenith coordinate, defined by the angle from the positive z axis (radians)

    Note: can take numpy array input and will output as numpy array.
    
    Citation: Based on cart2sph() in MATLAB, uses diffrent convention for the zenith.
    """
    r = np.sqrt(x**2+y**2+z**2)
    a = np.arctan2(y,x)
    ze = np.arctan2(np.sqrt(x**2+y**2),z)
    return r, a, ze

def Plot_Hydrogen_2D(n, l, m, rmax, res=1000, cont=True, scat=True, scatter_count=10000):
    """
    Plots a 2D image of the Hydrogen orbital as a heatmap, with either a contour or a Monte-Carlo based scatterplot.

    Takes input:
    n: the principal quantum number
    l: the angular momentum quantum number
    m: the magnetic quantum number
    rmax: the maximum radius in Bohr radii you would like to see

    (Optional)
    res: the resolution of the coordinate space, the number of subdivisions in the specified range [by default this value is set to 1000]
    cont: True or False of if a contour plot should be shown [by default this value is set to True]
    scat: True or False of if a scatter plot should be shown [by default this value is set to True]
    scatter_count: The number of points on the scatterplot [by default this value is set to 10000]
    """
    r = np.linspace(0, rmax, res)
    theta = np.linspace(0, 2 * np.pi, res)
    phi = np.linspace(0, 2 * np.pi, res)

    k = np.random.rand(scatter_count)
    t = np.random.rand(scatter_count)
    p = np.random.rand(scatter_count)
    Rs = np.interp(k, Norm_integrate(u, (n, l, r)), r)
    Ts = np.interp(t, Norm_integrate(uy, (l, m, theta, 0)), theta)
    Ps = 2 * np.pi * p

    r, theta = np.meshgrid(r, theta)
    fig1, ax1 = plt.subplots(subplot_kw={'projection':'polar'})
    ax1.set_facecolor('black')

    if cont==True:
        ax1.contourf(theta, r, np.abs(r * R(n, l, r) * Y(l, m, theta, 0))**2, levels=100, cmap='inferno')
        
    if scat==True:
        ax1.scatter(Ts, Rs, s=.1, alpha=0.7, c=np.abs(Rs * R(n, l, Rs) * Y(l, m, Ts, 0))**2, cmap='inferno')

def Plot_Hydrogen_3D(n, l, m, rmax, res=1000, scatter_count=1000, plot=True, plysave=False, cutout=False, 
                     heatmap=True, colormap='inferno', color='red'):
    """
    Creates a 3D scatter plot of the specified Hydrogen orbital with a heatmap applied.

    Takes input:
    n: the principal quantum number
    l: the angular momentum quantum number
    m: the magnetic quantum number
    rmax: the maximum radius in Bohr radii you would like to see

    (Optional)
    res: the resolution of the coordinate space, the number of subdivisions in the specified range [by default this value is set to 1000]
    scatter_count: The number of points on the scatterplot [by default this value is set to 10000]
    plot: True or False of whether the output should be plotted into a matplotlib 3D plot [by default this value is set to True]
    plysave: True or False of whether the plot is save to the file '...\Electron_Positions\Electron_Positions_(n, l, m).ply'[by default this value is set to False]
    cutout: True or False of whether the plot has a cut taken out of it [by default this value is set to False]
    heatmap: True or False [by default this value is set to True]
    colormap: string input [by default this value is set to 'red']
    color: string input [by default this value is set to 'inferno']

    Note: for the input res, the speed of the output is proportional to its value. i.e. doubling the value doubles the amount of time to output
    """

    r = np.linspace(0, rmax, res)
    theta = np.arccos(2 * np.linspace(0, 1, res) -1)
    phi = np.linspace(0, 2 * np.pi, res)

    k = np.random.rand(scatter_count)
    t = np.random.rand(scatter_count)
    p = np.random.rand(scatter_count)
    Rs = np.interp(k, Norm_integrate(u,  (n, l, r)         ),     r)
    Ts = np.interp(t, Norm_integrate(uy, (l, m, theta, 0)), theta)
    Ps = np.interp(p, Norm_integrate(up, (m, phi,)), phi)

    x, y, z = sph2crt(Rs, Ps, Ts)

    if cutout==True:
        x_elim, y_elim, z_elim = x, y, z
        Positions_elim = np.array([x_elim, y_elim, z_elim]).T
    
        i=0
        for x, y, z, in Positions_elim:
            if x>=0 and y<=0 and z>=0:
                # print(x, y, z)
                Positions_elim[i] = [0, 0, 0]
            i+=1
    
        Positions_elim = np.squeeze(Positions_elim)
        Positions = Positions_elim[~np.all(Positions_elim==0., axis=1)].T
        

    else:
        Positions = np.array([x, y, z])

    
    if heatmap==True:
        R2, P2, T2 = crt2sph(Positions[0], Positions[1], Positions[2])
        Y2 = np.abs(R2 * R(n, l, R2) * Y(l, m, T2, P2))**2
        color = Y2
        ColorMap = colormap

    else:
        ColorMap=None

    if plysave==True:
        # Pass xyz to Open3D.o3d.geometry.PointCloud and visualize
        pcd = o3d.geometry.PointCloud()
        pcd.points = o3d.utility.Vector3dVector(Positions.T)
        
        if heatmap==True:
            cmap = mpl.colormaps['inferno']
            pc_rgb = cmap(colors.Normalize(min(Y2), max(Y2))(Y2))[:,:3]
            pcd.colors = o3d.utility.Vector3dVector(pc_rgb)

        o3d.io.write_point_cloud(f"Electron_Positions/Electron_Positions_{n,l,m}.ply", pcd)
        
    if plot==True:
        fig2 = plt.figure() 
        ax2 = fig2.add_subplot(projection='3d')
        # if heatmap == True:
        #     kwargs = (c=np.abs(R2 * R(n, l, R2) * Y(l, m, T2, P2))**2, cmap='plasma')
        ax2.scatter(Positions[0], Positions[1], Positions[2], s=1, alpha=.5, c=color, cmap=ColorMap)
        ax2.set_xlabel('x axis')
        ax2.set_ylabel('y axis')
        ax2.set_zlabel('z axis')