'''
Calculates cubs (density and total ESP) and CP of total ECP for given wavefunction data (supports .wfx/.wfn/.fchk)
Then drawing them using Mayavi code

Input information:
python3.8 VisMap3.0.py *file with wfn data* -mode=new/old -vis=y/n
        -mode = old - tries to find generated WFNfile_ESP.cub WFNfile_Dens.cub and WFNfile_ESPCP.txt and use them.
                      if files cannot be found, offers user to generate them
              = new - (re)generates all cub and CP files
        -vis = y - launches Mayavi
             = n - doesn't launch visualization

'''


import os
import re
import sys
import numpy as np
import copy

Multiwfnpath='Multiwfn'

# dictionary of (nuclear_charge : Nucleus, vdW radii in Angstroms)
dnc2all = {1: ['H', 1.09, 0.23, 1.00794, 0.99609375, 0.9765625, 0.80078125],
           2: ['He', 1.40, 1.50, 4.002602, 0.72265625, 0.82421875, 0.9296875],
           3: ['Li', 1.82, 1.28, 6.941, 0.7421875, 0.7421875, 0.7421875],
           4: ['Be', 2.00, 0.96, 9.012182, 0.7421875, 0.7421875, 0.7421875],
           5: ['B', 2.00, 0.83, 10.811, 0.625, 0.234375, 0.234375],
           6: ['C', 1.70, 0.68, 12.0107, 0.328125, 0.328125, 0.328125],
           7: ['N', 1.55, 0.68, 14.0067, 0.1171875, 0.5625, 0.99609375],
           8: ['O', 1.52, 0.68, 15.9994, 0.99609375, 0.0, 0.0],
           9: ['F', 1.47, 0.64, 18.998403, 0.99609375, 0.99609375, 0.0],
           10: ['Ne', 1.54, 1.50, 20.1797, 0.72265625, 0.82421875, 0.9296875],
           11: ['Na', 2.27, 1.66, 22.98977, 0.7421875, 0.7421875, 0.7421875],
           12: ['Mg', 1.73, 1.41, 24.305, 0.7421875, 0.7421875, 0.7421875],
           13: ['Al', 2.00, 1.21, 26.981538, 0.7421875, 0.7421875, 0.7421875],
           14: ['Si', 2.10, 1.20, 28.0855, 0.82421875, 0.82421875, 0.82421875],
           15: ['P', 1.80, 1.05, 30.973761, 0.99609375, 0.546875, 0.0],
           16: ['S', 1.80, 1.02, 32.065, 0.99609375, 0.9609375, 0.55859375],
           17: ['Cl', 1.75, 0.99, 35.453, 0.0, 0.99609375, 0.0],
           18: ['Ar', 1.88, 1.51, 39.948, 0.72265625, 0.82421875, 0.9296875],
           19: ['K', 2.75, 2.03, 39.0983, 0.7421875, 0.7421875, 0.7421875],
           20: ['Ca', 2.00, 1.76, 40.078, 0.7421875, 0.7421875, 0.7421875],
           21: ['Sc', 2.00, 1.70, 44.95591, 0.7421875, 0.7421875, 0.7421875],
           22: ['Ti', 2.00, 1.60, 47.867, 0.7421875, 0.7421875, 0.7421875],
           23: ['V', 2.00, 1.53, 50.9415, 0.7421875, 0.7421875, 0.7421875],
           24: ['Cr', 2.00, 1.39, 51.9961, 0.7421875, 0.7421875, 0.7421875],
           25: ['Mn', 2.50, 1.61, 54.938049, 0.7421875, 0.7421875, 0.7421875],
           26: ['Fe', 2.00, 1.52, 55.845, 0.7421875, 0.7421875, 0.7421875],
           27: ['Co', 2.00, 1.26, 58.9332, 0.7421875, 0.7421875, 0.7421875],
           28: ['Ni', 1.63, 1.24, 58.6934, 0.7421875, 0.7421875, 0.7421875],
           29: ['Cu', 1.40, 1.32, 63.546, 0.99609375, 0.5078125, 0.27734375],
           30: ['Zn', 1.39, 1.22, 65.409, 0.7421875, 0.7421875, 0.7421875],
           31: ['Ga', 1.87, 1.22, 69.723, 0.7421875, 0.7421875, 0.7421875],
           32: ['Ge', 2.00, 1.17, 72.64, 0.7421875, 0.7421875, 0.7421875],
           33: ['As', 1.85, 1.21, 74.9216, 0.7421875, 0.7421875, 0.7421875],
           34: ['Se', 1.90, 1.22, 78.96, 0.7421875, 0.7421875, 0.7421875],
           35: ['Br', 1.85, 1.21, 79.904, 0.7421875, 0.5078125, 0.234375],
           36: ['Kr', 2.02, 1.50, 83.798, 0.72265625, 0.82421875, 0.9296875],
           37: ['Rb', 2.00, 2.20, 85.4678, 0.7421875, 0.7421875, 0.7421875],
           38: ['Sr', 2.00, 1.95, 87.62, 0.7421875, 0.7421875, 0.7421875],
           39: ['Y', 2.00, 1.90, 88.90585, 0.7421875, 0.7421875, 0.7421875],
           40: ['Zr', 2.00, 1.75, 91.224, 0.7421875, 0.7421875, 0.7421875],
           41: ['Nb', 2.00, 1.64, 92.90638, 0.7421875, 0.7421875, 0.7421875],
           42: ['Mo', 2.00, 1.54, 95.94, 0.7421875, 0.7421875, 0.7421875],
           43: ['Tc', 2.00, 1.47, 98.0, 0.7421875, 0.7421875, 0.7421875],
           44: ['Ru', 2.00, 1.46, 101.07, 0.7421875, 0.7421875, 0.7421875],
           45: ['Rh', 2.00, 1.45, 102.9055, 0.7421875, 0.7421875, 0.7421875],
           46: ['Pd', 1.63, 1.39, 106.42, 0.7421875, 0.7421875, 0.7421875],
           47: ['Ag', 1.72, 1.45, 107.8682, 0.99609375, 0.99609375, 0.99609375],
           48: ['Cd', 1.58, 1.44, 112.411, 0.7421875, 0.7421875, 0.7421875],
           49: ['In', 1.93, 1.42, 114.818, 0.7421875, 0.7421875, 0.7421875],
           50: ['Sn', 2.17, 1.39, 118.71, 0.7421875, 0.7421875, 0.7421875],
           51: ['Sb', 2.00, 1.39, 121.76, 0.7421875, 0.7421875, 0.7421875],
           52: ['Te', 2.06, 1.47, 127.6, 0.7421875, 0.7421875, 0.7421875],
           53: ['I', 2.58, 1.40, 126.90447, 0.625, 0.125, 0.9375],
           54: ['Xe', 2.16, 1.50, 131.293, 0.72265625, 0.82421875, 0.9296875],
           55: ['Cs', 2.00, 2.44, 132.90545, 0.7421875, 0.7421875, 0.7421875],
           56: ['Ba', 2.00, 2.15, 137.327, 0.7421875, 0.7421875, 0.7421875],
           57: ['La', 2.00, 2.07, 138.9055, 0.7421875, 0.7421875, 0.7421875],
           58: ['Ce', 2.00, 2.04, 140.116, 0.7421875, 0.7421875, 0.7421875],
           59: ['Pr', 2.00, 2.03, 140.90765, 0.7421875, 0.7421875, 0.7421875],
           60: ['Nd', 2.00, 2.01, 144.24, 0.7421875, 0.7421875, 0.7421875],
           61: ['Pm', 2.00, 1.99, 145.0, 0.7421875, 0.7421875, 0.7421875],
           62: ['Sm', 2.00, 1.98, 150.36, 0.7421875, 0.7421875, 0.7421875],
           63: ['Eu', 2.00, 1.98, 151.964, 0.7421875, 0.7421875, 0.7421875],
           64: ['Gd', 2.00, 1.96, 157.25, 0.7421875, 0.7421875, 0.7421875],
           65: ['Tb', 2.00, 1.94, 158.92534, 0.7421875, 0.7421875, 0.7421875],
           66: ['Dy', 2.00, 1.92, 162.5, 0.7421875, 0.7421875, 0.7421875],
           67: ['Ho', 2.00, 1.92, 164.93032, 0.7421875, 0.7421875, 0.7421875],
           68: ['Er', 2.00, 1.89, 167.259, 0.7421875, 0.7421875, 0.7421875],
           69: ['Tm', 2.00, 1.90, 168.93421, 0.7421875, 0.7421875, 0.7421875],
           70: ['Yb', 2.00, 1.87, 173.04, 0.7421875, 0.7421875, 0.7421875],
           71: ['Lu', 2.00, 1.87, 174.967, 0.7421875, 0.7421875, 0.7421875],
           72: ['Hf', 2.00, 1.75, 178.49, 0.7421875, 0.7421875, 0.7421875],
           73: ['Ta', 2.00, 1.70, 180.9479, 0.7421875, 0.7421875, 0.7421875],
           74: ['W', 2.00, 1.62, 183.84, 0.7421875, 0.7421875, 0.7421875],
           75: ['Re', 2.00, 1.51, 186.207, 0.7421875, 0.7421875, 0.7421875],
           76: ['Os', 2.00, 1.44, 190.23, 0.7421875, 0.7421875, 0.7421875],
           77: ['Ir', 2.00, 1.41, 192.217, 0.7421875, 0.7421875, 0.7421875],
           78: ['Pt', 1.72, 1.36, 195.078, 0.7421875, 0.7421875, 0.7421875],
           79: ['Au', 1.66, 1.50, 196.96655, 0.99609375, 0.83984375, 0.0],
           80: ['Hg', 1.55, 1.32, 200.59, 0.7421875, 0.7421875, 0.7421875],
           81: ['Tl', 1.96, 1.45, 204.3833, 0.7421875, 0.7421875, 0.7421875],
           82: ['Pb', 2.02, 1.46, 207.2, 0.7421875, 0.7421875, 0.7421875],
           83: ['Bi', 2.00, 1.48, 208.98038, 0.7421875, 0.7421875, 0.7421875],
           84: ['Po', 2.00, 1.40, 290.0, 0.7421875, 0.7421875, 0.7421875],
           85: ['At', 2.00, 1.21, 210.0, 0.7421875, 0.7421875, 0.7421875],
           86: ['Rn', 2.00, 1.50, 222.0, 0.72265625, 0.82421875, 0.9296875],
           87: ['Fr', 2.00, 2.60, 223.0, 0.7421875, 0.7421875, 0.7421875],
           88: ['Ra', 2.00, 2.21, 226.0, 0.7421875, 0.7421875, 0.7421875],
           89: ['Ac', 2.00, 2.15, 227.0, 0.7421875, 0.7421875, 0.7421875],
           90: ['Th', 2.00, 2.06, 232.0381, 0.7421875, 0.7421875, 0.7421875],
           91: ['Pa', 2.00, 2.00, 231.03588, 0.7421875, 0.7421875, 0.7421875],
           92: ['U', 1.86, 1.96, 238.02891, 0.7421875, 0.7421875, 0.7421875],
           93: ['Np', 2.00, 1.90, 237.0, 0.7421875, 0.7421875, 0.7421875],
           94: ['Pu', 2.00, 1.87, 244.0, 0.7421875, 0.7421875, 0.7421875],
           95: ['Am', 2.00, 1.80, 243.0, 0.7421875, 0.7421875, 0.7421875],
           96: ['Cm', 2.00, 1.69, 247.0, 0.7421875, 0.7421875, 0.7421875],
           97: ['Bk', 2.00, 1.54, 247.0, 0.7421875, 0.7421875, 0.7421875],
           98: ['Cf', 2.00, 1.83, 251.0, 0.7421875, 0.7421875, 0.7421875],
           99: ['Es', 2.00, 1.50, 252.0, 0.7421875, 0.7421875, 0.7421875],
           100: ['Fm', 2.00, 1.50, 257.0, 0.7421875, 0.7421875, 0.7421875],
           101: ['Md', 2.00, 1.50, 258.0, 0.7421875, 0.7421875, 0.7421875],
           102: ['No', 2.00, 1.50, 259.0, 0.7421875, 0.7421875, 0.7421875],
           103: ['Lr', 2.00, 1.50, 262.0, 0.7421875, 0.7421875, 0.7421875]
           }

dcpCol={'(3,-3)':(1,1,1),
        '(3,-1)':(0,1,0),
        '(3,+1)':(1,0,0),
        '(3,+3)':(0,0,1),
        }


#### Run Multiwfn
def Run_MWFN(mytext, needout=False):
    with open('myprog.inp', 'w') as inp:
        inp.write('\n'.join(mytext) + '\n')
    if needout == True:
        os.system('{0} {1} < {2} > {3}'.format(Multiwfnpath, inputfile, 'myprog.inp', 'myprog.out'))
    else:
        os.system('{0} {1} < {2}'.format(Multiwfnpath, inputfile, 'myprog.inp'))
    try:
        os.remove('myprog.inp')
    except FileNotFoundError:
        pass


def ReadCUB(inputfile):
    CENTERS=[]
    Scalars=[]
    with open(inputfile) as cube:
        lines=cube.readlines()[2:]
        nat=int(lines[0].split()[0])
        pointsv1 = int(lines[1].split()[0])
        pointsv2 = int(lines[2].split()[0])
        pointsv3 = int(lines[3].split()[0])

        origin=np.array([float(x) for x in lines[0].split()[-3:]])
        v1= np.array([float(x) for x in lines[1].split()[1:]])
        v2 = np.array([float(x) for x in lines[2].split()[1:]])
        v3 = np.array([float(x) for x in lines[3].split()[1:]])

        for i in range(nat):
            line=lines[4+i].split()
            CENTERS.append([int(line[0])]+[float(x) for x in line[2:]]+[dnc2all[int(line[0])][0] + str(i + 1)])
        for line in lines[4+nat:]:
            Scalars += [float(x) for x in line.split()]
        print('Found',len(Scalars), 'Scalars', 'in ', inputfile)
    XYZS_Data=[origin,v1,v2,v3,pointsv1,pointsv2,pointsv3,Scalars]
    return XYZS_Data,CENTERS


def CalcCub(IsCalced,fname):
    if IsCalced[0]==False:
        print('Calcualting Density cube')
        if IsCalced[1]:
            text = ['1000', '10', nproc, '5', '1', '8', fname+'_ESP.cub','2']
            Run_MWFN(text,False)
        else:
            text = ['1000', '10', nproc, '5', '1', '2', '2']
            Run_MWFN(text, False)
        os.rename('density.cub',fname + '_Dens.cub')
        IsCalced[0]=True
    if IsCalced[1]==False:
        print('Calcualting ESP cube')
        if IsCalced[0]:
            text = ['1000', '10', nproc, '5', '12', '8', fname+'_Dens.cub','2']
            Run_MWFN(text,False)
        else:
            text = ['1000', '10', nproc, '5', '12', '2', '2']
            Run_MWFN(text, False)
        os.rename('totesp.cub',fname + '_ESP.cub')
        IsCalced[1] = True


def CalcPoints(isoval):
    print('Searching surfanalysis.txt file with points')
    if not os.path.exists(fname+'_sa_' + str(isoval) + '.txt'):
        print('File not found. Calling Multiwfn for min/max locating')
        text = ['1000', '10', nproc, '12', '1', '1', str(isoval), '0', '1']
        Run_MWFN(text, False)
        os.rename('surfanalysis.txt', fname+'_sa_' + str(isoval) + '.txt')

    MAXMIN = []
    with open(fname + '_sa_' + str(isoval) + '.txt') as out:
        for line in out:
            line = line.split()
            if len(line) > 5:
                if '*' not in line and 'eV' not in line:
                    MAXMIN.append([float(line[3])] + [float(x) / 0.529 for x in line[4:]])
                elif '*' in line and 'eV' not in line:
                    MAXMIN.append([float(line[4])] + [float(x) / 0.529 for x in line[5:]])
    print('Located', len(MAXMIN), 'extremum points on the surface')
    return MAXMIN


def VisualizeData(CENTERS,CUBdat,CUBdatESP,xx,yy,zz):
    from mayavi import mlab
    # start mayavi
    mlab.figure(1, fgcolor=(1, 1, 1), bgcolor=(0, 0, 0))

    # Draw Atoms
    for Atom in CENTERS:
        color = tuple(dnc2all[Atom[0]][-3:])
        mlab.points3d(Atom[1], Atom[2], Atom[3], color=color, scale_factor=0.8, resolution=20, scale_mode='none')
        #mlab.text3d(Atom[1], Atom[2], Atom[3], Atom[-1], scale=(0.8,0.8,0.8))
        for Atom2 in CENTERS:
            dist = np.linalg.norm(np.array(Atom[1:-1]) - np.array(Atom2[1:-1]))
            if dist < (dnc2all[Atom[0]][1] + dnc2all[Atom2[0]][1]) and dist != 0:
                bond = np.stack([Atom[1:-1], Atom2[1:-1]], axis=-1)
                mlab.plot3d(bond[0], bond[1], bond[2])

    # Addind the scalar field
    src = mlab.pipeline.scalar_field(xx, yy, zz, CUBdat)

    # Addind new point_data to old scalar field
    src.image_data.point_data.add_array(CUBdatESP.T.ravel())
    src.image_data.point_data.get_array(1).name = 'ESP'
    src.update()

    # Selecting the first active attribute - scalar field and making contouros
    src2 = mlab.pipeline.set_active_attribute(src, point_scalars='scalar')
    contour = mlab.pipeline.contour(src2)
    contour.filter.contours = [0.001]
    isoval=0.001

    # Selecting the second attribute - ESP
    contour2 = mlab.pipeline.set_active_attribute(contour, point_scalars='ESP')
    surface=mlab.pipeline.surface(contour2, colormap='gist_rainbow', opacity=0.5)

    mlab.colorbar(object=surface, title='ESP, kkal/mol', orientation='vertical', nb_labels=3)

    while True:
        TODO=input('Choose other isosurface or generate min/max on the current plane?(isosurf value/gen/scan/kill/exit)\t')

        if TODO=='gen':
            MAXMIN=CalcPoints(isoval)

            if 'points' in locals():
                for i in range(len(points)):
                    try:
                        points[i].remove()
                        pointstext[i].remove()
                    except ValueError:
                        pass


            points=[]
            pointstext=[]
            for point in MAXMIN:
                points.append(mlab.points3d(point[1],point[2],point[3],color=(1,1,1), scale_factor=0.1, resolution=8, scale_mode='none'))
                pointstext.append(mlab.text3d(point[1],point[2],point[3],'{:10.2f}'.format(point[0]),scale=0.5))


        elif TODO=='scan':
            MAXMIN = CalcPoints(isoval)

            if 'points' in locals():
                for i in range(len(points)):
                    try:
                        points[i].remove()
                        pointstext[i].remove()
                    except ValueError:
                        pass
            points = []

            CPsnearUA = []
            for atom in CENTERS:
                if atom[0] in [1,8,9,16,17,34,35,52,53,84,85]:
                    for CP in MAXMIN:
                        if np.linalg.norm(np.array(CP[1:]) - atom[1:-1]) < 5.0 and CP not in CPsnearUA:
                            CPsnearUA.append(CP)
            MAXMIN=CPsnearUA

            points = []
            pointstext = []
            for point in MAXMIN:
                points.append(
                    mlab.points3d(point[1], point[2], point[3], color=(1, 1, 1), scale_factor=0.1, resolution=8,
                                  scale_mode='none'))
                pointstext.append(mlab.text3d(point[1], point[2], point[3], '{:10.2f}'.format(point[0]), scale=0.5))


        elif 'kill' in TODO:
            try:
                killval = float(TODO.split()[1])
                killpm = float(TODO.split()[2])
            except ValueError:
                print("Could't convert to float given number. Please, try again.")
                pass

            if 'points' in locals() and 'MAXMIN' in locals():
                for i in range(len(points)):
                    if MAXMIN[i][0]<killval+killpm and MAXMIN[i][0]>killval-killpm:
                        try:
                            print('removed', i, MAXMIN[i])
                            points[i].remove()
                            pointstext[i].remove()
                        except ValueError:
                            pass

        elif TODO=='exit':
            break

        else:
            if 'points' in locals():
                for i in range(len(points)):
                    try:
                        points[i].remove()
                        pointstext[i].remove()
                    except ValueError:
                        pass

            try:
                isoval=float(TODO)
                contour.filter.contours = [isoval]
            except ValueError:
                print("Could't convert to float given number. Please, try again.")
                pass

    mlab.show()


##### Main Code #####
### Input part
exe_path = os.path.abspath(os.path.dirname(sys.argv[0]))
systeminput = sys.argv
inputfile = systeminput[1]
s=inputfile.rfind('.')
fname=inputfile[:s]

for i in systeminput:
    if '-mode=' in i:
        mode = i[6:]
    if '-vis=' in i:
        vis = i[5]
    if '-nproc=' in i:
        nproc = i[7:]


if 'mode' not in globals():
    mode = 'old'
if 'vis' not in globals():
    vis = 'y'
if 'nproc' not in globals():
    nproc = '4'

if mode not in ['new','old']:
    mode = 'old'
if vis not in ['y','n']:
    vis = 'y'
try:
   nproc=str(int(nproc))
except ValueError:
   print('Could not convert given nproc to integer! Using 4 proc')
   nproc = '4'
   pass


if mode == 'old':
    IsCalced = [True, True]
    for i,x in enumerate([fname + '_Dens.cub',fname + '_ESP.cub']):
        if x in os.listdir(path='.'):
            print('Located ', x)
        else:
            IsCalced[i]=False
            print('Could not locate', x, " - It'll be (re)calculated")
elif mode == 'new':
    IsCalced = [False, False]

### Multiwfn part
CalcCub(IsCalced,fname)

### Visualization part
# reading all data, generating mashgrid, scalar matrices
CUB, CENTERS = ReadCUB(fname+'_Dens.cub')

origin, v1, v2, v3, pointsv1, pointsv2, pointsv3, Scalars = CUB
CUBdat=np.empty((pointsv1,pointsv2,pointsv3))
for x in range(pointsv1):
    for y in range(pointsv2):
        for z in range(pointsv3):
            CUBdat[x][y][z]=Scalars[x*pointsv2*pointsv3+y*pointsv3+z]
xx,yy,zz = np.mgrid[origin[0]:origin[0]+v1[0]*pointsv1:v1[0], origin[1]:origin[1]+v2[1]*pointsv2:v2[1], origin[2]:origin[2]+v3[2]*pointsv3:v3[2]]


CUB, CENTERS = ReadCUB(fname + '_ESP.cub')
TotESP = CUB[-1]
CUBdatESP = np.empty((pointsv1, pointsv2, pointsv3))
for x in range(pointsv1):
    for y in range(pointsv2):
        for z in range(pointsv3):
            CUBdatESP[x][y][z] = TotESP[x * pointsv2 * pointsv3 + y * pointsv3 + z]*627.0

# vis
if vis=='y':
    VisualizeData(CENTERS, CUBdat, CUBdatESP, xx, yy, zz)

