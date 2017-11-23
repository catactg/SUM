import os
import numpy as np
from basefunctions import *
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy
import seedselector
import csv


def compute_rmse(mesh_file1, mesh_file2, arrayname):
    """Compute root mean square error."""
    # read sum disk
    mesh1 = readvtp(mesh_file1)
    mesh2 = readvtp(mesh_file2)

    # compute RMSD
    # extract arrays
    array1 = vtk_to_numpy(mesh1.GetPointData().GetArray(arrayname))
    array2 = vtk_to_numpy(mesh2.GetPointData().GetArray(arrayname))
    # becuase there maybe nans in lat disks
    diffsquare = []
    for x in range(len(array1)):
        # they are both same type (ie either number or nan)
        if np.isnan(array1[x]) ==  np.isnan(array2[x]):
            # if not nan, compute rmse
            if not np.isnan(array1[x]):
                diffsquare.append((array1[x] - array2[x]) ** 2)

    rmse = np.sqrt( np.array(diffsquare) / float(len(diffsquare)) )
    # check what happens to nans
    return np.mean(rmse)

def get_mean_array_per_region(polydata, inarrayname):
    """Extract each region computate mean value per region"""
    labels = [  'ant','lat','laa_around','roof','post','isthmus','floor','septum',
                'lpv_sup_q1','lpv_sup_q2','lpv_sup_q3','lpv_sup_q4',
                'lpv_inf_q1','lpv_inf_q2','lpv_inf_q3','lpv_inf_q4',
                'rpv_sup_q1','rpv_sup_q2','rpv_sup_q3','rpv_sup_q4',
                'rpv_inf_q1','rpv_inf_q2','rpv_inf_q3','rpv_inf_q4']

    flatlabels = getflatregionslabels()

    outarray=[]
    for label in labels:
        # threshold each region
        region = cellthreshold(polydata,'sumlabels',flatlabels[label],flatlabels[label])

        # region should already has a single value
        array = vtk_to_numpy(region.GetCellData().GetArray(inarrayname))

        per = np.unique(array)
        per_nonan = per[ ~np.isnan(per)]
        if len( per_nonan ) > 1:
            print "WARNING: found multiple", inarrayname
            print "Region",label,"has mean value of",per_nonan
        else:
            print "Region",label,"has mean value of",per_nonan
            outarray.append(per[0])

    return outarray


def writesumarrays2csv(array1, array2, array3,
                       array1label, array2label, array3label,
                       ofile):
    """Write sum arrays to csv."""

    labels = [  'ant','lat','laa_around','roof','post','isthmus','floor','septum',
        'lpv_sup_q1','lpv_sup_q2','lpv_sup_q3','lpv_sup_q4',
        'lpv_inf_q1','lpv_inf_q2','lpv_inf_q3','lpv_inf_q4',
        'rpv_sup_q1','rpv_sup_q2','rpv_sup_q3','rpv_sup_q4',
        'rpv_inf_q1','rpv_inf_q2','rpv_inf_q3','rpv_inf_q4']

    flatlabels = getflatregionslabels()

    f = open(ofile, 'wb')
    # save headers
    line = 'regionid,sumlabel,' + array1label + ',' + array2label + ',' + array3label +'\n'
    f.write(line)
    for label,a1,a2,a3 in zip(labels,array1,array2,array3):
        line = str(flatlabels[label]) + ', ' + label + ', '
        line = line + str(a1) + ', ' + str(a2) + ', ' + str(a3)+'\n'
        f.write(line)
    f.close()

def mean_value_per_region(polydata, inarrayname, outarrayname):
    """Extract each region computate mean value per region"""
    labels = [  'ant','lat','laa_around','roof','post','isthmus','floor','septum',
                'lpv_sup_q1','lpv_sup_q2','lpv_sup_q3','lpv_sup_q4',
                'lpv_inf_q1','lpv_inf_q2','lpv_inf_q3','lpv_inf_q4',
                'rpv_sup_q1','rpv_sup_q2','rpv_sup_q3','rpv_sup_q4',
                'rpv_inf_q1','rpv_inf_q2','rpv_inf_q3','rpv_inf_q4']

    flatlabels = getflatregionslabels()

    for label in labels:
        # extract each region by thresholding label
        region = cellthreshold(polydata,'sumlabels',flatlabels[label],flatlabels[label])

        # compute mean value discarding zeros
        per = compute_nonzeromean(region, inarrayname)
        print "Region",label,"has mean value of",per

        polydata = set_cell_array_value_per_label(polydata,outarrayname,per,
            'sumlabels',flatlabels[label])

    return polydata


def area_coverage(polydata, inarrayname, outarrayname):
    """For the rest of the regions, extract each region computate area threshold / area total"""
    labels = ['ant','lat','laa_around','isthmus','floor','septum']
    # extract labels
    flatlabels = getflatregionslabels()

    for label in labels:
        # threshold each region
        region = cellthreshold(polydata,'sumlabels',flatlabels[label],flatlabels[label])
        region_t = triangulate(region) # for mass properties?
        area_total = surfacearea(region_t)

        # threshold each region
        subtarget = pointthreshold(region_t,inarrayname,1,1)
        if subtarget.GetNumberOfPoints()>0:
            area_subtarget = surfacearea(subtarget)
        else:
            area_subtarget = 0.
        # add value to array
        per = 100 * area_subtarget / area_total
        print "area_coverage",flatlabels[label],per
        polydata = set_cell_array_value_per_label(polydata,outarrayname,per,
            'sumlabels',flatlabels[label])
    return polydata

def circumferential_coverage(polydata, inarrayname, outarrayname):
    """From template generation we know the coordinates of the PV centers and radius"""
    centers =   [[-0.3, 0.10, 0],
                [-0.3, -0.15, 0],
                [0.25, 0.10, 0],
                [0.25, -0.15,0]]
    rad = 1 # far far away
    ref_axis = [0, 0, 1]
    # becuase of zflip, quadrant numbers are opposite to cartesian quadrants
    start_angles = [90, 0, 270, 180]

    # extract each PV's quadrants
    flatlabels = getflatregionslabels()

    # 4 veins
    pv_labels = ['lpv_sup','lpv_inf','rpv_sup','rpv_inf']
    for ind in range(0,len(pv_labels)):
        pv_label = pv_labels[ind]
        for q in range(1,5):
            # threshold each quadrant
            quadrant = cellthreshold(polydata,'sumlabels',
                flatlabels[pv_label + '_q' + str(q)],flatlabels[pv_label + '_q' + str(q)])

            n_total = 0
            n_count = 0
            center = centers[ind]
            start_angle = start_angles[q-1]
            tetha_array = np.linspace(0,90,90) # 0.5 degrees
            for tetha in tetha_array:
                # change the normal by a step
                x, y = polar2cart(rad,tetha + start_angle ,center)
                vector = subtractvectors(center,[x,y,0])
                vectorn = normalizevector(vector)
                normal = cross(vectorn,ref_axis)
                # make cuts
                ray = cutdataset(quadrant,center,normal)
                if ray.GetPointData().GetArray(inarrayname).GetNumberOfTuples()>0:
                    rayarray = vtk_to_numpy(ray.GetPointData().GetArray(inarrayname))

                    if np.sum(rayarray) > 0:
                        n_count += 1
                    if ray.GetNumberOfPoints() > 0:
                        n_total += 1
            # % coverage
            per = 100 * n_count / n_total
            print "circumferential_coverage",pv_label + '_q' + str(q),per
            polydata = set_cell_array_value_per_label(polydata,outarrayname,per,
                'sumlabels',flatlabels[pv_label + '_q' + str(q)])
    return polydata


def longitudinal_coverage(polydata, inarrayname, outarrayname):
    """From template generation we know the coordinates of posterior wall."""
    startpoints =   [[-0.12,-0.025,0],
                    [-0.236,0.2,0]]
    endpoints   =   [[0.072,-0.025,0],
                    [0.252,0.2,0]]
    normal = [1, 0, 0]
    step = 0.0001
    # extract roof and posterior wall
    flatlabels = getflatregionslabels()
    labels = ['post','roof']

    for ind in range(0,len(labels)):
        subtarget = cellthreshold(polydata,'sumlabels',
            flatlabels[labels[ind]],flatlabels[labels[ind]])

        # change the clippoint by a step
        currentpoint = startpoints[ind]

        n_total = 0
        n_count = 0
        x_array = np.linspace(currentpoint[0],endpoints[ind][0],100)
        for x in x_array:
            currentpoint[0] = x  # only on X-axis
            # make cuts
            ray = cutdataset(subtarget,currentpoint,normal)
            if ray.GetPointData().GetArray(inarrayname).GetNumberOfTuples()>0:
                rayarray = vtk_to_numpy(ray.GetPointData().GetArray(inarrayname))

                if np.sum(rayarray) > 0:
                    n_count += 1
                if ray.GetNumberOfPoints() > 0:
                    n_total += 1
        # % coverage
        per = 100 * n_count / n_total
        print "longitudinal_coverage",per
        polydata = set_cell_array_value_per_label(polydata,outarrayname,per,
            'sumlabels',flatlabels[labels[ind]])
    return polydata


def histogram_normalisation(polydata, arrayname):
    """Histogram normalisation. lge array to normalised histogram.
    Trying to reproduce approach of Udupa.
    Histogram normalisation targets, taken from C. Tobon-Gomez le simulation paper (EMBC)."""

    # extract array
    origarray = vtk_to_numpy(polydata.GetPointData().GetArray(arrayname))

    # make new array and initialise
    newarray = vtk.vtkDoubleArray()
    newarray.SetName(arrayname + '_norm')
    newarray.SetNumberOfTuples(polydata.GetNumberOfPoints())
    newarray.FillComponent(0, 0)

    # compute normalisation values
    s1 = min(origarray)
    s2 = max(origarray)
    pc1 = np.percentile(origarray,0.01)
    pc2 = np.percentile(origarray,99)

    if pc1>0 :
        f1_m = (0.00002)/(pc1-s1)
    else:
        f1_m = (0.00002)

    # interpolation functions
    f1_b = 0.00002-(f1_m*pc1)

    f5_m = ( 0.8-0.00002)/(pc2-pc1)
    f5_b = 0.8-(f5_m*pc2)

    f6_m = (1-0.8)/(s2-pc2)
    f6_b = 1-(f6_m*s2)

    # visit each value and remap
    for ind in range(len(origarray)):
        val=origarray[ind]

        if val<0:
            newval = 0
        else:
            if val<pc1:
                if pc1>0:
                    newval =(f1_m*val)+f1_b
                else:
                    newval =0
            else:
                if val<pc2:
                    newval =(f5_m*val)+f5_b
                else:
                    newval =(f6_m*val)+f6_b
        newarray.SetValue(ind,newval)

    polydata.GetPointData().AddArray(newarray)
    return polydata

def add_threshold_array(surface, arrayname, th_val):
    """Add array in which all values in  arrayname > th_val are set to 1"""
    # extract array
    array = surface.GetPointData().GetArray(arrayname)

    # make new array and initialise
    newarray = vtk.vtkDoubleArray()
    newarray.SetName(arrayname + '_th')
    newarray.SetNumberOfTuples(surface.GetNumberOfPoints())
    newarray.FillComponent(0, 0)

    # set to one if above threshold
    for p in range(surface.GetNumberOfPoints()):
        if array.GetValue(p) >= th_val:
            newarray.SetValue(p,1)

    surface.GetPointData().AddArray(newarray)
    return surface

def vein_clip_with_seeds(surface, surfacefileout):
    """Choosing 3 seeds for veing clipping to fix complicated meshes
    (i.e. with LAA- LSPV joint).
    Choose three seeds: a source, a target and the position of clipping.
    Suggestion: pick a source between right PVs. A target in the LAA.
    and a distance point in the LSPV (as far as possible from the LAA to get the
        distance of the clipping cylinder right.
        If it is better for the centerline to go on the LSPV,
        then click the target on the LSPV and the distance seed on the LAA."""

    nseeds = 3
    seeds = seed_interactor(surface)

    if not seeds.GetNumberOfIds() == nseeds:
        print 'You should select extactly',nseeds,' seeds. Try again!'
        seeds = seed_interactor(surface)


    # get point coordinates
    newpoints = vtk.vtkPoints()
    newvertices = vtk.vtkCellArray()
    points = []
    for s in range(seeds.GetNumberOfIds()):
        point = surface.GetPoint(seeds.GetId(s))
        points.append(point)
        pid = newpoints.InsertNextPoint(point)
        newvertices.InsertNextCell(1)
        newvertices.InsertCellPoint(pid)
    pointspd = vtk.vtkPolyData()
    pointspd.SetPoints(newpoints)
    pointspd.SetVerts(newvertices)
    writevtp(pointspd,surfacefileout)

    # centerline with source and target
    cl = vmtkcenterlines(surface, points[0], points[1])
    cl = vmtkcenterlinesmoothing(cl)

    # then find closest point to third seed
    # pointlocator
    locator = vtk.vtkPointLocator()
    locator.SetDataSet(cl)
    locator.BuildLocator()

    clippointid = locator.FindClosestPoint(points[2])

    # compute normal (inverted)
    closest_point = cl.GetPoint(clippointid)
    clipnormal = ( np.array(cl.GetPoint(clippointid+1)) -
                np.array(cl.GetPoint(clippointid)) )

    # find second point "upstream"
    vectorup = multiplyvector(clipnormal,5) # 5mm
    pointup = acumvectors(vectorup,closest_point)

    # radius = twice distance from seed point and centerline
    radius = 1.75 * euclideandistance(closest_point,points[2])
    body = cylinderclip(surface, closest_point, pointup, clipnormal, radius)

    surfacefilled = vmtksurfacecapper(body,nrings=4)

    return surfacefilled

def manual_widget_clip(mesh, runclipper=1):
    """Manual clipping to fix complicated meshes (i.e. with LAA- LSPV joint) """
    if runclipper ==1:
        # sphere widget doesn't work well
        clippled_surface= vmtksurfaceclipper(mesh)
    else:
        clippled_surface = mesh
    largest_surface = extractlargestregion(clippled_surface)
    # close clipped surface
    surfacefilled = vmtksurfacecapper(largest_surface,nrings=2)

    return surfacefilled

def select_seeds(surface, labels, surfacefileout, vis=False, avg=False):
    """Select 4 seeds, one per vein. A 5th seed may be selected as LAA.
        assuming order of seeds is consistent
        RSpv, RIpv, LIpv, LSpv
        4   1
        3   2"""

    # for avg mesh select LAA seed
    if avg:
        labelsrange = [76.0, 77.0, 79.0, 78.0, 36.0]
        nseeds= 5
    else:
        # for each PV
        labelsrange = [76.0, 77.0, 79.0, 78.0]
        nseeds =4

    seeds = seed_interactor(surface)

    # create the pointset
    newpoints = vtk.vtkPoints()
    newvertices = vtk.vtkCellArray()

    # create array on seeds with ground truth (GT) labels
    gtlabels_array = vtk.vtkDoubleArray()
    gtlabels_array.SetName(labels)

    if not seeds.GetNumberOfIds() == nseeds:
        print 'You should select extactly',nseeds,' seeds. Try again!'
        seeds = seed_interactor(surface)

    for s in range(seeds.GetNumberOfIds()):
        branchlabel = labelsrange[s]
        point = surface.GetPoint(seeds.GetId(s))
        pid = newpoints.InsertNextPoint(point)
        gtlabels_array.InsertNextValue(branchlabel)
        # Create the topology of the point (a vertex)
        newvertices.InsertNextCell(1)
        newvertices.InsertCellPoint(pid)

    pointspd = vtk.vtkPolyData()
    pointspd.SetPoints(newpoints)
    pointspd.SetVerts(newvertices)
    pointspd.GetPointData().AddArray(gtlabels_array)

    if vis:
        pointsgplyh = generateglyph(pointspd)
        visualise_default(pointsgplyh,surface,'seeds',labels,36,79)
    writevtp(pointspd,surfacefileout)

def seed_interactor(surface):
    """Interactor for seed selection."""
    computer = seedselector.vmtkPickPointSeedSelector()
    computer.SetSurface(surface)
    computer.Execute()

    return  computer.GetSourceSeedIds()

def multiple_seeds_to_csv_no_mitral(seedsfile, arrayname, labels, outfile):
    """Combine multiple seeds into a single seed per label."""

    f = open(outfile, 'wb')
    # write values of row by row
    # compute average s(s) per label and save
    allseeds = readvtp(seedsfile)

    # write values row by row
    for l in labels:

        currentseeds = pointthreshold(allseeds, arrayname,l,l,0)
        currentpoint = pointset_centreofmass(currentseeds)
        line = str(currentpoint[0]) + ',' + str(currentpoint[1]) + ',' + str(currentpoint[2]) + '\n'
        f.write(line)

    f.close()


def pv_centerlines_no_mitral(inputfile, seedsfile, outfile, pvends=1):
    """Create four pairs of centerlines, one for each vein (numbered 1 to 4)self.
    The first centerline of each pair runs from the vein's end to oppposite vein.
    The second centerline of each pair runs from the vein's end to the end
    of the nearby vein (vein 1 is near vein 2; vein 3 is near vein 4).
    Input:
    * surface mesh of atrium
    * point coordinates of veins' ends (point 0 to point 3)
    Ouput:
    * 4 pairs of centerlines"""


    surface = vmtksurfacereader(inputfile)
    points = np.loadtxt(seedsfile, delimiter=',').tolist()

    cl1 = vmtkcenterlines(surface, points[0], points[2] + points[3],pvends)
    cl2 = vmtkcenterlines(surface, points[1], points[2] + points[3],pvends)
    cl3 = vmtkcenterlines(surface, points[2], points[0] + points[1],pvends)
    cl4 = vmtkcenterlines(surface, points[3], points[0] + points[1],pvends)

    vmtksurfacewriter(cl1, outfile + 'clraw21.vtp')
    vmtksurfacewriter(cl2, outfile + 'clraw22.vtp')
    vmtksurfacewriter(cl3, outfile + 'clraw23.vtp')
    vmtksurfacewriter(cl4, outfile + 'clraw24.vtp')

def pv_centerlines_no_mitral_laa(inputfile, seedsfile, outfile, pvends=1):
    """Create four pairs of centerlines, one for each vein (numbered 1 to 5).
    The first centerline of each pair runs from the vein's end to oppposite vein.
    The second centerline of each pair runs from the vein's end to the end
    of the nearby vein (vein 1 is near vein 2; vein 3 is near vein 4).
    Input:
    * surface mesh of atrium
    * point coordinates of veins' ends (point 0 to point 4)
    Ouput:
    * 5 pairs of centerlines"""


    surface = vmtksurfacereader(inputfile)
    points = np.loadtxt(seedsfile, delimiter=',').tolist()

    cl1 = vmtkcenterlines(surface, points[0], points[2] + points[3],pvends)
    cl2 = vmtkcenterlines(surface, points[1], points[2] + points[3],pvends)
    cl3 = vmtkcenterlines(surface, points[2], points[0] + points[1],pvends)
    cl4 = vmtkcenterlines(surface, points[3], points[0] + points[1],pvends)
    cl5 = vmtkcenterlines(surface, points[4], points[0] + points[1],pvends)

    vmtksurfacewriter(cl1, outfile + 'clraw21.vtp')
    vmtksurfacewriter(cl2, outfile + 'clraw22.vtp')
    vmtksurfacewriter(cl3, outfile + 'clraw23.vtp')
    vmtksurfacewriter(cl4, outfile + 'clraw24.vtp')
    vmtksurfacewriter(cl5, outfile + 'clraw25.vtp')


def clip_veins_sections(inputfile, sufixfile, clspacing, maxslope,
                        skippointsfactor, highslope, bumpcriterion,
                        laa_seedon=0, visualise=False):
    """ We wish to clip the vein as close to the body as possible without
    including parts of the body or other veins. 'Trial' clips are
    obtained using vmtkcenterlinesections,  which creates for each point
    on the centerline a section perpendicular to the centerline and
    provides measures of the section such as maximum diameter.
    When the series of sections enter the atrium body, the maximum
    diameter increases significantly. To quantify the change in max
    diameter between one section and the next in terms of centerline
    spacing, we define 'slope'. When this slope exceeds a certain
    threshold, we assume to have entered the body. The clippoint is
    defined as the centerline point corresponding to the last section of
    the vein before entering the body. """

    surface = vmtksurfacereader(inputfile)  # atrium surface mesh

    # creating array to hold new autolabels
    if laa_seedon == 1:
        branchlabel= [0,77,76,78,79,37]
        nseeds = len(branchlabel)
    else:
        branchlabel= [0,77,76,78,79]
        nseeds = len(branchlabel)

    branch_array = vtk.vtkDoubleArray()
    branch_array.SetName('autolabels')
    branch_array.SetNumberOfTuples(surface.GetNumberOfPoints())
    surface.GetPointData().AddArray(branch_array)

    # initialize with bodylabel
    for i in range(surface.GetNumberOfPoints()):
        branch_array.SetValue(i, round(36))
    surface.GetPoints().Modified()

    for k in range(1, nseeds):
        print "Branchlabel",branchlabel[k]
        #-----------------------------------------------------------------------
        # Create vein's centerline
        #
        # Input is a pair of centerlines running from vein into the atrium body,
        # which bifurcate upon entering the atrium body. From this pair of
        # centerlines, we derive a single centerline from the vein's end until
        # the bifurcation using vmtkcenterlinemerge. Resampling and smoothing
        # operations are needed for later used vmtkcenterlinesections.
        #-----------------------------------------------------------------------

        cl = vmtksurfacereader(sufixfile + 'clraw2' + str(k) + '.vtp') # 2 means with endpoints
        cl = vmtkcenterlineresampling(cl, clspacing)
        cl = vmtkcenterlinesmoothing(cl)
        cl = vmtkbranchextractor(cl)
        vmtksurfacewriter(cl, sufixfile +'clbranch' + str(k) + '.vtp')

        #-----------------------------------------------------------------------
        # BUG FIX
        # for some anatomies, the branch extractor gives different
        # number of branches for both centerlines.
        # Calling cellthreshold seems to fix this, instead of previous
        # fix:
        # groupids = [0, 1, 2, 0, 1, 3]
        # for i in range(cl.GetNumberOfCells()):
        #     cl.GetCellData().GetArray('GroupIds').SetValue(i, groupids[i])
        #-----------------------------------------------------------------------
        cl = cellthreshold(cl,'GroupIds',0,0)


        cl = vmtkcenterlinemerge(cl)
        cl = extractcells(cl, [0])  # vein's centerline
        # original cl for clipping
        cl = vmtkcenterlineresampling(cl, clspacing)
        cl = vmtkcenterlineattributes(cl)
        vmtksurfacewriter(cl, sufixfile +'clvein' + str(k) + '.vtp')

        # removing endpoints that cause failure in vmtkcenterlinesections
        # cl only used in sections
        nskippoints = round(skippointsfactor*cl.GetNumberOfPoints())
        print "Skipping", nskippoints, "points from", cl.GetNumberOfPoints()
        cl = skippoints(cl,int(nskippoints))
        cl = vmtkcenterlineresampling(cl, clspacing)
        cl = vmtkcenterlineattributes(cl)


        #-----------------------------------------------------------------------
        # Trial clips
        # The algorithm is constructed such that the clippoint can not be the
        # first or last point of the centerline. This property is required for
        # calculating clipnormal as defined below.
        #-----------------------------------------------------------------------
        sections = vmtkcenterlinesections(surface, cl)

        closedarray = sections.GetCellData().GetArray('CenterlineSectionClosed')
        maxsizearray = (sections.GetCellData().
                        GetArray('CenterlineSectionMaxSize'))

        highcount = 0
        nbumpcriterion  = round(bumpcriterion*cl.GetNumberOfPoints())

        print "Anything below",nbumpcriterion, "points will be considered a bump."
        for i in range(1, sections.GetNumberOfCells()):

            # Skip sections that are preceded by open sections and skip first
            # 5 centerline points. This to avoid complications near the end of a
            # vein (holes in the surface mesh or sudden changes in centerline
            # direction might otherwise lead to exceeding the threshold far
            # from the atrium body).
            if closedarray.GetValue(i-1):
                # changed to "signed" difference to account for veins w
                # multiple outlets. This created a thin to wide vein change

                slope = (maxsizearray.GetValue(i) -
                            maxsizearray.GetValue(i-1))/ clspacing

                if slope > highslope:
                    highcount += 1
                else:
                    highcount = 0

                print i, slope, highcount

                if slope > maxslope:
                    break
                elif slope > highslope and highcount == nbumpcriterion:
                    break
                else:
                    pass

        if highcount == 0:
            clippointid = i - 1
        else:
            clippointid = i - (highcount)

        np.savetxt(sufixfile +'clippointid' + str(k) + '.csv',
                   np.array([clippointid])+int(nskippoints),
                   fmt='%i')

        #-----------------------------------------------------------------------
        # Prepare output
        #
        # Tried clipping to make nice ostia clips, but it's not so trivial.
        # Back to good old transfer labels
        #-----------------------------------------------------------------------

        vein = clip_vein(surface,cl,clippointid)
        surface = transfer_labels(surface,vein,'autolabels',round(branchlabel[k]))

        if visualise:
            visualise_color(vein,surface,'vein' + str(k))


    vmtksurfacewriter(surface, sufixfile + 'autolabels.vtp')


def clip_vein(surface, cl, clippointid):
    """Clip the vein at clippoint."""

    clippoint0 = cl.GetPoint(clippointid)
    clipnormal = (np.array(cl.GetPoint(clippointid+1)) -
                  np.array(cl.GetPoint(clippointid)))

    possvein = planeclip(surface, clippoint0, clipnormal)
    vein = extractclosestpointregion(possvein,clippoint0)

    return vein


def clip_vein_endpoint(surface, ifile_sufix, targetdistance,
                       laa_seedon = 0, specialvein=0, specialdist=0):
    """Clip vein the targetdistance away from the body."""
    regionslabels = getregionslabels()

    # extract the body from the surface
    # including all points (alloff=1) to avoid holes after appending
    body = pointthreshold(surface, 'autolabels',
                          regionslabels['body'], regionslabels['body'], 1)
    body = extractlargestregion(body)

    # initialize appender with the body
    appender = vtk.vtkAppendPolyData()
    appender.AddInputData(body)
    originaldist = targetdistance

    # if laa seed
    if laa_seedon == 1:
        indeces = ['pv1','pv2','pv3','pv4','laa']
        nseeds = len(indeces)
    else:
        indeces = ['pv1','pv2','pv3','pv4']
        nseeds = len(indeces)

    for k in range(1,nseeds+1):
        index = indeces[k-1]

        # extract vein
        # excluding some points (alloff=0)
        # to avoid overlapping edges after appending
        vein = pointthreshold(surface, 'autolabels', regionslabels[index],
                               regionslabels[index], 0)

        # load the centreline and the clipoint
        cl = readvtp(ifile_sufix + 'clvein' + str(k) + '.vtp')
        clippointid = int(np.loadtxt(ifile_sufix +
                                'clippointid' + str(k) + '.csv'))

        clippoint0 = cl.GetPoint(clippointid)
        clipnormal = (np.array(cl.GetPoint(clippointid + 1)) -
                      np.array(cl.GetPoint(clippointid )))

        abscissasarray = cl.GetPointData().GetArray('Abscissas')
        startabscissa = abscissasarray.GetValue(clippointid)
        currentabscissa = 0
        currentid = clippointid

        # if different distance for 1 vein
        if specialvein > 0:
            if regionslabels[index] == specialvein:
                targetdistance = specialdist
            else:
                targetdistance = originaldist

        # find clip point
        while ((currentabscissa < targetdistance) and
               (currentabscissa >= 0) and
               (currentid >= 0)):
            currentid -= 1
            currentabscissa = startabscissa - abscissasarray.GetValue(currentid)

        if currentid > 0:
            currentid = currentid + 1
        else:
            # vein ended before target distance
            # then clip 2 mm before end of centreline from end point
            currentid = 4

        # clip and append
        clippoint1 = cl.GetPoint(currentid)

        clippedvein = planeclip(vein, clippoint1, clipnormal, 0)

        # keep region closest to ostium point
        clippedvein = extractclosestpointregion(clippedvein, clippoint0)

        # clip generates new points to make a flat cut. The values may be interpolated.
        # we want all values to rounded to a certain label value.
        clippedvein = roundpointarray(clippedvein, 'autolabels')
        appender.AddInputData(clippedvein)

    # collect body + veins
    appender.Update()
    clippedsurface = appender.GetOutput()
    clippedsurface = cleanpolydata(clippedsurface)
    return clippedsurface


def getregionslabels():
    """Return dictionary linking regionids to anatomical locations."""
    regionslabels = {'body': 36,
                     'laa': 37,
                     'pv2': 76,
                     'pv1': 77,
                     'pv3': 78,
                     'pv4': 79}
    return regionslabels

def roundpointarray(polydata, name):
    """Round values in point array."""
    # get original array
    array = polydata.GetPointData().GetArray(name)

    # round labels
    for i in range(polydata.GetNumberOfPoints()):
        value = array.GetValue(i)
        array.SetValue(i, round(value))
    polydata.GetPoints().Modified()
    return polydata

def transfer_mitral_clip(surface, surfaceorig):
    """Extract mitral edge (largest) and use boolean operations to remove area
    below the mitral edge."""

    edges = extractboundaryedge(surfaceorig)
    mitraledge = extractlargestregion(edges)
    cover = vtk.vtkPolyData()
    generatecover(mitraledge, cover)
    mitraledge = append(cover,mitraledge)
    mitraledge = cleanpolydata(mitraledge)

    mitraledgedel = delaunay2D(mitraledge)
    mitraledgedel = cleanpolydata(mitraledgedel)
    mitraledgescaled = scalepolydata(mitraledgedel,1.25)

    mitralsurface = boolean_subtract(surface,mitraledgescaled)
    mitralsurface = extractlargestregion(mitralsurface)
    surface = add_point_array(surface,'mitral',10)
    surface = transfer_labels(surface,mitralsurface,'mitral',5)
    surface = round_labels_array(surface,'mitral',[5,10])
    surfaceclipped = pointthreshold(surface,'mitral',10,10,1) #all off

    surfaceclipped = remove_point_array(surfaceclipped,'mitral')

    return surfaceclipped

def cylinder_mitral_clip(surface, surfaceorig, vis=False, dist=10, arrayname='', arrayvalue=''):
    """Extract mitral edge (largest) and find smallest radius.
    Use cylinder to clip the mitral ring, instead of a plane,
    to avoid clip part of the LAA and/or PVS."""

    if arrayname == '':
        # extract edge
        edges = extractboundaryedge(surfaceorig)
        mitraledge = extractlargestregion(edges)
    else:
        mitraledge= pointthreshold(surfaceorig,
                                   arrayname,
                                   arrayvalue-0.5,
                                   arrayname+0.5)

    # clip with sphere becuase bodies are quite small
    center = pointset_centreofmass(mitraledge)
    normal = pointset_normal(mitraledge)
    closest_point = furthest_point_to_polydata(mitraledge,center)
    radius = euclideandistance(closest_point,center)


    vectordown = multiplyvector(normal,-dist)
    pointdown = acumvectors(vectordown,center)
    clip = cylinderclip(surface,pointdown,center,normal,radius)

    if vis:
        # Visualisation of coordinate system
        centerpd = point2vertexglyph(center)
        centerpdgl  = generateglyph(centerpd)
        axis = linesource(pointdown,center)
        allaxis = append(centerpdgl,axis)
        visualise_color(allaxis,surface,'clip')

    return clip


def edge_mitral_clip(surface, surfaceorig):
    """Extract mitral edge (largest) and fit a plane to it.
    Use plane to clip the mitral annulus."""

    # extrac mitral edge (largest)
    edges = extractboundaryedge(surfaceorig)
    mitraledge = extractlargestregion(edges)

    # fit plane to edge
    center = pointset_centreofmass(mitraledge)
    normal = pointset_normal(mitraledge)

    clip = planeclip(surface,center,normal,0) # insideout off

    return clip


def find_mitral_cylinder_pvs(surface, arrayname, outfile, scale=0.4, w=[0.7, 0.15, 0.15], vis=False):
    """Compute local coordinate system based on the body centroid and PVs centroid.
    The 3 axes are weighted as in w. The resulting vector is used to clip surface
    scale * radius away from the body centroid."""

    # extract body
    startedges = extractboundaryedge(surface)
    if startedges.GetNumberOfPoints()>0:
        surfacefilled = fillholes(surface,1000)
    else:
        surfacefilled = surface

    body = pointthreshold(surface,arrayname,36.0,36.0)
    bodycom = pointset_centreofmass(body)

    # average of left ostia to average of right ostia
    ostia = pointthreshold(surfacefilled,arrayname,78.0,79.0)
    edges = extractboundaryedge(ostia)

    leftcentroid = centroidofcentroids(edges)

    ostia = pointthreshold(surfacefilled,arrayname,76.0,77.0)
    edges = extractboundaryedge(ostia)

    rightcentroid = centroidofcentroids(edges)

    # final pvscom average of left and right
    pvscom = acumvectors(leftcentroid,rightcentroid)
    pvscom = dividevector(pvscom,2)

    # NOW AXES
    # Axis 1: Pvs com to body com
    pvdir = subtractvectors(bodycom,pvscom)
    pvdirn = normalizevector(pvdir)

    # Axis 2: normal to Pvs axis
    ostiadir1 = subtractvectors(leftcentroid,rightcentroid)
    ostiadirn = normalizevector(ostiadir1)

    ostiacross = cross(pvdirn,ostiadirn)
    ostiacrossn = normalizevector(ostiacross)

    # Axis 3: normal to axis 1 and 2
    pvcross = cross(ostiacrossn,pvdirn)
    pvcrossn = normalizevector(pvcross)

    # thought of using for weighting but defualt values seem all right
    bodylength= computelengthalongvector(body,bodycom,pvdirn)
    measurepoint = sumvectors(bodycom,scale*bodylength,pvdirn)
    bodythick = computelengthalongvector(body,measurepoint,ostiacrossn)
    bodywidth = computelengthalongvector(body,measurepoint,pvcrossn)


    pvdirnw = multiplyvector(pvdirn,w[0])
    ostiadirnw = multiplyvector(pvcrossn,w[1])
    ostiacrossnw = multiplyvector(ostiacrossn,w[2])

    plusvector = acumvectors(pvdirnw,ostiacrossnw)
    plusvector = acumvectors(plusvector,ostiadirnw)
    plusvectorn = normalizevector(plusvector)

    # clippoint with length vector
    # in very small bodies, modify scale
    if bodylength/bodythick<1.5:
        scale = 0.45
    clippoint = sumvectors(bodycom,scale*bodylength,plusvectorn)

    # current cut
    slicepv = cutdataset(surface,clippoint,plusvectorn)
    # if > 1 region --> clipping LAA as well.
    nr = countregions(slicepv)
    if nr > 1:
        # keep cut closest to point
        slicepv = extractclosestpointregion(slicepv,clippoint)
        # recompute center
        clippoint = pointset_centreofmass(slicepv)

    vectordown = multiplyvector(plusvectorn,10)
    pointdown = acumvectors(vectordown,clippoint)

    finalbody = cylinderclip(surface,clippoint,pointdown,plusvectorn,bodythick)


    if vis:
        # Visualisation of coordinate system
        pvscompd = point2vertexglyph(pvscom)
        pvscompdg  = generateglyph(pvscompd)

        plotpoint = sumvectors(pvscom,w[0]*bodylength,pvdirn)
        bodyaxis = linesource(plotpoint,pvscom)
        plotpoint = sumvectors(pvscom,w[1]*bodylength,pvcrossn)
        ostiaaxis = linesource(plotpoint,pvscom)
        plotpoint = sumvectors(pvscom,w[2]*bodylength,ostiacrossn)
        crossaxis = linesource(plotpoint,pvscom)

        allaxis = append(pvscompdg,ostiaaxis)
        allaxis = append(allaxis,bodyaxis)
        allaxis = append(allaxis,crossaxis)
        allaxis = append(allaxis,slicepv)

        writevtp(allaxis,outfile+'_axes.vtp')
        visualise_color(allaxis,surfacefilled,'plus')

    writevtp(finalbody,outfile+'.vtp')

    return finalbody

def find_mitral_plane_pvs(surface, arrayname, outfile, scale=0.4, w=[0.7,0.15,0.15], vis=False):
    """Compute local coordinate system based on the body centroid and PVs centroid.
    The 3 axes are weighted as in w. The resulting vector is used to clip surface
    scale * radius away from the body centroid."""

    # extract body
    startedges = extractboundaryedge(surface)
    if startedges.GetNumberOfPoints()>0:
        surfacefilled = fillholes(surface,1000)
    else:
        surfacefilled = surface

    body = pointthreshold(surface,arrayname,36.0,36.0)
    bodycom = pointset_centreofmass(body)

    # average of left ostia to average of right ostia
    ostia = pointthreshold(surfacefilled,arrayname,78.0,79.0)
    edges = extractboundaryedge(ostia)

    leftcentroid = centroidofcentroids(edges)

    ostia = pointthreshold(surfacefilled,arrayname,76.0,77.0)
    edges = extractboundaryedge(ostia)

    rightcentroid = centroidofcentroids(edges)

    # final pvscom average of left and right
    pvscom = acumvectors(leftcentroid,rightcentroid)
    pvscom = dividevector(pvscom,2)

    # NOW AXES
    # Axis 1: Pvs com to body com
    pvdir = subtractvectors(bodycom,pvscom)
    pvdirn = normalizevector(pvdir)

    # Axis 2: normal to Pvs axis
    ostiadir1 = subtractvectors(leftcentroid,rightcentroid)
    ostiadirn = normalizevector(ostiadir1)

    ostiacross = cross(pvdirn,ostiadirn)
    ostiacrossn = normalizevector(ostiacross)

    # Axis 3: normal to axis 1 and 2
    pvcross = cross(ostiacrossn,pvdirn)
    pvcrossn = normalizevector(pvcross)

    # thought of using for weighting but defualt values seem all right
    bodylength= computelengthalongvector(body,bodycom,pvdirn)
    measurepoint = sumvectors(bodycom,scale*bodylength,pvdirn)
    bodythick = computelengthalongvector(body,measurepoint,ostiacrossn)
    bodywidth = computelengthalongvector(body,measurepoint,pvcrossn)

    pvdirnw = multiplyvector(pvdirn,w[0])
    ostiadirnw = multiplyvector(pvcrossn,w[1])
    ostiacrossnw = multiplyvector(ostiacrossn,w[2])

    plusvector = acumvectors(pvdirnw,ostiacrossnw)
    plusvector = acumvectors(plusvector,ostiadirnw)
    plusvectorn = normalizevector(plusvector)

    #  clippoint with length vector
    # in very small bodies, modify scale
    if bodylength/bodythick<1.5:
        scale = 0.45
    clippoint = sumvectors(bodycom,scale*bodylength,plusvectorn)
    slicepv=cutdataset(surface,clippoint,plusvectorn)

    # adjust clipoint if slice has more than 1 edge
    # that means we are clipping through LAA and/or veins
    # make it further from centroid
    nr = countregions(slicepv)
    while nr > 1:

        scale = scale + 0.005
        clippoint = sumvectors(bodycom,scale*bodylength,plusvectorn)
        slicepv=cutdataset(surface,clippoint,plusvectorn)
        nr = countregions(slicepv)

    if vis:
        # Visualisation of coordinate system
        pvscompd = point2vertexglyph(pvscom)
        pvscompdg  = generateglyph(pvscompd)

        plotpoint = sumvectors(pvscom,w[0]*bodylength,pvdirn)
        bodyaxis = linesource(plotpoint,pvscom)
        plotpoint = sumvectors(pvscom,w[1]*bodylength,pvcrossn)
        ostiaaxis = linesource(plotpoint,pvscom)
        plotpoint = sumvectors(pvscom,w[2]*bodylength,ostiacrossn)
        crossaxis = linesource(plotpoint,pvscom)

        allaxis = append(pvscompdg,ostiaaxis)
        allaxis = append(allaxis,bodyaxis)
        allaxis = append(allaxis,crossaxis)
        allaxis = append(allaxis,slicepv)

        writevtp(allaxis,outfile+'_axes.vtp')
        visualise_color(allaxis,surfacefilled,'plus')

    # save mesh and txt
    insideout = 1
    finalbody = planeclip(surface,clippoint,plusvectorn,insideout) # insideout OFF
    finalbody = extractsurface(finalbody)

    writevtp(finalbody,outfile+'.vtp')
    # save CSV
    plane_to_csv(plusvectorn,clippoint,insideout,outfile + '.csv')

    return finalbody

def plane_to_csv(normal, point, onoff, outfile):
    """Save plane normal and center to csv."""

    f = open(outfile, 'wb')
    # write values of row by row

    line = str(normal[0]) + ',' + str(normal[1]) + ',' + str(normal[2]) + '\n'
    f.write(line)

    line = str(point[0]) + ',' + str(point[1]) + ',' + str(point[2]) + '\n'
    f.write(line)

    line = str(onoff) + ',' + ',' + '\n'
    f.write(line)

    f.close()

def initial_transform_pvends(source, target, arrayname, similarityon=0, laa_seedon=False):
    """Initial registration of each subject to the atlas mesh."""

    regionslabels = getregionslabels()

    # find the lmks
    sourcepoints = vtk.vtkPoints()
    targetpoints = vtk.vtkPoints()

    # extract edges
    edgessource = extractboundaryedge(source)
    edgestarget = extractboundaryedge(target)

    # add mitral seed
    mitraledge = pointthreshold(edgessource, arrayname, regionslabels['body'],
                regionslabels['body'],0)
    point = pointset_centreofmass(mitraledge)
    sourcepoints.InsertNextPoint(point)

    mitraledge = extractlargestregion(edgestarget)
    # find center of mass
    point = pointset_centreofmass(mitraledge)
    targetpoints.InsertNextPoint(point)

    # threshold each vein
    if laa_seedon:
        indexes = ['pv1','pv2','pv3','pv4','laa']
    else:
        indexes = ['pv1','pv2','pv3','pv4']

    for index in indexes:
        #sourcepoints
        veinedge = pointthreshold(edgessource, arrayname, regionslabels[index],
                regionslabels[index],0)
        # keep lasrgest only region
        largestedge = extractlargestregion(veinedge)
        # find centroid of each edge
        point = pointset_centreofmass(largestedge)
        sourcepoints.InsertNextPoint(point)

        #targetpoints
        veinedge = pointthreshold(edgestarget, arrayname, regionslabels[index],
                regionslabels[index],0)
        # keep largest only region
        largestedge = extractlargestregion(veinedge)
        # find centroid of each edge
        point = pointset_centreofmass(largestedge)
        targetpoints.InsertNextPoint(point)

    sourceintarget = transform_lmk(sourcepoints,targetpoints,source,similarityon)

    return sourceintarget

def zero_truncate_array(surface, arrayname, value=0):
    """Any value below zero is set to value."""
    array = surface.GetPointData().GetArray(arrayname)

    for i in range(surface.GetNumberOfPoints()):
        if array.GetValue(i) <= 0:
            array.SetValue(i, value)
    surface.GetPoints().Modified()
    return surface

def restore_nan_values(surface, arrayname, realmin):
    """Restore nan values."""
    array = surface.GetPointData().GetArray(arrayname)
    for p in range(array.GetNumberOfTuples()):
        val = array.GetValue(p)
        if val < realmin :
            array.SetValue(p,vtk.vtkMath.Nan())
    surface.GetPoints().Modified()
    return surface


def project_nan_array(source,
                      target,
                      arrayname,
                      point_map='',
                      mapstring=''):

    """Initialise with nan. Tranfer non-nan values to closest point."""
    sourcearray = source.GetPointData().GetArray(arrayname)

    # initialise nan
    newarray = vtk.vtkDoubleArray()
    newarray.SetName(arrayname)
    newarray.SetNumberOfTuples(target.GetNumberOfPoints())
    newarray.FillComponent(0, vtk.vtkMath.Nan())

    # pointlocator
    locator = vtk.vtkPointLocator()
    locator.SetDataSet(target)
    locator.BuildLocator()

    # visit each source value
    # find non-nan values and transfer
    nop = source.GetNumberOfPoints()
    for p in range(nop):
        val = sourcearray.GetValue(p)
        point = source.GetPoint(p)
        if not vtk.vtkMath.IsNan(val):
            closestpointid = locator.FindClosestPoint(point)
            newarray.SetValue(closestpointid,val)
            if point_map != '':
                point_map[ mapstring ].append( (p,closestpointid) )
    target.GetPointData().AddArray(newarray)

    return target

def original2unfold_point_map(point_map):
    """computing a single mapping from original to unfold"""

    labels = [
    'original2std',
    'std2average' ,
    'average2unfold']

    # todo: check what happens if we don't have
    # one to one mapping

    # turnign list of tuples in two list
    std2average = zip(*point_map['std2average'])
    average2unfold = zip(*point_map['average2unfold'])

    for orig,std in point_map['original2std']:

        std_index = std2average[0].index(std)
        avg_point = std2average[1][std_index]

        avg_index = average2unfold[0].index(avg_point)
        unfold_point = average2unfold[1][avg_index]

        point_map[ 'original2unfold' ].append( (orig,unfold_point) )

def point_map_dictionary(labels=''):
    if labels == '':
        labels = point_map_labels()
    dictionary={}
    for label in labels:
        dictionary[label] = []
    return dictionary

def point_map_labels():
    return [
    'original2std',
    'std2average' ,
    'average2unfold',
    'original2unfold']

def writepointmap2csv(point_map,
                       ofile,
                       labels=''):
    """Write point_map to csv."""
    # subset of pointids, in this order
    if labels == '':
        labels = point_map_labels()

    npoints = point_map_dictionary()
    maxpoints = 0
    for label in labels:
        npoints[label]=len(point_map[label])
        if len(point_map[label]) > maxpoints:
            maxpoints = len(point_map[label])

    with open(ofile, 'w') as csvfile:
        writer = csv.DictWriter(csvfile,labels)
        writer.writeheader()
        for row in range(maxpoints):
            minidict = point_map_dictionary(labels)
            for label in labels:
                if row < npoints[label]:
                    minidict[label] = point_map[label][row]
            writer.writerow(minidict)
        csvfile.close()

def find_closest_lmk(subjectmesh, atlasmesh):
    """with a predefined list of landmarks in the atlas (pointID).
    Find closest point in subject mesh."""

    lmks = [2774,1990,8774]
    subjectlmks = ['','','']

    # pointlocator
    locator = vtk.vtkPointLocator()
    locator.SetDataSet(subjectmesh)
    locator.BuildLocator()

    # visit each lmk
    # find points in subject mesh
    for p in len(lmks):
        point = atlasmesh.GetPoint(lmks[p])
        closestpointid = locator.FindClosestPoint(point)
        subjectlmks[p] = closestpointid

    return subjectlmks

def make_nan_glyph(source, arrayname, glyphsize=0.04):
    """Make glyph at on-nan values."""
    sourcearray = source.GetPointData().GetArray(arrayname)

    # create the pointset
    newpoints = vtk.vtkPoints()
    newvertices = vtk.vtkCellArray()

    # create array on non nan values
    newarray = vtk.vtkDoubleArray()
    newarray.SetName(arrayname)

    # visit each array value
    # find non-nan values
    for p in range(source.GetNumberOfPoints()):
        val = sourcearray.GetValue(p)
        point = source.GetPoint(p)

        if not vtk.vtkMath.IsNan(val):
            # save in pointset
            pid = newpoints.InsertNextPoint(point)
            newarray.InsertNextValue(val)
            # Create the topology of the point (a vertex)
            newvertices.InsertNextCell(1)
            newvertices.InsertCellPoint(pid)

    pointspd = vtk.vtkPolyData()
    pointspd.SetPoints(newpoints)
    pointspd.SetVerts(newvertices)
    pointspd.GetPointData().AddArray(newarray)

    glyph = generateglyph(pointspd,glyphsize)

    return glyph

