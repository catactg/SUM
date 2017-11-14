# from __future__ import division
import math
import vtk
import numpy as np
from colormaps import *
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy

#------------------------------------------------------------------------------
# Linear algebra
#------------------------------------------------------------------------------

def angle(v1, v2):
    return math.acos(dot(v1, v2) / (normvector(v1) * normvector(v2)))

def acumvectors(point1, point2):
    return [point1[0] + point2[0],
            point1[1] + point2[1],
            point1[2] + point2[2]]

def subtractvectors(point1, point2):
    return [point1[0] - point2[0],
            point1[1] - point2[1],
            point1[2] - point2[2]]

def dividevector(point, n):
    nr = float(n)
    return [point[0]/nr, point[1]/nr, point[2]/nr]

def multiplyvector(point, n):
    nr = float(n)
    return [nr*point[0], nr*point[1], nr*point[2]]

def sumvectors(vect1, scalar, vect2):
    return [vect1[0] + scalar*vect2[0],
            vect1[1] + scalar*vect2[1],
            vect1[2] + scalar*vect2[2]]

def cross(v1, v2):
    return [v1[1]*v2[2] - v1[2]*v2[1],
            v1[2]*v2[0] - v1[0]*v2[2],
            v1[0]*v2[1] - v1[1]*v2[0]]


def dot(v1, v2):
    return sum((a*b) for a, b in zip(v1, v2))


def euclideandistance(point1, point2):
    return math.sqrt((point1[0] - point2[0])**2 +
                     (point1[1] - point2[1])**2 +
                     (point1[2] - point2[2])**2)

def normvector(v):
    return math.sqrt(dot(v, v))


def normalizevector(v):
    norm = normvector(v)
    return [v[0] / norm,
            v[1] / norm,
            v[2] / norm]

def polar2cart(r, theta, center):
    theta_r = np.deg2rad(theta)
    x = r  * np.cos(theta_r) + center[0]
    y = r  * np.sin(theta_r) + center[1]
    return x, y


#------------------------------------------------------------------------------
# VTK
#------------------------------------------------------------------------------

def compute_nonzeromean(polydata,inarrayname,pointdata=True):

    if polydata.GetNumberOfPoints() > 0 :
        if pointdata:
            array = vtk_to_numpy(polydata.GetPointData().GetArray(inarrayname))
        else:
            array = vtk_to_numpy(polydata.GetCellData().GetArray(inarrayname))

        array_masked = np.ma.array( array, mask = (array <= 0) )
        array_nonan = array_masked[~np.isnan(array_masked)]
        return array_nonan.mean()
    else:
        return 0.

def transfer_labels(surface,ref,arrayname,value):

    # initiate point locator
    locator = vtk.vtkPointLocator()
    locator.SetDataSet(surface)
    locator.BuildLocator()

    # get array from surface
    array = surface.GetPointData().GetArray(arrayname)

    # go through each point of ref surface, determine closest point on surface,
    for i in range(ref.GetNumberOfPoints()):
        point = ref.GetPoint(i)
        closestpoint_id = locator.FindClosestPoint(point)
        array.SetValue(closestpoint_id, value)
    surface.GetPoints().Modified()
    return surface

def generatecover(edges, cover, arrayname=''):
    """Create caps for capping a surface with holes."""
    # create the building blocks of polydata.
    polys = vtk.vtkCellArray()
    points = vtk.vtkPoints()

    surfilt = vtk.vtkCleanPolyData()
    surfilt.SetInputData( edges )
    surfilt.Update()

    points.DeepCopy(surfilt.GetOutput().GetPoints())
    npoints = points.GetNumberOfPoints()

    if arrayname:
        # keep pre existing array
        array = surfilt.GetOutput().GetPointData().GetArray(arrayname)
        arraynp = vtk_to_numpy(array)
        array.InsertNextValue(np.mean(arraynp))

    # add centroid
    centr = np.zeros(3)
    for i in range( npoints ):
        pt = np.zeros(3)
        points.GetPoint(i,pt)
        centr = centr + pt

    centr = centr / npoints
    cntpt = points.InsertNextPoint(centr)

    # add cells
    for i in range(surfilt.GetOutput().GetNumberOfCells()):
        cell = surfilt.GetOutput().GetCell(i)
        polys.InsertNextCell(3)
        polys.InsertCellPoint(cell.GetPointId(0))
        polys.InsertCellPoint(cell.GetPointId(1))
        polys.InsertCellPoint(cntpt)

    # assign the pieces to the polydata
    cover.SetPoints(points)
    cover.SetPolys(polys)
    if arrayname:
        cover.GetPointData().AddArray(array)

def pointset_centreofmass(polydata):
    centre = [0, 0, 0]
    for i in range(polydata.GetNumberOfPoints()):
        point = [polydata.GetPoints().GetPoint(i)[0],
          polydata.GetPoints().GetPoint(i)[1],
          polydata.GetPoints().GetPoint(i)[2]]
        centre = acumvectors(centre,point)
    return dividevector(centre, polydata.GetNumberOfPoints())


def pointset_normal(polydata):
    # estimate the normal for a given set of points which are supposedly in a plane
    com = pointset_centreofmass(polydata)
    normal = [0, 0, 0]
    n = 0
    for i in range(polydata.GetNumberOfPoints()):
        point = [polydata.GetPoints().GetPoint(i)[0],
          polydata.GetPoints().GetPoint(i)[1],
          polydata.GetPoints().GetPoint(i)[2]]
        if n==0:
            vect = subtractvectors(point, com)
        else:
            vect2 = subtractvectors(point, com)
            crossprod = cross(vect, vect2)
            crossprod = normalizevector(crossprod)
            # make sure each normal is oriented coherently ...
            if n==1:
                normal2 = crossprod
            else:
                if dot(crossprod,normal2) < 0:
                    crossprod = [-crossprod[0],-crossprod[1],-crossprod[2]]
            normal = acumvectors(normal, crossprod)
        n += 1
    return normalizevector(normal)

def delaunay2D(polydata):
    delny = vtk.vtkDelaunay2D()
    delny.SetInputData(polydata)
    delny.SetTolerance(0.1)
    delny.SetAlpha(0.0)
    delny.BoundingTriangulationOff()
    delny.SetProjectionPlaneMode(vtk.VTK_BEST_FITTING_PLANE)
    delny.Update()

    return delny.GetOutput()


def smooth(polydata,iterations,factor):
    smoother = vtk.vtkSmoothPolyDataFilter()
    smoother.SetInputData(polydata)
    smoother.SetNumberOfIterations(iterations)
    smoother.FeatureEdgeSmoothingOn()
    smoother.SetRelaxationFactor(factor)
    smoother.Update()
    return smoother.GetOutput()

def add_cell_array(polydata,name,value):
    # create array and add a label
    array = vtk.vtkDoubleArray()
    array.SetName(name)
    array.SetNumberOfTuples(polydata.GetNumberOfCells())

    for i in range(polydata.GetNumberOfCells()):
        array.SetValue(i,value)
    polydata.GetCellData().AddArray(array)

    return polydata


def add_point_array(polydata,name,value):
    # create array and add a label
    array = vtk.vtkDoubleArray()
    array.SetName(name)
    array.SetNumberOfTuples(polydata.GetNumberOfPoints())

    for i in range(polydata.GetNumberOfPoints()):
        array.SetValue(i,value)
    polydata.GetPointData().AddArray(array)

    return polydata

def append(polydata1,polydata2):
    appender = vtk.vtkAppendPolyData()
    appender.AddInputData(polydata1)
    appender.AddInputData(polydata2)
    appender.Update()
    return appender.GetOutput()

def cellthreshold(polydata, arrayname, start=0, end=1):
    threshold = vtk.vtkThreshold()
    threshold.SetInputData(polydata)
    threshold.SetInputArrayToProcess(0,0,0,vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS,arrayname)
    threshold.ThresholdBetween(start,end)
    threshold.Update()

    surfer = vtk.vtkDataSetSurfaceFilter()
    surfer.SetInputConnection(threshold.GetOutputPort())
    surfer.Update()
    return surfer.GetOutput()

def centroidofcentroids(edges):
    # compute centroids of each edge
    # find average point
    acumvector = [0,0,0]
    rn = countregions(edges)
    for r in range(rn):
        oneedge = extractconnectedregion(edges,r)
        onecentroid = pointset_centreofmass(oneedge)
        acumvector = acumvectors(acumvector,onecentroid)
    finalcentroid = dividevector(acumvector,rn)
    return finalcentroid

def cleanpolydata(polydata):
    cleaner = vtk.vtkCleanPolyData()
    cleaner.SetInputData(polydata)
    cleaner.Update()
    return cleaner.GetOutput()

def computelengthalongvector(polydata,refpoint,vector):
    # polydata should be a closed surface

    # intersect with line
    point1 = refpoint
    point2 = sumvectors(refpoint,1000,vector) # far away point
    intersectpoints = intersectwithline(polydata,point1,point2)
    furthestpoint1 = furthest_point_to_polydata(intersectpoints,refpoint)

    # intersect with line the other way
    point1 = refpoint
    point2 = sumvectors(refpoint,-1000,vector) # far away point
    intersectpoints = intersectwithline(polydata,point1,point2)
    furthestpoint2 = furthest_point_to_polydata(intersectpoints,furthestpoint1)


    length = euclideandistance(furthestpoint1,furthestpoint2)
    return length

def countregions(polydata):
    # NOTE: preventive measures: clean before connectivity filter
    # to avoid artificial regionIds
    # It slices the surface down the middle
    surfer = vtk.vtkDataSetSurfaceFilter()
    surfer.SetInputData(polydata)
    surfer.Update()

    cleaner = vtk.vtkCleanPolyData()
    cleaner.SetInputConnection(surfer.GetOutputPort())
    cleaner.Update()

    connect = vtk.vtkPolyDataConnectivityFilter()
    connect.SetInputConnection(cleaner.GetOutputPort())
    connect.Update()
    return connect.GetNumberOfExtractedRegions()


def cutdataset(dataset, point, normal):
    cutplane = vtk.vtkPlane()
    cutplane.SetOrigin(point)
    cutplane.SetNormal(normal)
    cutter = vtk.vtkCutter()
    cutter.SetInputData(dataset)
    cutter.SetCutFunction(cutplane)
    cutter.Update()
    return cutter.GetOutput()

def cylinderclip(dataset, point0, point1,normal,radius):
    """Define cylinder. The cylinder is infinite in extent. We therefore have
    to truncate the cylinder using vtkImplicitBoolean in combination with
    2 clipping planes located at point0 and point1. The radius of the
    cylinder is set to be slightly larger than 'maxradius'."""

    rotationaxis = cross([0, 1, 0], normal)
    rotationangle = (180 / math.pi) * angle([0, 1, 0], normal)

    transform = vtk.vtkTransform()
    transform.Translate(point0)
    transform.RotateWXYZ(rotationangle, rotationaxis)
    transform.Inverse()

    cylinder = vtk.vtkCylinder()
    cylinder.SetRadius(radius)
    cylinder.SetTransform(transform)

    plane0 = vtk.vtkPlane()
    plane0.SetOrigin(point0)
    plane0.SetNormal([-x for x in normal])
    plane1 = vtk.vtkPlane()
    plane1.SetOrigin(point1)
    plane1.SetNormal(normal)

    clipfunction = vtk.vtkImplicitBoolean()
    clipfunction.SetOperationTypeToIntersection()
    clipfunction.AddFunction(cylinder)
    clipfunction.AddFunction(plane0)
    clipfunction.AddFunction(plane1)

    clipper = vtk.vtkClipPolyData()
    clipper.SetInputData(dataset)
    clipper.SetClipFunction(clipfunction)
    clipper.Update()

    return extractlargestregion(clipper.GetOutput())

def furthest_point_to_polydata(pointset,refpoint):
    # visist each point in pointset
    # selecte point furthest from reference point
    refdist = 0
    for i in range(pointset.GetNumberOfPoints()):

        dist = euclideandistance(pointset.GetPoint(i),refpoint)
        if dist > refdist:
            refdist = dist
            selectedpointid = i
    return pointset.GetPoint(selectedpointid)

def transfer_array_by_pointid(ref,target,arrayname,targetarrayname):
    # get array from reference
    refarray = ref.GetPointData().GetArray(arrayname)

    # create new array
    numberofpoints = target.GetNumberOfPoints()
    newarray = vtk.vtkDoubleArray()
    newarray.SetName(targetarrayname)
    newarray.SetNumberOfTuples(numberofpoints)

    # go through each point of target surface,
    for i in range(target.GetNumberOfPoints()):
        value = refarray.GetValue(i)
        newarray.SetValue(i, value)
    target.GetPointData().AddArray(newarray)

    return target

def extractboundaryedge(polydata):
    edge = vtk.vtkFeatureEdges()
    edge.SetInputData(polydata)
    edge.FeatureEdgesOff()
    edge.NonManifoldEdgesOff()
    edge.Update()
    return edge.GetOutput()


def extractcells(polydata, idlist):
    """Extract cells from polydata whose cellid is in idlist."""
    cellids = vtk.vtkIdList()  # specify cellids
    cellids.Initialize()
    for i in idlist:
        cellids.InsertNextId(i)

    extract = vtk.vtkExtractCells()  # extract cells with specified cellids
    extract.SetInputData(polydata)
    extract.AddCellList(cellids)

    geometry = vtk.vtkGeometryFilter()  # unstructured grid to polydata
    geometry.SetInputConnection(extract.GetOutputPort())
    geometry.Update()
    return geometry.GetOutput()

def skippoints(polydata,nskippoints):
    """Generate a single cell line from points in idlist."""

    # derive number of nodes
    numberofnodes = polydata.GetNumberOfPoints() - nskippoints

    # define points and line
    points = vtk.vtkPoints()
    polyline = vtk.vtkPolyLine()
    polyline.GetPointIds().SetNumberOfIds(numberofnodes)

    # assign id and x,y,z coordinates
    for i in range(nskippoints,polydata.GetNumberOfPoints()):
        pointid = i - nskippoints
        polyline.GetPointIds().SetId(pointid,pointid)
        point = polydata.GetPoint(i)
        points.InsertNextPoint(point)


    # define cell
    cells = vtk.vtkCellArray()
    cells.InsertNextCell(polyline)

    # add to polydata
    polyout = vtk.vtkPolyData()
    polyout.SetPoints(points)
    polyout.SetLines(cells)

    return polyout

def extractclosestpointregion(polydata, point=[0, 0, 0]):
    # NOTE: preventive measures: clean before connectivity filter
    # to avoid artificial regionIds
    # It slices the surface down the middle
    surfer = vtk.vtkDataSetSurfaceFilter()
    surfer.SetInputData(polydata)
    surfer.Update()

    cleaner = vtk.vtkCleanPolyData()
    cleaner.SetInputConnection(surfer.GetOutputPort())
    cleaner.Update()

    connect = vtk.vtkPolyDataConnectivityFilter()
    connect.SetInputConnection(cleaner.GetOutputPort())
    connect.SetExtractionModeToClosestPointRegion()
    connect.SetClosestPoint(point)
    connect.Update()
    return connect.GetOutput()

def extractconnectedregion(polydata, regionid):
    # NOTE: preventive measures: clean before connectivity filter
    # to avoid artificial regionIds
    # It slices the surface down the middle
    surfer = vtk.vtkDataSetSurfaceFilter()
    surfer.SetInputData(polydata)
    surfer.Update()

    cleaner = vtk.vtkCleanPolyData()
    cleaner.SetInputConnection(surfer.GetOutputPort())
    cleaner.Update()

    connect = vtk.vtkPolyDataConnectivityFilter()
    connect.SetInputConnection(cleaner.GetOutputPort())
    connect.SetExtractionModeToAllRegions()
    connect.ColorRegionsOn()
    connect.Update()
    surface = pointthreshold(connect.GetOutput(),'RegionId',float(regionid),float(regionid))
    return surface

def extractlargestregion(polydata):
    # NOTE: preventive measures: clean before connectivity filter
    # to avoid artificial regionIds
    # It slices the surface down the middle
    surfer = vtk.vtkDataSetSurfaceFilter()
    surfer.SetInputData(polydata)
    surfer.Update()

    cleaner = vtk.vtkCleanPolyData()
    cleaner.SetInputConnection(surfer.GetOutputPort())
    cleaner.Update()

    connect = vtk.vtkPolyDataConnectivityFilter()
    connect.SetInputConnection(cleaner.GetOutputPort())
    connect.SetExtractionModeToLargestRegion()
    connect.Update()

    # leaves phantom points ....
    cleaner = vtk.vtkCleanPolyData()
    cleaner.SetInputConnection(connect.GetOutputPort())
    cleaner.Update()
    return cleaner.GetOutput()

def extractsurface(polydata):
    surfer = vtk.vtkDataSetSurfaceFilter()
    surfer.SetInputData(polydata)
    surfer.Update()

    return surfer.GetOutput()

def fillholes(polydata,size):
    filler = vtk.vtkFillHolesFilter()
    filler.SetInputData(polydata)
    filler.SetHoleSize(size)
    filler.Update()

    return filler.GetOutput()

def getregionslabels():
    """Return dictionary linking regionids to anatomical locations."""
    regionslabels = {'body': 36,
                     'laa': 37,
                     'pv2': 76,
                     'pv1': 77,
                     'pv3': 78,
                     'pv4': 79}
    return regionslabels

def getflatregionslabels():
    """Return dictionary linking SUM regionids to anatomical locations."""
    regionslabels = {'ant': 1,
                     'lat': 2,
                     'lattop': 2,
                     'laa_around': 3,
                     'laa_bridge_left': 3,
                     'laa_bridge_right': 3,
                     'roof': 4,
                     'roof_addon': 4,
                     'post': 5,
                     'isthmus': 6,
                     'floor': 7,
                     'floor_addon': 7,
                     'septum': 8,
                     'lpv_sup_q1': 9,
                     'lpv_sup_q2': 10,
                     'lpv_sup_q3': 11,
                     'lpv_sup_q4': 12,
                     'lpv_inf_q1': 13,
                     'lpv_inf_q2': 14,
                     'lpv_inf_q3': 15,
                     'lpv_inf_q4': 16,
                     'rpv_sup_q1': 17,
                     'rpv_sup_q2': 18,
                     'rpv_sup_q3': 19,
                     'rpv_sup_q4': 20,
                     'rpv_inf_q1': 21,
                     'rpv_inf_q2': 22,
                     'rpv_inf_q3': 23,
                     'rpv_inf_q4': 24}


    return regionslabels

def generateglyph(polyIn,scalefactor=2):
    vertexGlyphFilter = vtk.vtkGlyph3D()
    sphereSource = vtk.vtkSphereSource()
    vertexGlyphFilter.SetSourceData(sphereSource.GetOutput())
    vertexGlyphFilter.SetInputData(polyIn)
    vertexGlyphFilter.SetColorModeToColorByScalar()
    vertexGlyphFilter.SetSourceConnection(sphereSource.GetOutputPort())
    vertexGlyphFilter.ScalingOn()
    vertexGlyphFilter.SetScaleFactor(scalefactor)
    vertexGlyphFilter.Update()
    return vertexGlyphFilter.GetOutput()

def intersectwithline(surface,p1,p2):

   # Create the locator
    tree = vtk.vtkOBBTree()
    tree.SetDataSet(surface)
    tree.BuildLocator()

    intersectPoints = vtk.vtkPoints()
    intersectCells = vtk.vtkIdList()

    tolerance=1.e-3
    tree.SetTolerance(tolerance)
    tree.IntersectWithLine(p1,p2,intersectPoints,intersectCells)

    return intersectPoints

def linesource(p1,p2):
    source = vtk.vtkLineSource()
    source.SetPoint1(p1[0],p1[1],p1[2])
    source.SetPoint2(p2[0],p2[1],p2[2])

    return source.GetOutput()

def planeclip(surface, point, normal, insideout=1):
    clipplane = vtk.vtkPlane()
    clipplane.SetOrigin(point)
    clipplane.SetNormal(normal)
    clipper = vtk.vtkClipPolyData()
    clipper.SetInputData(surface)
    clipper.SetClipFunction(clipplane)

    if insideout == 1:
        clipper.InsideOutOn()
    else:
        clipper.InsideOutOff()
    clipper.Update()
    return clipper.GetOutput()

def point2vertexglyph(point):
    points = vtk.vtkPoints()
    points.InsertNextPoint(point[0],point[1],point[2])

    poly = vtk.vtkPolyData()
    poly.SetPoints(points)

    glyph = vtk.vtkVertexGlyphFilter()
    glyph.SetInputData(poly)
    glyph.Update()
    return glyph.GetOutput()

def pointthreshold(polydata, arrayname, start=0, end=1,alloff=0):
    threshold = vtk.vtkThreshold()
    threshold.SetInputData(polydata)
    threshold.SetInputArrayToProcess(0,0,0,vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS,arrayname)
    threshold.ThresholdBetween(start,end)
    if (alloff):
        threshold.AllScalarsOff()
    threshold.Update()

    surfer = vtk.vtkDataSetSurfaceFilter()
    surfer.SetInputConnection(threshold.GetOutputPort())
    surfer.Update()
    return surfer.GetOutput()

def scalepolydata(polydata, scalefactor, inorigin=False):

    surfacecenter = polydata.GetCenter()
    toorigin = [0,0,0]
    toorigin[0] = -1*surfacecenter[0]
    toorigin[1] = -1*surfacecenter[1]
    toorigin[2] = -1*surfacecenter[2]

    # bring to origin + rotate + bring back
    transform = vtk.vtkTransform()
    transform.PostMultiply()
    transform.Translate(toorigin)
    transform.Scale(scalefactor,scalefactor,scalefactor)
    if not inorigin:
        transform.Translate(surfacecenter)

    transformfilter = vtk.vtkTransformFilter()
    transformfilter.SetTransform(transform)
    transformfilter.SetInputData(polydata)
    transformfilter.Update()

    return transformfilter.GetOutput()

def surfacearea(polydata):
    properties = vtk.vtkMassProperties()
    properties.SetInputData(polydata)
    properties.Update()
    return properties.GetSurfaceArea()

def triangulate(polydata):
    trianglefilter = vtk.vtkTriangleFilter()
    trianglefilter.SetInputData(polydata)
    trianglefilter.Update()
    return trianglefilter.GetOutput()

def transform_lmk(sourcepoints,targetpoints,surface,similarityon=False):
    lmktransform = vtk.vtkLandmarkTransform()
    lmktransform.SetSourceLandmarks(sourcepoints)

    lmktransform.SetTargetLandmarks(targetpoints)
    if similarityon:
        lmktransform.SetModeToSimilarity()
    else:
        lmktransform.SetModeToAffine()

    lmktransform.Update()

    transformfilter = vtk.vtkTransformPolyDataFilter()
    transformfilter.SetInputData(surface)
    transformfilter.SetTransform(lmktransform)
    transformfilter.Update()

    return transformfilter.GetOutput()

def round_labels_array(surface,arrayname,labels):
    """Any value that is not part of the labels is rounded to minvalue."""
    # threshold range step = 1
    minval = min(labels)
    maxval = max(labels)
    dif = np.zeros(len(labels))
    for val in range(minval,maxval+1):
        mindif = 10000
        closestlabel = 0
        patch = pointthreshold(surface,arrayname,val-0.5,val+0.5,1) # all off
        if patch.GetNumberOfPoints()>0:
            for l in range(0,len(labels)):
                dif[l] = val - labels[l]
            mindif = min(abs(dif))
            if (mindif > 0.01):
                # found points to round
                transfer_labels(surface,patch,arrayname,minval)
    return surface

def visualise_default(surface,ref,case,arrayname,mini,maxi):
    """Visualise surface with a default parameters."""

    #Create a lookup table to map cell data to colors
    lut = vtk.vtkLookupTable()
    lut.SetNumberOfTableValues(255)
    lut.SetValueRange(0, 255)

    # qualitative data from colorbrewer  --> matching qualitative colormap of Paraview
    lut.SetTableValue(0     , 0     , 0     , 0, 1)  #Black
    lut.SetTableValue(mini, 1,1,1, 1) # white
    lut.SetTableValue(mini+1, 77/255.,175/255., 74/255. ,     1)  # green
    lut.SetTableValue(maxi-3, 152/255.,78/255.,163/255., 1) # purple
    lut.SetTableValue(maxi-2, 255/255.,127/255., 0., 1) # orange
    lut.SetTableValue(maxi-1, 55/255., 126/255., 184/255., 1) # blue
    lut.SetTableValue(maxi, 166/255.,86/255.,40/255., 1) # brown
    lut.Build()

    # create a text actor
    txt = vtk.vtkTextActor()
    txt.SetInput(case)
    txtprop=txt.GetTextProperty()
    txtprop.SetFontFamilyToArial()
    txtprop.SetFontSize(18)
    txtprop.SetColor(0, 0, 0)
    txt.SetDisplayPosition(20, 30)

    # create a rendering window, renderer, and renderwindowinteractor
    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    iren = vtk.vtkRenderWindowInteractor()
    style = vtk.vtkInteractorStyleTrackballCamera()
    iren.SetInteractorStyle(style)
    iren.SetRenderWindow(renWin)

    # surface mapper and actor
    surfacemapper = vtk.vtkPolyDataMapper()
    surfacemapper.SetInputData(surface)
    surfacemapper.SetScalarModeToUsePointFieldData()
    surfacemapper.SelectColorArray(arrayname)
    surfacemapper.SetLookupTable(lut)
    surfacemapper.SetScalarRange(0,255)
    surfaceactor = vtk.vtkActor()
    surfaceactor.SetMapper(surfacemapper)

    # refsurface mapper and actor
    refmapper = vtk.vtkPolyDataMapper()
    refmapper.SetInputData(ref)
    refmapper.SetScalarModeToUsePointFieldData()
    refmapper.SelectColorArray(arrayname)
    refmapper.SetLookupTable(lut)
    refmapper.SetScalarRange(0,255)
    refactor = vtk.vtkActor()
    refactor.GetProperty().SetOpacity(0.7)
    refactor.SetMapper(refmapper)


    # Remove existing lights and add lightkit lights
    ren.RemoveAllLights()
    lightkit = vtk.vtkLightKit()
    lightkit.AddLightsToRenderer(ren)

    # assign actors to the renderer
    ren.AddActor(refactor)
    ren.AddActor(surfaceactor)
    ren.AddActor(txt)

    # set the background and size; zoom in; and render
    ren.SetBackground(1, 1, 1)
    renWin.SetSize(1280, 960)
    ren.ResetCamera()
    ren.GetActiveCamera().Zoom(1)

    # enable user interface interactor
    iren.Initialize()
    renWin.Render()
    iren.Start()

    outcam = ren.GetActiveCamera()

def visualise_default_continuous(surface,overlay,case,arrayname,edges=0,flip=0,colormap='BlYl',
    interact=1,filename='./screenshot.png',LegendTitle='',mini='',maxi='',mag=2):
    """Visualise surface with a continuos colormap according to 'arrayname'."""

    # surface mapper and actor
    surfacemapper = vtk.vtkPolyDataMapper()
    surfacemapper.SetInputData(surface)
    surfacemapper.SelectColorArray(arrayname)

    #Create a lookup table to map cell data to colors
    if surface.GetPointData().GetArray(arrayname):
        array = vtk_to_numpy(surface.GetPointData().GetArray(arrayname))
        surfacemapper.SetScalarModeToUsePointFieldData()
    else:
        array = vtk_to_numpy(surface.GetCellData().GetArray(arrayname))
        surfacemapper.SetScalarModeToUseCellFieldData()


    if not maxi:
        maxi = np.nanmax(array)

    if not mini:
        # mini might be zero, and is true
        if mini == 0:
            mini = mini
        else:
            mini = np.nanmin(array)

    colors = getcolors(colormap)
    numcolors = int(len(colors))

    if colormap == '24_regions':
        lut = vtk.vtkLookupTable()
        lut.SetNumberOfTableValues(24)

        # don't interpolate
        for c in range(len(colors)):
            this = colors[c]
            lut.SetTableValue(this[0], this[1], this[2], this[3])
        lut.Build()
        surfacemapper.SetLookupTable(lut)
        surfacemapper.SetScalarRange(1,24)
        surfacemapper.InterpolateScalarsBeforeMappingOn()
    else:
        lut = vtk.vtkColorTransferFunction()
        lut.SetColorSpaceToHSV()
        lut.SetNanColor(0.5,0.5,0.5)
        for c in range(len(colors)):
                cmin = colors[0][0]
                cmax = colors[numcolors-1][0]
                this = colors[c]
                rat = (this[0] - cmin) / (cmax - cmin)
                t = (maxi - mini) * rat + mini
                lut.AddRGBPoint(t, this[1], this[2], this[3])
        lut.Build()
        surfacemapper.SetLookupTable(lut)
        surfacemapper.SetScalarRange(mini,maxi)
        surfacemapper.InterpolateScalarsBeforeMappingOn()

    # create a text actor
    txt = vtk.vtkTextActor()
    txt.SetInput(case)
    txtprop=txt.GetTextProperty()
    txtprop.SetFontFamilyToArial()
    txtprop.SetFontSize(30)
    txtprop.SetColor(0, 0, 0)
    txt.SetDisplayPosition(20, 30)

    # create a rendering window, renderer, and renderwindowinteractor
    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    iren = vtk.vtkRenderWindowInteractor()
    style = vtk.vtkInteractorStyleTrackballCamera()
    iren.SetInteractorStyle(style)
    iren.SetRenderWindow(renWin)

    surfaceactor = vtk.vtkActor()
    surfaceactor.SetMapper(surfacemapper)

    # overlay mapper and actor
    refmapper = vtk.vtkPolyDataMapper()
    refmapper.SetInputData(overlay)
    refmapper.SetScalarModeToUsePointFieldData()

    if edges == 1:
        refactor = vtk.vtkActor()
        refactor.GetProperty().SetOpacity(1)
        if colormap == 'erdc_rainbow_grey':
            refactor.GetProperty().SetColor(1,1,1)
        else:
            refactor.GetProperty().SetColor(0,0,0)
        refactor.GetProperty().SetRepresentationToWireframe()
        refactor.GetProperty().SetLineWidth(6.)
        refactor.SetMapper(refmapper)
    else:
        refactor = vtk.vtkActor()
        refactor.GetProperty().SetOpacity(0.5)
        refactor.GetProperty().SetColor(1, 0, 0)
        refactor.SetMapper(refmapper)

    ren.AddActor(surfaceactor)
    ren.AddActor(refactor)

    if LegendTitle:
        ScalarBarActor = vtk.vtkScalarBarActor()
        # colormap
        ScalarBarActor.SetLookupTable(surfaceactor.GetMapper().GetLookupTable())

        #labels format
        ScalarBarActor.GetLabelTextProperty().ItalicOff()
        ScalarBarActor.GetLabelTextProperty().BoldOn()
        ScalarBarActor.GetLabelTextProperty().ShadowOff()
        ScalarBarActor.GetLabelTextProperty().SetFontFamilyToArial()
        ScalarBarActor.GetLabelTextProperty().SetFontSize(100)
        ScalarBarActor.GetLabelTextProperty().SetColor(0.,0.,0.)
        ScalarBarActor.SetLabelFormat('%.2f')

        if colormap == '24_regions':
            ScalarBarActor.GetLabelTextProperty().BoldOff()
            ScalarBarActor.GetLabelTextProperty().SetFontSize(100)
            ScalarBarActor.SetNumberOfLabels(24)
            ScalarBarActor.SetLabelFormat('%.0f')

        # orientation
        ScalarBarActor.SetMaximumWidthInPixels(175)
        ScalarBarActor.SetPosition(0.90,0.15)
        ren.AddActor(ScalarBarActor)
    else:
        ren.AddActor(txt)


    # set the background and size; zoom in; and render
    ren.SetBackground(1, 1, 1)
    renWin.SetSize(875, 800)
    ren.ResetCamera()


    if flip ==1:
        # flip to foot to head position
        # default values from paraview
        aCamera = vtk.vtkCamera()
        aCamera.SetViewUp(0, 1, 0)
        aCamera.SetPosition(0,0,-3.17)
        aCamera.SetFocalPoint(0,0,0)
        aCamera.SetClippingRange(3.14,3.22)
        aCamera.SetParallelScale(0.83)
        ren.SetActiveCamera(aCamera)

    ren.GetActiveCamera().Zoom(1.4)
    iren.Initialize()
    renWin.Render()

    # Remove existing lights and add lightkit lights
    ren.RemoveAllLights()
    lightkit = vtk.vtkLightKit()
    lightkit.AddLightsToRenderer(ren)

    # enable user interface interactor
    if interact ==1:

        iren.Start()
    else:
        # save as png
        ## Screenshot
        windowToImageFilter = vtk.vtkWindowToImageFilter()

        windowToImageFilter.SetInput(renWin)
        windowToImageFilter.SetMagnification(mag) #set the resolution of the output image (3 times the current resolution of vtk render window)
        windowToImageFilter.SetInputBufferTypeToRGBA() #also record the alpha (transparency) channel
        windowToImageFilter.Update()

        pngwriter = vtk.vtkPNGWriter()
        pngwriter.SetFileName(filename)
        pngwriter.SetInputConnection(windowToImageFilter.GetOutputPort())
        pngwriter.Write()


def visualise_nan_glyph(surface,glyph,overlay,case,arrayname,edges=0,flip=0,colormap='BlYl',
    interact=1,filename='./screenshot.png',LegendTitle='',mini='',maxi=''):
    """Visualise surface with glyphs colormapped according to 'arrayname'."""

    # glyph mapper
    glyphmapper = vtk.vtkPolyDataMapper()
    glyphmapper.SetInputData(glyph)
    glyphmapper.SelectColorArray(arrayname)


    # only pointdata
    array = vtk_to_numpy(surface.GetPointData().GetArray(arrayname))
    glyphmapper.SetScalarModeToUsePointFieldData()

    if not maxi:
        maxi = np.nanmax(array)

    if not mini:
        # mini might be zero, and is true
        if mini == 0:
            mini = mini
        else:
            mini = np.nanmin(array)

    colors = getcolors(colormap)
    numcolors = int(len(colors))

    if colormap == '24_regions':
        lut = vtk.vtkLookupTable()
        lut.SetNumberOfTableValues(24)

        # don't interpolate
        for c in range(len(colors)):
            this = colors[c]
            lut.SetTableValue(this[0], this[1], this[2], this[3])
        lut.Build()
        glyphmapper.SetLookupTable(lut)
        glyphmapper.SetScalarRange(1,24)
        glyphmapper.InterpolateScalarsBeforeMappingOn()
    else:
        lut = vtk.vtkColorTransferFunction()
        lut.SetColorSpaceToHSV()
        lut.SetNanColor(0.5,0.5,0.5)
        for c in range(len(colors)):
                cmin = colors[0][0]
                cmax = colors[numcolors-1][0]
                this = colors[c]
                rat = (this[0] - cmin) / (cmax - cmin)
                t = (maxi - mini) * rat + mini
                lut.AddRGBPoint(t, this[1], this[2], this[3])
        lut.Build()
        glyphmapper.SetLookupTable(lut)
        glyphmapper.SetScalarRange(mini,maxi)
        glyphmapper.InterpolateScalarsBeforeMappingOn()

    # create a text actor
    txt = vtk.vtkTextActor()
    txt.SetInput(case)
    txtprop=txt.GetTextProperty()
    txtprop.SetFontFamilyToArial()
    txtprop.SetFontSize(30)
    txtprop.SetColor(0, 0, 0)
    txt.SetDisplayPosition(20, 30)

    # create a rendering window, renderer, and renderwindowinteractor
    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    iren = vtk.vtkRenderWindowInteractor()
    style = vtk.vtkInteractorStyleTrackballCamera()
    iren.SetInteractorStyle(style)
    iren.SetRenderWindow(renWin)

    # surface mapper
    surfacemapper = vtk.vtkPolyDataMapper()
    surfacemapper.SetInputData(surface)
    surfaceactor = vtk.vtkActor()
    surfaceactor.SetMapper(surfacemapper)
    surfaceactor.GetProperty().SetOpacity(0.5)
    surfaceactor.GetProperty().SetColor(0.5,0.5,0.5)

    # glyph in scalar colors
    glyphactor = vtk.vtkActor()
    glyphactor.SetMapper(glyphmapper)

    # overlay mapper and actor
    refmapper = vtk.vtkPolyDataMapper()
    refmapper.SetInputData(overlay)
    refmapper.SetScalarModeToUsePointFieldData()

    if edges == 1:
        refactor = vtk.vtkActor()
        refactor.GetProperty().SetOpacity(1)
        refactor.GetProperty().SetColor(0, 0, 0)
        refactor.GetProperty().SetRepresentationToWireframe()
        refactor.GetProperty().SetLineWidth(6.)
        refactor.SetMapper(refmapper)
    else:
        refactor = vtk.vtkActor()
        refactor.GetProperty().SetOpacity(0.5)
        refactor.GetProperty().SetColor(1, 0, 0)
        refactor.SetMapper(refmapper)

    ren.AddActor(surfaceactor)
    ren.AddActor(refactor)
    ren.AddActor(glyphactor)

    if LegendTitle:
        ScalarBarActor = vtk.vtkScalarBarActor()
        # colormap
        ScalarBarActor.SetLookupTable(glyphactor.GetMapper().GetLookupTable())

        #labels format
        ScalarBarActor.GetLabelTextProperty().ItalicOff()
        ScalarBarActor.GetLabelTextProperty().BoldOn()
        ScalarBarActor.GetLabelTextProperty().ShadowOff()
        ScalarBarActor.GetLabelTextProperty().SetFontFamilyToArial()
        ScalarBarActor.GetLabelTextProperty().SetFontSize(100)
        ScalarBarActor.GetLabelTextProperty().SetColor(0.,0.,0.)
        ScalarBarActor.SetLabelFormat('%.2f')

        if colormap == '24_regions':
            ScalarBarActor.GetLabelTextProperty().BoldOff()
            ScalarBarActor.GetLabelTextProperty().SetFontSize(100)
            ScalarBarActor.SetNumberOfLabels(24)
            ScalarBarActor.SetLabelFormat('%.0f')

        # orientation
        ScalarBarActor.SetMaximumWidthInPixels(200)
        ScalarBarActor.SetPosition(0.90,0.15)
        ren.AddActor(ScalarBarActor)
    else:
        ren.AddActor(txt)


    # set the background and size; zoom in; and render
    ren.SetBackground(1, 1, 1)
    renWin.SetSize(875, 800)
    ren.ResetCamera()


    if flip ==1:
        # default values from paraview
        aCamera = vtk.vtkCamera()
        aCamera.SetViewUp(0, 1, 0)
        aCamera.SetPosition(0,0,-3.17)
        aCamera.SetFocalPoint(0,0,0)
        aCamera.SetClippingRange(3.14,3.22)
        aCamera.SetParallelScale(0.83)
        ren.SetActiveCamera(aCamera)

    ren.GetActiveCamera().Zoom(1.4)
    iren.Initialize()
    renWin.Render()

    # Remove existing lights and add lightkit lights
    ren.RemoveAllLights()
    lightkit = vtk.vtkLightKit()
    lightkit.AddLightsToRenderer(ren)

    # enable user interface interactor
    if interact ==1:

        iren.Start()
    else:
        # save as png
        ## Screenshot
        windowToImageFilter = vtk.vtkWindowToImageFilter()

        windowToImageFilter.SetInput(renWin)
        windowToImageFilter.SetMagnification(3) #set the resolution of the output image (3 times the current resolution of vtk render window)
        windowToImageFilter.SetInputBufferTypeToRGBA() #also record the alpha (transparency) channel
        windowToImageFilter.Update()

        pngwriter = vtk.vtkPNGWriter()
        pngwriter.SetFileName(filename)
        pngwriter.SetInputConnection(windowToImageFilter.GetOutputPort())
        pngwriter.Write()

def visualise_color(surface,ref,case):
    """Visualise surface in solid color and 'ref' in trasparent."""
    # create a text actor
    txt = vtk.vtkTextActor()
    txt.SetInput(case)
    txtprop=txt.GetTextProperty()
    txtprop.SetFontFamilyToArial()
    txtprop.SetFontSize(18)
    txtprop.SetColor(0, 0, 0)
    txt.SetDisplayPosition(20, 30)

    # create a rendering window, renderer, and renderwindowinteractor
    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    iren = vtk.vtkRenderWindowInteractor()
    style = vtk.vtkInteractorStyleTrackballCamera()
    iren.SetInteractorStyle(style)
    iren.SetRenderWindow(renWin)

    # surface mapper and actor
    surfacemapper = vtk.vtkPolyDataMapper()
    surfacemapper.SetInputData(surface)
    surfacemapper.SetScalarModeToUsePointFieldData()
    surfaceactor = vtk.vtkActor()
    surfaceactor.GetProperty().SetColor(288/255, 26/255, 28/255)
    surfaceactor.SetMapper(surfacemapper)

    # refsurface mapper and actor
    refmapper = vtk.vtkPolyDataMapper()
    refmapper.SetInputData(ref)
    refmapper.SetScalarModeToUsePointFieldData()

    refactor = vtk.vtkActor()
    refactor.GetProperty().SetOpacity(0.5)
    refactor.GetProperty().SetColor(1, 1, 1)
    refactor.SetMapper(refmapper)


    # assign actors to the renderer
    ren.AddActor(surfaceactor)
    ren.AddActor(refactor)
    ren.AddActor(txt)

    # set the background and size; zoom in; and render
    ren.SetBackground(1, 1, 1)
    renWin.SetSize(800 , 800)
    ren.ResetCamera()
    ren.GetActiveCamera().Zoom(1)

    # enable user interface interactor
    iren.Initialize()
    renWin.Render()
    iren.Start()

#------------------------------------------------------------------------------
# Input/Output
#------------------------------------------------------------------------------

def readvtp(filename, dataarrays=True):
    """Read polydata in vtp format."""
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(filename)
    reader.Update()
    if dataarrays == False:
        for i in range(reader.GetNumberOfPointArrays()):
            arrayname = reader.GetPointArrayName(i)
            reader.SetPointArrayStatus(arrayname, 0)
        for i in range(reader.GetNumberOfCellArrays()):
            arrayname = reader.GetCellArrayName(i)
            reader.SetPointArrayStatus(arrayname, 0)
        reader.Update()
    return reader.GetOutput()

def readpolydatavtk(filename, dataarrays=True):
    """Read polydata in vtk format."""
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(filename)
    reader.Update()
    return reader.GetOutput()

def set_cell_array_value(polydata,name,value):
    """Set all components of data array to 'value'."""
    # get original array
    array = polydata.GetCellData().GetArray(name)

    # round labels
    for i in range(polydata.GetNumberOfCells()):
        array.SetValue(i,value)
    polydata.GetPoints().Modified()

    return polydata

def set_cell_array_value_per_label(polydata,name,value,labelname,label):
    """Set all components of data array to 'value' if they belong to a label."""
    # get original array
    array = polydata.GetCellData().GetArray(name)
    labelarray = polydata.GetCellData().GetArray(labelname)

    # round labels
    for i in range(polydata.GetNumberOfCells()):
        if labelarray.GetValue(i) == label:
            array.SetValue(i,value)
    polydata.GetPoints().Modified()

    return polydata


def vtk2vtp(inputfile, outputfile):
    """Read a vtk polydata and save as vtp."""
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(inputfile)
    reader.Update()
    surface = reader.GetOutput()
    writevtp(surface,outputfile)

def ply2vtk(filename1,filename2):
    """Read a ply file and save as vtk ascii."""
    reader = vtk.vtkPLYReader()
    reader.SetFileName(filename1)
    reader.Update()

    writer = vtk.vtkPolyDataWriter()
    writer.SetInputConnection(reader.GetOutputPort())
    writer.SetFileTypeToASCII()
    writer.SetFileName(filename2)
    writer.Write()

def writeply(surface,filename):
    """Write mesh as ply file."""
    writer = vtk.vtkPLYWriter()
    writer.SetInputData(surface)
    writer.SetFileTypeToASCII()
    writer.SetFileName(filename)
    writer.Write()

def writevtk(surface,filename):
    """Write vtk polydata file."""
    writer = vtk.vtkPolyDataWriter()
    writer.SetInputData(surface)
    writer.SetFileTypeToASCII()
    writer.SetFileName(filename)
    writer.Write()

def writevtp(surface, filename):
    """Write vtp polydata file."""
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetInputData(surface)
    writer.SetFileName(filename)
    writer.Write()


#------------------------------------------------------------------------------
# VMTK functions
#------------------------------------------------------------------------------

def vmtksurfacereader(filename):
    reader = vmtkscripts.vmtkSurfaceReader()
    reader.InputFileName = filename
    reader.Execute()
    return reader.Surface


def vmtkcenterlineresampling(centerline, length=.1):
    resampler = vmtkscripts.vmtkCenterlineResampling()
    resampler.Centerlines = centerline
    resampler.Length = length
    resampler.Execute()
    return resampler.Centerlines


def vmtkcenterlinesmoothing(centerline, iterations=100, factor=0.1):
    smoother = vmtkscripts.vmtkCenterlineSmoothing()
    smoother.Centerlines = centerline
    smoother.NumberOfSmoothingIterations = iterations
    smoother.SmoothingFactor = factor
    smoother.Execute()
    return smoother.Centerlines

def vmtkbranchextractor(centerline):
    extractor = vmtkscripts.vmtkBranchExtractor()
    extractor.Centerlines = centerline
    extractor.RadiusArrayName = 'MaximumInscribedSphereRadius'
    extractor.Execute()
    return extractor.Centerlines

def vmtksurfacewriter(polydata, filename):
    writer = vmtkscripts.vmtkSurfaceWriter()
    writer.Surface = polydata
    writer.OutputFileName = filename
    writer.Execute()

def vmtkcenterlinemerge(centerline, length=.1):
    merger = vmtkscripts.vmtkCenterlineMerge()
    merger.Centerlines = centerline
    merger.Length = length
    merger.RadiusArrayName = 'MaximumInscribedSphereRadius'
    merger.GroupIdsArrayName = 'GroupIds'
    merger.CenterlineIdsArrayName = 'CenterlineIds'
    merger.BlankingArrayName = 'Blanking'
    merger.TractIdsArrayName = 'TractIds'
    merger.Execute()
    return merger.Centerlines

def vmtkcenterlineattributes(centerline):
    computer = vmtkscripts.vmtkCenterlineAttributes()
    computer.Centerlines = centerline
    computer.Execute()
    return computer.Centerlines

def vmtkcenterlinesections(surface, centerline):
    sectioner = vmtkscripts.vmtkCenterlineSections()
    sectioner.Surface = surface
    sectioner.Centerlines = centerline
    sectioner.Execute()
    return sectioner.CenterlineSections

def vmtkcenterlines(surface, sourcepoints, targetpoints,endpoints=0):
    computer = vmtkscripts.vmtkCenterlines()
    computer.Surface = surface
    computer.SeedSelectorName = 'pointlist'
    computer.SourcePoints = sourcepoints
    computer.TargetPoints = targetpoints
    computer.AppendEndPoints = endpoints

    computer.Execute()
    return computer.Centerlines

def vmtksurfacecapper(surface,method='centerpoint',nrings=4,const=0.1, interactive=0):
    capper = vmtkscripts.vmtkSurfaceCapper()
    capper.Surface = surface
    capper.Method = method
    capper.NumberOfRings = nrings
    capper.ConstraintFactor = const
    capper.Interactive = interactive
    capper.Execute()
    return capper.Surface

def vmtksurfaceclipper(surface):
    computer = vmtkscripts.vmtkSurfaceClipper()
    computer.Surface = surface
    computer.Execute()
    return computer.Surface
