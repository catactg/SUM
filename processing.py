import os
import sys
import shutil
from SUM_utils import *
from constants import MATLAB_BIN_PATH, CURRENTS_BUILD_PATH, MESHLABSERVER_PATH

def run_standardization(meshfile, atlaspath, dist, maxslope,
            clspacing, skippointsfactor, highslope, bumpcriterion, pvends,
            use_seed_selector=False, use_laa_seed=False, visualise=True):
    """Runs standardization on original mesh.
    Saves output at the location of the input meshfile.

    Arguments:
        datapath            =   full path to mesh file
        atlaspath           =   root path of average mesh and unfold disk
        dist                =   length of PVs to keep
        maxslope            =   anything above this is ostium
        clspacing           =   resample the centerline with this spacing
        skippointsfactor    =   percentage of points to ignore at begining of centerline
        highslope           =   above this slope we start counting as potential ostium location
        bumpcriterion       =   ostium if slope higher than highslope and above bump criterion
        pvends              =   enabe PV ends in vmtkcenterlines
        use_laa_seed        =   enable seed for appendage (True|False)
        use_seed_selector   =   enable selection of seeds (True|False)
        visualise           =   enable visualisation of intermin results (True|False)
    """

    fileroot = os.path.dirname(meshfile)
    filenameroot = os.path.splitext(os.path.basename(meshfile))[0]

    inputfile = os.path.join(fileroot, filenameroot + '.vtk')
    outfile = os.path.join(fileroot, filenameroot + '.vtp')
    vtk2vtp(inputfile, outfile)

    # basic VMTK cleaning
    inputfile = os.path.join(fileroot, filenameroot + '.vtp')
    surface = vmtksurfacereader(inputfile)
    outfile = os.path.join(fileroot, filenameroot + 'kite.vtp')
    surface = vmtksurfacekiteremoval(surface)
    vmtksurfacewriter(surface, outfile)
    surface = cleanpolydata(surface)

    if visualise:
        visualise_color(surface, surface, 'original')

    # save surface as ply
    outfile = os.path.join(fileroot, filenameroot + '.ply')
    writeply(surface, outfile)

    inputfile = os.path.join(fileroot, filenameroot + '.ply')
    outfile = os.path.join(fileroot, filenameroot + 'poisson.ply')
    currentscript = './poissonOT8_QEC.mlx'
    print currentscript
    os.system(MESHLABSERVER_PATH + ' -i ' + inputfile + ' -o ' + outfile +
            ' -s ' + currentscript)

    inputfile = os.path.join(fileroot, filenameroot + 'poisson.ply')
    outfile = os.path.join(fileroot, filenameroot + 'poisson.vtk')
    ply2vtk(inputfile, outfile)

    # prepare files
    inputfile = os.path.join(fileroot, filenameroot + 'poisson.vtk')
    inputsurface = readpolydatavtk(inputfile)
    inputfile = os.path.join(fileroot, filenameroot + 'poisson.vtp')
    writevtp(inputsurface, inputfile)

    # compute seeds
    inputfile = os.path.join(fileroot, filenameroot + 'poisson.vtk')
    outputfile = os.path.join(fileroot, filenameroot + 'seeds.vtp')
    surface = readpolydatavtk(inputfile)

    if not os.path.exists(outputfile) or use_seed_selector:
        select_seeds(surface,
                     'GTLabels',
                     outputfile,
                     visualise,
                     use_laa_seed)

    # prepare seeds
    seedsfile = os.path.join(fileroot, filenameroot + 'seeds.vtp')
    outseedsfile = os.path.join(fileroot, filenameroot + 'seeds.csv')

    if use_laa_seed:
        # 36 is laa
        multiple_seeds_to_csv_no_mitral(seedsfile,
                                        'GTLabels',
                                        [77, 76, 78, 79, 36],
                                        outseedsfile)
        # new pv centerlines no mitral
        seedsfile = os.path.join(fileroot,filenameroot + 'seeds.csv')
        outfile = os.path.join(fileroot,filenameroot + '_')
        pv_centerlines_no_mitral_laa(inputfile,
                                     seedsfile,
                                     outfile,
                                     pvends)

        # other settings for clipping
        specialvein = 37
        specialdistance = 2.0

    else:
        # no 36 because LAA is not used for centerlines
        multiple_seeds_to_csv_no_mitral(seedsfile,
                                        'GTLabels',
                                        [77, 76, 78, 79],
                                        outseedsfile)
        # new pv centerlines no mitral
        seedsfile = os.path.join(fileroot, filenameroot + 'seeds.csv')
        outfile = os.path.join(fileroot, filenameroot + '_')
        pv_centerlines_no_mitral(inputfile,
                                 seedsfile,
                                 outfile,
                                 pvends)

        # other settings for clipping
        specialvein = 0 # disabled
        specialdistance = 0 # disabled


    # label PVs automatically
    inputfile = os.path.join(fileroot, filenameroot + 'poisson.vtp')
    outfile = os.path.join(fileroot, filenameroot + '_')
    clip_veins_sections(inputfile,
                        outfile,
                        clspacing,
                        maxslope,
                        skippointsfactor,
                        highslope,
                        bumpcriterion,
                        use_laa_seed,
                        visualise)

    # clip PV end points
    sufixfile = os.path.join(fileroot, filenameroot + '_')
    inputfile = os.path.join(fileroot, filenameroot + '_autolabels.vtp')
    inputsurface = readvtp(inputfile)
    stdmesh = clip_vein_endpoint(inputsurface,
                                 sufixfile,
                                 dist,
                                 use_laa_seed,
                                 specialvein,
                                 specialdistance)


    # save
    o_file = os.path.join(fileroot, filenameroot + '_clipped.vtk')
    writevtk(stdmesh, o_file)

    if visualise:
        origmesh = readpolydatavtk(os.path.join(fileroot, filenameroot +  '.vtk'))
        visualise_default(stdmesh, origmesh, 'standard mesh', 'autolabels', 36, 79)

def run_currents(meshfile, atlaspath, mitralcliptype, pvcliptype,
                 use_similarity = False, use_laa_seed = False, visualise=True):
    """
    Performs elastic registration between two surfaces

    Arguments:
        meshfile            =   full path to mesh file
        atlaspath           =   root path of average mesh and unfold disk
        mitralcliptype      =   defines type of mitral clip (manual|auto)
        pvcliptype          =   type of PV clip (short|long)
        use_similarity      =   use similarity transform to initialise mesh registration (True|False)
        use_laa_seed        =   enable seed for appendage (True|False)
        visualise           =   enable visualisation of intermin results (True|False)

    """

    exe_root= os.path.join(MATLAB_BIN_PATH + ' -nodesktop -nosplash -r ')
    exe_matlab = exe_root + '"cd ' + CURRENTS_BUILD_PATH + '; match2vtks('

    fileroot = os.path.dirname(meshfile)
    filenameroot = os.path.splitext(os.path.basename(meshfile))[0]

    # mitral clip
    i_file = os.path.join(fileroot, filenameroot + '.vtk')
    surfaceorig = readpolydatavtk(i_file)
    s_file = os.path.join(fileroot, filenameroot + '_clipped.vtk')
    surface = readpolydatavtk(s_file)

    # type of mitral clip
    if mitralcliptype == 'manual':
        print "Previously defined mitral edge"
        surfaceclipped = cylinder_mitral_clip(surface, surfaceorig, visualise)

    if mitralcliptype == 'auto':
        print "Auto clipping plane"
        w=[0.95,0.05,0.0]
        o_file = os.path.join(fileroot, filenameroot + '_clipped_mitral')

        if use_laa_seed:
            # if we labeled the LAA, we can use a plane for mitral clipping
            surfaceclipped = find_mitral_plane_pvs(surface, 'autolabels', o_file, 0.35, w, 0)
        else:
            # otherwise, a cylinder which is contained and will preserve the LAA structure
            surfaceclipped = find_mitral_cylinder_pvs(surface, 'autolabels', o_file, 0.35, w, 0)

    o_file = os.path.join(fileroot, filenameroot + '_clipped_mitral.vtk')
    writevtk(surfaceclipped, o_file)

    if visualise:
        visualise_color(surfaceclipped, surface, 'mitral clip')

    # using clipped for lmk selection
    i_file = os.path.join(fileroot, filenameroot + '_clipped_mitral.vtk')
    source = readpolydatavtk(i_file)

    if pvcliptype == 'short':
        avg_file = os.path.join(atlaspath, 'affine_average_clipped_shortpvs.vtp')
    else:
        avg_file = os.path.join(atlaspath, 'affine_average_clipped_longpvs.vtp')
    target = readvtp(avg_file)

    # intialise
    deformed = initial_transform_pvends(source,
                                        target,
                                        'autolabels',
                                        use_similarity,
                                        use_laa_seed)
    o_file = os.path.join(fileroot, filenameroot +'_pvends_mitral_init.vtk')
    writevtk(deformed, o_file)

    if visualise:
        visualise_color(deformed, target, 'initialisation')

    if pvcliptype == 'short':
        avg_file = os.path.abspath(os.path.join(atlaspath, 'affine_average_clipped_shortpvs.vtk'))
    else:
        avg_file = os.path.abspath(os.path.join(atlaspath, 'affine_average_clipped_longpvs.vtk'))

    i_file = os.path.abspath(os.path.join(fileroot, filenameroot +'_pvends_mitral_init.vtk'))
    o_file = os.path.abspath(os.path.join(fileroot, filenameroot +'_pvends_mitral_init_currents'))

    # run currents
    line = exe_matlab + '\'' + i_file + '\''+ ','
    line = line + '\'' + avg_file + '\''+ ',' + '\'' + o_file + '\''+ ',' + '\''
    line = line + '1' + '\''+ ',' + '\'0.0001\'); quit;"'
    print line
    os.system(line)

def run_sum(meshfile, atlaspath, pvcliptype, glyphon, restorenans,
            arraysource, arraytarget, colormap, value_range, visualise= True):
    """
    Computes standardized unfold map
    Arguments:
            meshfile    =   full path to mesh file
            atlaspath   =   root path of average mesh and unfold disk
            pvcliptype  =   type of PV clip (short|long)
            glyphon     =   enalbe glyph for point_by_point EP measurements (True|False)
            restorenans =   restore nan values if previously set to a different value (True|False)
            arraysource =   original scalar arrayname in scalarfile
            arraytarget =   new scalar arrayname for output file
            colormap    =   colormap to use for PNG screen shot
            value_range =   pair with min and max values for colormap scaling
            visualise   =   enables visualition of intermin results (True|False)
    """

    fileroot = os.path.dirname(meshfile)
    filenameroot = os.path.splitext(os.path.basename(meshfile))[0]
    outprefix = os.path.join(fileroot,filenameroot + '_' + arraytarget)

    if not glyphon:
        # Reload scalar values on stdmeshes
        i_file = os.path.join(fileroot, filenameroot + '_clipped_mitral.vtk')
        o_file = os.path.join(fileroot, filenameroot + '_clipped_mitral_scalars.vtp')
        s_file = meshfile

        target = readpolydatavtk(i_file)
        source = readpolydatavtk(s_file)

        surfproj = vmtksurfaceprojection(target, source)
        writevtp(surfproj, o_file)

        # Transfer scalar values on registered meshes by pointid
        i_file = os.path.join(fileroot, filenameroot + '_pvends_mitral_init_currents.vtk')
        o_file = os.path.join(fileroot, filenameroot + '_pvends_mitral_init_scalars.vtp')
        s_file = os.path.join(fileroot, filenameroot + '_clipped_mitral_scalars.vtp')

        surface = readpolydatavtk(i_file)
        source = readvtp(s_file)

        targetle = transfer_array_by_pointid(source, surface, arraysource, arraytarget)
        writevtp(targetle, o_file)

        # Transfer scalar values from registered meshes to average mesh
        if pvcliptype == 'short':
            i_file = os.path.join(atlaspath, 'affine_average_clipped_shortpvs.vtp')
        else:
            i_file = os.path.join(atlaspath, 'affine_average_clipped_longpvs.vtp')

        o_file = os.path.join(fileroot, filenameroot + '_pvends_mitral_init_on_average.vtp')
        s_file = os.path.join(fileroot, filenameroot + '_pvends_mitral_init_scalars.vtp')

        surface = readvtp(i_file)
        source = readvtp(s_file)

        # any value < 0 will be truncated
        sourcenolb = zero_truncate_array(source, arraytarget, 0.)
        surfproj = vmtksurfaceprojection(surface, sourcenolb)
        writevtp(surfproj, o_file)

        if visualise:
            visualise_color(source, surface, 'currents')

        # Transfer scalar values from average mesh to "unfold" mesh
        i_file = os.path.join(atlaspath, 'affine_average_clipped_shortpvs.vtp')
        o_file = os.path.join(fileroot, filenameroot + '_pvends_mitral_init_on_average_unfold.vtp')
        s_file = os.path.join(fileroot, filenameroot + '_pvends_mitral_init_on_average.vtp')

        surface = readvtp(i_file)
        source = readvtp(s_file)

        surfproj = vmtksurfaceprojection(surface, source)
        writevtp(surfproj, o_file)

        # Transfer scalar values from "unfold" mesh to disk by pointid
        i_file = os.path.join(atlaspath, 'disk_basedon3d_regions_new_origin.vtp')
        o_file = (outprefix + '_disk_basedon3d.vtp')
        s_file = os.path.join(fileroot, filenameroot + '_pvends_mitral_init_on_average_unfold.vtp')

        surface = readvtp(i_file)
        source = readvtp(s_file)
        targetle = transfer_array_by_pointid(source, surface, arraytarget, arraytarget)
        writevtp(targetle, o_file)

        i_file = os.path.join(atlaspath, 'disk_uniform_labels_pvs.vtp')
        o_file = (outprefix + '_disk_uniform.vtp')
        s_file = (outprefix + '_disk_basedon3d.vtp')

        surface = readvtp(i_file)
        source = readvtp(s_file)
        targetle = vmtksurfaceprojection(surface, source)
        writevtp(targetle, o_file)

        o_file = (outprefix + '_disk_uniform.png')

        i_file = os.path.join(atlaspath, 'disk_uniform_edges.vtp')
        overlay = readvtp(i_file)

        visualise_default_continuous(targetle,
                                     overlay,
                                     arraytarget,
                                     arraytarget,
                                     1,
                                     1,
                                     colormap,
                                     0,
                                     o_file,
                                     arraytarget,
                                     value_range[0],
                                     value_range[1])

    if glyphon:
        # non continous scalar values
        i_file = os.path.join(fileroot, filenameroot + '_clipped_mitral.vtk')
        o_file = os.path.join(fileroot, filenameroot + '_clipped_mitral_scalars.vtp')
        s_file = meshfile

        surface = readpolydatavtk(i_file)
        source = readpolydatavtk(s_file)

        # since the values are single points, they are easily lost in the vmtkprojection
        # solution: visit every non nan value and transfer it to the closest point
        if visualise:
            visualise_default(surface, surface, arraytarget, 'autolabels', 36, 79)

        if restorenans:
            source = restore_nan_values(source,  arraytarget,  0.0001)

        point_map = point_map_dictionary()

        surface = project_nan_array(
                                    source,
                                    surface,
                                    arraytarget,
                                    point_map,
                                    'original2std'
                                    )
        writevtp(surface, o_file)

        # Transfer scalar values on registered meshes by pointid
        i_file = os.path.join(fileroot, filenameroot + '_pvends_mitral_init_currents.vtk')
        o_file = os.path.join(fileroot, filenameroot + '_pvends_mitral_init_scalars.vtp')
        s_file = os.path.join(fileroot, filenameroot + '_clipped_mitral_scalars.vtp')

        surface = readpolydatavtk(i_file)
        source = readvtp(s_file)

        if visualise:
            visualise_default(surface, surface, arraytarget, 'autolabels', 36, 79)

        surface = transfer_array_by_pointid(source, surface, arraytarget, arraytarget)
        writevtp(surface, o_file)

        # Transfer scalar values from registered meshes to average mesh
        if pvcliptype == 'short':
            i_file = os.path.join(atlaspath, 'affine_average_clipped_shortpvs.vtp')
        else:
            i_file = os.path.join(atlaspath, 'affine_average_clipped_longpvs.vtp')
        o_file = os.path.join(fileroot, filenameroot + '_pvends_mitral_init_on_average.vtp')
        s_file = os.path.join(fileroot, filenameroot + '_pvends_mitral_init_scalars.vtp')

        surface = readvtp(i_file)
        source = readvtp(s_file)
        surface = project_nan_array(
                                    source,
                                    surface,
                                    arraytarget,
                                    point_map,
                                    'std2average'
                                    )
        writevtp(surface, o_file)


        # Transfer scalar values from average mesh to "unfold" mesh
        i_file = os.path.join(atlaspath, 'affine_average_clipped_shortpvs.vtp')
        o_file = os.path.join(fileroot, filenameroot + '_pvends_mitral_init_on_average_unfold.vtp')
        s_file = os.path.join(fileroot, filenameroot + '_pvends_mitral_init_on_average.vtp')

        surface = readvtp(i_file)
        source = readvtp(s_file)

        surface = project_nan_array(
                                source,
                                surface,
                                arraytarget,
                                point_map,
                                'average2unfold'
                                )
        writevtp(surface, o_file)


        # Transfer scalar values from "unfold" mesh to disk by pointid
        i_file = os.path.join(atlaspath, 'disk_basedon3d_regions_new_origin.vtp')
        o_file = (outprefix + '_disk_basedon3d.vtp')
        s_file = os.path.join(fileroot, filenameroot + '_pvends_mitral_init_on_average_unfold.vtp')

        surface = readvtp(i_file)
        source = readvtp(s_file)

        surface = transfer_array_by_pointid(source, surface, arraytarget, arraytarget)
        writevtp(surface, o_file)

        i_file = os.path.join(atlaspath, 'disk_uniform_labels_pvs.vtp')
        s_file = (outprefix + '_disk_basedon3d.vtp')

        surface = readvtp(i_file)
        source = readvtp(s_file)

        # region edges
        i_file = os.path.join(atlaspath, 'disk_uniform_edges.vtp')
        overlay = readvtp(i_file)

        surface = project_nan_array(
                                source,
                                surface,
                                arraytarget
                                )

        # make a glyph for every point a non "nan" value
        glyphsize = 0.015
        glyph = make_nan_glyph(surface, arraytarget,glyphsize)
        o_file = (outprefix + '_disk_uniform_glyph.png')
        writevtp(glyph, outprefix + '_disk_uniform_glyph.vtp')

        visualise_nan_glyph(surface,
                            glyph,
                            overlay,
                            arraytarget,
                            arraytarget,
                            1,
                            1,
                            colormap,
                            0,
                            o_file,
                            arraytarget,
                            value_range[0],
                            value_range[1])

        o_file = (outprefix + '_disk_uniform.vtp')
        writevtp(surface, o_file)

        # find original2unfold point_map
        writepointmap2csv(point_map,  outprefix + '_intermediate_point_map.csv')
        original2unfold_point_map(point_map)
        writepointmap2csv(point_map,  outprefix + '_point_map.csv', ['original2unfold'])

def run_quantification(atlaspath, unfolddisktarget, datatype, arraytarget,
                       colormaptarget, threshold_value, value_range,
                       paired_unfold_disk='', paired_array='',
                       paired_colormap='', paired_range=''):
    """
    Computes metrics per region (e.g. lesion extent (lge and force)
                                       and average force per segment)

    Arguments:
        atlaspath           =   root path of average mesh and unfold disk
        unfolddisktarget    =   full path to target unfold disk (either lge or carto)
        datatype            =   data type to process: force | lge
        arraytarget         =   array to process on target mesh
        colormaptarget      =   colormap for arrayref visualisation
        threshold_value     =   above this value data will be measured for coverage
                                (to disable use 'NaN').
        value_range         =   pair with min and max values for visulization
        paired_unfold_disk    =   full path to a paired unfold disk (e.g. force).
                                The edges of the current unfold will be ovarlayed
                                on the paired unfold disk (NaN | path_to_mesh).
        paired_array        =   array to process on paired unfold
        paired_colormap     =   colormap for paired_array visualisation
        paired_range         =   pair with min and max values for visulization
                                """

    fileroot = os.path.dirname(unfolddisktarget)
    filenameroot = os.path.splitext(os.path.basename(unfolddisktarget))[0]
    outtarget = os.path.join (fileroot,filenameroot)

    # if a threshold value is provided, any region above that value will be
    # used to quantify lesion extent
    if threshold_value:
        print "Computing coverage of", arraytarget, "above: ", threshold_value
        target = readvtp(unfolddisktarget)

        if datatype == 'lge':
            print "Running histogram normalisation..."
            target = histogram_normalisation(target, arraytarget)
            arraysufix='_norm'
        else:
            print "No histogram normalisation"
            arraysufix=''

        targetnorm_th= add_threshold_array(target, arraytarget + arraysufix, float(threshold_value))
        o_file = os.path.join (fileroot, filenameroot + arraysufix + '.vtp')
        writevtp(target, o_file)

        # compute edges
        thregion = pointthreshold(targetnorm_th, arraytarget + arraysufix + '_th', 1, 1)
        overlay = extractboundaryedge(thregion)

        o_file = outtarget + arraysufix + '_th_edges.png'
        visualise_default_continuous(target,
                                     overlay,
                                     arraytarget + arraysufix,
                                     arraytarget+ arraysufix,
                                     1,
                                     1,
                                     colormaptarget,
                                     0,
                                     o_file, colormaptarget)

        if paired_unfold_disk:
            # get file root
            paired_root = os.path.dirname(paired_unfold_disk)
            paired_filename = os.path.splitext(os.path.basename(paired_unfold_disk))[0]
            outref = os.path.join (paired_root,paired_filename)

            # add current edges onto paired unfold disk
            paired = readvtp(paired_unfold_disk)
            o_file = outref + '_' + arraytarget + '_edges.png'
            print "Saving paired overlay to",o_file

            visualise_default_continuous(paired,
                                         overlay,
                                         paired_array,
                                         paired_array,
                                         1,
                                         1,
                                         paired_colormap,
                                         0,
                                         o_file,
                                         paired_colormap,
                                         paired_range[0],
                                         paired_range[1])

        # compute % coverage per region
        # add array for results.
        target_array = add_cell_array(target, arraytarget + '_per', 0)

        # then populate each region with value
        target_array = area_coverage(target_array, arraytarget + arraysufix + '_th', arraytarget + '_per')
        target_array = circumferential_coverage(target_array, arraytarget + arraysufix + '_th', arraytarget + '_per')
        target_array = longitudinal_coverage(target_array, arraytarget + arraysufix + '_th', arraytarget + '_per')

        o_file = outtarget + arraysufix + '_per.vtp'
        writevtp(target_array, o_file)

        # percentage plus regions averlay
        i_file = os.path.join(atlaspath, 'disk_uniform_edges.vtp')
        edges = readvtp(i_file)
        o_file = outtarget + arraysufix + '_per_0_100.png'
        visualise_default_continuous(target_array,
                                     edges,
                                     arraytarget + '_per',
                                     arraytarget + '_per',
                                     1,
                                     1,
                                     colormaptarget,
                                     0,
                                     o_file,
                                     'per',
                                     0,
                                     100)

    else:
        # compute mean value per region
        print "Averaging", arraytarget, "per region"

        target = readvtp(unfolddisktarget)
        target_array = add_cell_array(target, arraytarget + '_mean', 0.)

        target_array = mean_value_per_region(target_array,
                                           arraytarget,
                                           arraytarget + '_mean')

        o_file = outtarget + '_mean.vtp'
        writevtp(target_array, o_file)

        # percentage plus regions overlay
        i_file = os.path.join(atlaspath, 'disk_uniform_edges.vtp')
        edges = readvtp(i_file)
        o_file = outtarget + '_mean.png'
        visualise_default_continuous(target_array,
                                     edges,
                                     arraytarget + '_mean',
                                     arraytarget + '_mean',
                                     1,
                                     1,
                                     colormaptarget,
                                     0,
                                     o_file,
                                     'mean',
                                     value_range[0],
                                     value_range[1])
