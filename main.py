import os
import sys
import shutil
import argparse
from processing import *

def get_parameters(datatype):
    """"Parameters specific to each tested data type:
    * contact-catheter force (force)
    * late gadolinium enhancement (lge)
    * local activation time (lat)

    These parameters can be modified to fit new data."""

    if datatype == 'force':
        return  dict(restorenans = False,
                     sourcearrays = ['FTI'],
                     targetarrays = ['FTI'],
                     colormaps = ['erdc_rainbow_grey'],
                     value_range = [[0,800]], # for colormap

                     # parameters used if a paired_unfold_disk is provided
                     paired_arrays = [''],
                     paired_colormaps = [''],
                     paired_range = [['','']],

                     pvends = 1,
                     skippointsfactor = 0.05, # percentage of points to ignore at beginning of centerline
                     clspacing = 0.4, # resample the centerline with this spacing
                     maxslope = 5, # anything above this is ostium
                     highslope = 1.2, # above this slope we start counting
                     bumpcriterion = 0.05, # ostium if slope higher than highslope and above bump criterion
                     threshold_value = 100, # threshold value for FTI average
                     quantification_range=[0,800])

    elif datatype == 'lge':
        return  dict(restorenans = False,
                     sourcearrays = ['scalars'],
                     targetarrays = ['lge'],
                     colormaps = ['YlOrRd'],
                     value_range = [['','']], # with empty values, the min and max will be computed from array

                     # parameters used if a paired_unfold_disk is provided
                     paired_arrays = ['FTI'],
                     paired_colormaps = ['erdc_rainbow_grey'],
                     paired_range = [[0,800]],

                     pvends = 1,
                     skippointsfactor = 0.1, # percentage of points to ignore at beginning of centerline
                     clspacing = 0.4, # resample the centerline with this spacing
                     maxslope = 5, # anything above this is ostium
                     highslope = 1.2, # above this slope we start counting
                     bumpcriterion = 0.05, # ostium if slope higher than highslope and above bump criterion
                     threshold_value = 0.5, # threshold value for lge segmentation
                     quantification_range=['',''])

    elif datatype == 'lat':
            return  dict(restorenans = True,
                        sourcearrays = ['CS_LAT','HRA_LAT'],
                        targetarrays = ['CS_LAT','HRA_LAT'],
                        colormaps = ['erdc_rainbow_inv',
                                  'erdc_rainbow_inv'],
                        value_range = [ [0,160],[0,160]],

                        # parameters used if a paired_unfold_disk is provided
                        paired_arrays = ['',''],
                        paired_colormaps = ['',''],
                        paired_range = [['',''],['','']],

                        pvends = 0, # Enforce the centerline to reach the end boundary of the surface.
                                    # Appends a segment at the end of the centerline to reach the seed point.
                                    # Useful for a blunt vein. Turn on or off (1|0)
                        skippointsfactor = 0.1, # percentage of points to ignore at beginning of centerline
                        clspacing = 0.4, # resample the centerline with this spacing
                        maxslope = 5, # anything above this is ostium
                        highslope = 1.2, # above this slope we start counting
                        bumpcriterion = 0.05, # ostium if slope higher than highslope and above bump criterion
                        threshold_value = '', # disables thresholding
                        quantification_range=[0,160])
    else:
        print "Unrecognised data type"


def extend_arguments_dictionary(args_dict):

    """Default variables """
    FULL_PATH = os.getcwd()
    args_dict['atlaspath']= os.path.join(FULL_PATH,'./atlas')
    args_dict['fileroot']= os.path.dirname(args_dict['meshfile'])
    args_dict['filenameroot'] = os.path.splitext(
                                    os.path.basename(args_dict['meshfile']))[0]

    if args_dict['use_laa_seed'] and args_dict['pvcliptype'] == 'long':
        print "WARNING: LAA seed only available with short PVs. Switching configuration."
        args_dict['pvcliptype'] = 'short'

    if args_dict['pvcliptype'] == 'short':
        args_dict['pvdist'] = 2 # length of PVs to keep
    elif args_dict['pvcliptype'] == 'long':
        args_dict['pvdist'] = 10 # length of PVs to keep
    else:
        print "Unrecognised pvcliptype"

    args_dict.update( get_parameters(args_dict['datatype']) )
    return args_dict

def main(args):

    # Appends other variables to arguments
    args_dict = extend_arguments_dictionary(vars(args))

    # Run the pipeline
    # Steps can be skipped if input data already exists
    if not args_dict['skip_standardization']:
        run_standardization(args_dict['meshfile'],
                args_dict['atlaspath'],
                args_dict['pvdist'],
                args_dict['maxslope'],
                args_dict['clspacing'],
                args_dict['skippointsfactor'],
                args_dict['highslope'],
                args_dict['bumpcriterion'],
                args_dict['pvends'],
                args_dict['use_seed_selector'],
                args_dict['use_laa_seed'],
                args_dict['visualize'])

    if not args_dict['skip_currents']:
        run_currents(args_dict['meshfile'],
                     args_dict['atlaspath'],
                     args_dict['mitral_clip_type'],
                     args_dict['pvcliptype'],
                     args_dict['use_similarity'],
                     args_dict['use_laa_seed'],
                     args_dict['visualize'])

    if not args_dict['skip_sum']:
        for arrayind in range(len(args_dict['targetarrays'])):
            run_sum(args_dict['meshfile'],
                    args_dict['atlaspath'],
                    args_dict['pvcliptype'],
                    args_dict['use_glyphs'],
                    args_dict['restorenans'],
                    args_dict['sourcearrays'][arrayind],
                    args_dict['targetarrays'][arrayind],
                    args_dict['colormaps'][arrayind],
                    args_dict['value_range'][arrayind],
                    args_dict['visualize'])

    if not args_dict['skip_quantification']:
        for arrayind in range(len(args_dict['targetarrays'])):
            unfoldiskfile = os.path.join(
                                args_dict['fileroot'],
                                args_dict['filenameroot'] +
                                '_' +
                                args_dict['targetarrays'][arrayind] +
                                '_disk_uniform.vtp')
            # quantifies extent
            run_quantification(args_dict['atlaspath'],
                               unfoldiskfile,
                               args_dict['datatype'],
                               args_dict['targetarrays'][arrayind],
                               args_dict['colormaps'][arrayind],
                               args_dict['threshold_value'],
                               args_dict['quantification_range'],
                               args_dict['paired_unfold_disk'],
                               args_dict['paired_arrays'][arrayind],
                               args_dict['paired_colormaps'][arrayind],
                               args_dict['paired_range'][arrayind])

            #quantifies average values
            run_quantification(args_dict['atlaspath'],
                               unfoldiskfile,
                               args_dict['datatype'],
                               args_dict['targetarrays'][arrayind],
                               args_dict['colormaps'][arrayind],
                               '', # disables thresholding
                               args_dict['quantification_range'])

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = 'Standardized unfold map')

    parser.add_argument("--meshfile", type=str,
                    help="Full path to mesh file in VTK format", required=True)

    parser.add_argument("--datatype", type=str,
        help="Data type to process: force | lge | lat", required=True)

    # optional inputs
    parser.add_argument("--pvcliptype", type=str, default='short',
        help="How to clip the pulmonary veins: short | long" )

    parser.add_argument("--paired_unfold_disk", type=str, default='',
        help="Full path to a paired unfold disk. The edges of the current unfold " +
        "will be overlaid on the paired unfold disk" )

    parser.add_argument("--mitral_clip_type", type=str, default='auto',
        help="The algorithm will compute a mitral plane based on the body and PVs centroids. " +
        "If the meshfile already contains a mitral plane clip, use the 'manual' option to preserve it" )

    parser.add_argument("--use_seed_selector", action='store_true', default=False,
                        help="Seed selection will always run for a new case. " +
                        "To select new seeds, activate this flag")

    parser.add_argument("--use_laa_seed", action='store_true', default=False,
                        help="Activate seed for appendage")

    parser.add_argument("--use_similarity", action='store_true', default=False,
                        help="The method uses an affine transform " +
                        "to initialize mesh registration. " +
                        "Activate this flag to use a similarity transform")

    parser.add_argument("--use_glyphs", action='store_true', default=False,
                        help="Activate glyph on visualization. " +
                        "Recommended for spare measurements (e.g. some lat meshes)")

    parser.add_argument('--visualize', action='store_true', default=False,
                        help="Activate visualization")

    parser.add_argument('--skip_standardization', action='store_true', default=False,
                        help="Skip run_standardization step")

    parser.add_argument('--skip_currents', action='store_true', default=False,
                        help="Skip run_currents step")

    parser.add_argument('--skip_sum',action='store_true', default=False,
                        help="Skip run_sum step")

    parser.add_argument('--skip_quantification', action='store_true', default=False,
                        help="Skip run_quantification step")

    args = parser.parse_args()

    main(args)


