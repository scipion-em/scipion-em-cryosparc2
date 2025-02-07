import argparse
import sys
import json
import os
import time


def cs2Star(args):
    from glob import glob
    import pandas as pd
    import numpy as np
    import logging
    from pyem import metadata
    from pyem import star
    # see https://numpy.org/doc/stable/reference/generated/numpy.load.html
    # for an explanation of MAX_HEADER_SIZE
    MAX_HEADER_SIZE = 50000

    log = logging.getLogger('root')
    hdlr = logging.StreamHandler(sys.stdout)
    log.addHandler(hdlr)
    log.setLevel(logging.getLevelName(logging.INFO))

    if args.swapxy:
        log.warning(
            "Axis swapping is now the default and --swapxy has no effect. "
            "Use --noswapxy if unswapping is needed (unlikely).")

    if args.input[0].endswith(".cs"):
        log.debug("Detected CryoSPARC 2+ .cs file")
        cs = np.load(args.input[0], max_header_size=MAX_HEADER_SIZE)
        if args.first10k:
            cs = cs[:10000]

        if args.movies:
            if not os.path.isdir(args.output):
                log.error("%s is not a directory" % args.output)
                return 1
            log.info("Writing per-movie star files into %s" % args.output)
            trajdir = os.path.dirname(os.path.dirname(args.input[0]))
            if len(args.input) > 1:
                pts = [a for a in args.input[1:] if a.endswith(".cs")]
            else:
                pts = None
            data_general = metadata.cryosparc_2_cs_movie_parameters(cs,
                                                                    passthroughs=pts,
                                                                    trajdir=trajdir,
                                                                    path=args.micrograph_path)
            data_general[star.Relion.MICROGRAPHMETADATA] = data_general[
                star.Relion.MICROGRAPH_NAME].apply(
                lambda x: os.path.join(args.output, os.path.basename(
                    x.rstrip(".mrc")) + ".star"))
            for mic in metadata.cryosparc_2_cs_motion_parameters(cs,
                                                                 data_general,
                                                                 trajdir=trajdir):
                fn = mic[star.Relion.GENERALDATA][
                    star.Relion.MICROGRAPHMETADATA]
                log.debug("Writing %s" % fn)
                star.write_star_tables(fn, mic)
            fields = [star.Relion.VOLTAGE, star.Relion.CS, star.Relion.AC,
                      star.Relion.MICROGRAPHORIGINALPIXELSIZE,
                      star.Relion.MICROGRAPHPIXELSIZE,
                      star.Relion.MICROGRAPH_NAME,
                      star.Relion.MICROGRAPHMETADATA,
                      star.Relion.MICROGRAPHBINNING, star.Relion.OPTICSGROUP]
            if len(args.input) > 1 and args.input[-1].endswith(".star"):
                mic_star = args.input[-1]
                star.write_star(mic_star, data_general[
                    [f for f in fields if f in data_general]])
            return 0

        try:
            log.info("Parsing .cs file...")
            df = metadata.parse_cryosparc_2_cs(cs, passthroughs=args.input[1:],
                                               minphic=args.minphic,
                                               boxsize=args.boxsize,
                                               swapxy=args.noswapxy,
                                               invertx=args.invertx,
                                               inverty=args.inverty)
        except Exception as e:
            log.error(e, exc_info=True)
            log.error("Required fields could not be mapped. Are you using the "
                      "right input file(s)?")
            return 1
    else:
        log.debug("Detected CryoSPARC 0.6.5 .csv file")
        if len(args.input) > 1:
            log.error("Only one file at a time supported for "
                      "CryoSPARC 0.6.5 .csv format")
            return 1
        meta = metadata.parse_cryosparc_065_csv(
            args.input[0])  # Read cryosparc metadata file.
        df = metadata.cryosparc_065_csv2star(meta, args.minphic)

    if args.cls is not None:
        df = star.select_classes(df, args.cls)

    if args.flipy:
        log.info("Flipping refined shifts in Y")
        df[star.Relion.ORIGINY] = -df[star.Relion.ORIGINY]
        log.info("Flipping particle orientation through XZ plane")
        df = star.transform_star(df,
                                 np.array([[1, 0, 0], [0, -1, 0], [0, 0, -1]]),
                                 leftmult=True)

    if args.strip_uid is not None:
        df = star.strip_path_uids(df, inplace=True, count=args.strip_uid)

    if args.copy_micrograph_coordinates is not None:
        df = star.augment_star_ucsf(df, inplace=True)
        coord_star = pd.concat(
            (star.parse_star(inp, keep_index=False, augment=True) for inp in
             glob(args.copy_micrograph_coordinates)), join="inner")
        key = star.merge_key(df, coord_star, threshold=0)
        if key is None and not args.strip_uid:
            log.debug("Merge key not found, removing leading UIDs")
            df = star.strip_path_uids(df, inplace=True)
            key = star.merge_key(df, coord_star)
        log.debug("Coordinates merge key: %s" % key)
        if args.cached or key == star.Relion.IMAGE_NAME:
            fields = star.Relion.MICROGRAPH_COORDS
        else:
            fields = star.Relion.MICROGRAPH_COORDS + [star.UCSF.IMAGE_INDEX,
                                                      star.UCSF.IMAGE_PATH]
        n = df.shape[0]
        df = star.smart_merge(df, coord_star, fields=fields, key=key)
        star.simplify_star_ucsf(df)
        if df.shape[0] != n:
            log.warning("%d / %d particles remain after coordinate merge" % (
            df.shape[0], n))

    if args.micrograph_path is not None:
        df = star.replace_micrograph_path(df, args.micrograph_path,
                                          inplace=True)

    if args.transform is not None:
        r = np.array(json.loads(args.transform))
        df = star.transform_star(df, r, inplace=True)

    try:
        df = star.check_defaults(df, inplace=True)
    except Exception as e:
        print(e)

    if args.relion2:
        df = star.remove_new_relion31(df, inplace=True)
        star.write_star(args.output, df, resort_records=True, optics=False)
    else:
        df = star.remove_deprecated_relion2(df, inplace=True)
        # Changing NaN values. These values denote erroneous coordinates
        if hasattr(df, 'rlnAnglePsi'):
            nanValues = len(df.rlnAnglePsi.values[np.isnan(df.rlnAnglePsi.values)])
            if nanValues:
                df.rlnAnglePsi.values[np.isnan(df.rlnAnglePsi.values)] = 0
                log.warning("WARNING: %d dataframes contains erroneous "
                            "coordinates. These coordinates are removed" % nanValues)
        star.write_star(args.output, df, resort_records=True, optics=True)

    log.info("Output fields: %s" % ", ".join(df.columns))
    return 0


def defineArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("input",
                        help="Cryosparc metadata .csv (v0.6.5) or .cs (v2+) files",
                        nargs="*")
    parser.add_argument("output", help="Output .star file")
    parser.add_argument("--movies",
                        help="Write per-movie star files into output directory",
                        action="store_true")
    parser.add_argument("--boxsize",
                        help="Cryosparc refinement box size (if different from particles)",
                        type=float)
    # parser.add_argument("--passthrough", "-p", help="List file required for some Cryosparc 2+ job types")
    parser.add_argument("--class",
                        help="Keep this class in output, may be passed multiple times",
                        action="append", type=int, dest="cls")
    parser.add_argument("--minphic",
                        help="Minimum posterior probability for class assignment",
                        type=float, default=0)
    parser.add_argument("--stack-path", help="Path to single particle stack",
                        type=str)
    parser.add_argument("--micrograph-path",
                        help="Replacement path for micrographs or movies")
    parser.add_argument("--copy-micrograph-coordinates",
                        help="Source for micrograph paths and particle coordinates (file or quoted glob)",
                        type=str)
    parser.add_argument("--swapxy",
                        help="Swap X and Y axes when converting particle coordinates from normalized to absolute",
                        action="store_true")
    parser.add_argument("--noswapxy",
                        help="Do not swap X and Y axes when converting particle coordinates",
                        action="store_false")
    parser.add_argument("--invertx", help="Invert particle coordinate X axis",
                        action="store_true")
    parser.add_argument("--inverty", help="Invert particle coordinate Y axis",
                        action="store_false")
    parser.add_argument("--flipy", help="Invert refined particle Y shifts",
                        action="store_true")
    parser.add_argument("--cached",
                        help="Keep paths from the Cryosparc 2+ cache when merging coordinates",
                        action="store_true")
    parser.add_argument("--transform",
                        help="Apply rotation matrix or 3x4 rotation plus translation matrix to particles (Numpy format)",
                        type=str)
    parser.add_argument("--relion2", "-r2", help="Relion 2 compatible outputs",
                        action="store_true")
    parser.add_argument("--strip-uid",
                        help="Strip all leading UIDs from file names",
                        nargs="?", default=None, const=-1,
                        type=int)
    parser.add_argument("--10k",
                        help="Only read first 10,000 particles for rapid testing.",
                        action="store_true",
                        dest="first10k")
    parser.add_argument("--loglevel", "-l", type=str, default="WARNING",
                        help="Logging level and debug output")
    return parser


if __name__ == "__main__":
    parser = defineArgs()
    argsList = [sys.argv[1]]
    if len(sys.argv) > 3:
        argsList.append(sys.argv[3])
    argsList.append(sys.argv[2])
    args = parser.parse_args(argsList)
    sys.exit(cs2Star(args))
