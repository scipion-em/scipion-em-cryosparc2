import argparse
import sys
import re
import json

def cs2Star(args):
    from glob import glob
    import pandas as pd
    import numpy as np
    import logging
    from pyem import metadata
    from pyem import star

    log = logging.getLogger('root')
    hdlr = logging.StreamHandler(sys.stdout)
    log.addHandler(hdlr)
    log.setLevel(logging.getLevelName(args.loglevel.upper()))
    if args.input[0].endswith(".cs"):
        log.debug("Detected CryoSPARC 2+ .cs file")
        cs = np.load(args.input[0])
        try:
            df = metadata.parse_cryosparc_2_cs(cs, passthroughs=args.input[1:],
                                               minphic=args.minphic,
                                               boxsize=args.boxsize,
                                               swapxy=args.swapxy,
                                               invertx=args.invertx,
                                               inverty=args.inverty)
        except (KeyError, ValueError) as e:
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

    def strip_path_uids1(df, inplace=False, count=-1):
        df = df if inplace else df.copy()
        pat = re.compile("[0-9]{21}_")
        if star.UCSF.IMAGE_PATH in df:
            df[star.UCSF.IMAGE_PATH] = df[star.UCSF.IMAGE_PATH].str.replace(pat, "",
                                                                  regex=True,
                                                                  n=count)
        return df

    if args.strip_uid is not None:
        df = star.strip_path_uids(df, inplace=True, count=args.strip_uid)
        df = strip_path_uids1(df, inplace=True, count=args.strip_uid)

    if args.copy_micrograph_coordinates is not None:
        df = star.augment_star_ucsf(df, inplace=True)
        coord_star = pd.concat(
            (star.parse_star(inp, keep_index=False, augment=True) for inp in
             glob(args.copy_micrograph_coordinates)), join="inner")
        key = star.merge_key(df, coord_star)
        if key is None:
            log.debug("Merge key not found, removing leading UIDs")
            df = star.strip_path_uids(df, inplace=True)
            key = star.merge_key(df, coord_star)
        log.debug("Coordinates merge key: %s" % key)
        if args.cached or key == star.Relion.IMAGE_NAME:
            fields = star.Relion.MICROGRAPH_COORDS
        else:
            fields = star.Relion.MICROGRAPH_COORDS + [star.UCSF.IMAGE_INDEX,
                                                      star.UCSF.IMAGE_PATH]
        df = star.smart_merge(df, coord_star, fields=fields, key=key)
        star.simplify_star_ucsf(df)

    if args.micrograph_path is not None:
        df = star.replace_micrograph_path(df, args.micrograph_path,
                                          inplace=True)

    if args.transform is not None:
        r = np.array(json.loads(args.transform))
        df = star.transform_star(df, r, inplace=True)

    if args.relion2:
        df = star.remove_new_relion31(df, inplace=True)
        star.write_star(args.output, df, resort_records=True, optics=False)
    else:
        # df = star.remove_deprecated_relion2(df, inplace=True)
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
                        help="Replacement path for micrographs")
    parser.add_argument("--copy-micrograph-coordinates",
                        help="Source for micrograph paths and particle coordinates (file or quoted glob)",
                        type=str)
    parser.add_argument("--swapxy",
                        help="Swap X and Y axes when converting particle coordinates from normalized to absolute",
                        action="store_true")
    parser.add_argument("--invertx", help="Invert particle coordinate X axis",
                        action="store_true")
    parser.add_argument("--inverty", help="Invert particle coordinate Y axis",
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
                        nargs="?", default=0, type=int)
    parser.add_argument("--loglevel", "-l", type=str, default="WARNING",
                        help="Logging level and debug output")
    return parser

if __name__ == "__main__":
    parser = defineArgs()
    argsList = [sys.argv[1], sys.argv[2]]
    args = parser.parse_args(argsList)
    sys.exit(cs2Star(args))