# **************************************************************************
# *
# * Authors: Yunior C. Fonseca Reyna    (cfonseca@cnb.csic.es)
# *
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
import  subprocess

import emtable
import numpy as np
import os
import argparse
import sys
import logging
logger = logging.getLogger(__name__)

from emtable.metadata import _guessType
from pwem.emlib.image import ImageHandler
import pwem.emlib.metadata as md
from pwem.objects import (String, Integer, Transform, Particle,
                          Coordinate, Acquisition, CTFModel)
from pyworkflow.object import ObjectWrap, Float
import pyworkflow.utils as pwutils
from pwem.constants import *

from ..constants import *
from .. import Plugin


def convertCs2Star(argsList):
    input = os.path.abspath(argsList[0])
    output = os.path.abspath(argsList[1])
    args = ' %s %s' % (input, output)
    if len(argsList) > 2:
        passthrough = os.path.abspath(argsList[2])
        args += ' %s' % passthrough
    cryosparcScriptPath = os.path.join(os.path.dirname(__file__),
                                       CRYOSPARC_CS2STAR_SCRIPT)
    cmd = (Plugin.getCondaActivationCmd() + Plugin.getPyemEnvActivation() + ' && python3 ' +
           cryosparcScriptPath + args)

    logger.info("convertCs2Star: %s" % cmd)

    process = subprocess.Popen(cmd, shell=True, cwd=os.getcwd(), stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    out, error = process.communicate()
    logger.info(out.decode())
    logger.error(error.decode())


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
    # parser.add_argument("--passthrough", "-p",
    #                     help="List file required for some Cryosparc 2+ job types")
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


def addRandomSubset(img, imgRow):
    halve = 1 + (img.getObjId() % 2)
    imgRow.setValue(md.RLN_PARTICLE_RANDOM_SUBSET, int(halve))


def cryosparcToLocation(filename):
    """ Return a location (index, filename) given
    a cryoSPARC filename with the index@filename structure. """
    if '@' in filename:
        indexStr, fn = filename.split('@')
        return int(indexStr), str(fn)
    else:
        return NO_INDEX, str(filename)


def setOfImagesToMd(imgSet, imgMd, imgToFunc, **kwargs):
    """ This function will fill Relion metadata from a SetOfMicrographs
    Params:
        imgSet: the set of images to be converted to metadata
        md: metadata to be filled
        rowFunc: this function can be used to setup the row before
            adding to meta
    """

    if 'alignType' not in kwargs:
        kwargs['alignType'] = imgSet.getAlignment()

    for img in imgSet:
        objId = imgMd.addObject()
        imgRow = md.Row()
        imgToFunc(img, imgRow, **kwargs)
        imgRow.writeToMd(imgMd, objId)


def particleToRow(part, partRow, **kwargs):
    """ Set labels values from Particle to md row. """
    coord = part.getCoordinate()
    if coord is not None:
        coordinateToRow(coord, partRow, copyId=False)
    if part.hasMicId():
        partRow.setValue(md.RLN_MICROGRAPH_ID, int(part.getMicId()))
        # If the row does not contain the micrograph name
        # use a fake micrograph name using id to relion
        # could at least group for CTF using that
        if not partRow.hasLabel(md.RLN_MICROGRAPH_NAME):
            partRow.setValue(md.RLN_MICROGRAPH_NAME,
                             'fake_micrograph_%06d.mrc' % part.getMicId())
    if part.hasAttribute('_rlnParticleId'):
        partRow.setValue(md.RLN_PARTICLE_ID, int(part._rlnParticleId.get()))

    if kwargs.get('fillRandomSubset') and part.hasAttribute('_rlnRandomSubset'):
        partRow.setValue(md.RLN_PARTICLE_RANDOM_SUBSET,
                         int(part._rlnRandomSubset.get()))
        if part.hasAttribute('_rlnBeamTiltX'):
            partRow.setValue('rlnBeamTiltX',
                             float(part._rlnBeamTiltX.get()))
            partRow.setValue('rlnBeamTiltY',
                             float(part._rlnBeamTiltY.get()))

    imageToRow(part, partRow, md.RLN_IMAGE_NAME, **kwargs)


def imageToRow(img, imgRow, imgLabel=RELIONCOLUMNS.rlnImageName.value,
               **kwargs):
    # Provide a hook to be used if something is needed to be
    # done for special cases before converting image to row
    preprocessImageRow = kwargs.get('preprocessImageRow', None)
    if preprocessImageRow:
        preprocessImageRow(img, imgRow)

    setRowId(imgRow, img)  # Set the id in the metadata as MDL_ITEM_ID
    index, fn = img.getLocation()
    # check if the is a file mapping
    filesDict = kwargs.get('filesDict', {})
    filename = filesDict.get(fn, fn)

    imgRow.set(imgLabel, locationToCryosparc(index, filename))

    if kwargs.get('writeCtf', True) and img.hasCTF():
        ctfModelToRow(img.getCTF(), imgRow)

    # alignment is mandatory at this point, it should be check
    # and detected defaults if not passed at readSetOf.. level
    alignType = kwargs.get('alignType')

    if alignType != ALIGN_NONE and img.hasTransform():
        alignmentToRow(img.getTransform(), imgRow, alignType)

    if kwargs.get('writeAcquisition', True) and img.hasAcquisition():
        acquisitionToRow(img.getAcquisition(), imgRow)

    # Write all extra labels to the row
    dictLabel = {}
    objectToRow(img, imgRow, dictLabel,
                extraLabels=IMAGE_EXTRA_LABELS + kwargs.get('extraLabels', []))

    # Provide a hook to be used if something is needed to be
    # done for special cases before converting image to row
    postprocessImageRow = kwargs.get('postprocessImageRow', None)
    if postprocessImageRow:
        postprocessImageRow(img, imgRow)


def acquisitionToRow(acquisition, ctfRow):
    """ Set labels values from acquisition to md row. """
    objectToRow(acquisition, ctfRow, ACQUISITION_DICT)


def ctfModelToRow(ctfModel, ctfRow):
    """ Set labels values from ctfModel to md row. """
    # Refresh phase shift!
    phaseShift = ctfModel.getPhaseShift()

    if phaseShift is not None:
        ctfRow.set(RELIONCOLUMNS.rlnPhaseShift.value, phaseShift)

    objectToRow(ctfModel, ctfRow, CTF_DICT, extraLabels=CTF_EXTRA_LABELS)


def locationToCryosparc(index, filename):
    """ Convert an index and filename location
    to a string with @ as expected in cryoSPARC.
    """
    if index != NO_INDEX:
        return "%06d@%s" % (index, filename)

    return filename


def alignmentToRow(alignment, alignmentRow, alignType):
    """
    is2D == True-> matrix is 2D (2D images alignment)
            otherwise matrix is 3D (3D volume alignment or projection)
    invTransform == True  -> for xmipp implies projection
                          -> for xmipp implies alignment
    """
    is2D = alignType == ALIGN_2D
    is3D = alignType == ALIGN_3D
    inverseTransform = alignType == ALIGN_PROJ
    matrix = alignment.getMatrix()
    shifts, angles = geometryFromMatrix(matrix, inverseTransform)

    alignmentRow.set(RELIONCOLUMNS.rlnOriginX.value, shifts[0])
    alignmentRow.set(RELIONCOLUMNS.rlnOriginY.value, shifts[1])

    if is2D:
        angle = angles[0] + angles[2]
        alignmentRow.set(RELIONCOLUMNS.rlnAnglePsi.value, -angle)

        flip = bool(np.linalg.det(matrix[0:2, 0:2]) < 0)
        if flip:
            logger.debug("FLIP in 2D not implemented")
    elif is3D:
        raise Exception("3D alignment conversion for Relion not implemented. "
                        "It seems the particles were generated with an "
                        "incorrect alignment type. You may either re-launch "
                        "the protocol that generates the particles "
                        "with angles or set 'Consider previous alignment?' "
                        "to No")
    else:
        alignmentRow.set(RELIONCOLUMNS.rlnOriginZ.value, shifts[2])
        alignmentRow.set(RELIONCOLUMNS.rlnAngleRot.value, angles[0])
        alignmentRow.set(RELIONCOLUMNS.rlnAngleTilt.value, angles[1])
        alignmentRow.set(RELIONCOLUMNS.rlnAnglePsi.value, angles[2])


def geometryFromMatrix(matrix, inverseTransform):
    from pwem.convert.transformations import (translation_from_matrix,
                                              euler_from_matrix)

    if inverseTransform:
        from numpy.linalg import inv
        matrix = inv(matrix)
        shifts = -translation_from_matrix(matrix)
    else:
        shifts = translation_from_matrix(matrix)
    angles = -np.rad2deg(euler_from_matrix(matrix, axes='szyz'))
    return shifts, angles


def coordinateToRow(coord, coordRow, copyId=True):
    """ Set labels values from Coordinate coord to md row. """
    if copyId:
        setRowId(coordRow, coord)
    objectToRow(coord, coordRow, COOR_DICT, extraLabels=COOR_EXTRA_LABELS)
    if coord.getMicName():
        micName = coord.getMicName()
        coordRow.setValue(md.RLN_MICROGRAPH_NAME, str(micName.replace(" ", "")))
    else:
        if coord.getMicId():
            coordRow.setValue(md.RLN_MICROGRAPH_NAME, str(coord.getMicId()))


def objectToRow(obj, row, attrDict, extraLabels=[]):
    """ This function will convert an EMObject into a XmippMdRow.
    Params:
        obj: the EMObject instance (input)
        row: the XmippMdRow instance (output)
        attrDict: dictionary with the map between obj attributes(keys) and
            row MDLabels in Xmipp (values).
        extraLabels: a list with extra labels that could be included
            as _xmipp_labelName
    """
    row.set(RELIONCOLUMNS.rlnEnabled.value, obj.isEnabled())

    for attr, label in attrDict.items():
        if hasattr(obj, attr):
            valueType = md.label2Python(label)
            row.set(label, valueType(getattr(obj, attr).get()))

    attrLabels = attrDict.values()

    for label in extraLabels:
        attrName = '_' + label
        if label not in attrLabels and hasattr(obj, attrName):
            value = obj.getAttributeValue(attrName)
            row.setValue(label, value)


def setRowId(mdRow, obj, label=RELIONCOLUMNS.rlnImageId.value):
    mdRow.set(label, int(obj.getObjId()))


def convertBinaryVol(vol, outputDir):
    """ Convert binary volume to a format read by Cryosparc.
    Params:
        vol: input volume object to be converted.
        outputDir: where to put the converted file(s)
    Return:
        new file name of the volume (converted or not).
    """

    ih = ImageHandler()

    # This approach can be extended when
    # converting from a binary file format that
    # is not read from Relion
    def convertToMrc(fn):
        """ Convert from a format that is not read by Relion
        to mrc format.
        """
        from os.path import join
        newFn = join(outputDir, pwutils.replaceBaseExt(fn, 'mrc'))
        ih.convert(fn, newFn)
        return newFn

    ext = vol.getFileName()

    if not ext.endswith('.mrc'):
        fn = convertToMrc(vol.getFileName())
    else:
        fn = vol.getFileName()
    return fn


def createItemMatrix(item, row, align):
    item.setCTF(rowToCtfModel(row))
    pixelSize = item.getSamplingRate()
    item.setTransform(rowToAlignment(row, align, pixelSize))


def rowToAlignment(alignmentRow, alignType, pixelSize=1.0):
    """
    is2D == True-> matrix is 2D (2D images alignment)
            otherwise matrix is 3D (3D volume alignment or projection)
    invTransform == True  -> for xmipp implies projection
    """
    if alignType == ALIGN_3D:
        raise Exception("3D alignment conversion for Relion not implemented.")

    is2D = alignType == ALIGN_2D
    inverseTransform = alignType == ALIGN_PROJ
    if alignmentRow.hasAnyColumn(ALIGNMENT_DICT.values()):
        alignment = Transform()
        angles = np.zeros(3)
        shifts = np.zeros(3)
        shifts[0] = alignmentRow.get(RELIONCOLUMNS.rlnOriginXAngst.value, default=0.)/pixelSize
        shifts[1] = alignmentRow.get(RELIONCOLUMNS.rlnOriginYAngst.value, default=0.)/pixelSize
        if not is2D:
            angles[0] = alignmentRow.get(RELIONCOLUMNS.rlnAngleRot.value, default=0.)
            angles[1] = alignmentRow.get(RELIONCOLUMNS.rlnAngleTilt.value, default=0.)
            angles[2] = alignmentRow.get(RELIONCOLUMNS.rlnAnglePsi.value, default=0.)
            shifts[2] = alignmentRow.get(RELIONCOLUMNS.rlnOriginZAngst.value, default=0.)/pixelSize
        else:
            angles[2] = - alignmentRow.get(RELIONCOLUMNS.rlnAnglePsi.value, default=0.)
        M = matrixFromGeometry(shifts, angles, inverseTransform)
        alignment.setMatrix(M)
    else:
        alignment = None

    return alignment


def setCryosparcAttributes(obj, objRow, *labels):
    """ Set an attribute to obj from a label that is not
    basic ones. The new attribute will be named _rlnLabelName
    and the datatype will be set correctly.
    """
    for label in labels:
        value = objRow.get(label)
        valueType = _guessType(value)
        if valueType is int:
            value = Integer(value)
        elif valueType is float:
            value = Float(value)
        else:
            value = String(value)
        setattr(obj, '_%s' % label, value)


def matrixFromGeometry(shifts, angles, inverseTransform):
    """ Create the transformation matrix from a given
    2D shifts in X and Y...and the 3 euler angles.
    """
    from pwem.convert.transformations import euler_matrix
    radAngles = -np.deg2rad(angles)
    M = euler_matrix(radAngles[0], radAngles[1], radAngles[2], 'szyz')
    if inverseTransform:
        from numpy.linalg import inv
        M[:3, 3] = -shifts[:3]
        M = inv(M)
    else:
        M[:3, 3] = shifts[:3]

    return M


def convertBinaryFiles(imgSet, outputDir, extension='mrcs', **kwargs):
    """ Convert binary images files to a format read by Cryosparc.
    Params:
        imgSet: input image set to be converted.
        outputDir: where to put the converted file(s)
    Return:
        A dictionary with old-file as key and new-file as value
        If empty, not conversion was done.
    """
    filesDict = {}
    ih = ImageHandler()
    outputRoot = os.path.join(outputDir, 'input')
    # Get the extension without the dot
    stackFiles = imgSet.getFiles()
    ext = pwutils.getExt(next(iter(stackFiles)))[1:]
    rootDir = pwutils.commonPath(list(stackFiles))

    def getUniqueFileName(fn, extension):
        """ Get an unique file for either link or convert files.
        It is possible that the base name overlap if they come
        from different runs. (like particles.mrcs after relion preprocess)
        """
        from os.path import join
        newFn = join(outputRoot, pwutils.replaceBaseExt(fn, extension))
        newRoot = pwutils.removeExt(newFn)

        values = filesDict.values()
        counter = 1

        while newFn in values:
            counter += 1
            newFn = '%s_%05d.%s' % (newRoot, counter, extension)

        return newFn

    def createBinaryLink(fn):
        """ Just create a link named .mrcs to cryoSPARC understand
        that it is a binary stack file and not a volume.
        """
        newFn = getUniqueFileName(fn, extension)
        if not os.path.exists(newFn):
            pwutils.createAbsLink(os.path.abspath(fn), newFn)
            logger.debug("   %s -> %s" % (newFn, fn))
        return newFn

    def convertStack(fn):
        """ Convert from a format that is not read by Cryosparc
        to an spider stack.
        """
        newFn = getUniqueFileName(fn, 'mrc')
        ih.convertStack(fn, newFn)
        logger.debug("   %s -> %s" % (fn, newFn))
        return newFn

    def replaceRoot(fn):
        """ Link create to the root folder, so just replace that
        in the name, no need to do anything else.
        """
        return fn.replace(rootDir, outputRoot)

    if ext == extension:
        logger.debug("convertBinaryFiles: creating soft links.")
        logger.debug("   Root: %s -> %s" % (outputRoot, rootDir))
        mapFunc = replaceRoot
        pwutils.createAbsLink(os.path.abspath(rootDir), outputRoot)
    elif ext == 'mrc' and extension == 'mrcs':
        logger.debug("convertBinaryFiles: creating soft links (mrcs -> mrc).")
        mapFunc = createBinaryLink
    elif ext.endswith('hdf') or ext.endswith(
            'stk'):  # assume eman .hdf format or .stk format
        logger.debug("convertBinaryFiles: converting stacks. (%s -> %s)"
              % (ext, extension))
        mapFunc = convertStack
    else:
        mapFunc = None

    if mapFunc is not None:
        pwutils.makePath(outputRoot)
        for fn in stackFiles:
            newFn = mapFunc(fn)  # convert or link
            filesDict[fn] = newFn  # map new filename

    return filesDict


def writeSetOfParticles(imgSet, fileName, extraPath):
    args = {'outputDir': extraPath,
            'fillMagnification': True,
            'fillRandomSubset': True}
    # try:
    #     logger.info('Trying to generate the star file with Relion convert...')
    #     from relion import convert
    #     alignType = ALIGN_PROJ if imgSet.hasAlignmentProj() else ALIGN_NONE
    #     args['alignType'] = alignType
    #     args['incompatibleExtensions'] = ['hdf', 'stk']
    #     convert.writeSetOfParticles(imgSet, fileName, **args)
    #     logger.info('The star file was generate successfully ...')
    # except Exception:
    if imgSet.hasAlignmentProj() and imgSet.getAttributeValue("_rlnRandomSubset") is None:
        args['postprocessImageRow'] = addRandomSubset
    logger.info('The star file generation with Relion convert failed ...')
    logger.info('Trying to generate the star file with cryoSPARC convert ...')
    cryosPARCwriteSetOfParticles(imgSet, fileName, **args)
    logger.info('The star file was generate successfully ...')


def cryosPARCwriteSetOfParticles(imgSet, starFile, outputDir, **kwargs):
    if outputDir is not None:
        filesDict = convertBinaryFiles(imgSet, outputDir)
        kwargs['filesDict'] = filesDict
    partMd = md.MetaData()
    setOfImagesToMd(imgSet, partMd, particleToRow, **kwargs)

    if kwargs.get('fillMagnification', False):
        pixelSize = imgSet.getSamplingRate()
        mag = imgSet.getAcquisition().getMagnification()
        detectorPxSize = mag * pixelSize / 10000

        partMd.fillConstant(md.RLN_CTF_MAGNIFICATION, mag)
        partMd.fillConstant(md.RLN_CTF_DETECTOR_PIXEL_SIZE, detectorPxSize)
    else:
        # Remove Magnification from metadata to avoid wrong values of pixel size.
        # In Relion if Magnification and DetectorPixelSize are in metadata,
        # pixel size is ignored in the command line.
        partMd.removeLabel(md.RLN_CTF_MAGNIFICATION)

    blockName = kwargs.get('blockName', 'particles')
    partMd.write('%s@%s' % (blockName, starFile))


def rowToCtfModel(ctfRow):
    """ Create a CTFModel from a row of a meta """
    if ctfRow.hasAllColumns(CTF_DICT.values()):
        ctfModel = CTFModel()

        rowToObject(ctfRow, ctfModel, CTF_DICT, extraLabels=CTF_EXTRA_LABELS)
        if ctfRow.hasColumn(RELIONCOLUMNS.rlnPhaseShift.value):
            ctfModel.setPhaseShift(
                ctfRow.get(RELIONCOLUMNS.rlnPhaseShift.value, 0))
        ctfModel.standardize()
        setPsdFiles(ctfModel, ctfRow)
    else:
        ctfModel = None

    return ctfModel


def setPsdFiles(ctfModel, ctfRow):
    """ Set the PSD files of CTF estimation related
    to this ctfModel. The values will be read from
    the ctfRow if present.
    """
    for attr, label in CTF_PSD_DICT.items():
        if ctfRow.hasColumn(label):
            setattr(ctfModel, attr, String(ctfRow.get(label)))


def rowToObject(row, obj, attrDict, extraLabels=[]):
    """ This function will convert from a XmippMdRow to an EMObject.
    Params:
        row: the XmippMdRow instance (input)
        obj: the EMObject instance (output)
        attrDict: dictionary with the map between obj attributes(keys) and
            row MDLabels in Xmipp (values).
        extraLabels: a list with extra labels that could be included
            as properties with the label name such as: _rlnSomeThing
    """
    obj.setEnabled(row.get(RELIONCOLUMNS.rlnEnabled.value, 1) > 0)

    for attr, label in attrDict.items():
        value = row.get(label)
        if not hasattr(obj, attr):
            setattr(obj, attr, ObjectWrap(value))
        else:
            getattr(obj, attr).set(value)

    attrLabels = attrDict.values()

    for label in extraLabels:
        if label not in attrLabels and row.hasColumn(label):
            setattr(obj, '_' + label, ObjectWrap(row.get(label)))


def setObjId(obj, mdRow, label=RELIONCOLUMNS.rlnImageId.value):
    obj.setObjId(mdRow.get(label, None))


def rowToParticle(partRow, particleClass=Particle, **kwargs):
    """ Create a Particle from a row of a meta """
    img = particleClass()

    # Provide a hook to be used if something is needed to be
    # done for special cases before converting image to row
    samplingRate = kwargs.get('samplingRate', None)
    if samplingRate:
        img.setSamplingRate(samplingRate)

    preprocessImageRow = kwargs.get('preprocessImageRow', None)
    if preprocessImageRow:
        preprocessImageRow(img, partRow)

    # Decompose Relion filename
    index, filename = cryosparcToLocation(
        partRow.get(RELIONCOLUMNS.rlnImageName.value))
    img.setLocation(index, filename)

    if partRow.hasColumn(RELIONCOLUMNS.rlnClassNumber.value):
        img.setClassId(partRow.get(RELIONCOLUMNS.rlnClassNumber.value))

    if kwargs.get('readCtf', True):
        img.setCTF(rowToCtfModel(partRow))

    # alignment is mandatory at this point, it should be check
    # and detected defaults if not passed at readSetOf.. level
    alignType = kwargs.get('alignType')

    if alignType != ALIGN_NONE:
        samplingRate = kwargs.get('samplingRate')
        img.setTransform(rowToAlignment(partRow, alignType, samplingRate))

    if kwargs.get('readAcquisition', True):
        img.setAcquisition(rowToAcquisition(partRow))

    if kwargs.get('magnification', None):
        img.getAcquisition().setMagnification(kwargs.get("magnification"))

    setObjId(img, partRow)
    # Read some extra labels
    extraLabel = {}
    rowToObject(partRow, img, extraLabel,
                extraLabels=IMAGE_EXTRA_LABELS + kwargs.get('extraLabels', []))

    img.setCoordinate(rowToCoordinate(partRow))

    # copy micId if available from row to particle
    if partRow.hasColumn(RELIONCOLUMNS.rlnMicrographId.value):
        img.setMicId(partRow.get(RELIONCOLUMNS.rlnMicrographId.value))

    # copy particleId if available from row to particle
    if partRow.hasColumn(RELIONCOLUMNS.rlnParticleId.value):
        img._rlnParticleId = Integer(
            partRow.get(RELIONCOLUMNS.rlnParticleId.value))

    # Provide a hook to be used if something is needed to be
    # done for special cases before converting image to row
    postprocessImageRow = kwargs.get('postprocessImageRow', None)
    if postprocessImageRow:
        postprocessImageRow(img, partRow)
    return img


def rowToCoordinate(coordRow):
    """ Create a Coordinate from a row of a meta """
    # Check that all required labels are present in the row
    if coordRow.hasAllColumns(COOR_DICT):
        coord = Coordinate()
        rowToObject(coordRow, coord, COOR_DICT, extraLabels=COOR_EXTRA_LABELS)

        micName = None

        if coordRow.hasColumn(RELIONCOLUMNS.rlnMicrographId.value):
            micId = int(coordRow.get(RELIONCOLUMNS.rlnMicrographId.value))
            coord.setMicId(micId)
            # If RLN_MICROGRAPH_NAME is not present, use the id as a name
            micName = micId

        if coordRow.hasLabel(RELIONCOLUMNS.rlnMicrographName.value):
            micName = coordRow.get(RELIONCOLUMNS.rlnMicrographName.value)

        coord.setMicName(micName)

    else:
        coord = None

    return coord


def rowToAcquisition(acquisitionRow):
    """ Create an acquisition from a row of a meta """
    if acquisitionRow.hasAllColumns(ACQUISITION_DICT):
        acquisition = Acquisition()
        rowToObject(acquisitionRow, acquisition, ACQUISITION_DICT)
    else:
        acquisition = None

    return acquisition


def readSetOfParticles(filename, partSet, **kwargs):
    """read from Relion image meta
        filename: The metadata filename where the images are.
        imgSet: the SetOfParticles that will be populated.
        rowToParticle: this function will be used to convert the row to Object
    """
    for imgRow in emtable.Table.iterRows(filename):
        img = rowToParticle(imgRow, **kwargs)
        partSet.append(img)

    partSet.setHasCTF(img.hasCTF())
    partSet.setAlignment(kwargs['alignType'])


if __name__ == "__main__":
    parser = defineArgs()
    sys.exit(convertCs2Star(parser.parse_args()))
