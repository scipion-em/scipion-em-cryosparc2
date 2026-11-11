"""Reading a classification's particles only from the job's own alignments: some cryoSPARC versions copy into them
the alignments of earlier jobs (a 2D classification of refined particles also carries their 3D refinement), and
the conversion would give every particle the refinement's single class."""
import os
import tempfile
import unittest

import numpy as np

from cryosparc2.convert import (keepJobAlignments, convertCs2Star, ALIGNMENTS_2D,
                                ALIGNMENTS_3D_CLASSES)

N = 12
CLASSES_2D = np.arange(N, dtype=np.uint32) % 4  # 4 classes, 3 particles each


def _particles(with2D=True, with3D=True, withMulti=False):
    """A .cs like a 2D classification of refined particles: the 2D classes and the refinement's single class."""
    fields = [("uid", "<u8"), ("blob/path", "S40"), ("blob/idx", "<u4"), ("blob/shape", "<u4", (2,)),
              ("blob/psize_A", "<f4")]
    if with3D:
        fields += [("alignments3D/split", "<u4"), ("alignments3D/shift", "<f4", (2,)),
                   ("alignments3D/pose", "<f4", (3,)), ("alignments3D/class_posterior", "<f4"),
                   ("alignments3D/class", "<u4")]
    if with2D:
        fields += [("alignments2D/split", "<u4"), ("alignments2D/shift", "<f4", (2,)),
                   ("alignments2D/pose", "<f4"), ("alignments2D/class_posterior", "<f4"),
                   ("alignments2D/class", "<u4")]
    if withMulti:
        fields += [("alignments3D_multi/class_posterior", "<f4", (2,)), ("alignments3D_multi/class", "<u4", (2,))]
    cs = np.zeros(N, dtype=fields)
    cs["uid"] = np.arange(1, N + 1)
    cs["blob/path"] = b"J1/imported/particles.mrcs"
    cs["blob/idx"] = np.arange(N)
    cs["blob/shape"] = 64
    cs["blob/psize_A"] = 1.5
    if with3D:
        cs["alignments3D/class_posterior"] = 1
        cs["alignments3D/class"] = 0  # a refinement: one class
        cs["alignments3D/pose"] = 0.3
    if with2D:
        cs["alignments2D/class_posterior"] = 1
        cs["alignments2D/class"] = CLASSES_2D
        cs["alignments2D/pose"] = 0.5
    return cs


class TestKeepJobAlignments(unittest.TestCase):

    def test_A2DClassificationLeavesOutTheEarlier3DAlignments(self):
        kept = keepJobAlignments(_particles(), ALIGNMENTS_2D)
        self.assertFalse([n for n in kept.dtype.names if n.startswith("alignments3D")])
        self.assertEqual(kept["alignments2D/class"].tolist(), CLASSES_2D.tolist())
        self.assertEqual(kept["uid"].tolist(), list(range(1, N + 1)))

    def test_A3DClassificationLeavesOutTheEarlierRefinementAnd2D(self):
        kept = keepJobAlignments(_particles(withMulti=True), ALIGNMENTS_3D_CLASSES)
        names = kept.dtype.names
        self.assertFalse([n for n in names if n.startswith(("alignments3D/", "alignments2D/"))])
        self.assertIn("alignments3D_multi/class_posterior", names)

    def test_TheCsStaysAsItIsWithoutTheJobsOwnAlignments(self):
        refinedOnly = _particles(with2D=False)
        self.assertIs(keepJobAlignments(refinedOnly, ALIGNMENTS_2D), refinedOnly)
        self.assertIs(keepJobAlignments(refinedOnly, ALIGNMENTS_3D_CLASSES), refinedOnly)
        plain2D = _particles(with3D=False)
        self.assertIs(keepJobAlignments(plain2D, ALIGNMENTS_2D), plain2D)


class TestConvert2DClassesOfRefinedParticles(unittest.TestCase):
    """The whole conversion (pyem included): every 2D class reaches the star, not the refinement's single one."""

    def test_TheStarKeepsThe2DClasses(self):
        from emtable import Table
        with tempfile.TemporaryDirectory() as folder:
            csFile = os.path.join(folder, "J9_020_particles.cs")
            with open(csFile, "wb") as f:
                np.save(f, _particles())
            starFile = os.path.join(folder, "particles.star")

            convertCs2Star([csFile, starFile], alignments=ALIGNMENTS_2D)
            classes = [int(row.get("rlnClassNumber")) for row in Table.iterRows("particles@" + starFile)]
            self.assertEqual(sorted(classes), sorted((CLASSES_2D + 1).tolist()))

            convertCs2Star([csFile, starFile])  # as before: the refinement's class wins
            classes = [int(row.get("rlnClassNumber")) for row in Table.iterRows("particles@" + starFile)]
            self.assertEqual(set(classes), {1})


if __name__ == "__main__":
    unittest.main()
