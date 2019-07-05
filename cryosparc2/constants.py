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
from pyworkflow.em.constants import (
    SYM_CYCLIC, SYM_TETRAHEDRAL, SYM_OCTAHEDRAL, SYM_I222,
    SYM_I222r)

# Not yet release in Scipion, once is released this try /catch import can be removed
try:
    from pyworkflow.em.constants import SYM_DIHEDRAL_Y
except:
    from pyworkflow.em import SCIPION_SYM_NAME
    SYM_DIHEDRAL_Y = 12
    SCIPION_SYM_NAME[SYM_DIHEDRAL_Y] = "Dyn"

CRYOSPARC_HOME = 'CRYOSPARC_HOME'

# Supported versions:
V2_5 = 'V2_5'
V2_8 = 'V2_8'
V2_9 = 'V2_9'

# Suffix for all project
suffix = 'ProjCryoSparc_'



CS_SYM_NAME = {}
CS_SYM_NAME[SYM_CYCLIC] = 'Cn'
CS_SYM_NAME[SYM_DIHEDRAL_Y] = 'Dn'
CS_SYM_NAME[SYM_TETRAHEDRAL] = 'T'
CS_SYM_NAME[SYM_OCTAHEDRAL] = 'O'
CS_SYM_NAME[SYM_I222] = 'I1'
CS_SYM_NAME[SYM_I222r] = 'I2'