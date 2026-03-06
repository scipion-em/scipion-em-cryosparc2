# **************************************************************************
# *
# * Authors:     David Herreros (dherreros@cnb.csic.es)     [2]
# *
# * [1] MRC Laboratory of Molecular Biology, MRC-LMB
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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


import os
import shutil
import numpy as np

import cryosparc2
from cryosparc2.utils import generateFlexVolumes


class HeterogeneityProgramInterface:
    def __init__(self, _path_template: str, _program_loading_params: dict):
        self.model = self.prepare_heterogeneity_program(**_program_loading_params)
        self.path_template = _path_template

    def prepare_heterogeneity_program(self, **kwargs) -> object:
        model = {"projectId": kwargs.pop("projectId"), "workSpaceId": kwargs.pop("workSpaceId"),
                 "projectPath": kwargs.pop("projectPath"),
                 "trainJobId": kwargs.pop("trainJobId"), "csGPU": kwargs.pop("csGPU")}
        return model

    def decode_state_from_latent(self, latent: np.array) -> None:
        cryosparc2.Plugin._defineVariables()
        flexGeneratorJob = generateFlexVolumes(latent, self.model["projectId"],
                                               self.model["workSpaceId"],
                                               self.model["trainJobId"],
                                               gpu=self.model["csGPU"])
        flexGeneratorJob = str(flexGeneratorJob.get())
        for idx in range(latent.shape[0]):
            volume_path = os.path.join(self.model["projectPath"], flexGeneratorJob,
                                       flexGeneratorJob + "_series_000",
                                       flexGeneratorJob + "_series_000_frame_{:03d}.mrc".format(idx))
            shutil.copyfile(volume_path, self.path_template.format(idx + 1))
