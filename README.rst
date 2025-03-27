=================
cryoSPARC2 plugin
=================

This plugin allows to use cryoSPARC2 programs within the Scipion framework

`CryoSPARC <https://cryosparc.com/>`_ is a backend and frontend software system
that provides data processing and image analysis capabilities for single particle
cryo-EM, along with a browser based user interface and command line tools.

You will need to use `3.0.0 <https://scipion-em.github.io/docs/release-3.0.0/docs/scipion-modes/how-to-install.html>`_ version of Scipion to run these protocols.

* **2D Classification**: Classify particles into multiple 2D classes to facilitate stack cleaning and removal of junk particles.  Also useful as a sanity check to investigate particle quality.
* **3D Ab-Initio Reconstruction**:  Reconstruct a single (homogeneous) or multiple (heterogeneous) 3-D maps from a set of particles, without any initial models or starting structures required.
* **Particle Subtraction**: Subtract projections of a masked volume from particles.
* **Heterogeneous Refinement**: Heterogeneous Refinement simultaneously classifies particles and refines structures from n initial structures, usually obtained following an Ab-Initio Reconstruction
* **Local CTF Refinement (per-particle defocus)**: Local CTF Refinement performs per-particle defocus estimation for each particle in a dataset, against a given 3D reference structure.
* **Global CTF Refinement (per-group beam tilt, trefoil, spherical aberration, tetrafoil)**: Global CTF Refinement performs per-exposure-group CTF parameter refinement of higher-order aberrations, against a given 3D reference.
* **Sharppening**: Sharpen a volume following refinement.
* **Helical 3D Refinement**: Reconstruct and refine a homogeneous helical assembly, with or without imposition and refinement of symmetry parameters.
* **3D Homogeneous Refinement(new)**: Rapidly refine a single homogeneous structure to high-resolution and validate using the gold-standard FSC. Using new faster GPU code, and support for higher-order aberration (beam tilt, spherical aberration, trefoil, tetrafoil) correction and per-particle defocus refinement on the fly.
* **3D Non uniform Refinement(new)**: Apply non-uniform refinement to achieve higher resolution and map quality. Specially designed for small proteins and membrane proteins.
* **3D Local Refinement(new)**  Refine a masked region within a consensus structure by allowing particle alignments to vary only slightly.
* **Symmetry Expansion**: Duplicate particles around a point-group symmetry.
* **Homogeneous Reconstruction**: Reconstruct half-maps from input particles with alignments
* **3D Classification**: Classify particles into multiple 3D classes and optimize 3D class densities (currently, without re-aligning particle pose or shift).
* **3D Variability Analysis**: Protocol to compute the principle modes of variability with a dataset of aligned particles
* **3D Variability Display**: Protocol to create various versions of a 3D variability result that can be used for display
* **Blob Picker**: Automatically picks particles by searching for Gaussian signals.
* **Patch CTF Estimation**:  Patch-based CTF estimation automatically estimates defocus variation for tilted, bent, deformed samples and is accurate for all particle sizes and types including flexible and membrane proteins.
* **3D Flex Data Prep**: Prepares particles for use in 3DFlex training and reconstruction. At the same  way,  Takes in a consensus (rigid) refinement density map, plus optionally a segmentation and generates a tetrahedral mesh for 3DFlex.
* **3D Flex Mesh Prep**: Takes in a consensus (rigid) refinement density map, plus optionally a segmentation and generates a tetrahedral mesh for 3DFlex. See Mesh Generation below.
* **3D Flex Training**: Uses a mesh and prepared particles (at a downsampled resolution) to train a 3DFlex model. Parameters control the number of latent dimensions, size of the model, and training hyperparameters. This job outputs checkpoints during training.
* **3D Flex Reconstruction**: Takes in a checkpoint from training as well as prepared high-resolution particles and performs high-resolution refinement using L-BFGS under the 3DFlex model. This is the stage at which improvements to density in high-res regions are computed. Outputs two half-maps that can be used for FSC validation, sharpening, and other downstream tasks.
* **3D Flex Generator**: Takes in a checkpoint from training and generates volume series from it, to show what the model is learning about the motion of the particle. This job can be run while training is ongoing to see progress along the way. This job can also optionally take in a high-resolution density map (e.g., from 3D Flex Reconstruction) and  will upscale the deformation model and apply deformations to the high resolution map.


**Latest plugin version**
==========================

**v4.2.1**
-----------

* **fixed**      Fixing errors related the classification protocols

**v4.2.0**
-----------

* **new**        Compatibility with cryoSPARC v4.6.2
* **new**        Adding 3D Flex Generator Protocol
* **updated**    Generating the particles and the consensus volume in the 3D Flex Data Prepare protocol
* **updated**    Validating the mask in the 3d Flex Mesh protocol


**v4.1.6**
-----------

* **fixed**      Fix parameter mask help as described in https://github.com/scipion-em/scipion-em-cryosparc2/issues/170

**v4.1.5**
-----------

* **new**        Compatibility with cryoSPARC v4.6


**v4.1.4**
-----------

* **new**        Compatibility with cryoSPARC v4.5.3
* **new**        Flexibility protocols can launch CS viewer

**v4.1.3**
-----------

* **new**        Compatibility with cryoSPARC v4.5.1
* **new**        Registering flex particles
* **new**        Integration with FlexUtils plugin

**v4.1.2**
-----------
* **fixed**       Tolerating deletion of projects within CS as well as their folders in the file system

* **new**         Add new protocols:
                    * **3D Flex Data Prep**: Prepares particles for use in 3DFlex training and reconstruction. At the same  way,  Takes in a consensus (rigid) refinement density map, plus optionally a segmentation and generates a tetrahedral mesh for 3DFlex.
                    * **3D Flex Mesh Prep**: Takes in a consensus (rigid) refinement density map, plus optionally a segmentation and generates a tetrahedral mesh for 3DFlex. See Mesh Generation below.
                    * **3D Flex Training**: Uses a mesh and prepared particles (at a downsampled resolution) to train a 3DFlex model. Parameters control the number of latent dimensions, size of the model, and training hyperparameters. This job outputs checkpoints during training.
                    * **3D Flex Reconstruction**: Takes in a checkpoint from training as well as prepared high-resolution particles and performs high-resolution refinement using L-BFGS under the 3DFlex model. This is the stage at which improvements to density in high-res regions are computed. Outputs two half-maps that can be used for FSC validation, sharpening, and other downstream tasks.

* **new**         Allowing Scipion to import coordinates


**Installing the plugin**
=========================

In order to install the plugin follow these instructions:

1. **Install the plugin**

.. code-block::

     scipion installp -p scipion-em-cryosparc2

or through the **plugin manager** by launching Scipion and following **Configuration** >> **Plugins**


2. Install **CryoSPARC software**

CryoSPARC v2 software will *NOT* be installed automatically with the plugin. The
independent installation of CryoSPARC software suite by the user is required
before running the programs.

To install CryoSPARC v2 software review the detailed system requirements and install
instructions available `here <https://cryosparc.com/docs/reference/install/>`_.
These cover workstation and cluster installs, file configuration and how to update
cryoSPARC v2 when new versions become available.

3. Add the following variables to the scipion config file (run scipion3 config --show to open it)

   .. code-block::

       # The root directory where cryoSPARC code and dependencies is installed.
       CRYOSPARC_HOME = <install_path>   (CRYOSPARC_DIR will work for legacy reasons)
       
       # full name of the initial admin account to be created
       CRYOSPARC_USER = <user_name>

       # Optional variables
       ---------------------

       # The password with which cryoSPARC was installed.
       # This is only required for the use of the Flexutils plugin and its
       # connection to the 3D flex training protocol.
       CRYOSPARC_PASSWORD = <password>

       #Folder (available to all workers) where scipion will create cryosparc projects
       CRYO_PROJECTS_DIR = <path> (default to <CRYOSPARC_HOME>/scipion_projects)

       # Specifies whether the CS installation is standalone or not. If False,
       # it is assumed that CS is installed in a cluster. If the variable is not
       # defined, by default assume that the installation is standalone and its
       # value would be True
       CRYOSPARC_STANDALONE_INSTALLATION = <True or False>

       # Name of the default lane where the protocols will be launched
       CRYOSPARC_DEFAULT_LANE = <lane name>



**To install in development mode**

- Clone or download the plugin repository

.. code-block::

          git clone https://github.com/scipion-em/scipion-em-cryosparc2.git

- Install the plugin in developer mode.

.. code-block::

  scipion installp -p local/path/to/scipion-em-cryosparc2 --devel

===============
Buildbot status
===============

Status devel version:

.. image:: http://scipion-test.cnb.csic.es:9980/badges/cryosparc2_devel.svg

Status production version:

.. image:: http://scipion-test.cnb.csic.es:9980/badges/cryosparc2_prod.svg

