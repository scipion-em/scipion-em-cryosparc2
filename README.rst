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
* **3D Homogeneous Refinement**: Rapidly refine a single homogeneous structure to high-resolution and validate using the gold-standard FSC.
* **3D Non uniform Refinement**: Apply non-uniform refinement to acheive higher resolution and map quality
* **Particle Subtraction**: Subtract projections of a masked volume from particles.
* **3D Local Refinement** Naive local refinement.
* **Heterogeneous Refinement (3D Classification)**: Heterogeneous Refinement simultaneously classifies particles and refines structures from n initial structures, usually obtained following an Ab-Initio Reconstruction
* **Local CTF Refinement (per-particle defocus)**: Local CTF Refinement performs per-particle defocus estimation for each particle in a dataset, against a given 3D reference structure.
* **Global CTF Refinement (per-group beam tilt, trefoil, spherical aberration, tetrafoil)**: Global CTF Refinement performs per-exposure-group CTF parameter refinement of higher-order aberrations, against a given 3D reference.
* **Sharppening**: Sharpen a volume following refinement.
* **Helical 3D Refinement**: Reconstruct and refine a homogeneous helical assembly, with or without imposition and refinement of symmetry parameters.
* **3D Homogeneous Refinement(new)**: Rapidly refine a single homogeneous structure to high-resolution and validate using the gold-standard FSC. Using new faster GPU code, and support for higher-order aberration (beam tilt, spherical aberration, trefoil, tetrafoil) correction and per-particle defocus refinement on the fly.
* **3D Non uniform Refinement(new)**: Apply non-uniform refinement to achieve higher resolution and map quality. Specially designed for small proteins and membrane proteins.
* **3D Local Refinement(new)**  Refine a masked region within a consensus structure by allowing particle alignments to vary only slightly.
* **Symmetry Expansion**: Duplicate particles around a point-group symmetry.
* **Homogeneous Reconstruction**: Reconstruct half-maps from input particles with alignments

**Latest plugin version**
=========================

**v3.3.1**
----------
* **new**      : Compatibility with cryoSPARC v3.3.1



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

3. Add the following variables bellow the PACKAGES section at ``~/.config/scipion/scipion.conf`` file:

   .. code-block::

       # The root directory where cryoSPARC code and dependencies is installed.
       CRYOSPARC_HOME = <install_path>   (CRYOSPARC_DIR will work for legacy reasons)
       
       # full name of the initial admin account to be created
       CRYOSPARC_USER = <user_name>

       # Optional variables

       #Folder (available to all workers) where scipion will create cryosparc projects
       CRYO_PROJECTS_DIR = <path> (default to <CRYOSPARC_HOME>/scipion_projects)


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

