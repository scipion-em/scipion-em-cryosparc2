=================
cryoSPARC2 plugin
=================

This plugin allows to use cryoSPARC2 programs within the Scipion framework

`CryoSPARC <https://cryosparc.com/>`_ is a backend and frontend software system
that provides data processing and image analysis capabilities for single particle
cryo-EM, along with a browser based user interface and command line tools.


**Install this plugin**
=======================

You will need to use `2.0.0 <https://scipion-em.github.io/docs/release-2.0.0/docs/scipion-modes/how-to-install.html>`_ version of Scipion to run these protocols.

* **2D Classification**: Classify particles into multiple 2D classes to facilitate stack cleaning and removal of junk particles.  Also useful as a sanity check to investigate particle quality.
* **3D Homogeneous Refinement**: Rapidly refine a single homogeneous structure to high-resolution and validate using the gold-standard FSC
* **3D Ab-Initio Reconstruction**:  Reconstruct a single (homogeneous) or multiple (heterogeneous) 3-D maps from a set of particles, without any initial models or starting structures required.

For now, the plugin must be installed in developer mode.
For that, follow these instructions:

1. **Install the plugin**


- **Stable version**

.. code-block::

      scipion installp -p scipion-em-cryosparc2

OR

  - through the **plugin manager GUI** by launching Scipion and following **Configuration** >> **Plugins**


- **Developer's version**

- Clone or download the plugin repository

.. code-block::

            git clone https://github.com/scipion-em/scipion-em-cryosparc2.git

- Install the plugin in developer mode.

.. code-block::

    scipion installp -p local/path/to/scipion-em-cryosparc2 --devel


2. Install **CryoSPARC software**

    CryoSPARC v2 software will *NOT* be installed automatically with the plugin. The
    independent installation of CryoSPARC software suite by the user is required
    before running the programs.

    To install CryoSPARC v2 software review the detailed system requirements and install
    instructions available `here <https://cryosparc.com/docs/reference/install/>`_.
    These cover workstation and cluster installs, file configuration and how to update
    cryoSPARC v2 when new versions become available.

3. Add the following variables in the ``scipion.conf`` file:

   .. code-block::

       # The root directory where cryoSPARC code and dependencies will be installed.
       CRYOSPARC_DIR = <install_path>   
       
       # full name of the initial admin account to be created
       CRYOSPARC_USER = <user_name>
       
       # path on the worker node to a writable directory residing on the local SSD
       CRYOSSD_DIR = <ssd_path>





