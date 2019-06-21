=================
cryoSPARC2 plugin
=================

This plugin allows to use cryoSPARC2 programs within the Scipion framework

`CryoSPARC <https://cryosparc.com/>`_ is a backend and frontend software system
that provides data processing and image analysis capabilities for single particle
cryo-EM, along with a browser based user interface and command line tools.


**Install this plugin**
-----------------------

You will need to use `2.0.0 <https://scipion-em.github.io/docs/release-2.0.0/docs/scipion-modes/how-to-install.html>`_ version of Scipion to run these protocols.

For now, is recommended that the plugin be installed in developer mode.
For that, follow these instructions:

1. Clone or download the plugin repository

.. code-block::

            git clone https://github.com/scipion-em/scipion-em-cryosparc2.git

2. Install the plugin in developer mode.

.. code-block::

    scipion installp -p local/path/to/scipion-em-cryosparc2 --devel


3. Install **CryoSPARC software**

    CryoSPARC v2 software will *NOT* be installed automatically with the plugin. The
    independent installation of CryoSPARC software suite by the user is required
    before running the programs.

    To install CryoSPARC v2 software review the detailed system requirements and install
    instructions available `here <https://cryosparc.com/docs/reference/install/>`_.
    These cover workstation and cluster installs, file configuration and how to update
    cryoSPARC v2 when new versions become available.

4. Add the following variables in the ``scipion.conf`` file:

   .. code-block::

       CRYOSPARC_DIR = <install_path>
       CRYOSPARC_USER = <user name with which cryosparc was installed>
       CRYOSSD_DIR = <path on the worker node to a writable directory residing on the local SSD>





