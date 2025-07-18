.. _SpackUberenv:

Spack and Uberenv
=================

GEOS is transitioning to a new `Uberenv <https://github.com/LLNL/uberenv>`_ and `Spack <https://github.com/spack/spack/>`_ system for building our dependencies. We refer the reader to the `Uberenv documentation <https://uberenv.readthedocs.io/en/latest/>`_ and `Spack documentation <https://spack.readthedocs.io/en/latest/index.html>`_, in particular the Spack documentation sections worth reading are:

* `Specs and dependencies <https://spack.readthedocs.io/en/latest/basic_usage.html#specs-dependencies>`_
* `Virtual dependencies <https://spack.readthedocs.io/en/latest/basic_usage.html#virtual-dependencies>`_
* `Manual compiler configuration <https://spack.readthedocs.io/en/latest/getting_started.html?highlight=compilers.yaml#manual-compiler-configuration>`_
* `External packages <https://spack.readthedocs.io/en/latest/packages_yaml.html#external-packages>`_


Building the dependencies with Uberenv involves running:

.. code-block:: console

    ./scripts/uberenv/uberenv.py --spack-env-file=/path/to/your/spack.yaml --spec="%clang@14.0.6" --prefix=/output/path/to/third-party/GEOS/libraries

The Uberenv command **requires** the ``--spack-env-file`` option to pass in a ``spack.yaml`` configuration file, also known in Spack as an environment file, to tell Spack where pre-installed system libraries and compilers are located. See :ref:`SpackYaml` for more on how to create a ``spack.yaml`` file.

.. warning::
  Without a ``spack.yaml`` configuration file, Spack will ignore any system libraries you have installed and will go down a rabbit hole building dependencies. Without a ``spack.yaml`` configuration file, Spack will not allow you to choose the compiler to build with. The Uberenv command **must** be used with a ``spack.yaml`` configuration file, it **should never be used without one**.

This Uberenv command uses the ``--spec`` option to determine what compiler you want to use, and where you can also specify how you want your dependencies built (see :ref:`BuildConfig`). In the the Uberenv command above, Clang 14.0.6 is chosen from the ``spack.yaml`` file.

This Uberenv command will create a directory specified by ``--prefix=directory-name`` (or a directory named ``uberenv_libs`` if the option is not provided) in the current working directory, clone Spack into ``directory-name/spack`` and install the dependencies into ``directory-name/system_dependent_path``. This Uberenv command will then spit out host-config files (see :ref:`HostConfig`) in the current directory which you can use to build GEOS or LvArray.

For more Uberenv command-line options, you can run the ``uberenv.py`` script with the ``--help`` option or consult the `Uberenv command line options documentation <https://uberenv.readthedocs.io/en/latest/#command-line-options>`_.

.. note::
  There is no requirement that your environment file be named ``spack.yaml`` when it is passed to Uberenv using the ``--spack-env-file`` command line option.

.. note::
  On LC systems only, there is not a requirement to specify ``--spack-env-file``. This is because Uberenv uses the environment variable ``SYS_TYPE`` in combination with the ``.uberenv_config.json`` Uberenv configuration file to determine the folder name that contains the required ``spack.yaml`` file (e.g. ``scripts/spack_configs/blueos_3_ppc64le_ib_p9/spack.yaml``). More information on Uberenv configuration behavior can be found in `Uberenv spack configurations documentation <https://uberenv.readthedocs.io/en/latest/#spack-configurations>`_.

.. _SpackYaml:

spack.yaml
----------

The ``spack.yaml`` configuration file tells Spack where it can find relevant packages and compilers to build GEOS third-party dependencies. Without ``spack.yaml``, building the dependencies will take significantly longer.

There are many examples and resources available for constructing a ``spack.yaml`` file:

* GEOS's LC configuration files for `toss_4_x86_64_ib <https://github.com/GEOS-DEV/thirdPartyLibs/blob/master/scripts/spack_configs/toss_4_x86_64_ib/spack.yaml>`_ and `blueos_3_ppc64le_ib_p9 <https://github.com/GEOS-DEV/thirdPartyLibs/blob/master/scripts/spack_configs/blueos_3_ppc64le_ib_p9/spack.yaml>`_. Additionally, the header of these configuration files include the Spack spec to pass to ``--spec`` for different compilers and package variants.
* LLNL's shared Spack configurations for RADIUSS projects: https://github.com/LLNL/radiuss-spack-configs/tree/main
* NERSC Spack Infrastructure: https://github.com/NERSC/spack-infrastructure/tree/main
* Shared Spack configuration files with other HPC sites: https://github.com/spack/spack-configs
* The documentation list mentioned above: :ref:`SpackUberenv`

spack.yaml from scratch
~~~~~~~~~~~~~~~~~~~~~~~

If the examples and resources listed in :ref:`SpackYaml` are not applicable to your system, or you would like to see what packages are already installed on your system, you can call Uberenv with the following option:

.. code-block:: console

    ./scripts/uberenv/uberenv.py --setup-and-env-only

This command will setup Spack and ask Spack create a ``spack.yaml`` environment file for you. Uberenv will invoke `spack compiler find <https://spack.readthedocs.io/en/latest/getting_started.html#compiler-configuration>`_ and `spack external find <https://spack.readthedocs.io/en/latest/packages_yaml.html#automatically-find-external-packages>`_ to find pre-installed compilers and packages on your system.

.. note::
  These spack commands called underneath are not perfect and are not guaranteed to find every compiler and every package available on your system. In fact, the commands can be non-deterministic, and you can get two different ``spack.yaml`` environment files from two separate invocations! Also, not all the packages found are relevant dependencies or sub-dependencies for GEOS. Of the packages that are relevant, they may not be the right version!

  This command should be used as a first approximation of your system environment, to determine the paths where more suitable compilers and packages are potentially located.

Required package versions in spack.yaml
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the LC configuration files mentioned in :ref:`SpackYaml` section, you will see a list of packages with the ``require`` keyword:

.. code-block:: console

  hypre:
    require: "@git.06da35b1a4b1066a093bc0c6c48aee12bee74cd4"
  ...

This tells Spack that GEOS always `requires <https://spack.readthedocs.io/en/latest/packages_yaml.html#requirements-syntax>`_ a specific commit of ``hypre``, a commit on the latest develop branch in this case. Ideally, package versions should be specified in the `GEOS Spack package file  <https://github.com/GEOS-DEV/thirdPartyLibs/blob/master/scripts/spack_packages/packages/geosx/package.py>`_. However, when a version of a package is newer than what Spack knows about or an unversioned commit is needed, the Spack package syntax cannot express that requirement. As a result:

.. warning::
  Every ``spack.yaml`` file **must have the list of packages with the "require" keyword** to ensure these packages are built with the right versions or commits. See LC configuration files mentioned in :ref:`SpackYaml` section for the list.



Uberenv configuration file
--------------------------

Uberenv needs a `.uberenv_config.json <https://github.com/GEOS-DEV/thirdPartyLibs/blob/master/.uberenv_config.json>`_ configuration file to function as a submodule.
Details on the various configuration options can be found in `Uberenv project configuration documentation <https://uberenv.readthedocs.io/en/latest/#project-configuration>`_. The most notable option for maintenance is ``spack_commit``, which is the Spack commit that Uberenv checkouts to build the dependencies.


pygeosx
-------

.. warning::
  The spack build system for ``pygeosx`` is a work in progress.

It is worth noting that GEOS has `two project json files <https://uberenv.readthedocs.io/en/latest/#project-configuration>`_ (``.uberenv_config.json`` and ``scripts/pygeosx_configs/pygeosx.json``) and two configuration directories for LC systems (``scripts/spack_configs`` and ``scripts/pygeosx_configs``). The ``.uberenv_config.json`` project json file and ``scripts/spack_configs`` directory is for building GEOS dependencies. The ``scripts/pygeosx_configs/pygeosx.json`` project json file and ``scripts/pygeosx_configs`` directory is for building ``pygeosx`` dependencies.This is because ``pygeosx`` has a separate list of required compilers and packages to build from GEOS (e.g. ``pygeosx``'s numpy dependency recommends building with gcc and using openblas for BLAS/LAPACK). When not building ``pygeosx``, other dependencies of GEOS still depend on python. An existing system version of python will work just fine, and can be put in GEOS's ``spack.yaml`` to prevent Spack from building its own verion of python. By default, Uberenv will find and use ``.uberenv_config.json`` to build GEOS, but you can use the ``--project-json`` command line option to target ``scripts/pygeosx_configs/pygeosx.json`` to build ``pygeosx``:

.. code-block:: console

    ./scripts/uberenv/uberenv.py --spack-config-dir=/path/to/your/config/directory/ --spec="%clang@14.0.6" --project-json="scripts/pygeosx_configs/pygeosx.json"

.. note::
    When building ``pygeosx``, Spack will build various python packages, however by default they are not installed in python. There are various ways of accomplishing `this <https://spack.readthedocs.io/en/latest/basic_usage.html#extensions-python-support>`_, but the recommended approach is to use spack environments. Once you build ``pygeosx`` using Uberenv, Spack will create a view that ensures the Spack-built python can find the built python packages. For example, with a default ``uberenv_libs`` directory of dependencies, the path to the view of python will be ``uberenv_libs/._view/*/bin/python3``. If you want to use your ``pygeosx`` python3 executable in GEOS, you will need to update your host-config's ``Python3_ROOT_DIR`` and ``Python3_EXECUTABLE`` to the path to Spack's view of python.


.. _BuildConfig:

Build Configuration
-------------------

.. warning::
	The spack build system is undergoing updates. The ``petsc`` variant and others are still a work in progress.

The GEOS Spack package has a lot of options, or what Spack calls variants, for controlling which dependencies you would like to build and how you'd like them built. The `GEOS Spack package file  <https://github.com/GEOS-DEV/thirdPartyLibs/blob/master/scripts/spack_packages/packages/geosx/package.py>`_ has variants that are marked with ``variant()`` in the file.

For example if you wanted to build with GCC 8.3.1, without Caliper and with Hypre as the Linear Algebra Interface, your spec would be ``%gcc@8.3.1 ~caliper lai=hypre``.

The GEOS Spack package lists out the libraries that GEOS depends ons. These dependencies are marked with ``depends_on()`` in the file.

Using the Spack spec syntax, you can inturn specify variants for each of the dependencies of GEOS. For example, you could modify the spec above to build RAJA in debug mode by using ``%gcc@8.3.1 ~caliper lai=hypre ^raja build_type=Debug``. When building with Uberenv, Spack should print out a table containing the full spec for every dependency it will build. If you would like to look at the variants for say RAJA in more detail, you can find the package file at ``uberenv_libs/spack/var/spack/repos/builtin/packages/raja/package.py``, by using `file finder <https://docs.github.com/en/get-started/accessibility/keyboard-shortcuts#source-code-browsing>`_ on the `Spack Github website <https://github.com/spack/spack>`_, or by searching for the package at https://packages.spack.io/.


.. _HostConfig:

Host-Config Generation
----------------------

The logic for generating the host-configs can be found in the `GEOS spack recipe <https://github.com/GEOS-DEV/thirdPartyLibs/blob/master/scripts/spack_packages/packages/geosx/package.py>`_. The GEOS host-config is generated by the ``geos_hostconfig()`` function, while the LvArray host-config is generated by the ``lvarray_hostconfig()`` function. After successfully building all the third-party dependencies, Spack will call these two functions to populate two host-configs based on information it knows about the dependencies.

.. note::
  The host-config generation is currently based on LC systems, and the generated host-config may be missing or have incorrect details for your system (e.g. choice of ``MPIEXEC_NUMPROC_FLAG``). Please modify the python functions and/or host-configs generated as needed.


LC TPL Build Script
-------------------

On LC systems, it is necessary to update the third-party library installations after a change to the configuration. The ``setupLC-TPL-uberenv.bash`` `script <https://github.com/GEOS-DEV/thirdPartyLibs/blob/master/scripts/setupLC-TPL-uberenv.bash>`_ is used to build the third-party libraries on multiple LC systems using uberenv:

.. code-block:: console

    ./setupLC-TPL-uberenv.bash /path/to/shared/installation/directory

This command will also generate a LvArray and GEOS host-config for each specified machine and compiler combination.

.. note::
  The terminal output from the command may fail to update. In that case, you can track the progress of the builds by looking at the generated \*.log file associated with each machine and compiler combination (e.g. ruby-gcc-12.log).


Adding a Dependency (Advanced)
------------------------------

Adding a dependency to GEOS is straight forward **if** the dependency already builds with Spack. If that is the case, then all you need to do is add a ``depends_on('cool-new-library')`` to the GEOS ``package.py`` file. If however the dependency doesn't have a Spack package, you will have to add one by creating a ``cool-new-library/package.py`` file in the ``scripts/spack_packages/packages`` directory and adding the logic to build it there. For instructions on how to create a package recipe from scratch, Spack has provided a `Spack Packing Guide <https://spack.readthedocs.io/en/latest/packaging_guide.html>`_.

Oftentimes (unfortunately), even when a package already exists in Spack, it might not work out of the box for your system. In this case copy over the existing ``package.py`` file from the Spack repository into ``scripts/spack_packages/packages/cool-new-library/package.py``, as if you were adding a new package, and perform your modifications there. Once you have the package working, copy the package back into the Spack repository and commit+push your changes to Spack.
