Qt5 Installation
================

At its core, MountainLab uses Qt5, a cross-platform C++ toolkit. You'll need it in order to compile MountainLab.

Using the package manager
-------------------------

If you are on a later version of Ubuntu (such as 16.04), you can get away with installing Qt5 using the package manager (which is great news, see :doc:`mountainlab_installation`). Otherwise, that method may not give a recent enough Qt5 version. In that case (for example if you are on 14.04, or using Mac or other Linux flavor), you should install Qt5 via qt.io website (see below).

If you've got Ubuntu 16.04 or later (good news). See :doc:`mountainlab_installation` 

Installing Qt5 from the www.qt.io
---------------------------------

As mentioned above you are not using a later version of Ubuntu, you probably won't get a recent enough version from the package manager. In that case follow these directions to install a recent version of Qt5:

Go to https://www.qt.io/download/ and click through that you want to install the open source version. Download and run the appropriate installer. Note that you do not need to set up a Qt account -- you can skip that step in the install wizard.

We recommend installing this in your home directory, which does not require admin privileges.

Once installed you will need to prepend the path to qmake to your PATH environment variable. On my system that is /home/magland/Qt/5.7/gcc_64/bin. You may instead do sudo ln -s /home/magland/Qt/5.7/gcc_64/bin/qmake /usr/local/bin/qmake.

Anaconda users may need to un-export anaconda/miniconda path in order to make qt5 from the operating system available rather than the one supplied with anaconda. To do this, edit your ~/.bashrc file, comment out the export command containing anaconda or miniconda path, and open a new terminal. Make sure you are using your OS' installation by running which qmake