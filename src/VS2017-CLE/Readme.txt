INSTRUCTIONS FOR COMPILING SWMM5.EXE USING MICROSOFT VISUAL STUDIO 2017
=======================================================================

1. Create a sub-directory named VS2017-CLE under the directory where
   the SWMM 5 Engine source code files are stored and copy
   VS2017-CLE.VCPROJ and swmm5.lib to it.

2. Launch Visual Studio 2017 and use the File >> Open command to open
   the VS2017-CLE.VCPROJ file.

3. Issue the Build >> Configuration Manager command and select the
   Release configuration.

4. Issue the Build >> Build VS2017-CLE command to build SWMM5.EXE
   (which will appear in the Release subdirectory underneath the
   VS2017-CLE directory).
