# vmd-tprreader
Molfileplugin for VMD that lets me read in GROMACS tpr files, with some limitations. These limitations include:

* Only single precision tpr files
* No position restraints or alchemical changes
* The files can't be made by a very old file generator. GROMACS versions prior to 4.0 may be too old to load.
* Zero guarantees it will work on non-linux systems

To use, you are required to recompile VMD (or at least the molfile plugins), having placed the source files under plugins/molfile_plugin/src, and editing the Makefile in plugins/molfile_plugin. Specifically, you will need to add "tprplugin" to STATICPLUGINS, "tprplugin.so" to PLUGINS, "${ARCHDIR}/tprplugin-s.o" to ARCHIVEOBJS, and add the targets required for building:

```
tprplugin.so: ${ARCHDIR}/tprplugin.o
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)
${ARCHDIR}/tprplugin.o: tprplugin.C ${PLUGINAPI}
	$(CXX) $(CXXFLAGS) $(SHLDFLAGS) $(INCDIR) -c $< $(COPTO)$@
${ARCHDIR}/tprplugin-s.o: tprplugin.C ${PLUGINAPI}
	${CXX} ${SCXXFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_tprplugin" -c $< $(COPTO)$@
```
	
Now I think this is a very useful feature, however it is not quite ready to be incorparated into VMD proper, as it depends on <rpc/rpc.h>, which is the standard XDR library on Linux systems, but does not have an equivalent in Windows. If you think this is important, help me out by eliminating the dependencies on xdr_read functions! Some of these are equivalent to utilities found in Gromacs.h in the VMD source tree.
