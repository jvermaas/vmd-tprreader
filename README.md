# vmd-tprreader
Molfileplugin for VMD that lets me read in GROMACS tpr files, with some limitations. These limitations include:

* Only single precision tpr files
* No alchemical changes
* The files can't be made by a very old file generator. GROMACS versions prior to 4.0 may be too old to load.
* Zero guarantees it will work on non-linux systems, although the current version works under windows.

To use, you are required to recompile VMD (or at least the molfile plugins), having placed the source files under plugins/molfile_plugin/src, and editing the Makefile in plugins/molfile_plugin. Specifically, you will need to add "tprplugin" to STATICPLUGINS, "tprplugin.so" to PLUGINS, "${ARCHDIR}/tprplugin-s.o" to ARCHIVEOBJS, and add the targets required for building:

```
tprplugin.so: ${ARCHDIR}/tprplugin.o
	$(SHLD) $(LOPTO)${ARCHDIR}/$@ $? $(LDFLAGS)
${ARCHDIR}/tprplugin.o: tprplugin.C ${PLUGINAPI}
	$(CXX) $(CXXFLAGS) $(SHLDFLAGS) $(INCDIR) -c $< $(COPTO)$@
${ARCHDIR}/tprplugin-s.o: tprplugin.C ${PLUGINAPI}
	${CXX} ${SCXXFLAGS} $(INCDIR) $(DEF)"VMDPLUGIN=molfile_tprplugin" -c $< $(COPTO)$@
```
	
Now I think this is a very useful feature, and is getting closer to being incorporated into VMD proper.
Currently, we use the GROMACS internal xdr headers, which while intended to work on Windows, also work on Linux to avoid depending on `<rpc.h>`.

## Testing

If you want to see how this works on a set of tprs, you can do the following:

```bash
make
./tprtest testtprs/*
```

That would tell you if your tpr files are readable or not, along with some diagnostic information.
