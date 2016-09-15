#define VMDPLUGIN_STATIC
#include "molfile_plugin.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <rpc/rpc.h>

/* I am making a very concious decision here. I'm only going to TRY supporting gromacs TPX files that were written by 
version 4.0 or later. This means TPX format version 58 and later.
*/
#define STRLEN 4096
#define SAVELEN 8
#include "xdrread.h"
#include "ffread.h"
//This holds the meat and potatoes, the results after extracting all the topology info from the tpr file.
struct tprdata {
	int version; //Version of the file format (fver)
	int wversion; //Version of the generation code that made the file. (Equivalent to fgen)
	int natoms, ngtc; //This is the TOTAL number of atoms.
	float boxdims[9];
	char *symtab;//Looks like internal gromacs names. Truncate to 8 charachters for storage
	//(longer ones aren't related to atom names, types or residues)
	int symtablen, nmoltypes, nmolblock;
	int *atomsinmol, *resinmol;
	int *molnames;
	float **charges, **masses;
	int **resids, **ptypes;
	unsigned short **types;
	int **atomnameids, **atomtypeids, **resnames, **atomicnumbers;
	int **interactionlist[F_NRE];
	int *nr[F_NRE];
	int *molbtype, *molbnmol, *molbnatoms;
	XDR* xdrptr;
	FILE *f;
	int readcoordinates;
};

template<typename real>
void readff(XDR* xdrs, int version) {
	int atnr, ntypes, i, j;
	int *functype;
	double reppow =12.0;
	float fudge;
	xdr_int(xdrs, &atnr);
	//printf("Atnr: %d\n", atnr);
	xdr_int(xdrs, &ntypes);
	//printf("Ntypes: %d\n", ntypes);
	functype = (int *) malloc(sizeof(int) * ntypes);
	for (i=0; i < ntypes; i++)
		xdr_int(xdrs, &functype[i]);
	if (version >= 66)
		xdr_double(xdrs, &reppow);
	xdr_float(xdrs, &fudge);
	for (i = 0; i < ntypes; i++) {
		for (j = 0; (j < NFTUPD); j++) {
			if (version < ftupd[j].fvnr && functype[i] >= ftupd[j].ftype)
				functype[i] += 1;
		}
		readparams<float>(xdrs, version, functype[i]);
	}
	free(functype);
}
//do_atomtypes
void read_atomtypes(XDR* xdrs, int version) {
	int i;
	int numtypes = readInt(xdrs);
	//printf("Num atom types:%d\n", numtypes);
	//Read a bunch of stuff that is set optionally.
	for (i = 0; i < numtypes; i++) {
		readReal<float>(xdrs);//Radius?
		readReal<float>(xdrs);//Volume?
		readReal<float>(xdrs);//Surface tension?
		if (version >= 40) {
			readReal<float>(xdrs);
		}
		if (version >= 60) {
			readReal<float>(xdrs);
			readReal<float>(xdrs);
		}
	}
}
//Reading in this case means discarding.
void read_cmap(XDR* xdrs) {
	int i, ngrid, gs;
	ngrid = readInt(xdrs);
	gs = readInt(xdrs);
	//printf("ngrid: %d, gs: %d\n", ngrid, gs);
	//printf("This many reals should be read: %d\n", ngrid * gs * gs);
	for (i = 0; i < ngrid * gs * gs; i++) {
		readReal<float>(xdrs);
		readReal<float>(xdrs);
		readReal<float>(xdrs);
		readReal<float>(xdrs);
	}
}
void read_groups(XDR* xdrs, int ngrp) {
	int g, i, j, tmp;
	//do_grps
	//printf("ngrp = %d\n", ngrp);
	for(i = 0; i < ngrp; i++) {
		tmp = readInt(xdrs);
		//printf("Number of reads = %d\n", tmp);
		for (j = 0; j < tmp; j++) {
			readInt(xdrs);
		}
	}
	//ngrpname
	tmp = readInt(xdrs);
	//printf("ngrpname = %d\n", tmp);
	for (i = 0; i < tmp; i++) {
		g = readInt(xdrs);
		//printf("symtabid = %d\n", g);
	}
	for (g = 0; g < ngrp; g++) {
		tmp = readInt(xdrs);
		//printf("This many chars: %d\n", tmp);
		if (tmp) {
			for (i = 0; i < tmp; i++) {
				readChar(xdrs);
			}
		}
	}
}
int readtprAfterPrecision (tprdata *tpr) {
	//Start at do_tpx, trace down to do_tpxheader
	XDR *xdrs = tpr->xdrptr;
	tpr->readcoordinates = 0;
	float dummy;
	int tmp, i, j,k;
	unsigned int ui;
	int hasIR, hasCoord, hasV, hasF, hasTop, hasDim, hasIntermoleculeBonds, hasGBSA;
	int numcmap;
	float boxoffset[9];
	float boxv[9];
	unsigned short sdummy;
	bool bClear;
	int sum = 0;
	tpr->version = readInt(xdrs);
	printf("File Format Version: %d\n", tpr->version);

	if (tpr->version >= 77 && tpr->version <= 79) {
		printString(xdrs);
	}

	tpr->wversion = readInt(xdrs);
	printf("Generator Version: %d\n", tpr->wversion);
	
	if (tpr->version >= 81) {
		printString(xdrs);
		//I dunno why this is wrong. Should say "release", but it totally doesn't.
		xdr_setpos(xdrs, 4 + xdr_getpos(xdrs));
	}
	//Bailouts if things are too new/we can't guarantee accurately reading them.
	if (tpr->wversion > 26 || tpr->version <= 57) {
		printf("Your file cannot be read, as it has version %d, but we can read from version 57 to at least 103.\n", tpr->version);
		printf("The generator version for your file is %d, but we can only read up to 26\n", tpr->wversion);
		return MOLFILE_ERROR;
	}

	//If we were dealing with TPA files, we'd need to do something here. Not a clue what.
	//do_section line in the gromacs source.
	tpr->natoms = readInt(xdrs);
	//printf("Natoms: %d\n", tpr->natoms);
	tpr->ngtc = readInt(xdrs);
	//printf("Ngtc: %d\n", tpr->ngtc);
	if (tpr->version < 62) {
		i = readInt(xdrs);
		dummy = readReal<float>(xdrs);
	}

	if (tpr->version >=79) {
		xdr_int(xdrs, &tmp);
		//printf("FEP state: %d\n", tmp);
	}
	dummy = readReal<float>(xdrs);
	//printf("Lambda: %f\n", dummy);
	hasIR = readInt(xdrs);
	//printf("IR state?: %d\n", hasIR);
	hasTop = readInt(xdrs);
	//printf("Contains topology: %d\n", hasTop);
	hasCoord = readInt(xdrs);
	//printf("Contains coordinates: %d\n", hasCoord);
	hasV = readInt(xdrs);
	//printf("Contains velocities: %d\n", hasV);
	hasF = readInt(xdrs);
	//printf("Contains forces: %d\n", hasF);
	hasDim = readInt(xdrs);
	//printf("Contains dimensions: %d\n", hasDim);
	//End of do_tpxheader
	if (hasDim) {
		readvector(xdrs, (tpr->boxdims), 9);
		//printf("Box dim: (%f %f %f), (%f %f %f), (%f %f %f)\n", tpr->boxdims[0],tpr->boxdims[1],tpr->boxdims[2],
		//	tpr->boxdims[3],tpr->boxdims[4],tpr->boxdims[5],tpr->boxdims[6],tpr->boxdims[7],tpr->boxdims[8]);
		if (tpr->version >= 51) {
			readvector(xdrs, boxoffset, 9);
		}
		else {
			boxoffset[0] = boxoffset[1] = boxoffset[2] = 0;
		}
		//printf("Box offset: (%f %f %f), (%f %f %f), (%f %f %f)\n", boxoffset[0],boxoffset[1],boxoffset[2],
		//	boxoffset[3],boxoffset[4],boxoffset[5],boxoffset[6],boxoffset[7],boxoffset[8]);
		readvector(xdrs, boxv, 9);
		//printf("Box vel: (%f %f %f), (%f %f %f), (%f %f %f)\n", boxv[0],boxv[1],boxv[2],
		//	boxv[3],boxv[4],boxv[5],boxv[6],boxv[7],boxv[8]);
	}
	//Dump the thermostat storage stuff. VMD wouldn't know what to do with it anyway.
	if (tpr->ngtc) {
		if (tpr->version < 69) {
			for (i = 0; i < tpr->ngtc; i++) {
				dummy = readReal<float>(xdrs);
				//printf("NGTC %d: %f\n", i, dummy);
			}
		}
		for (i = 0; i < tpr->ngtc; i++) {
			dummy = readReal<float>(xdrs);
			//printf("NGTC %d: %f\n", i, dummy);
		}
	}
	//do_mtop starts here, which starts by reading the symtab (do_symtab)
	tpr->symtablen = readInt(xdrs);
	//printf("Symtab length: %d\n", tpr->symtablen);
	tpr->symtab = (char*) malloc(SAVELEN * tpr->symtablen * sizeof(char));
	for (i = 0; i < tpr->symtablen; i++) {
		xdr_setpos(xdrs, xdr_getpos(xdrs)+4);
		saveString(xdrs, &(tpr->symtab[SAVELEN * i]));
	}
	/*for (i = 0; i < tpr->symtablen; i++) {
		for (j = 0; j < SAVELEN; j++)
			printf("%c", tpr->symtab[SAVELEN*i+j]);
		printf("\n");
	}*/
	tmp = readInt(xdrs);
	/*printf ("System name: ");
	for (j = 0; j < SAVELEN; j++)
		printf("%c", tpr->symtab[SAVELEN*tmp+j]);
	printf("\n");*/
	//Now "read" in the forcefield. This SHOULD be do_ffparams
	//printf("%d\n", xdr_getpos(xdrs));
	readff<float>(xdrs, tpr->version);
	//printf("%d\n", xdr_getpos(xdrs));

	//Read in the type of molecules.
	tpr->nmoltypes = readInt(xdrs);
	//printf("nmoltypes: %d\n", tpr->nmoltypes);
	tpr->atomsinmol = (int*) malloc(sizeof(int) * tpr->nmoltypes);
	//printf("Malloced a thing\n");
	tpr->molnames = (int*) malloc(sizeof(int) * tpr->nmoltypes);
	tpr->resinmol = (int*) malloc(sizeof(int) * tpr->nmoltypes);
	tpr->charges = (float**) malloc(sizeof(float*) * tpr->nmoltypes);
	tpr->masses = (float**) malloc(sizeof(float*) * tpr->nmoltypes);
	tpr->types = (unsigned short**) malloc(sizeof(unsigned short*) * tpr->nmoltypes);
	tpr->ptypes = (int**) malloc(sizeof(int*) * tpr->nmoltypes);
	tpr->resids = (int**) malloc(sizeof(int*) * tpr->nmoltypes);
	tpr->atomnameids = (int**) malloc(sizeof(int*) * tpr->nmoltypes);
	tpr->atomtypeids = (int**) malloc(sizeof(int*) * tpr->nmoltypes);
	tpr->atomicnumbers = (int**) malloc(sizeof(int*) * tpr->nmoltypes);
	tpr->resnames = (int**) malloc(sizeof(int*) * tpr->nmoltypes);
	for (j = 0; j < F_NRE; j++) {
		tpr->interactionlist[j] = (int**) malloc(sizeof(int*) * tpr->nmoltypes);
		tpr->nr[j] = (int*) malloc(sizeof(int) * tpr->nmoltypes);
	}
	
	//printf("Malloced things\n");
	//do_moltype
	for (i = 0; i < tpr->nmoltypes; i++) {
		tpr->molnames[i] = readInt(xdrs);
		/*printf ("Mol name %d: ", i);
		for (j = 0; j < SAVELEN; j++)
			printf("%c", tpr->symtab[SAVELEN*tpr->molnames[i]+j]);
		printf("\n");*/
		tpr->atomsinmol[i] = readInt(xdrs);
		tpr->resinmol[i] = readInt(xdrs);
		tpr->charges[i] = (float*) malloc(sizeof(float) * tpr->atomsinmol[i]);
		tpr->masses[i] = (float*) malloc(sizeof(float) * tpr->atomsinmol[i]);
		tpr->types[i] = (unsigned short*) malloc(sizeof(unsigned short) * tpr->atomsinmol[i]);
		tpr->ptypes[i] = (int*) malloc(sizeof(int) * tpr->atomsinmol[i]);
		tpr->resids[i] = (int*) malloc(sizeof(int) * tpr->atomsinmol[i]);
		tpr->atomicnumbers[i] = (int*) malloc(sizeof(int) * tpr->atomsinmol[i]);
		//printf("%d atoms, %d residues in molecule %d\n", tpr->atomsinmol[i], tpr->resinmol[i], i);
		for (j=0; j < tpr->atomsinmol[i]; j++) {
			tpr->masses[i][j] = readReal<float>(xdrs);
			tpr->charges[i][j] = readReal<float>(xdrs);
			readReal<float>(xdrs); //mB
			readReal<float>(xdrs); //cB
			xdr_u_short(xdrs, &(tpr->types[i][j]));
			xdr_u_short(xdrs, &sdummy);//typeB
			tpr->ptypes[i][j] = readInt(xdrs);
			tpr->resids[i][j] = readInt(xdrs);
			if (tpr->version >= 52) {
				tpr->atomicnumbers[i][j] = readInt(xdrs);
			}
			//if (i == 0)
			//	printf("%d: mass %f charge %f type %d ptype %d resid %d, periodic table number %d\n", j,
			//		masses[i][j], charges[i][j], types[i][j], ptypes[i][j],resids[i][j], tmp);
		}
		tpr->atomnameids[i] = (int*) malloc(sizeof(int) * tpr->atomsinmol[i]);
		readintvector(xdrs, tpr->atomnameids[i], tpr->atomsinmol[i]);
		tpr->atomtypeids[i] = (int*) malloc(sizeof(int) * tpr->atomsinmol[i]);
		readintvector(xdrs, tpr->atomtypeids[i], tpr->atomsinmol[i]);
		for (j = 0; j < tpr->atomsinmol[i]; j++) {
			readReal<float>(xdrs);
		}
		//Read residues
		tpr->resnames[i] = (int*) malloc(sizeof(int) * tpr->resinmol[i]);
		//printf("%d residues in mol %d\n", tpr->resinmol[i], i);
		for (j = 0; j < tpr->resinmol[i]; j++) {
			tpr->resnames[i][j] = readInt(xdrs);
			if (tpr->version >=63)
				xdr_setpos(xdrs, 8 + xdr_getpos(xdrs));
			else
				tpr->resnames[i][j] += 1;
			/*if ( i == 0) {
				for (k = 0; k < SAVELEN; k++)
					printf("%c", tpr->symtab[SAVELEN*tpr->resnames[i][j]+k]);
				printf("\n");
			}*/
		}
		//do_ilists
		//printf ("Atoms in mol %d: %d\n", i, tpr->atomsinmol[i]);
		//Here is where we would read the i-lists, however we only need a subset.
		//Based on the stuff in src/gmxlib/ifunc.c, we only need bonds (index 0),
		//angles (index 10), dihedrals (index 18 and 19), impropers (index 21 and 22), and CMAP cross terms (index 25)
		for (j = 0; j < F_NRE; j++) {
			bClear = FALSE;
			for (k = 0; k < NFTUPD; k++) {
				if (tpr->version < ftupd[k].fvnr && j == ftupd[k].ftype) {
					bClear = TRUE;
				}
			}
			if (bClear) {
				tpr->nr[j][i] = 0;
			}
			else {
				tpr->nr[j][i] = readInt(xdrs);
				//printf("j, k, interactions: %d, %d, %d\n", j, k, tpr->nr[j][i]);
			}
			if (tpr->nr[j][i] != 0) {
				tpr->interactionlist[j][i] = (int*) malloc(sizeof(int) * tpr->nr[j][i]);
				//printf("Interaction id: %d number of interactants %d\n", j, tpr->nr[j][i]);
				readintvector(xdrs, tpr->interactionlist[j][i], tpr->nr[j][i]);
				/*for (k=0; k < tpr->nr[j][i]; k++) {
					printf("%d ", tpr->interactionlist[j][i][k]);
				}
				printf("\n");*/
			}
			else {
				tpr->interactionlist[j][i] = NULL;
			}
		}
		//end do_ilists


		//do_block and do_blocka. Remove these. VMD doesn't know what they are.
		//They refer to atomgroups and exclusions. Useful for dynamics, not useful for visualization.
		k = readInt(xdrs);
		for (j = 0; j <= k; j++) {
			readInt(xdrs);
		}
		k = readInt(xdrs);
		tmp = readInt(xdrs);
		for (j = 0; j <= k; j++) {
			readInt(xdrs);
		}
		for (j = 0; j < tmp; j++) {
			readInt(xdrs);
		}
	}
	//end of do_moltype
	tpr->nmolblock = readInt(xdrs);
	//printf("nmolblock: %d\n", tpr->nmolblock);
	//do_molblock
	tpr->molbtype = (int*) malloc(tpr->nmolblock * sizeof(int));
	tpr->molbnmol = (int*) malloc(tpr->nmolblock * sizeof(int));
	tpr->molbnatoms = (int*) malloc(tpr->nmolblock * sizeof(int));
	for (i = 0; i < tpr->nmolblock; i++) {
		tpr->molbtype[i] = readInt(xdrs);
		tpr->molbnmol[i] = readInt(xdrs);
		tpr->molbnatoms[i] = readInt(xdrs);
		k = readInt(xdrs);
		//printf("posresXA: %d\n", k);
		for (j = 0; j < k; j++) {
			readReal<float>(xdrs);
		}
		k = readInt(xdrs);
		//printf("posresXB: %d\n", k);
		for (j = 0; j < k; j++) {
			readReal<float>(xdrs);
		}
		//printf("Segname: %d ", tpr->molbtype[i]);
		/*for (j = 0; j < SAVELEN; j++) {
			printf("%c", tpr->symtab[SAVELEN*tpr->molnames[tpr->molbtype[i]]+j]);
		}*/
		//printf("\n\tNumatoms: %d\n", tpr->molbnatoms[i]);
		//printf("\tNumcopies: %d\n", tpr->molbnmol[i]);
		//sum += tpr->molbnatoms[i] * tpr->molbnmol[i];
	}
	//Burn off the number of atoms after the do_molblock
	tmp = readInt(xdrs);
	//printf("What is this? %d\n", tmp);
	if (tpr->version >= 103) {
		hasIntermoleculeBonds = readInt(xdrs);
		//printf("intermolecularbondeds %d\n", hasIntermoleculeBonds);
		if (hasIntermoleculeBonds) {
			printf("Systems with intermolecular bonds are not supported. The file reports %d intermolecular bonds.\n", hasIntermoleculeBonds);
			return MOLFILE_ERROR;
		}
	}
	//do_atomtypes
	read_atomtypes(xdrs, tpr->version);
	//printf("Reading cmap terms\n");
	//do_cmap
	if (tpr->version >= 65) {
		read_cmap(xdrs);
	}
	//do_groups
	//printf("Reading groups\n");
	read_groups(xdrs, egcNR);
	//printf("Returning control\n");
	return MOLFILE_SUCCESS;
}



static int read_tpr_structure(void *mydata, int *optflags, molfile_atom_t *atoms) {
	tprdata *tpr = (tprdata *) mydata;
	//printf("Structure stuff\n");
	*optflags = MOLFILE_MASS | MOLFILE_CHARGE; // | MOLFILE_ATOMICNUMBER; //For the coarse grained residues I'm testing, the atomic number is -1. That is bad.
	int idx = 0;
	int i, j, k;
	for (i = 0; i < tpr->nmolblock; i++) {
		for (j = 0; j < tpr->molbnmol[i]; j++) {
			for (k = 0; k < tpr->molbnatoms[i]; k++) {
				//printf("%d %d %d\n", i, j, k);
				//printf("%s\n",&(tpr->symtab[SAVELEN*tpr->molnames[tpr->molbtype[i]]]));
				strcpy(atoms[idx].segid, &(tpr->symtab[SAVELEN*tpr->molnames[tpr->molbtype[i]]]));
				//printf("%s\n",&(tpr->symtab[SAVELEN*tpr->atomnameids[tpr->molbtype[i]][k]]));
				strcpy(atoms[idx].name, &(tpr->symtab[SAVELEN*tpr->atomnameids[tpr->molbtype[i]][k]]));
				//printf("%s\n",&(tpr->symtab[SAVELEN*tpr->atomtypeids[tpr->molbtype[i]][k]]));
				strcpy(atoms[idx].type, &(tpr->symtab[SAVELEN*tpr->atomtypeids[tpr->molbtype[i]][k]]));
				//printf("%d\n",tpr->resids[tpr->molbtype[i]][k]);
				atoms[idx].resid = tpr->resids[tpr->molbtype[i]][k];
				//printf("%s\n",&(tpr->symtab[SAVELEN*tpr->resnames[tpr->molbtype[i]][atoms[idx].resid]]));
				strcpy(atoms[idx].resname, &(tpr->symtab[SAVELEN*tpr->resnames[tpr->molbtype[i]][atoms[idx].resid]]));
				//printf("%f\n",tpr->masses[tpr->molbtype[i]][k]);
				atoms[idx].mass = tpr->masses[tpr->molbtype[i]][k];
				//printf("%f\n",tpr->charges[tpr->molbtype[i]][k]);
				atoms[idx].charge = tpr->charges[tpr->molbtype[i]][k];
				//atoms[idx].atomicnumber = tpr->atomicnumbers[tpr->molbtype[i]][k];
				idx++;
			}
		}
	}
	//printf("I'm done here\n");
	return MOLFILE_SUCCESS;
}

static void *open_tpr_read(const char *filename, const char *, 
	int *natoms) {
	tprdata *tprdat = NULL;
	FILE *fin;
	XDR *xdrs = new XDR;
	int i, j, tmp;
	int precision;
	if (!(fin = fopen(filename, "r"))) {
		fprintf(stderr, "tprplugin) Cannot open tpr file '%s'\n", filename);
		return NULL;
	}
	xdrstdio_create(xdrs, fin, XDR_DECODE);
	xdr_int(xdrs, &i);
	//printf("%d\n", i);
	printString(xdrs);
	precision = readInt(xdrs);
	//printf("%d\n", precision);
	//printf("I have the precision! %d\n", precision);
	if (precision == 4) {
		tprdat = new tprdata;
		tprdat->f = fin;
		tprdat->xdrptr = xdrs;
		if (readtprAfterPrecision(tprdat) != MOLFILE_SUCCESS) {
			delete tprdat;
			return NULL;
		}
		*natoms = tprdat->natoms;
		//printf("Total number of atoms: %d, %d\n", *natoms, tprdat->natoms);
		//printf("Finished initial reading\n");
		return tprdat;
	}
	else {
		printf("Illegal precision (requires single)\n");
		return NULL;
	}
	fclose(fin);
	return NULL;
}

static int read_tpr_timestep(void *v, int natoms, molfile_timestep_t *ts) {
	tprdata *tpr = (tprdata *)v;
	XDR *xdrs = tpr->xdrptr;
	int tmp;
	if (tpr->readcoordinates) {
		return MOLFILE_ERROR;
	}
	//printf("tpr->natoms %d\n", tpr->natoms);
	//Get the positions.
	//printf("Reading timestep. tsptr = %d, tscoords =%d\n", ts, ts->coords);
	if (ts != NULL) {
		//printf("I'm in here!\n");
		//tmp = readInt(xdrs);
		/*printf("%d, %d\n", tmp, natoms);
		if (tmp != natoms) {
			printf("This should not be possible. But it happened!\n");
			return MOLFILE_ERROR;
		}*/
		readvector(xdrs, ts->coords, 3 * tpr->natoms);
		for (int i = 0; i < 3 * natoms; i++) {
			ts->coords[i] = 10 * ts->coords[i];
		}
		/*for (int i = 0; i < natoms; i+=1000) {
			printf("Atom %d: %f %f %f\n", i, ts->coords[3*i+0], ts->coords[3*i+1],ts->coords[3*i+2]);
			printf("Atom %d: %f %f %f\n", i+1, ts->coords[3*i+3], ts->coords[3*i+4],ts->coords[3*i+5]);
		}*/
		/*for (int i = 0; i < 3*natoms; i++) {
			if (! isfinite(ts->coords[i])) {
				printf("The %d coordinate of atom %d is not finite!\n", (i % 3), i / 3);
			}
		}*/
		ts->A = tpr->boxdims[0] * 10;
		ts->B = tpr->boxdims[4] * 10;
		ts->C = tpr->boxdims[8] * 10;
		ts->alpha = 90 - asin(tpr->boxdims[3] * tpr->boxdims[6] + tpr->boxdims[4] * tpr->boxdims[7] + tpr->boxdims[5] * tpr->boxdims[8]) * 90.0 / M_PI_2;
		ts->beta  = 90 - asin(tpr->boxdims[0] * tpr->boxdims[6] + tpr->boxdims[1] * tpr->boxdims[7] + tpr->boxdims[2] * tpr->boxdims[8]) * 90.0 / M_PI_2;
		ts->gamma = 90 - asin(tpr->boxdims[0] * tpr->boxdims[3] + tpr->boxdims[1] * tpr->boxdims[4] + tpr->boxdims[2] * tpr->boxdims[5]) * 90.0 / M_PI_2;
		//printf("%f %f %f %f %f %f\n", ts->A, ts->B, ts->C, ts->alpha, ts->beta, ts->gamma);
	}
	tpr->readcoordinates = 1;
	return MOLFILE_SUCCESS;
}
static int read_tpr_bonds(void *v, int *nbonds, int **fromptr, int **toptr, 
						 float **bondorder, int **bondtype, 
						 int *nbondtypes, char ***bondtypename) {
	tprdata *tpr = (tprdata *)v;
	int i, j, k, l, itraction, mtype;
	const int bondinteractions = 4;
	int bondinteraction[bondinteractions] = {F_BONDS,F_G96BONDS,F_CONSTR,F_SETTLE};
	int aoffset = 0;
	int boffset = 0;
	*nbonds = 0;
	*fromptr = NULL;
	*toptr = NULL;
	*bondorder = NULL;
	*bondtype = NULL;
	*nbondtypes = 0;
	*bondtypename = NULL;
	//Gromacs made this complicated. I'm doing to try and do this the simple way, checking one at a time to see if the interaction types aren't empty.
	for (i = 0; i < tpr->nmolblock; i++) {
		mtype = tpr->molbtype[i];
		//printf("%d atoms in block %d\n", tpr->atomsinmol[i], i);
		//printf("Type %d\n", tpr->molbtype[i]);
		for (j = 0; j < tpr->molbnmol[i]; j++) {
			for (k = 0; k < bondinteractions; k++) {
				itraction = bondinteraction[k];
				if (k == 3) {
					//printf("%d elements in settle\n", tpr->nr[itraction][mtype]);
					*nbonds += 2 * tpr->nr[itraction][mtype] / 4;
				} else {
					*nbonds += tpr->nr[itraction][mtype] / 3;
				}
				//printf("nbonds = %d\n", *nbonds);
			}
		}
	}
	int *fromlist = (int *) malloc(*nbonds * sizeof(int));
	int *tolist = (int *) malloc(*nbonds * sizeof(int));
	int *bolist = (int *) malloc(*nbonds * sizeof(int));
	//printf("I Malloced\n");
	for (i = 0; i < tpr->nmolblock; i++) {
		mtype = tpr->molbtype[i];
		for (j = 0; j < tpr->molbnmol[i]; j++) {
			for (k = 0; k < bondinteractions; k++) {
				itraction = bondinteraction[k];
				if (k==3) {
					for (l = 0; l < 2 * tpr->nr[itraction][mtype] / 4; l++) {
						fromlist[boffset + l] = 1 + tpr->interactionlist[itraction][mtype][4*(l/2)+1] + aoffset;
						tolist[boffset + l] = 1 + tpr->interactionlist[itraction][mtype][4*(l/2)+2+l] + aoffset;
						bolist[boffset + l] = k;
					}
					boffset += 2 * tpr->nr[itraction][mtype] / 4;
					/*for (l = 0; l < tpr->nr[itraction][mtype]; l++) {
						printf("%d ", tpr->interactionlist[itraction][mtype][l]);
					}
					printf("\n");*/
				}
				else {
					for (l = 0; l < tpr->nr[itraction][mtype] / 3; l++) {
						//The 1+ comes because the from and to lists are 1-indexed, since psfs are 1 indexed.
						fromlist[boffset + l] = 1 + tpr->interactionlist[itraction][mtype][3*l+1] + aoffset;
						tolist[boffset + l] = 1 + tpr->interactionlist[itraction][mtype][3*l+2] + aoffset;
						bolist[boffset + l] = k;
					}
					boffset += tpr->nr[itraction][mtype] / 3;
				}
			}
			aoffset += tpr->atomsinmol[mtype];
		}
	}
	/*printf("Bondlist:\n");
	for (i = 0; i < *nbonds; i++) {
		printf("%d to %d\n", fromlist[i], tolist[i]);
	}*/
	*fromptr = fromlist;
	*toptr = tolist;
	return MOLFILE_SUCCESS;
}

static int read_tpr_angles(void *v, int *numangles, int **angles, 
						 int **angletypes, int *numangletypes, 
						 char ***angletypenames, int *numdihedrals,
						 int **dihedrals, int **dihedraltypes, 
						 int *numdihedraltypes, char ***dihedraltypenames,
						 int *numimpropers, int **impropers, 
						 int **impropertypes, int *numimpropertypes, 
						 char ***impropertypenames, int *numcterms, 
						 int **cterms, int *ctermcols, int *ctermrows) {
	tprdata *tpr = (tprdata *)v;
	int i, j, k, l, itraction, mtype;
	int aoffset = 0;
	int boffset = 0;
	/* initialize data to zero */
	*numangles         = 0;
	*angles            = NULL;
	*angletypes        = NULL;
	*numangletypes     = 0;
	*angletypenames    = NULL;
	*numdihedrals      = 0;
	*dihedrals         = NULL;
	*dihedraltypes     = NULL;
	*numdihedraltypes  = 0;
	*dihedraltypenames = NULL;
	*numimpropers      = 0;
	*impropers         = NULL;
	*impropertypes     = NULL;
	*numimpropertypes  = 0;
	*impropertypenames = NULL;
	*numcterms         = 0;
	*cterms            = NULL;
	*ctermrows         = 0;
	*ctermcols         = 0;
	const int angleinteractions = 6;
	int angleinteraction[angleinteractions] = {F_ANGLES, F_G96ANGLES, F_RESTRANGLES, F_LINEAR_ANGLES, F_UREY_BRADLEY, F_SETTLE};
	for (i = 0; i < tpr->nmolblock; i++) {
		mtype = tpr->molbtype[i];
		for (j = 0; j < tpr->molbnmol[i]; j++) {
			for (k = 0; k < angleinteractions; k++) {
				itraction = angleinteraction[k];
				*numangles += tpr->nr[itraction][mtype] / 4;
			}
		}
	}
	//printf("Numangles read in: %d\n", *numangles);
	int *anglelist = (int *) malloc(3 * (*numangles) * sizeof(int));
	//printf("I Malloced\n");
	for (i = 0; i < tpr->nmolblock; i++) {
		mtype = tpr->molbtype[i];
		//printf("%d atoms in block %d\n", tpr->atomsinmol[i], i);
		for (j = 0; j < tpr->molbnmol[i]; j++) {
			for (k = 0; k < angleinteractions; k++) {
				itraction = angleinteraction[k];
				for (l = 0; l < tpr->nr[itraction][mtype] / 4; l++) {
					//The 1+ comes because the angle lists are 1-indexed, since psfs are 1 indexed.
					if (k == 5) {
						anglelist[boffset + 3*l + 0] = 1 + tpr->interactionlist[itraction][mtype][4*l+2] + aoffset;
						anglelist[boffset + 3*l + 1] = 1 + tpr->interactionlist[itraction][mtype][4*l+1] + aoffset;
						anglelist[boffset + 3*l + 2] = 1 + tpr->interactionlist[itraction][mtype][4*l+3] + aoffset;
					}
					else {
						anglelist[boffset + 3*l + 0] = 1 + tpr->interactionlist[itraction][mtype][4*l+1] + aoffset;
						anglelist[boffset + 3*l + 1] = 1 + tpr->interactionlist[itraction][mtype][4*l+2] + aoffset;
						anglelist[boffset + 3*l + 2] = 1 + tpr->interactionlist[itraction][mtype][4*l+3] + aoffset;
					}
				}
				boffset += 3*(tpr->nr[itraction][mtype] / 4);
			}
			aoffset += tpr->atomsinmol[mtype];
		}
	}
	*angles = anglelist;

	aoffset = 0;
	boffset = 0;
	const int pdihsinteractions = 1;
	int pdihsinteraction[pdihsinteractions] = {F_PDIHS};
	for (i = 0; i < tpr->nmolblock; i++) {
		mtype = tpr->molbtype[i];
		for (j = 0; j < tpr->molbnmol[i]; j++) {
			for (k = 0; k < pdihsinteractions; k++) {
				itraction = pdihsinteraction[k];
				*numdihedrals += tpr->nr[itraction][mtype] / 5;
			}
		}
	}
	printf("Numdihedrals read in: %d\n", *numdihedrals);
	int *pdihlist = (int *) malloc(4 * (*numdihedrals) * sizeof(int));
	//printf("I Malloced\n");
	for (i = 0; i < tpr->nmolblock; i++) {
		mtype = tpr->molbtype[i];
		//printf("%d atoms in block %d\n", tpr->atomsinmol[i], i);
		for (j = 0; j < tpr->molbnmol[i]; j++) {
			for (k = 0; k < pdihsinteractions; k++) {
				itraction = pdihsinteraction[k];
				for (l = 0; l < tpr->nr[itraction][mtype] / 5; l++) {
					//The 1+ comes because the angle lists are 1-indexed, since psfs are 1 indexed.
					pdihlist[boffset + 4*l + 0] = 1 + tpr->interactionlist[itraction][mtype][5*l+1] + aoffset;
					pdihlist[boffset + 4*l + 1] = 1 + tpr->interactionlist[itraction][mtype][5*l+2] + aoffset;
					pdihlist[boffset + 4*l + 2] = 1 + tpr->interactionlist[itraction][mtype][5*l+3] + aoffset;
					pdihlist[boffset + 4*l + 3] = 1 + tpr->interactionlist[itraction][mtype][5*l+4] + aoffset;
				}
				boffset += 4*(tpr->nr[itraction][mtype] / 5);
			}
			aoffset += tpr->atomsinmol[mtype];
		}
	}
	*dihedrals = pdihlist;

	aoffset = 0;
	boffset = 0;
	const int idihsinteractions = 1;
	int idihsinteraction[idihsinteractions] = {F_IDIHS};
	for (i = 0; i < tpr->nmolblock; i++) {
		mtype = tpr->molbtype[i];
		for (j = 0; j < tpr->molbnmol[i]; j++) {
			for (k = 0; k < idihsinteractions; k++) {
				itraction = idihsinteraction[k];
				*numimpropers += tpr->nr[itraction][mtype] / 5;
			}
		}
	}
	printf("Numimpropers read in: %d\n", *numimpropers);
	int *idihlist = (int *) malloc(4 * (*numimpropers) * sizeof(int));
	//printf("I Malloced\n");
	for (i = 0; i < tpr->nmolblock; i++) {
		mtype = tpr->molbtype[i];
		//printf("%d atoms in block %d\n", tpr->atomsinmol[i], i);
		for (j = 0; j < tpr->molbnmol[i]; j++) {
			for (k = 0; k < idihsinteractions; k++) {
				itraction = idihsinteraction[k];
				for (l = 0; l < tpr->nr[itraction][mtype] / 5; l++) {
					//The 1+ comes because the angle lists are 1-indexed, since psfs are 1 indexed.
					idihlist[boffset + 4*l + 0] = 1 + tpr->interactionlist[itraction][mtype][5*l+1] + aoffset;
					idihlist[boffset + 4*l + 1] = 1 + tpr->interactionlist[itraction][mtype][5*l+2] + aoffset;
					idihlist[boffset + 4*l + 2] = 1 + tpr->interactionlist[itraction][mtype][5*l+3] + aoffset;
					idihlist[boffset + 4*l + 3] = 1 + tpr->interactionlist[itraction][mtype][5*l+4] + aoffset;
				}
				boffset += 4*(tpr->nr[itraction][mtype] / 5);
			}
			aoffset += tpr->atomsinmol[mtype];
		}
	}
	*impropers = idihlist;

	aoffset = 0;
	boffset = 0;
	const int cmapinteraction = 1;
	int cmapinteractions[cmapinteraction] = {F_CMAP};
	for (i = 0; i < tpr->nmolblock; i++) {
		mtype = tpr->molbtype[i];
		for (j = 0; j < tpr->molbnmol[i]; j++) {
			for (k = 0; k < cmapinteraction; k++) {
				itraction = cmapinteractions[k];
				*numcterms += tpr->nr[itraction][mtype] / 6;
			}
		}
	}
	//printf("Numcmap terms read in: %d\n", *numcterms);
	int *cmaplist = (int *) malloc(8 * (*numcterms) * sizeof(int));
	//printf("I Malloced\n");
	for (i = 0; i < tpr->nmolblock; i++) {
		mtype = tpr->molbtype[i];
		//printf("%d atoms in block %d\n", tpr->atomsinmol[i], i);
		for (j = 0; j < tpr->molbnmol[i]; j++) {
			for (k = 0; k < cmapinteraction; k++) {
				itraction = cmapinteractions[k];
				for (l = 0; l < tpr->nr[itraction][mtype] / 6; l++) {
					//The 1+ comes because the angle lists are 1-indexed, since psfs are 1 indexed.
					cmaplist[boffset + 5*l + 0] = 1 + tpr->interactionlist[itraction][mtype][6*l+1] + aoffset;
					cmaplist[boffset + 5*l + 1] = 1 + tpr->interactionlist[itraction][mtype][6*l+2] + aoffset;
					cmaplist[boffset + 5*l + 2] = 1 + tpr->interactionlist[itraction][mtype][6*l+3] + aoffset;
					cmaplist[boffset + 5*l + 3] = 1 + tpr->interactionlist[itraction][mtype][6*l+4] + aoffset;
					cmaplist[boffset + 5*l + 4] = 1 + tpr->interactionlist[itraction][mtype][6*l+2] + aoffset;
					cmaplist[boffset + 5*l + 5] = 1 + tpr->interactionlist[itraction][mtype][6*l+3] + aoffset;
					cmaplist[boffset + 5*l + 6] = 1 + tpr->interactionlist[itraction][mtype][6*l+4] + aoffset;
					cmaplist[boffset + 5*l + 7] = 1 + tpr->interactionlist[itraction][mtype][6*l+5] + aoffset;
				}
				boffset += 5*(tpr->nr[itraction][mtype] / 6);
			}
			aoffset += tpr->atomsinmol[mtype];
		}
	}
	*cterms = cmaplist;

	*ctermcols = 0;
	*ctermrows = 0;
	return MOLFILE_SUCCESS;
}

static void close_tpr_read(void *mydata) {
	int i, j;
	//printf("Closing TPR file\n");
	tprdata *tpr = (tprdata *)mydata;
	fclose(tpr->f);
	for (i = 0; i < tpr->nmoltypes; i++) {
		//printf("%d of %d (%d atoms)\n", i+1, tpr->nmoltypes, tpr->atomsinmol[i]);
		free(tpr->charges[i]);
		free(tpr->masses[i]);
		free(tpr->types[i]);
		free(tpr->ptypes[i]);
		free(tpr->resids[i]);
		free(tpr->atomnameids[i]);
		free(tpr->atomtypeids[i]);
		free(tpr->resnames[i]);
		free(tpr->atomicnumbers[i]);
	}
	for (i = 0; i < F_NRE; i++) {
		if (tpr->interactionlist[i] != NULL) {
			for (j = 0; j < tpr->nmoltypes; j++) {
				if (tpr->interactionlist[i][j] != NULL)
					free(tpr->interactionlist[i][j]);
			}
			free(tpr->interactionlist[i]);
		}
		if (tpr->nr[i] != NULL) {
			free(tpr->nr[i]);
		}
	}
	free(tpr->atomicnumbers);
	free(tpr->molnames);
	free(tpr->molbtype);
	free(tpr->molbnmol);
	free(tpr->molbnatoms);
	free(tpr->resnames);
	free(tpr->atomtypeids);
	free(tpr->atomnameids);
	free(tpr->types);
	free(tpr->ptypes);
	free(tpr->resids);
	free(tpr->charges);
	free(tpr->masses);
	free(tpr->symtab);
	free(tpr->atomsinmol);
	free(tpr->resinmol);
	xdr_destroy(tpr->xdrptr);
	delete tpr->xdrptr;
	delete tpr;
	//printf("TPR file read completely\n");
}


/*
 * Initialization stuff down here
 */

static molfile_plugin_t plugin;

VMDPLUGIN_API int VMDPLUGIN_init(void) {
	memset(&plugin, 0, sizeof(molfile_plugin_t));
	plugin.abiversion = vmdplugin_ABIVERSION;
	plugin.type = MOLFILE_PLUGIN_TYPE;
	plugin.name = "tpr";
	plugin.prettyname = "Gromacs Binary Topology";
	plugin.author = "Josh Vermaas";
	plugin.majorv = 5;
	plugin.minorv = 1;//Corresponds to the Gromacs version I was basing this on.
	plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
	plugin.filename_extension = "tpr";
	plugin.open_file_read = open_tpr_read;
	plugin.read_structure = read_tpr_structure;
	plugin.read_next_timestep = read_tpr_timestep;
	plugin.read_bonds = read_tpr_bonds;
	plugin.read_angles = read_tpr_angles;
	plugin.close_file_read = close_tpr_read;
	return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
	(*cb)(v, (vmdplugin_t *)&plugin);
	return 0;
}

VMDPLUGIN_API int VMDPLUGIN_fini(void) { return VMDPLUGIN_SUCCESS; }

/*

int main (int argc, char *argv[]) {
	FILE *fin;
	XDR xdrs;
	int i, j, tmp;
	int precision;
	//fin = fopen("topol.tpr", "r");
	//fin = fopen("topol100.tpr", "r");
	fin = fopen("3WVF/topol.tpr", "rb");
	xdrstdio_create(&xdrs, fin, XDR_DECODE);
	xdr_int(&xdrs, &i);
	printf("%d\n", i);
	printString(&xdrs);
	precision = readInt(&xdrs);
	printf("%d\n", precision);
	if (precision == 4) {
		tprdata *tprdat = new tprdata;
		readtprAfterPrecision(&xdrs, tprdat);
		//return 0;
	}
	/*else if (precision == 8) {
		tprdata<double> *tprdat = new tprdata;
		readtprAfterPrecision<double>(&xdrs, tprdat);
		return 0;
	}*//*
	else {
		printf("Illegal precision (neither single nor double)\n");
		return 1;
	}
	fclose(fin);
	return 0;
}*/