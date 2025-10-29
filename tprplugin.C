#define VMDPLUGIN_STATIC
#include "molfile_plugin.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


/* I am making a very concious decision here. I'm only going to TRY supporting gromacs TPX files that were written by
version 4.0 or later. This means TPX format version 58 and later.
*/

#include "Gromacs.h"

#define STRLEN 4096
#define SAVELEN 8
#define TRUE (1)
#define FALSE (0)

#include "tprplugin.h"

int readff(md_file *mf, int version) {
    int atnr, ntypes, i, j;
    int *functype;
    double reppow =12.0;
    float fudge;
    if (trx_int(mf, &atnr)) return MOLFILE_ERROR;
    #ifdef TPRDEBUG
    printf("Atnr: %d\n", atnr);
    #endif
    if (trx_int(mf, &ntypes)) return MOLFILE_ERROR;
    #ifdef TPRDEBUG
    printf("Ntypes: %d\n", ntypes);
    #endif
    functype = (int *) malloc(sizeof(int) * ntypes);
    for (i=0; i < ntypes; i++)
    	if (trx_int(mf, &functype[i])) return MOLFILE_ERROR;
    if (version >= 66)
    	if (trx_double(mf, &reppow)) return MOLFILE_ERROR;
   	if (trx_real(mf, &fudge)) return MOLFILE_ERROR;
    for (i = 0; i < ntypes; i++) {
        for (j = 0; (j < NFTUPD); j++) {
            if (version < ftupd[j].fvnr && functype[i] >= ftupd[j].ftype)
                functype[i] += 1;
        }
        //Equivalent to do_iparams
        readparams(mf, version, functype[i]);
    }
    free(functype);
    return MOLFILE_SUCCESS;
}
//do_atomtypes
int read_atomtypes(md_file *mf, int version) {
    int i;
    int err = 0;
    int numtypes;
    err |= trx_int(mf, &numtypes);
#ifdef TPRDEBUG
    printf("Num atom types:%d\n", numtypes);
#endif
    //Read a bunch of stuff that is set optionally.
    for (i = 0; i < numtypes; i++) {
        if (version < 113) {
            err |= trx_real(mf, NULL);//Radius?
            err |= trx_real(mf, NULL);//Volume?
            err |= trx_real(mf, NULL);//Surface tension?
        }
        if (version >= 40) {
            err |= trx_real(mf, NULL);
        }
        if (version >= 60 && version < 113) {
            err |= trx_real(mf, NULL);
            err |= trx_real(mf, NULL);
        }
    }
    return err;
}
//Reading in this case means discarding.
int read_cmap(md_file *mf) {
    int i, ngrid, gs;
    int err = 0;
    err |= trx_int(mf, &ngrid);
    err |= trx_int(mf, &gs);
    #ifdef TPRDEBUG
    printf("ngrid: %d, gs: %d\n", ngrid, gs);
    printf("This many reals should be read: %d\n", ngrid * gs * gs);
    #endif
    for (i = 0; i < ngrid * gs * gs; i++) {
    	err |= trx_real(mf, NULL);
    	err |= trx_real(mf, NULL);
    	err |= trx_real(mf, NULL);
    	err |= trx_real(mf, NULL);
    }
    return err;
}
//do_groups
int read_groups(md_file *mf, int ngrp, tprdata *tpr) {
	int g, i, j, tmp;
	int err = 0;
	//do_grps
	#ifdef TPRDEBUG
	printf("ngrp = %d\n", ngrp);
	#endif
	for(i = 0; i < ngrp; i++) {
		err |= trx_int(mf, &tmp);
		#ifdef TPRDEBUG
		printf("Number of reads = %d\n", tmp);
		#endif
		for (j = 0; j < tmp; j++) {
			trx_int(mf, NULL);
		}
	}
	//end do_grps
	//ngrpname
	err |= trx_int(mf, &tmp);
	#ifdef TPRDEBUG
	printf("ngrpname = %d\n", tmp);
	#endif
	for (i = 0; i < tmp; i++) {
		err |= trx_int(mf, &g);
		#ifdef TPRDEBUG
		printf("symtabid = %d\n", g);
		for (j = 0; j < SAVELEN; j++)
			printf("%c", tpr->symtab[SAVELEN*g+j]);
		printf("\n");
		#endif
	}
	for (g = 0; g < ngrp; g++) {
		err |= trx_int(mf, &tmp);
		#ifdef TPRDEBUG
		printf("This many chars: %d\n", tmp);
		#endif
		if (tmp) {
			//For reasons that aren't clear to me, writer version 27 has a different way of reading/writing the chars.
			if (tpr->wversion >=27) {
				// Advance the pointer by tmp bytes.
				fseek(mf->f, tmp, SEEK_CUR);
			}
			else {
				fseek(mf->f, 4*tmp, SEEK_CUR);
			}
		}
	}
	return err;
}

#define MIN(a, b) ((a)<(b)? (a):(b))

void tpr_save_string(md_file* mf, char* saveloc, int genversion) {
    #ifdef _WIN32
        long long int len = 0;
    #else
        long len = 0;
    #endif
    int i, j;
    fpos_t pos;
    char buf[STRLEN];
    trx_long(mf, &len);
    fread(buf, 1, int(len), mf->f);
    // GROMACS is weird. Before writer version 27, the reads were always aligned to 4 bytes.
    // In subsequent versions, they are not. So to maintain backwards compatability, add an
    // extra seek.
    if (genversion < 27 && len % 4) {
    	fseek(mf->f, 4 - (len % 4), SEEK_CUR);
    }

    for (i = 0; i < MIN(int(len), (SAVELEN-1)); i++) {
        saveloc[i] = buf[i];
    }
    saveloc[i] = '\0';

}

int readtprAfterPrecision (tprdata *tpr) {
	//Start at do_tpx, trace down to do_tpxheader
	md_file *mf = tpr->mf;
	char buf[STRLEN];
	tpr->readcoordinates = 0;
	long fsize;
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

	if (trx_int(mf, &(tpr->version))) return MOLFILE_ERROR;
	printf("File Format Version: %d\n", tpr->version);

	if (tpr->version >= 77 && tpr->version <= 79) {
		trx_string(mf, buf, STRLEN);
	}

	if (trx_int(mf, &(tpr->wversion))) return MOLFILE_ERROR;
	printf("Generator Version: %d\n", tpr->wversion);

	//Bailouts if things are too new/we can't guarantee accurately reading them.
	if (tpr->wversion > 28 || tpr->version <= 57) {
		printf("Your file cannot be read, as it has version %d, but we can read from version 57 to at least 128.\n", tpr->version);
		printf("The generator version for your file is %d, but we can only read up to 28\n", tpr->wversion);
		return MOLFILE_ERROR;
	}
	
	if (tpr->version >= 81) {
		j = tpr_string(mf, buf, STRLEN);
		fread(buf, 4, 1, mf->f);
		// //I dunno why this is wrong. Should say "release", but it totally doesn't.

	}

	//If we were dealing with TPA files, we'd need to do something here. Not a clue what.
	//do_section line in the gromacs source.
	if (trx_int(mf, &(tpr->natoms))) return MOLFILE_ERROR;
	#ifdef TPRDEBUG
	printf("Natoms: %d\n", tpr->natoms);
	#endif

	if (trx_int(mf, &(tpr->ngtc))) return MOLFILE_ERROR;
	#ifdef TPRDEBUG
	printf("Ngtc: %d\n", tpr->ngtc);
	#endif
	if (tpr->version < 62) {
		if (trx_int(mf, &tmp)) return MOLFILE_ERROR;
		if (trx_real(mf, &dummy)) return MOLFILE_ERROR;
	}

	if (tpr->version >=79) {
		if (trx_int(mf, &tmp)) return MOLFILE_ERROR;
		#ifdef TPRDEBUG
		printf("FEP state: %d\n", tmp);
		#endif
	}
	if (trx_real(mf, &dummy)) return MOLFILE_ERROR;
	#ifdef TPRDEBUG
	printf("Lambda: %f\n", dummy);
	#endif
	if (trx_int(mf, &hasIR)) return MOLFILE_ERROR;
	#ifdef TPRDEBUG
	printf("IR state?: %d\n", hasIR);
	#endif
	if (trx_int(mf, &hasTop)) return MOLFILE_ERROR;
	#ifdef TPRDEBUG
	printf("Contains topology: %d\n", hasTop);
	#endif
	if (trx_int(mf, &hasCoord)) return MOLFILE_ERROR;
	#ifdef TPRDEBUG
	printf("Contains coordinates: %d\n", hasCoord);
	#endif
	if (trx_int(mf, &hasV)) return MOLFILE_ERROR;
	#ifdef TPRDEBUG
	printf("Contains velocities: %d\n", hasV);
	#endif
	tpr->has_velocities = 0;
	if(hasV) tpr->has_velocities = 1;
	if (trx_int(mf, &hasF)) return MOLFILE_ERROR;
	#ifdef TPRDEBUG
	printf("Contains forces: %d\n", hasF);
	#endif
	if (trx_int(mf, &hasDim)) return MOLFILE_ERROR;
	#ifdef TPRDEBUG
	printf("Contains dimensions: %d\n", hasDim);
	#endif
	if (tpr->wversion >= 27) {
		if (trx_long(mf, &fsize)) return MOLFILE_ERROR;
		#ifdef TPRDEBUG
		printf("Filesize: %d\n", fsize);
		#endif
	}

	//End of do_tpxheader
	//Now start in do_tpx_body
	//Lead with do_tprx_state_first
	if (hasDim) {
		if (trx_rvector(mf, &(tpr->boxdims[0]))) return MOLFILE_ERROR;
		if (trx_rvector(mf, &(tpr->boxdims[3]))) return MOLFILE_ERROR;
		if (trx_rvector(mf, &(tpr->boxdims[6]))) return MOLFILE_ERROR;
		#ifdef TPRDEBUG
		printf("Box dim: (%f %f %f), (%f %f %f), (%f %f %f)\n", tpr->boxdims[0],tpr->boxdims[1],tpr->boxdims[2],
		 	tpr->boxdims[3],tpr->boxdims[4],tpr->boxdims[5],tpr->boxdims[6],tpr->boxdims[7],tpr->boxdims[8]);
		#endif
		if (tpr->version >= 51) {
			if (trx_rvector(mf, &(boxoffset[0]))) return MOLFILE_ERROR;
			if (trx_rvector(mf, &(boxoffset[3]))) return MOLFILE_ERROR;
			if (trx_rvector(mf, &(boxoffset[6]))) return MOLFILE_ERROR;
		}
		else {
			boxoffset[0] = boxoffset[1] = boxoffset[2] = 0;
		}
		#ifdef TPRDEBUG
		printf("Box offset: (%f %f %f), (%f %f %f), (%f %f %f)\n", boxoffset[0],boxoffset[1],boxoffset[2],
			boxoffset[3],boxoffset[4],boxoffset[5],boxoffset[6],boxoffset[7],boxoffset[8]);
		#endif
		if (trx_rvector(mf, &(boxv[0]))) return MOLFILE_ERROR;
		if (trx_rvector(mf, &(boxv[3]))) return MOLFILE_ERROR;
		if (trx_rvector(mf, &(boxv[6]))) return MOLFILE_ERROR;
		#ifdef TPRDEBUG
		printf("Box vel: (%f %f %f), (%f %f %f), (%f %f %f)\n", boxv[0],boxv[1],boxv[2],
			boxv[3],boxv[4],boxv[5],boxv[6],boxv[7],boxv[8]);
		#endif
	}
	//Dump the thermostat storage stuff. VMD wouldn't know what to do with it anyway.
	if (tpr->ngtc) {
		if (tpr->version < 69) {
			for (i = 0; i < tpr->ngtc; i++) {
				if (trx_real(mf, &dummy)) return MOLFILE_ERROR;
				//printf("NGTC %d: %f\n", i, dummy);
			}
		}
		for (i = 0; i < tpr->ngtc; i++) {
			if (trx_real(mf, &dummy)) return MOLFILE_ERROR;
			//printf("NGTC %d: %f\n", i, dummy);
		}
	}

	//do_tpx_mtop points to do_mtop, which starts here.
	//do_mtop starts by reading the symtab (do_symtab)
	if (trx_int(mf, &(tpr->symtablen))) return MOLFILE_ERROR;
	#ifdef TPRDEBUG
	printf("Symtab length: %d\n", tpr->symtablen);
	#endif
	tpr->symtab = (char*) malloc(SAVELEN * tpr->symtablen * sizeof(char));
	for (i = 0; i < tpr->symtablen; i++) {
		#ifdef TPREXTRADEBUG
		printf("i=%d\n", i);
		#endif
		tpr_save_string(mf, &(tpr->symtab[SAVELEN * i]),tpr->wversion);
	}
	#ifdef TPRDEBUG
	for (i = 0; i < tpr->symtablen; i++) {
		for (j = 0; j < SAVELEN; j++)
			printf("%c", tpr->symtab[SAVELEN*i+j]);
		printf("\n");
	}
	#endif

	if (trx_int(mf, &tmp)) return MOLFILE_ERROR;
	/*printf ("System name: ");
	for (j = 0; j < SAVELEN; j++)
		printf("%c", tpr->symtab[SAVELEN*tmp+j]);
	printf("\n");*/

	//Now "read" in the forcefield. This SHOULD be do_ffparams
	readff(mf, tpr->version);

	//Read in the type of molecules.
	if (trx_int(mf, &(tpr->nmoltypes))) return MOLFILE_ERROR;
	#ifdef TPRDEBUG
	printf("nmoltypes: %d\n", tpr->nmoltypes);
	#endif
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
		if (trx_int(mf, &(tpr->molnames[i]))) return MOLFILE_ERROR;
		#ifdef TPRDEBUG
		printf ("Mol name %d: ", i);
		for (j = 0; j < SAVELEN; j++)
			printf("%c", tpr->symtab[SAVELEN*tpr->molnames[i]+j]);
		printf("\n");
		#endif
		if (trx_int(mf, &(tpr->atomsinmol[i]))) return MOLFILE_ERROR;
		if (trx_int(mf, &(tpr->resinmol[i]))) return MOLFILE_ERROR;
		tpr->charges[i] = (float*) malloc(sizeof(float) * tpr->atomsinmol[i]);
		tpr->masses[i] = (float*) malloc(sizeof(float) * tpr->atomsinmol[i]);
		tpr->types[i] = (unsigned short*) malloc(sizeof(unsigned short) * tpr->atomsinmol[i]);
		tpr->ptypes[i] = (int*) malloc(sizeof(int) * tpr->atomsinmol[i]);
		tpr->resids[i] = (int*) malloc(sizeof(int) * tpr->atomsinmol[i]);
		tpr->atomicnumbers[i] = (int*) malloc(sizeof(int) * tpr->atomsinmol[i]);
		#ifdef TPRDEBUG
		printf("%d atoms, %d residues in molecule %d\n", tpr->atomsinmol[i], tpr->resinmol[i], i);
		#endif
		for (j=0; j < tpr->atomsinmol[i]; j++) {
			if (trx_real(mf, &(tpr->masses[i][j]))) return MOLFILE_ERROR;
			if (trx_real(mf, &(tpr->charges[i][j]))) return MOLFILE_ERROR;
			if (trx_real(mf, NULL)) return MOLFILE_ERROR; //mB
			if (trx_real(mf, NULL)) return MOLFILE_ERROR; //cB

			if (trx_int(mf, &tmp)) return MOLFILE_ERROR;
			tpr->types[i][j] = (unsigned short) tmp;
			if (trx_int(mf, &tmp)) return MOLFILE_ERROR;//typeB
			if (tpr->wversion >=27) {
				fseek(mf->f, -4, SEEK_CUR);
			}
			if (trx_int(mf, &(tpr->ptypes[i][j]))) return MOLFILE_ERROR;
			if (trx_int(mf, &(tpr->resids[i][j]))) return MOLFILE_ERROR;
			if (tpr->version >= 52) {
				if (trx_int(mf, &(tpr->atomicnumbers[i][j]))) return MOLFILE_ERROR;
			}
			#ifdef TPREXTRADEBUG
			if (i == 0)
				printf("%d: mass %f charge %f type %d ptype %d resid %d, periodic table number %d\n", j,
					tpr->masses[i][j], tpr->charges[i][j], tpr->types[i][j], tpr->ptypes[i][j],tpr->resids[i][j], tmp);
			#endif
		}
		tpr->atomnameids[i] = (int*) malloc(sizeof(int) * tpr->atomsinmol[i]);
		if (tpr_ivector(mf, tpr->atomnameids[i], tpr->atomsinmol[i])) return MOLFILE_ERROR;
		tpr->atomtypeids[i] = (int*) malloc(sizeof(int) * tpr->atomsinmol[i]);
		if (tpr_ivector(mf, tpr->atomtypeids[i], tpr->atomsinmol[i])) return MOLFILE_ERROR;
		for (j = 0; j < tpr->atomsinmol[i]; j++) {
			if (trx_real(mf, NULL)) return MOLFILE_ERROR;
		}

		//Read residues
		tpr->resnames[i] = (int*) malloc(sizeof(int) * tpr->resinmol[i]);
		#ifdef TPRDEBUG
		printf("%d residues in mol %d\n", tpr->resinmol[i], i);
		#endif
		for (j = 0; j < tpr->resinmol[i]; j++) {
			if (trx_int(mf, &(tpr->resnames[i][j]))) return MOLFILE_ERROR;
			if (tpr->version >=63) {
				if (tpr->wversion < 27)
					fseek(mf->f, 8, SEEK_CUR);
				else
					fseek(mf->f, 5, SEEK_CUR);
			}
			else
				tpr->resnames[i][j] += 1;
			#ifdef TPREXTRADEBUG
			printf("%d\n", tpr->resnames[i][j]);
			if ( i == 0) {
				for (k = 0; k < SAVELEN; k++)
					printf("%c", tpr->symtab[SAVELEN*tpr->resnames[i][j]+k]);
				printf("\n");
			}
			#endif
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
				if (trx_int(mf, &(tpr->nr[j][i]))) return MOLFILE_ERROR;
				#ifdef TPREXTRADEBUG
				printf("j, k, interactions: %d, %d, %d\n", j, k, tpr->nr[j][i]);
				#endif
			}
			if (tpr->nr[j][i] != 0) {
				tpr->interactionlist[j][i] = (int*) malloc(sizeof(int) * tpr->nr[j][i]);
				#ifdef TPRDEBUG
				printf("Interaction id: %d number of interactants %d\n", j, tpr->nr[j][i]);
				#endif
				tpr_ivector(mf, tpr->interactionlist[j][i], tpr->nr[j][i]);
				#ifdef TPREXTRADEBUG
				for (k=0; k < tpr->nr[j][i]; k++) {
					printf("%d ", tpr->interactionlist[j][i][k]);
				}
				printf("\n");
				#endif
			}
			else {
				tpr->interactionlist[j][i] = NULL;
			}
		}
		//end do_ilists


		//do_block and do_blocka. Remove these. VMD doesn't know what they are.
		//They refer to atomgroups and exclusions. Useful for dynamics, not useful for visualization.
		if (trx_int(mf, &k)) return MOLFILE_ERROR;
		for (j = 0; j <= k; j++) {
			if (trx_int(mf, NULL)) return MOLFILE_ERROR;
		}
		if (trx_int(mf, &k)) return MOLFILE_ERROR;
		if (trx_int(mf, &tmp)) return MOLFILE_ERROR;
		for (j = 0; j <= k; j++) {
			if (trx_int(mf, NULL)) return MOLFILE_ERROR;
		}
		for (j = 0; j < tmp; j++) {
			if (trx_int(mf, NULL)) return MOLFILE_ERROR;
		}
	}
	//end of do_moltype

	if (trx_int(mf, &(tpr->nmolblock))) return MOLFILE_ERROR;
	#ifdef TPRDEBUG
	printf("nmolblock: %d\n", tpr->nmolblock);
	#endif
	//do_molblock
	tpr->molbtype = (int*) malloc(tpr->nmolblock * sizeof(int));
	tpr->molbnmol = (int*) malloc(tpr->nmolblock * sizeof(int));
	tpr->molbnatoms = (int*) malloc(tpr->nmolblock * sizeof(int));
	for (i = 0; i < tpr->nmolblock; i++) {
		if (trx_int(mf, &(tpr->molbtype[i]))) return MOLFILE_ERROR;
		if (trx_int(mf, &(tpr->molbnmol[i]))) return MOLFILE_ERROR;
		if (trx_int(mf, &(tpr->molbnatoms[i]))) return MOLFILE_ERROR;
		if (trx_int(mf, &k)) return MOLFILE_ERROR;
		#ifdef TPRDEBUG
		printf("posresXA: %d\n", k);
		#endif
		//Position restraints have a float for every coordinate, hence the multiplier by 3.
		for (j = 0; j < 3 * k; j++) {
			if (trx_real(mf, NULL)) return MOLFILE_ERROR;
		}
		if (trx_int(mf, &k)) return MOLFILE_ERROR;
		#ifdef TPRDEBUG
		printf("posresXB: %d\n", k);
		#endif
		for (j = 0; j < 3 * k; j++) {
			if (trx_real(mf, NULL)) return MOLFILE_ERROR;
		}
		#ifdef TPRDEBUG
		printf("Segname: %d ", tpr->molbtype[i]);
		for (j = 0; j < SAVELEN; j++) {
			printf("%c", tpr->symtab[SAVELEN*tpr->molnames[tpr->molbtype[i]]+j]);
		}
		printf("\n\tNumatoms: %d\n", tpr->molbnatoms[i]);
		printf("\tNumcopies: %d\n", tpr->molbnmol[i]);
		#endif
		//sum += tpr->molbnatoms[i] * tpr->molbnmol[i];
	}
	//Burn off the number of atoms after the do_molblock
	if (trx_int(mf, &tmp)) return MOLFILE_ERROR;
#ifdef TPRDEBUG
    printf("What is this? %d. It should be the number of atoms.\n", tmp);
#endif
    if (tpr->version >= 103) { //103 is tpxv_IntermolecularBondeds
    	if (trx_int(mf, &hasIntermoleculeBonds)) return MOLFILE_ERROR;
#ifdef TPRDEBUG
        printf("intermolecularbondeds %d\n", hasIntermoleculeBonds);
#endif
        if (hasIntermoleculeBonds) {
            printf("Systems with intermolecular bonds are not supported. The file reports %d intermolecular bonds.\n", hasIntermoleculeBonds);
            return MOLFILE_ERROR;
        }
        if (tpr->wversion >= 27) {
        	fseek(mf->f, -3, SEEK_CUR);
        }
    }

    if (tpr->version < 128) //128 is tpxv_RemoveAtomtypes
    {
	    //do_atomtypes
	    if (read_atomtypes(mf, tpr->version)) return MOLFILE_ERROR;
	}

#ifdef TPRDEBUG
    printf("Reading cmap terms\n");
#endif
    //do_cmap
    if (tpr->version >= 65) {
        if (read_cmap(mf)) return MOLFILE_ERROR;
    }
    //do_groups
#ifdef TPRDEBUG
    printf("Reading groups\n");
#endif
    read_groups(mf, egcNR, tpr);
    if (tpr->version >= 120) {
    	long len;
    	if (trx_long(mf, &len)) return MOLFILE_ERROR;
        int* jj = new int[len];
        #ifdef TPRDEBUG
        printf("Intermolecular Exclusions: %d\n", len);
        printf("%d\n", ftell(mf->f));
        #endif
        tpr_ivector(mf, jj, len);
    }
#ifdef TPRDEBUG
    printf("This is my file position: %d\n", ftell(mf->f));
    printf("Returning control\n");
#endif
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
    char buf[STRLEN];
    md_file *mf = new md_file;

    int i, j;
    if (!(fin = fopen(filename, "rb"))) {
        fprintf(stderr, "tprplugin) Cannot open tpr file '%s'\n", filename);
        return NULL;
    }
    mf->f = fin;
    if (trx_int(mf, &i)) {
    	fprintf(stderr, "tprplugin) Could not read initial integer from file.\n");
        return NULL;
    }
	if (i > STRLEN) {//If i value is large, everything in the file should be endian swapped.
		mf->rev = 1;
	}
	j = tpr_string(mf, buf, STRLEN);
	if (trx_int(mf, &(mf->prec))) {
		fprintf(stderr, "tprplugin) Could not read precision from file.\n");
        return NULL;
	}
    if (mf->prec == 4) {
        tprdat = new tprdata;
        memset(tprdat, 0, sizeof(tprdata));
        tprdat->mf = mf;
        if (readtprAfterPrecision(tprdat) != MOLFILE_SUCCESS) {
            delete tprdat;
            return NULL;
        }
        *natoms = tprdat->natoms;
        return tprdat;
    }
    else {
    	fprintf(stderr, "tprplugin) Illegal precision (requires single)\n");
        return NULL;
    }
    fclose(fin);
    return NULL;
}

static int read_tpr_timestep(void *v, int natoms, molfile_timestep_t *ts) {
    tprdata *tpr = (tprdata *)v;
    md_file *mf = tpr->mf;
    if (tpr->readcoordinates) {
        return MOLFILE_ERROR;
    }
    //Get the positions.
    if (ts != NULL) {
        tpr_rvector(mf, ts->coords, 3 * tpr->natoms);
        for (int i = 0; i < 3 * natoms; i++) {
            ts->coords[i] = 10 * ts->coords[i]; //A
        }
        #ifdef TPRDEBUG
        printf("\nAtoms Coordinates: (A)\n");
        for (int i = 0; i < natoms; i++) {
            printf("x[%d]: %f %f %f\n", i,   ts->coords[3*i+0], ts->coords[3*i+1],ts->coords[3*i+2]);
        }
        printf("coordinate end position: %d\n", ftell(mf->f));

        for (int i = 0; i < 3*natoms; i++) {
            if (! isfinite(ts->coords[i])) {
                printf("The %d coordinate of atom %d is not finite!\n", (i % 3), i / 3);
            }
        }
        #endif

        if(tpr->has_velocities)
        {
            tpr_rvector(mf, ts->velocities, 3 * tpr->natoms);
            for (int i = 0; i < 3 * natoms; i++)
            {
                ts->velocities[i] *= 10; //A/ps
                //fprintf(stderr, "%f\n", ts->velocities[i]);
            }
        }

        #ifdef TPRDEBUG
        if(tpr->has_velocities)
        {
            printf("\nAtoms Velocities: (A/ps)\n");
            for (int i = 0; i < natoms; i++)
            {
                printf("v[%d]: %f %f %f\n", i,   ts->velocities[3*i+0], ts->velocities[3*i+1],ts->velocities[3*i+2]);
            }
            printf("velocity end position: %d\n", ftell(mf->f));
        }
        #endif

        ts->A = sqrt(tpr->boxdims[0]*tpr->boxdims[0] + tpr->boxdims[1]*tpr->boxdims[1] + tpr->boxdims[2] * tpr->boxdims[2]) * 10;
        ts->B = sqrt(tpr->boxdims[3]*tpr->boxdims[3] + tpr->boxdims[4]*tpr->boxdims[4] + tpr->boxdims[5] * tpr->boxdims[5]) * 10;
        ts->C = sqrt(tpr->boxdims[6]*tpr->boxdims[6] + tpr->boxdims[7]*tpr->boxdims[7] + tpr->boxdims[8] * tpr->boxdims[8]) * 10;

        if(ts->A <= 0 || ts->B <= 0 || ts->C <=0)
        {
            ts->A     = ts->B    = ts->C     = 0;
            ts->alpha = ts->beta = ts->gamma = 0;
        }
        else
        {
        	ts->alpha = acos((tpr->boxdims[3] * tpr->boxdims[6] + tpr->boxdims[4] * tpr->boxdims[7] + tpr->boxdims[5] * tpr->boxdims[8])*100/(ts->A*ts->B)) * 90.0 / M_PI_2;
        	ts->beta  = acos((tpr->boxdims[0] * tpr->boxdims[6] + tpr->boxdims[1] * tpr->boxdims[7] + tpr->boxdims[2] * tpr->boxdims[8])*100/(ts->A*ts->C)) * 90.0 / M_PI_2;
        	ts->gamma = acos((tpr->boxdims[0] * tpr->boxdims[3] + tpr->boxdims[1] * tpr->boxdims[4] + tpr->boxdims[2] * tpr->boxdims[5])*100/(ts->B*ts->C)) * 90.0 / M_PI_2;
        }
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
    //printf("Numdihedrals read in: %d\n", *numdihedrals);
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
    //printf("Numimpropers read in: %d\n", *numimpropers);
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
                    cmaplist[boffset + 8*l + 0] = 1 + tpr->interactionlist[itraction][mtype][6*l+1] + aoffset;
                    cmaplist[boffset + 8*l + 1] = 1 + tpr->interactionlist[itraction][mtype][6*l+2] + aoffset;
                    cmaplist[boffset + 8*l + 2] = 1 + tpr->interactionlist[itraction][mtype][6*l+3] + aoffset;
                    cmaplist[boffset + 8*l + 3] = 1 + tpr->interactionlist[itraction][mtype][6*l+4] + aoffset;
                    cmaplist[boffset + 8*l + 4] = 1 + tpr->interactionlist[itraction][mtype][6*l+2] + aoffset;
                    cmaplist[boffset + 8*l + 5] = 1 + tpr->interactionlist[itraction][mtype][6*l+3] + aoffset;
                    cmaplist[boffset + 8*l + 6] = 1 + tpr->interactionlist[itraction][mtype][6*l+4] + aoffset;
                    cmaplist[boffset + 8*l + 7] = 1 + tpr->interactionlist[itraction][mtype][6*l+5] + aoffset;
                }
                boffset += 8*(tpr->nr[itraction][mtype] / 6);
            }
            aoffset += tpr->atomsinmol[mtype];
        }
    }
    *cterms = cmaplist;

    *ctermcols = 0;
    *ctermrows = 0;
    return MOLFILE_SUCCESS;
}

static int read_tpr_timestep_metadata(void *v, molfile_timestep_metadata_t *metadata)
{
    tprdata *tpr = (tprdata *)v;
    if(tpr->has_velocities == 1)
    {
        metadata->has_velocities = 1;
    }
    else
    {
        metadata->has_velocities = 0;
    }
    return MOLFILE_SUCCESS;
}

static void close_tpr_read(void *mydata) {
	int i, j;
	//printf("Closing TPR file\n");
	tprdata *tpr = (tprdata *)mydata;
	fclose(tpr->mf->f);
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
	delete tpr->mf;
	delete tpr;
	//printf("TPR file read completely\n");
}


/*
 * Initialization stuff down here
 */

#ifndef TPRTEST

static molfile_plugin_t tpr_plugin;

VMDPLUGIN_API int VMDPLUGIN_init(void) {
    memset(&tpr_plugin, 0, sizeof(molfile_plugin_t));
    tpr_plugin.abiversion = vmdplugin_ABIVERSION;
    tpr_plugin.type = MOLFILE_PLUGIN_TYPE;
    tpr_plugin.name = "tpr";
    tpr_plugin.prettyname = "Gromacs Binary Topology";
    tpr_plugin.author = "Josh Vermaas";
    tpr_plugin.majorv = 2023;
    tpr_plugin.minorv = 0;//Corresponds to the Gromacs version I was basing this on.
    tpr_plugin.is_reentrant = VMDPLUGIN_THREADUNSAFE;
    tpr_plugin.filename_extension = "tpr";
    tpr_plugin.open_file_read = open_tpr_read;
    tpr_plugin.read_timestep_metadata = read_tpr_timestep_metadata;
    tpr_plugin.read_structure = read_tpr_structure;
    tpr_plugin.read_next_timestep = read_tpr_timestep;
    tpr_plugin.read_bonds = read_tpr_bonds;
    tpr_plugin.read_angles = read_tpr_angles;
    tpr_plugin.close_file_read = close_tpr_read;
    return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
    (*cb)(v, (vmdplugin_t *)&tpr_plugin);
    return 0;
}

VMDPLUGIN_API int VMDPLUGIN_fini(void) { return VMDPLUGIN_SUCCESS; }

#else

int main (int argc, char *argv[]) {
	int natoms;
	tprdata *tprdat = NULL;
	for (int fcount=1; fcount < argc; fcount++) {
		printf("Attempting to read %s\n", argv[fcount]);
		tprdat = (tprdata*)open_tpr_read(argv[fcount], NULL, &natoms);
		if (tprdat != NULL){
	        printf("Total number of atoms: %d, %d\n", natoms, tprdat->natoms);
			printf("Finished initial reading\n");
			molfile_timestep_t *ts = new molfile_timestep_t;
			ts->coords = new float[3*tprdat->natoms];
			ts->velocities = new float[3*tprdat->natoms];
			read_tpr_timestep(tprdat, tprdat->natoms, ts);
	    }
	    else {
	    	fprintf(stderr, "open_tpr_read failed\n");
	        return NULL;
	    }
	    FILE *fin = tprdat->mf->f;
		fseek(fin, 0L, SEEK_END);
		long length = ftell(fin);
        printf("END : %ld\n", length);
		fclose(fin);
	}
	return 0;
}

#endif

