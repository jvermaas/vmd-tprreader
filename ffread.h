//From topology/topology.h
enum {
    egcTC,    egcENER,   egcACC, egcFREEZE,
    egcUser1, egcUser2,  egcVCM, egcCompressedX,
    egcORFIT, egcQMMM,
    egcNR
};

//From topology/ifunc.h
enum
{
    F_BONDS,
    F_G96BONDS,
    F_MORSE,
    F_CUBICBONDS,
    F_CONNBONDS,
    F_HARMONIC,
    F_FENEBONDS,
    F_TABBONDS,
    F_TABBONDSNC,
    F_RESTRBONDS,
    F_ANGLES,
    F_G96ANGLES,
    F_RESTRANGLES,
    F_LINEAR_ANGLES,
    F_CROSS_BOND_BONDS,
    F_CROSS_BOND_ANGLES,
    F_UREY_BRADLEY,
    F_QUARTIC_ANGLES,
    F_TABANGLES,
    F_PDIHS,
    F_RBDIHS,
    F_RESTRDIHS,
    F_CBTDIHS,
    F_FOURDIHS,
    F_IDIHS,
    F_PIDIHS,
    F_TABDIHS,
    F_CMAP,
    F_GB12,
    F_GB13,
    F_GB14,
    F_GBPOL,
    F_NPSOLVATION,
    F_LJ14,
    F_COUL14,
    F_LJC14_Q,
    F_LJC_PAIRS_NB,
    F_LJ,
    F_BHAM,
    F_LJ_LR,
    F_BHAM_LR,
    F_DISPCORR,
    F_COUL_SR,
    F_COUL_LR,
    F_RF_EXCL,
    F_COUL_RECIP,
    F_LJ_RECIP,
    F_DPD,
    F_POLARIZATION,
    F_WATER_POL,
    F_THOLE_POL,
    F_ANHARM_POL,
    F_POSRES,
    F_FBPOSRES,
    F_DISRES,
    F_DISRESVIOL,
    F_ORIRES,
    F_ORIRESDEV,
    F_ANGRES,
    F_ANGRESZ,
    F_DIHRES,
    F_DIHRESVIOL,
    F_CONSTR,
    F_CONSTRNC,
    F_SETTLE,
    F_VSITE2,
    F_VSITE2FD,
    F_VSITE3,
    F_VSITE3FD,
    F_VSITE3FAD,
    F_VSITE3OUT,
    F_VSITE4FD,
    F_VSITE4FDN,
    F_VSITEN,
    F_COM_PULL,
    F_DENSITYFITTING,
    F_EQM,
    F_EPOT,
    F_EKIN,
    F_ETOT,
    F_ECONSERVED,
    F_TEMP,
    F_VTEMP,
    F_PDISPCORR,
    F_PRES,
    F_DVDL_CONSTR,
    F_DVDL,
    F_DKDL,
    F_DVDL_COUL,
    F_DVDL_VDW,
    F_DVDL_BONDED,
    F_DVDL_RESTRAINT,
    F_DVDL_TEMPERATURE, /* not calculated for now, but should just be the energy (NVT) or enthalpy (NPT), or 0 (NVE) */
    F_NRE /* This number is for the total number of energies      */
};

//From tpxio.c
typedef struct {
    int fvnr;  /* file version number in which the function type first appeared */
    int ftype; /* function type */
} t_ftupd;
static const t_ftupd ftupd[] = {
    { 20, F_CUBICBONDS        },
    { 20, F_CONNBONDS         },
    { 20, F_HARMONIC          },
    { 34, F_FENEBONDS         },
    { 43, F_TABBONDS          },
    { 43, F_TABBONDSNC        },
    { 70, F_RESTRBONDS        },
    { 98, F_RESTRANGLES },
    { 76, F_LINEAR_ANGLES     },
    { 30, F_CROSS_BOND_BONDS  },
    { 30, F_CROSS_BOND_ANGLES },
    { 30, F_UREY_BRADLEY      },
    { 34, F_QUARTIC_ANGLES    },
    { 43, F_TABANGLES         },
    { 98, F_RESTRDIHS },
    { 98, F_CBTDIHS },
    { 26, F_FOURDIHS          },
    { 26, F_PIDIHS            },
    { 43, F_TABDIHS           },
    { 65, F_CMAP              },
    { 60, F_GB12              },
    { 61, F_GB13              },
    { 61, F_GB14              },
    { 72, F_GBPOL             },
    { 72, F_NPSOLVATION       },
    { 41, F_LJC14_Q           },
    { 41, F_LJC_PAIRS_NB      },
    { 32, F_BHAM_LR           },
    { 32, F_RF_EXCL           },
    { 32, F_COUL_RECIP        },
    { 93, F_LJ_RECIP          },
    { 46, F_DPD               },
    { 30, F_POLARIZATION      },
    { 36, F_THOLE_POL         },
    { 90, F_FBPOSRES          },
    { 22, F_DISRESVIOL        },
    { 22, F_ORIRES            },
    { 22, F_ORIRESDEV         },
    { 26, F_DIHRES            },
    { 26, F_DIHRESVIOL        },
    { 49, F_VSITE4FDN         },
    { 50, F_VSITEN            },
    { 46, F_COM_PULL          },
    { 20, F_EQM               },
    { 46, F_ECONSERVED        },
    { 69, F_VTEMP},
    { 66, F_PDISPCORR         },
    { 54, F_DVDL_CONSTR       },
    { 76, F_ANHARM_POL        },
    { 79, F_DVDL_COUL         },
    { 79, F_DVDL_VDW         },
    { 79, F_DVDL_BONDED      },
    { 79, F_DVDL_RESTRAINT    },
    { 79, F_DVDL_TEMPERATURE  },
    { 117, F_DENSITYFITTING },
    { 118, F_VSITE2FD }
};
#define asize(a) ((int)(sizeof(a)/sizeof((a)[0])))
#define NFTUPD asize(ftupd)


void printString(XDR *xdrs) {
	int len, i;
	char buf[STRLEN];
	xdr_int(xdrs, &len);
	xdr_opaque(xdrs, buf, len);
	/*for (i = 0; i < len; i++) {
		printf("%c", buf[i]);
	}
	printf("\n");*/
}
// void saveString(XDR* xdrs, char* saveloc) {
// 	int len = 0;
//     int i;
// 	char buf[STRLEN];
//     xdr_int(xdrs, &len);
//     printf("Length : %ld ", len);
// 	xdr_opaque(xdrs, buf, len);
// 	for (i = 0; i < MIN(len, SAVELEN-1); i++) {
// 		printf("%c", buf[i]);
// 		saveloc[i] = buf[i];
// 	}
//     saveloc[i] = '\0';
// 	printf("\n");
// }

void saveString(XDR* xdrs, char* saveloc, int genversion) {
    long int len = 0;
    int i;
    char buf[STRLEN];
    xdr_int64_t(xdrs, &len);    
    xdr_opaque(xdrs, buf, len);
    for (i = 0; i < MIN(int(len), (SAVELEN-1)); i++) {
        saveloc[i] = buf[i];
    }
    if (genversion >= 27) {
        int j = len % 4;
        //If it wasn't divisible by 4 in the first place, shift the file pointer back.
        //XDR reads are always byte aligned, but starting in version 27, GROMACS stopped
        //using the XDR library exclusively.
        if (j) {
            xdr_setpos(xdrs, xdr_getpos(xdrs) - (4-j));
        }  
    }
    saveloc[i] = '\0';
}

template<typename real>
void readparams (XDR* xdrs, int file_version, int ftype) {
	switch (ftype)
    {
        case F_ANGLES:
        case F_G96ANGLES:
        case F_BONDS:
        case F_G96BONDS:
        case F_HARMONIC:
        case F_IDIHS:
        	//Read the 4 required harmonics.
        	readReal<real>(xdrs);
        	readReal<real>(xdrs);
        	readReal<real>(xdrs);
        	readReal<real>(xdrs);
            break;
        case F_RESTRANGLES:
            readReal<real>(xdrs);
            readReal<real>(xdrs);
        case F_LINEAR_ANGLES:
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            break;
        case F_FENEBONDS:
        	readReal<real>(xdrs);
            readReal<real>(xdrs);
            break;
        case F_RESTRBONDS:
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            break;
        case F_TABBONDS:
        case F_TABBONDSNC:
        case F_TABANGLES:
        case F_TABDIHS:
            readReal<real>(xdrs);
            readInt(xdrs);
            readReal<real>(xdrs);
            break;
        case F_CROSS_BOND_BONDS:
        	readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            break;
        case F_CROSS_BOND_ANGLES:
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            break;
        case F_UREY_BRADLEY:
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            if (file_version >= 79)
            {
                readReal<real>(xdrs);
	            readReal<real>(xdrs);
	            readReal<real>(xdrs);
	            readReal<real>(xdrs);
            }
            break;
        case F_QUARTIC_ANGLES:
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            break;
        case F_BHAM:
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            break;
        case F_MORSE:
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            if (file_version >= 79)
            {
                readReal<real>(xdrs);
	            readReal<real>(xdrs);
	            readReal<real>(xdrs);
            }
            break;
        case F_CUBICBONDS:
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            break;
        case F_CONNBONDS:
            break;
        case F_POLARIZATION:
            readReal<real>(xdrs);
            break;
        case F_ANHARM_POL:
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            break;
        case F_WATER_POL:
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            break;
        case F_THOLE_POL:
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            break;
        case F_LJ:
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            break;
        case F_LJ14:
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            break;
        case F_LJC14_Q:
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            break;
        case F_LJC_PAIRS_NB:
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            break;
        case F_PDIHS:
        case F_PIDIHS:
        case F_ANGRES:
        case F_ANGRESZ:
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readInt(xdrs);
            break;
        case F_RESTRDIHS:
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            break;
        case F_DISRES:
            readInt(xdrs);
            readInt(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            break;
        case F_ORIRES:
            readInt(xdrs);
            readInt(xdrs);
            readInt(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            break;
        case F_DIHRES:
            if (file_version < 82)
            {
                readInt(xdrs);
            	readInt(xdrs);
            }
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            if (file_version >= 82)
            {
                readReal<real>(xdrs);
	            readReal<real>(xdrs);
	            readReal<real>(xdrs);
            }
            break;
        case F_POSRES:
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            break;
        case F_FBPOSRES:
            readInt(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            break;
        case F_CBTDIHS:
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            break;
        case F_RBDIHS:
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            break;
        case F_FOURDIHS:
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            break;
        case F_CONSTR:
        case F_CONSTRNC:
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            break;
        case F_SETTLE:
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            break;
        case F_VSITE2FD:
        case F_VSITE2:
            readReal<real>(xdrs);
            break;
        case F_VSITE3:
        case F_VSITE3FD:
        case F_VSITE3FAD:
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            break;
        case F_VSITE3OUT:
        case F_VSITE4FD:
        case F_VSITE4FDN:
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            readReal<real>(xdrs);
            break;
        case F_VSITEN:
            readInt(xdrs);
            readReal<real>(xdrs);
            break;
        case F_GB12:
        case F_GB13:
        case F_GB14:
            /* We got rid of some parameters in version 68 */
            if (file_version < 68)
            {
                readReal<real>(xdrs);
	            readReal<real>(xdrs);
	            readReal<real>(xdrs);
	            readReal<real>(xdrs);

            }
            if (file_version < 113)
            {
                readReal<real>(xdrs);
                readReal<real>(xdrs);
                readReal<real>(xdrs);
                readReal<real>(xdrs);
                readReal<real>(xdrs);
            }
            break;
        case F_CMAP:
            readInt(xdrs);
            readInt(xdrs);
            break;
    }
}