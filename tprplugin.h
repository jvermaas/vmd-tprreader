/*
Header that describes what is all *in* a tpr file, defining energies,
force terms, and other stuff that we need to be interoperable with a tpr.
*/

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
    F_VSITE1,
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
    F_ENNPOT,
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

//This holds the meat and potatoes, the results after extracting all the topology info from the tpr file.
struct tprdata {
    int precision;
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
    md_file *mf;
    int has_velocities;
    int readcoordinates;
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
    { 118, F_VSITE2FD },
    { 121, F_VSITE1 },
    { 137, F_ENNPOT}
};
#define asize(a) ((int)(sizeof(a)/sizeof((a)[0])))
#define NFTUPD asize(ftupd)

//Equivalent to do_iparams
void readparams (md_file *mf, int file_version, int ftype) {
    switch (ftype)
    {
        case F_ANGLES:
        case F_G96ANGLES:
        case F_BONDS:
        case F_G96BONDS:
        case F_HARMONIC:
        case F_IDIHS:
            //Read the 4 required harmonics.
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            break;
        case F_RESTRANGLES:
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            break;
        case F_LINEAR_ANGLES:
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            break;
        case F_FENEBONDS:
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            break;
        case F_RESTRBONDS:
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            break;
        case F_TABBONDS:
        case F_TABBONDSNC:
        case F_TABANGLES:
        case F_TABDIHS:
            trx_real(mf, NULL);
            trx_int(mf, NULL);
            trx_real(mf, NULL);
            break;
        case F_CROSS_BOND_BONDS:
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            break;
        case F_CROSS_BOND_ANGLES:
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            break;
        case F_UREY_BRADLEY:
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            if (file_version >= 79)
            {
                trx_real(mf, NULL);
                trx_real(mf, NULL);
                trx_real(mf, NULL);
                trx_real(mf, NULL);
            }
            break;
        case F_QUARTIC_ANGLES:
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            break;
        case F_BHAM:
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            break;
        case F_MORSE:
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            if (file_version >= 79)
            {
                trx_real(mf, NULL);
                trx_real(mf, NULL);
                trx_real(mf, NULL);
            }
            break;
        case F_CUBICBONDS:
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            break;
        case F_CONNBONDS:
            break;
        case F_POLARIZATION:
            trx_real(mf, NULL);
            break;
        case F_ANHARM_POL:
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            break;
        case F_WATER_POL:
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            break;
        case F_THOLE_POL:
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            break;
        case F_LJ:
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            break;
        case F_LJ14:
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            break;
        case F_LJC14_Q:
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            break;
        case F_LJC_PAIRS_NB:
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            break;
        case F_PDIHS:
        case F_PIDIHS:
        case F_ANGRES:
        case F_ANGRESZ:
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_int(mf, NULL);
            break;
        case F_RESTRDIHS:
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            break;
        case F_DISRES:
            trx_int(mf, NULL);
            trx_int(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            break;
        case F_ORIRES:
            trx_int(mf, NULL);
            trx_int(mf, NULL);
            trx_int(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            break;
        case F_DIHRES:
            if (file_version < 82)
            {
                trx_int(mf, NULL);
                trx_int(mf, NULL);
            }
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            if (file_version >= 82)
            {
                trx_real(mf, NULL);
                trx_real(mf, NULL);
                trx_real(mf, NULL);
            }
            break;
        case F_POSRES:
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            break;
        case F_FBPOSRES:
            trx_int(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            break;
        case F_CBTDIHS:
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            break;
        case F_RBDIHS:
        case F_FOURDIHS:
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            break;
        case F_CONSTR:
        case F_CONSTRNC:
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            break;
        case F_SETTLE:
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            break;
        case F_VSITE1: break; //VSite1 has no parameters
        case F_VSITE2FD:
        case F_VSITE2:
            trx_real(mf, NULL);
            break;
        case F_VSITE3:
        case F_VSITE3FD:
        case F_VSITE3FAD:
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            break;
        case F_VSITE3OUT:
        case F_VSITE4FD:
        case F_VSITE4FDN:
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            trx_real(mf, NULL);
            break;
        case F_VSITEN:
            trx_int(mf, NULL);
            trx_real(mf, NULL);
            break;
        case F_GB12:
        case F_GB13:
        case F_GB14:
            /* We got rid of some parameters in version 68 */
            if (file_version < 68)
            {
                trx_real(mf, NULL);
                trx_real(mf, NULL);
                trx_real(mf, NULL);
                trx_real(mf, NULL);

            }
            if (file_version < 113)
            {
                trx_real(mf, NULL);
                trx_real(mf, NULL);
                trx_real(mf, NULL);
                trx_real(mf, NULL);
                trx_real(mf, NULL);
            }
            break;
        case F_CMAP:
            trx_int(mf, NULL);
            trx_int(mf, NULL);
            break;
    }
}
