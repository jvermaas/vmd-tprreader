#ifdef _WIN32
#define MIN(a, b) ((a)<(b)? (a):(b))
#endif

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
void printStringTPR(tprdata* tpr) {
    int len, i;
    char buf[STRLEN];
    len = readIntTPR(tpr);
    fread(buf, len, 1, tpr->f);
    for (i = 0; i < len; i++) {
        printf("%c", buf[i]);
    }
    printf("\n");
}
// void saveString(XDR* xdrs, char* saveloc) {
//  int len = 0;
//     int i;
//  char buf[STRLEN];
//     xdr_int(xdrs, &len);
//     printf("Length : %ld ", len);
//  xdr_opaque(xdrs, buf, len);
//  for (i = 0; i < MIN(len, SAVELEN-1); i++) {
//      printf("%c", buf[i]);
//      saveloc[i] = buf[i];
//  }
//     saveloc[i] = '\0';
//  printf("\n");
// }

void saveString(XDR* xdrs, char* saveloc, int genversion) {
#ifdef _WIN32
    long long int len = 0;
#else
    long len = 0;
#endif
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
            break;
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
        case F_VSITE1: break; //VSite1 has no parameters
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
