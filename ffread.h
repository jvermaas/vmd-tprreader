/*
This is the header to read in the force field parameters encoded in a tpr file.
VMD does not need these per se, or know what to do with them, but because
these are included *before* the positions that we do want, we need a way of reading
past them.
*/
// #define MIN(a, b) ((a)<(b)? (a):(b))

// void printString(XDR *xdrs) {
//     int len, i;
//     char buf[STRLEN];
//     xdr_int(xdrs, &len);
//     xdr_opaque(xdrs, buf, len);
//     /*for (i = 0; i < len; i++) {
//         printf("%c", buf[i]);
//     }
//     printf("\n");*/
// }
// void printStringTPR(tprdata* tpr) {
//     int len, i;
//     char buf[STRLEN];
//     len = readIntTPR(tpr);
//     fread(buf, len, 1, tpr->f);
//     for (i = 0; i < len; i++) {
//         printf("%c", buf[i]);
//     }
//     printf("\n");
// }
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

// void saveString(XDR* xdrs, char* saveloc, int genversion) {
// #ifdef _WIN32
//     long long int len = 0;
// #else
//     long len = 0;
// #endif
//     int i;
//     char buf[STRLEN];
//     xdr_int64_t(xdrs, &len);
//     xdr_opaque(xdrs, buf, len);
//     for (i = 0; i < MIN(int(len), (SAVELEN-1)); i++) {
//         saveloc[i] = buf[i];
//     }
//     if (genversion >= 27) {
//         int j = len % 4;
//         //If it wasn't divisible by 4 in the first place, shift the file pointer back.
//         //XDR reads are always byte aligned, but starting in version 27, GROMACS stopped
//         //using the XDR library exclusively.
//         if (j) {
//             xdr_setpos(xdrs, xdr_getpos(xdrs) - (4-j));
//         }
//     }
//     saveloc[i] = '\0';
// }

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
