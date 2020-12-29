#ifdef _WIN32
#include "gmx_internal_xdr.h"
#endif

#include "tprplugin.C"
int main (int argc, char *argv[]) {
	for (int fcount=1; fcount < argc; fcount++) {
		printf("Attempting to read %s\n", argv[fcount]);
		FILE *fin;
		int i, j, tmp;
		fin = fopen(argv[fcount], "rb");
		tprdata *tprdat = new tprdata;
		tprdat->f = fin;
		i = readIntTPR(tprdat);
		printf("%d\n", i);
		i = readIntTPR(tprdat);
		printf("%d\n", i);
		i = readIntTPR(tprdat);
		printf("%d\n", i);
		fseek(fin, 0, SEEK_SET);
		XDR *xdrs = new XDR;
		xdrstdio_create(xdrs, fin, XDR_DECODE);
		xdr_int(xdrs, &i);
		//printf("%d\n", i);
		printString(xdrs);
		int precision = readInt(xdrs);
		if (precision == 4) {
			tprdat = new tprdata;
			tprdat->f = fin;
			tprdat->xdrptr = xdrs;
			if (readtprAfterPrecision(tprdat) != MOLFILE_SUCCESS) {
				delete tprdat;
				return NULL;
			}
			int natoms = tprdat->natoms;
			//printf("Total number of atoms: %d, %d\n", *natoms, tprdat->natoms);
			//printf("Finished initial reading\n");
			molfile_timestep_t *ts = new molfile_timestep_t;
			ts->coords = new float[3*tprdat->natoms];
			ts->velocities = new float[3*tprdat->natoms];
			read_tpr_timestep(tprdat, tprdat->natoms, ts);
		}
		else {
			printf("Illegal precision (requires single)\n");
			return NULL;
		}
		fseek(fin, 0L, SEEK_END);
		long length = ftell(fin);
        printf("END : %ld\n", length);
		fclose(fin);
	}
	return 0;
}
