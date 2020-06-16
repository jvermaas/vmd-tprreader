#include "tprplugin.C"
int main (int argc, char *argv[]) {
	for (int fcount=1; fcount < argc; fcount++) {
		printf("Attempting to read %s\n", argv[fcount]);
		FILE *fin;
		XDR *xdrs = new XDR;
		int i, j, tmp;
		int precision;
		fin = fopen(argv[fcount], "rb");
		xdrstdio_create(xdrs, fin, XDR_DECODE);
		xdr_int(xdrs, &i);
		printf("%d\n", i);
		printString(xdrs);
		precision = readInt(xdrs);
		printf("%d\n", precision);
		if (precision == 4) {
			tprdata *tprdat = new tprdata;
			tprdat->f = fin;
			tprdat->xdrptr = xdrs;
			readtprAfterPrecision(tprdat);

			molfile_timestep_t *ts = new molfile_timestep_t;
			ts->coords = new float[3*tprdat->natoms];
			read_tpr_timestep(tprdat, tprdat->natoms, ts);
			//return 0;
		}
		// else if (precision == 8) {
		// 	tprdata<double> *tprdat = new tprdata;
		// 	readtprAfterPrecision<double>(&xdrs, tprdat);
		// 	return 0;
		// }
		else {
			printf("Illegal precision (requires single)\n");
			return 1;
		}
		fclose(fin);
		delete xdrs;
	}
	return 0;
}