#include "tprplugin.C"
int main (int argc, char *argv[]) {
	for (int fcount=1; fcount < argc; fcount++) {
		printf("Attempting to read %s\n", argv[fcount]);
		tprdata *tprdat = NULL;
    	FILE *fin;
    	char buf[STRLEN];
    	md_file *mf = new md_file;
    	int i, j;
		fin = fopen(argv[fcount], "rb");
		mf->f = fin;
		if (trx_int(mf, &i)) {
	    	fprintf(stderr, "tprplugin) Could not read initial integer from file.\n");
	        return NULL;
	    }
		if (i > STRLEN) {//If i value is large, everything in the file should be endian swapped.
			mf->rev = 1;
			printf("Reverse endian\n");
		}
		j = tpr_string(mf, buf, STRLEN);
		for (i = 0; i < j; i++) {
	        printf("%c", buf[i]);
	    }
	    printf("\n");
		if (trx_int(mf, &(mf->prec))) {
			fprintf(stderr, "tprplugin) Could not read precision from file.\n");
	        return NULL;
		}
	    printf("I have the precision! %d\n", mf->prec);
	    if (mf->prec == 4) {
	        tprdat = new tprdata;
	        memset(tprdat, 0, sizeof(tprdata));
	        tprdat->mf = mf;
	        if (readtprAfterPrecision(tprdat) != MOLFILE_SUCCESS) {
	        	printf("I wasn't successful\n");
	            delete tprdat;
	            return NULL;
	        }
	        int natoms = tprdat->natoms;
	        printf("Total number of atoms: %d, %d\n", natoms, tprdat->natoms);
			printf("Finished initial reading\n");
			molfile_timestep_t *ts = new molfile_timestep_t;
			ts->coords = new float[3*tprdat->natoms];
			ts->velocities = new float[3*tprdat->natoms];
			read_tpr_timestep(tprdat, tprdat->natoms, ts);
	    }
	    else {
	    	fprintf(stderr, "tprplugin) Illegal precision (requires single)\n");
	        return NULL;
	    }
		fseek(fin, 0L, SEEK_END);
		long length = ftell(fin);
        printf("END : %ld\n", length);
		fclose(fin);
	}
	return 0;
}
