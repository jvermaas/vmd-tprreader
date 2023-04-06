/*
These define our little intermediate procedures to read reals, ints, chars, and
other things based on the XDR spec. Unfortunately, XDR has become super niche,
and so even things like python are deprecating the libraries that handle XDRs.
*/
float readReal (md_file* mf) {
    float tmp;
    trx_real(mf, &tmp);
    return tmp;
}
inline int readInt (XDR* xdrs) {
    int tmp;
    trx_int(xdrs, &tmp);
    return tmp;
}
inline long long int readInt64 (XDR* xdrs) {
#ifdef _WIN32
    long long int tmp;
#else
    long tmp;
#endif
    xdr_int64_t (xdrs, &tmp);
    return tmp;
}
inline unsigned char readChar (XDR* xdrs) {
    unsigned char tmp;
    xdr_u_char (xdrs, &tmp);
    return tmp;
}
void readvector (XDR *xdrs, float *arr, int len) {
    int i;
    for (i = 0; i < len; i++) {
        arr[i] = readReal<float>(xdrs);
    }
}

void readintvector (XDR *xdrs, int *arr, int len) {
    int i;
    for (i = 0; i < len; i++) {
        arr[i] = readInt(xdrs);
    }
}
