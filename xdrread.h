template<typename real>
float readReal (XDR* xdrs) {
    real dummy;
    switch (sizeof(real)) {
        case 4:
        xdr_float(xdrs, (float*)&dummy); break;
        case 8:
        xdr_double(xdrs, (double*)&dummy); break;
    }
    return (float)dummy;
}
inline int readInt (XDR* xdrs) {
    int tmp;
    xdr_int (xdrs, &tmp);
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
