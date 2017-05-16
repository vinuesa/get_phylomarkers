s/^#undef  *\([ABCDEFGHIJKLMNOPQRSTUVWXYZ_]\)/#undef CLUSTAL_OMEGA_\1/
s/^#undef  *\([abcdefghijklmnopqrstuvwxyz]\)/#undef _clustal_omega_\1/
s/^#define  *\([ABCDEFGHIJKLMNOPQRSTUVWXYZ_][abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_]*\)\(.*\)/#ifndef CLUSTAL_OMEGA_\1\
#define CLUSTAL_OMEGA_\1\2\
#endif/
s/^#define  *\([abcdefghijklmnopqrstuvwxyz][abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_]*\)\(.*\)/#ifndef _clustal_omega_\1\
#define _clustal_omega_\1\2\
#endif/
