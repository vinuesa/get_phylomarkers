bool ppma_check_data ( int xsize, int ysize, int maxrgb, int *rarray,
       int *garray, int *barray );
bool ppma_example ( int xsize, int ysize, int *rarray, int *garray, int *barray );

bool ppma_read ( char *file_in_name, int *xsize, int *ysize, int *maxrgb,
       int **rarrary, int **garray, int **barray );
bool ppma_read_data ( ifstream &file_in, int xsize, int ysize, int *rarray,
       int *garray, int *barray );
bool ppma_read_header ( ifstream &file_in, int *xsize, int *ysize, int *maxrgb );
bool ppma_read_test ( char *file_in_name );

bool ppma_write ( char *file_out_name, int xsize, int ysize, int *rarray, 
      int *garray, int *barray );
bool ppma_write_data ( ofstream &file_out, int xsize, int ysize, int *rarray,
       int *garray, int *barray );
bool ppma_write_header ( ofstream &file_out, char *file_out_name, int xsize, 
       int ysize, int maxrgb );
bool ppma_write_test ( char *file_out_name );
