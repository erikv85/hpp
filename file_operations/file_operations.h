/* read_doubles_from_file: this function reads n double numbers from
   the given file into the buffer pointer to by p, which must have
   been allocated before calling this function. The function returns 0
   if successful, or -1 if it failed (for example if the file could
   not be found or if the file had wrong size). */
int read_doubles_from_file(int n, double* p, const char* fileName);

/* write_doubles_to_file: this function writes n double numbers to the
   given file from the buffer pointer to by p. The function returns 0
   if successful, or -1 if it failed. */
int write_doubles_to_file(int n, const double* p, const char* fileName);

