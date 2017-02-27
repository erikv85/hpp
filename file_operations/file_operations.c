#include "file_operations.h"
#include <stdio.h>

int read_doubles_from_file(int n, double* p, const char* fileName) {
  /* Open input file and determine its size. */
  FILE* input_file = fopen(fileName, "rb");
  if(!input_file) {
    printf("read_doubles_from_file error: failed to open input file '%s'.\n", fileName);
    return -1;
  }
  /* Get filesize using fseek() and ftell(). */
  fseek(input_file, 0L, SEEK_END);
  size_t fileSize = ftell(input_file);
  /* Now use fseek() again to set file position back to beginning of the file. */
  fseek(input_file, 0L, SEEK_SET);
  if(fileSize != n * sizeof(double)) {
    printf("read_doubles_from_file error: size of input file '%s' does not match the given n.\n", fileName);
    printf("For n = %d the file size is expected to be (n * sizeof(double)) = %lu but the actual file size is %lu.\n",
	   n, n * sizeof(double), fileSize);
    return -1;
  }
  /* Read contents of input_file into buffer. */
  size_t noOfItemsRead = fread(p, sizeof(char), fileSize, input_file);
  if(noOfItemsRead != fileSize) {
    printf("read_doubles_from_file error: failed to read file contents into buffer.\n");
    return -1;
  }
  /* OK, now we have the file contents in the buffer.
     Since we are finished with the input file, we can close it now. */
  if(fclose(input_file) != 0) {
    printf("read_doubles_from_file error: error closing input file.\n");
    return -1;
  }
  /* Return 0 to indicate success. */
  return 0;
}

int write_doubles_to_file(int n, const double* p, const char* fileName) {
  if(n <= 0) {
    printf("write_doubles_to_file error: (n <= 0).\n");
    return -1;
  }
  FILE* output_file = fopen(fileName, "wb");
  if(!output_file) {
    printf("write_doubles_to_file error: failed to open output file '%s'.\n", fileName);
    return -1;
  }
  /* Now write to the file using fwrite(). */
  size_t noOfItemsWritten = fwrite(p, sizeof(double), n, output_file);
  if(noOfItemsWritten != (size_t)n) {
    printf("write_doubles_to_file error: failed to write buffer contents to output file.\n");
    return -1;
  }
  /* OK, we have successfully written to the output file. Now close it. */
  if(fclose(output_file) != 0) {
    printf("write_doubles_to_file error: error closing output file.\n");
    return -1;
  }
  /* Return 0 to indicate success. */
  return 0;
}
