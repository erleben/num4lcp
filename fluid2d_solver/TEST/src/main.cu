#include <cusp/gallery/poisson.h>
#include <cusp/print.h>
#include <cusp/coo_matrix.h>
#include <cusp/csr_matrix.h>


int main(int argc, char **argv)
{
  std::vector<int>   row_indices(13);
  std::vector<int>   column_indices(13);
  std::vector<float> values(13);

  row_indices[0]   = 0; column_indices[0]  = 0; values[0]  = 2.0;
  row_indices[1]   = 1; column_indices[1]  = 1; values[1]  = 2.0;
  row_indices[2]   = 2; column_indices[2]  = 2; values[2]  = 2.0;
  row_indices[3]   = 3; column_indices[3]  = 3; values[3]  = 2.0;
  row_indices[4]   = 4; column_indices[4]  = 4; values[4]  = 2.0;

  row_indices[8]   = 4; column_indices[8]  = 3; values[8]  = 1.0;
  row_indices[5]   = 1; column_indices[5]  = 0; values[5]  = 1.0;
  row_indices[6]   = 2; column_indices[6]  = 1; values[6]  = 1.0;
  row_indices[7]   = 3; column_indices[7]  = 2; values[7]  = 1.0;

  row_indices[12]  = 3; column_indices[12] = 4; values[12] = 3.0;
  row_indices[9]   = 0; column_indices[9]  = 1; values[9]  = 3.0;
  row_indices[10]  = 1; column_indices[10] = 2; values[10] = 3.0;
  row_indices[11]  = 2; column_indices[11] = 3; values[11] = 3.0;

  cusp::coo_matrix<int, float, cusp::host_memory> Acoo(5,5,13);

  Acoo.row_indices    = row_indices;
  Acoo.column_indices = column_indices;
  Acoo.values         = values;

  std::cout << "Acoo = " << std::endl;
  cusp::print(Acoo);

  if(Acoo.is_sorted_by_row())
    std::cout << "SORTED BY ROW" << std::endl;

  if(Acoo.is_sorted_by_row_and_column())
    std::cout << "SORTED BY ROW AND COLUMN" << std::endl;

  Acoo.sort_by_row_and_column();


  if(Acoo.is_sorted_by_row())
    std::cout << "SORTED BY ROW" << std::endl;

  if(Acoo.is_sorted_by_row_and_column())
    std::cout << "SORTED BY ROW AND COLUMN" << std::endl;

  std::cout << "Acoo = " << std::endl;
  cusp::print(Acoo);

  cusp::coo_matrix<int, float, cusp::host_memory> Acsr = Acoo;

  std::cout << "Acsr = " << std::endl;
  cusp::print(Acsr);

  std::cout << "Acoo = " << std::endl;
  cusp::print(Acoo);

	return 0;
}
