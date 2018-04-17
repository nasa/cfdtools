//
// A demonstration of using the KDTREE2 C++ routines, and timing.
// This file is in the public domain.
//

#include "kdtree2.hpp"

#include <boost/multi_array.hpp>
#include <boost/random.hpp>

using namespace boost; 


static minstd_rand generator(42u); 
static uniform_real<> uni_dist(0,1); 
variate_generator<minstd_rand&,uniform_real<> > uni(generator,uni_dist); 

float random_variate() {
  // between [0,1)
  return(uni()); 
}

//
// define, for convenience a 2d array of floats. 
//
typedef multi_array<float,2> array2dfloat;


#include <ctime>

float time_a_search(kdtree2* tree, int nn, int nsearch) {
  int dim = tree->dim;
  vector<float> query(dim);
  kdtree2_result_vector result; 

  clock_t t0, t1; 

  t0 = clock();

  for (int i=0; i<nsearch;i++) {
      for (int j=0; j<dim; j++) query[j] = random_variate(); 
      tree->n_nearest(query,nn,result);
  }

  t1 = clock(); 

  return(static_cast<float> 
	 (static_cast<double> (t1-t0) / static_cast<double> (CLOCKS_PER_SEC) ));
}

void time_random_searches(kdtree2* tree, int nn) {
  // emit the number of searches per second.
  int nsearch;

  nsearch = 50;
  for(;;) {
    float t = time_a_search(tree,nn,nsearch); 
    if (t < 1.0) {
      nsearch *= 5; 
      continue; 
    } else {
      float sps = float(nsearch) / t ;
      cout << "C++ impl, for nn=" << nn << " searches/sec = " << sps << "\n";
      return;
    }
  }
}

int main() {
  array2dfloat data(extents[10][3]);  // declare a 10000 x 3 array.
  array2dfloat realdata; 

  // notice it is in C-standard layout. 
  kdtree2* tree;
  kdtree2_result_vector res; 
  int N, dim;

  if (false) {
    for (int i=0; i<10; i++) {
      for (int j=0; j<3; j++)
	data[i][j] = static_cast<float> (3*i+j);
    }
    tree = new kdtree2(data,true); 
    
    //    tree->dump_data(); 
    //data[0][0]=666.0;  // mutate it underneath.  DO NOT DO THIS IN REAL USE
    //tree->dump_data(); // test to see it change there.
    
    tree->n_nearest_around_point(1,1,1,res); 
    for (unsigned int i=0; i<res.size(); i++) {
      printf("result[%d]= (%d,%f)\n",i,res[i].idx,res[i].dis);
    }
    
    delete tree;
  }
  printf("Give me N, and dim (e.g. '1000 3').  No commas!");
  scanf("%d %d",&N,&dim);
  printf("I found N=%d,dim=%d\n",N,dim);
  realdata.resize(extents[N][dim]); 
    
  for (int i=0; i<N; i++) {
    for (int j=0; j<dim; j++) 
      realdata[i][j] = random_variate();
  }
  
  tree = new kdtree2(realdata,true);
  tree->sort_results = true;
  cout << "Tree created, now testing against brute force..."; 
  {
    vector<float> query(dim); 
    kdtree2_result_vector result, resultbrute;
    int nn = 10; 

    for (int i=0; i<50; i++) {
      for (int j=0; j<dim; j++) query[j] = random_variate(); 

      tree->n_nearest_brute_force(query,nn,resultbrute);
      tree->n_nearest(query,nn,result); // search for 10 of them.

      for (int k=0; k<nn; k++) {
	if ((resultbrute[k].dis != result[k].dis) ||
	    (resultbrute[k].idx != result[k].idx)) {
	  cout << "Mismatch! nn=" << k << " brute=[" << 
	    resultbrute[k].dis << "," << resultbrute[k].idx << 
	    "] tree=[" << result[k].dis << "," << result[k].idx << "]\n"; 
	}
      }
    }
  }
  cout << "\nTesting complete.  Now testing timing...\n";
  tree->sort_results = false;

  {

    int nnarray[] = {1,5,10,25,500} ;
    for (int i=0; i< 5; i++) {
      time_random_searches(tree,nnarray[i]); 
    }
    
  }


}

