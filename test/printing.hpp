
template <class Type>
inline
Type print(matrix<Type> m){
  int rows=m.rows();
  int cols=m.cols();

    for(int x=0;x<rows;x++)
    {
        for(int y=0;y<cols;y++) 
        {
            std::cout<<m(x,y) << " ";
        }
    std::cout<<std::endl; 
    }

  return 0;
}


template <class Type>
inline
Type print(array<Type> m){
  int rows=m.dim[0];
  int cols=m.dim[1];

    for(int x=0;x<rows;x++)
    {
        for(int y=0;y<cols;y++) 
        {
            std::cout<<m(x,y) << " ";
        }
    std::cout<<std::endl; 
    }

  return 0;
}


template <class Type>
inline
Type print(vector<Type> v){
  int len=v.size();

    for(int x=0;x<len;x++)
    {
      std::cout<< v(x) << " ";
    }
  std::cout<<std::endl; 
  return 0;
}