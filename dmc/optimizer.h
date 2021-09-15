class stochasticasticReconfigurationOptimizer
{
private:
  using vector_t=std::vector<real_t>;
  
  stochasticasticReconfigurationOptimizer(int nP);

  
  vector_t deltaStep();
  
private:

  void buildSMatrix(vector_t & estimators);
  void buildForce(vector_t & estimators);

  
  Eigen::MatrixXd Smatrix;
  
  vector_t delta;

}
