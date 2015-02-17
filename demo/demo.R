qdemo1 = function()
{
  if(require(Rcpp) != 1)
    install.packages("Rcpp");
  require(QUBIC);
  data(data);
  results = qubic(data);
  return (results);
}

qdemo2 = function()
{
  if(require(Rcpp) != 1)
    install.packages("Rcpp");
  require(QUBIC);
  data(data);
  results = qubic(data, file = "test");
  return (results);
}

qdemo3 = function()
{
  if(require(Rcpp) != 1)
    install.packages("Rcpp");
  require(QUBIC);
  data(data);
  results = qubic(data, file = "Example", q = 0.06, c = 0.95, f = 1, k = 2, r = 1, o = 100, d = "F");
  return (results);
}

qdemo4 = function()
{
  if(require(Rcpp) != 1)
    install.packages("Rcpp");
  require(QUBIC);
  data(toy);
  results = qubic(toy, file = "toy", d = "T");
  return (results);
}

qdemo5 = function()
{
  if(require(Rcpp) != 1)
    install.packages("Rcpp");
  require(QUBIC);
  data(data);
  results = qubic(data);
  thisBC = getBC(results, numBC = 3);
  return (thisBC);
}

qdemo6 = function()
{
  if(require(Rcpp) != 1)
    install.packages("Rcpp");
  require(QUBIC);
  data(toy);
  results = qubic(toy, file = "toy", d = "T");
  thisBC = getBC(results);
  return (thisBC);
}

cat("Use the function :\n\n");
cat("results = qdemo1();\n");
cat("results = qdemo2();\n");
cat("results = qdemo3();\n");
cat("results = qdemo4();\n");
cat("thisBC = qdemo5();\n");
cat("thisBC = qdemo6();\n\n");

demo = function()
{
  cat("Use the function :\n\n");
  cat("results = qdemo1();\n");
  cat("results = qdemo2();\n");
  cat("results = qdemo3();\n");
  cat("results = qdemo4();\n");
  cat("thisBC = qdemo5();\n");
  cat("thisBC = qdemo6();\n\n");
}