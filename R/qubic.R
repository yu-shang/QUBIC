.onLoad <- function(lib, pkg)
{
  library.dynam("QUBIC", pkg, lib);
}

qubic <- function(data = data(data), file = "rQUBIC", q = 0.06, c = 0.95, f = 1, k = 2, r = 1, o = 100, d = "F")
{
  data = as.matrix(data);
  rows = nrow(data);
  cols = ncol(data);
  rowsnames = rownames(data);
  colsnames = colnames(data);
  Bre = 1;
  r_d = 0;
  if(length(rowsnames) == 0)
    rowsnames = paste("gene", 1:rows, sep = "")
  if(length(colsnames) == 0)
    colsnames = paste("cond", 1:cols, sep = "")
  message = paste("Matrix contains", rows, "genes by", cols, "conditions");  #File %s 
  file = c(file,file);
  if(rows < 3 || cols < 3)
  {
    cat("neither rows number nor cols number can be too small\n");
    Bre = 0;
  }
  if( q > 0.5 || q <= 0 )
  {
    cat("-q quantile discretization should be (0,.5]\n");
    Bre = 0;
  }
  if( f > 1 || f < 0 )
  {
    cat("-f overlapping filtering should be [0,1.]\n");
    Bre = 0;
  }
  if( c > 1 || c <= 0.5 )
  {
    cat("-c noise ratio should be (0.5,1]\n");
    Bre = 0;
  }
  if( k < 2 && k != -1 )
  {
    cat("-k minimum column width should be >= 2\n");
    Bre = 0;
  }
  if( o <= 0 )
  {
    cat("-n number of blocks to report should be > 0\n");
    Bre = 0;
  }
  if (Bre == 0)
  {
    .C("r_puts");
    cat("\n[Usage]\nresults = qubic(data, [argument list]);\nor\nresults = qubic(data, file = 'rQUBIC', q = 0.06, c = 0.95, f = 1, k = 2, r = 1, o = 100, d = 'F');")
    return (NULL);
  }
  cat(message);cat("\n");
  if( (d == "t") || (d == "T") || (d == "True") || (d == "TRUE") || (d == "true") )
  {
    r_d = 1;
    d = "T";
  }
  else
  {
    d = "F";
  }
  BC = matrix(-1, o, 3);
  BC[,1] = 1:o;
  rownames(BC) = paste("BC", 1:o, sep = "");
  colnames(BC) = c("BC", "Genes", "Conds");
  results = .C( "r_main", data = as.double(data), rowsnames = as.character(rowsnames), colsnames = as.character(colsnames), rows = as.integer(rows), cols = as.integer(cols), file = as.character(file), q = as.double(q), c = as.double(c), f = as.double(f), k = as.integer(k), r = as.integer(r), o = as.integer(o), r_d = as.integer(r_d));
  message = paste("\n\nDiscretization rules are written to", paste(results$file[1], ".rules", sep = ""), "\nFormatted data are written to", paste(results$file[1], ".chars", sep = ""), "\nThe clusters are written to", paste(results$file[1], ".blocks", sep = ""), "\n"); 
  cat(message);
  results$data = matrix(results$data, rows, cols);
  rownames(results$data) = rowsnames;
  colnames(results$data) = colsnames;
  temp = matrix(c(results$q, results$c, results$f, results$k, results$r, results$o, d), , 7);
  colnames(temp) = c("-q", "-c", "-f", "-k", "-r", "-o", "-d");
  rownames(temp) = c("Argument List");
  rac = matrix(c(results$rows, results$cols), 1, 2);
  colnames(rac) = c("ROWS", "COLS");
  rownames(rac) = c("");
  filename = paste(results$file[1], ".blocks", sep = "");
  if(file.exists(filename) == "TRUE")
  {
    filename = c(filename,filename);
    BCresults = .C("cgetbc", BC = as.double(BC), o = as.integer(o), filename = as.character(filename), file = as.character(results$file));
    BCresults$BC = matrix(BCresults$BC, o, 3);
    rownames(BCresults$BC) = paste("BC", 1:o, sep = "");
    colnames(BCresults$BC) = c("BC", "Genes", "Conds");
    BCresults$BC = BCresults$BC[,2:3];
    for(times in 1:o)
    {
      if(BCresults$BC[times,1] == -1)
      {
        times = times - 1;
        BCresults$BC = BCresults$BC[1:times,];
        break;
      }
    }
    results = list(Discrete = results$data, ROWS_COLS = rac, BC = BCresults$BC, WriteToFile = results$file[1], ArgumentList = temp);
    message = paste("\nThe Biggest ");
    cat(message);
    thisBC = getBC(results);
    return (results);
  }
  else
  {
    results = list(Discrete = results$data, ROWS_COLS = rac, WriteToFile = results$file[1], ArgumentList = temp);
    return (results);
  }
}

getBC <- function(results, numBC = 1)
{
  filename = paste(results$WriteToFile, ".bc", sep = "");
  if( (length(results$BC) == 0) || (file.exists(filename) == "FALSE") )
  {
    cat("The first argument must be the result of QUBIC\n");
    return (NULL);
  }
  if(nrow(results$BC) < numBC)
  {
    cat("The BC you want to get does not exist\n");
    return (NULL);
  }
  nr = results$BC[numBC,1];
  nc = results$BC[numBC,2];
  if( (nr < 1) || (nc < 1) )
  {
    cat("QUBIC should have to run first\n");
    return (NULL);
  }
  allBC = readLines(filename, n = numBC*3);
  thisBC = list(genes = as.character(strsplit(allBC[numBC*3-1], "\t")[[1]]), conds = as.character(strsplit(allBC[numBC*3], "\t")[[1]]), num_gene = length(as.character(strsplit(allBC[numBC*3-1], "\t")[[1]])), num_cond = length(as.character(strsplit(allBC[numBC*3], "\t")[[1]])) )
  message = paste(" Bicluster contains", thisBC$num_gene, "genes and", thisBC$num_cond, "conditions\n"); 
  cat(message);
  return (thisBC);
}
